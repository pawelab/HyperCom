// screen_grids.cpp: bounds of the sender and receiver scores of cell x cell grids without ranking them, for the background pairs of
// prioritize_lr's p-values (R/screen.R)
//
// For each grid this gives every cell a lower and an upper bound of the sender and receiver scores that score_grids.cpp returns. Each
// entry Q = (a * b) / F of a within-sample block (the kernel's expression, so the same doubles, and the same skipped entries) is placed
// in a value bucket of its order-preserving key. With S_d the number of entries in the buckets below bucket d, an entry of bucket d has
// min rank - 1 in [S_d, S_{d+1} - 1] whatever the order inside the bucket. A score is the mean of its entries' (min rank - 1) / (K - 1),
// so the sums of these integer bounds over the same entries (exact in uint64), divided as the kernel divides, bound it up to the
// rounding of the divisions (~1e-16 relative; R/screen.R widens the bounds by 1e-12). Any partition of the keys into ordered buckets
// gives valid bounds; the buckets only set how tight they are.
//
// Buckets: from a sample of the entries (every SAMPLE-th entry of each block, at least SAMPLE_MIN entries in all, at most about
// 2^SAMPLE_BITS), the sampled key range
// is cut into 2^COARSE_BITS coarse cells, and each cell into the power of two of equal-width sub-buckets nearest above its sampled share
// of NB_TARGET buckets; at most NB_CAP buckets, so their ids fit uint16, and keys outside the sampled range go to two catch-alls.
// Locating a key reads one 8-byte coarse entry. Two passes per grid: the first counts the entries of each bucket and keeps each
// entry's bucket id (2 B), the second adds the bounds of the ids in position order. Teams of threads screen grids from one queue, one
// grid per team at a time; the members of a team share its ids and split its columns, with per-member counts and row sums merged
// (integers: the bounds do not depend on the threads or teams). Teams are as small as the ids of all teams at once fit id_budget
// bytes: one thread each for small grids, all threads in one team for big ones. Memory: 2 B per within-sample entry per team, beside
// the 8 B of the blocks, and 24 B per cell per thread for the row sums. The ids can be kept from one call to the next (keep_ids):
// made once, never zero-filled (the count pass writes every id before the sum pass reads it), freed by screen_release_cpp; new
// memory is slow to first touch, on Windows a batch of 60,000 cells screened 8-35% faster with kept ids. A thread that cannot start ends the run with an error (the
// barriers let the others go). Worker threads take the calling thread's SSE control word, as in score_grids.cpp. Do not compile
// with -ffast-math.

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include "grid_layout.h"
#include "r_alloc.h"
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#define HYPERCOM_X86 1
#endif

namespace {

// Tuning knobs: speed and tightness only, never a result
const size_t NB_TARGET = 32768;         // buckets aimed at (about 2x after rounding to powers of two, at most NB_CAP)
const int COARSE_BITS = 12;             // coarse cells of the sampled key range
const size_t SAMPLE = 61;               // bucket bounds from every 61st entry of a block
const size_t SAMPLE_MIN = (size_t) 1 << 16;  // or more often, so that a grid gives at least this many,
const int SAMPLE_BITS = 20;             // or less often, every (entries >> 20)-th (odd), when that is sparser: big grids sample faster
const size_t CHUNK = 256;               // rows of Q made at once
const size_t COLS = 8;                  // columns per unit of work
const unsigned SPIN = 100;              // busy-wait rounds before a waiting thread yields
// Hard limits
const size_t NB_CAP = 65533;            // buckets 1..nb and the catch-alls 0, nb + 1 stay below SKIP
const uint16_t SKIP = 0xFFFF;           // id of a skipped entry

typedef double v2d __attribute__((vector_size(16)));  // two doubles, one SSE2 instruction

// Bucket ids kept across calls: one buffer per team, grown when a call needs more, freed by screen_release_cpp
std::vector<std::unique_ptr<uint16_t[]>> kept_ids;
std::vector<size_t> kept_size;

unsigned get_mxcsr(){
#ifdef HYPERCOM_X86
  return _mm_getcsr();
#else
  return 0;
#endif
}

void set_mxcsr(unsigned m){
#ifdef HYPERCOM_X86
  _mm_setcsr(m);
#else
  (void) m;
#endif
}

inline uint64_t sortable_key(double v){
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | ((uint64_t) 1 << 63));
}

// sortable_key on the bits u of a value, the sign flip as a mask: ~u if negative, else u | sign bit
inline uint64_t flip_key(uint64_t u){
  return u ^ ((uint64_t) ((int64_t) u >> 63) | ((uint64_t) 1 << 63));
}

// The value with bits u is 0 or NaN: shifted left by one, 0 wraps to the top and NaN lies above infinity's 0xFFE0...0
inline bool zero_or_nan(uint64_t u){
  return (u << 1) - 1 >= (uint64_t) 0xFFE0000000000000ULL;
}

inline int bit_width(uint64_t x){
  return x ? 64 - __builtin_clzll(x) : 0;
}

inline bool skipped(size_t ii, size_t jj, double v){
  return ii == jj || v == 0 || std::isnan(v);
}

// Q[k] = (a[k] * bj) / f[k] for k < m, two at a time
inline void make_q(const double *a, double bj, const double *f, size_t m, double *q){
  const v2d b2 = {bj, bj};
  size_t k = 0;
  for(; k + 1 < m; k += 2){
    v2d av, fv;
    std::memcpy(&av, a + k, sizeof av);
    std::memcpy(&fv, f + k, sizeof fv);
    const v2d qv = (av * b2) / fv;
    std::memcpy(q + k, &qv, sizeof qv);
  }
  for(; k < m; ++k) q[k] = (a[k] * bj) / f[k];
}

// The blocks: cells of each, offsets, and units of work (COLS columns of one block)
struct Layout {
  std::vector<std::vector<int>> members;
  std::vector<uint64_t> size_of;        // the block size of each cell
  std::vector<size_t> off;
  size_t total = 0;
  struct Unit {
    size_t blk, j0, j1;
  };
  std::vector<Unit> units;
};

Layout make_layout(const int *sample, size_t n){
  Layout L;
  L.members = sample_blocks(sample, n);
  L.size_of.assign(n, 0);
  for(const auto &m : L.members){
    for(int k : m) L.size_of[k] = m.size();
  }
  for(size_t blk = 0; blk < L.members.size(); ++blk){
    const size_t s = L.members[blk].size();
    L.off.push_back(L.total);
    L.total += s * s;
    for(size_t j0 = 0; j0 < s; j0 += COLS) L.units.push_back({blk, j0, std::min(s, j0 + COLS)});
  }
  return L;
}

// a and b of each block, contiguous
struct GridData {
  std::vector<std::vector<double>> a, b;
};

void gather(const Layout &L, const double *a, const double *b, GridData &gd){
  gd.a.resize(L.members.size());
  gd.b.resize(L.members.size());
  for(size_t blk = 0; blk < L.members.size(); ++blk){
    const std::vector<int> &m = L.members[blk];
    gd.a[blk].resize(m.size());
    gd.b[blk].resize(m.size());
    for(size_t k = 0; k < m.size(); ++k){
      gd.a[blk][k] = a[m[k]];
      gd.b[blk][k] = b[m[k]];
    }
  }
}

struct Coarse {
  uint32_t base;                        // first bucket of the cell
  uint32_t shift;                       // sub-bucket = (key offset within the cell) >> shift
};

struct Buckets {
  uint64_t lo = 0, mask = 0;
  int cshift = 0;
  size_t ncoarse = 0, nb = 0;           // buckets 1..nb; 0 below the sampled range, nb + 1 above
  std::vector<Coarse> cell;
  inline size_t locate(uint64_t k) const {
    if(k < lo) return 0;
    const uint64_t x = k - lo, c = x >> cshift;
    if(c >= ncoarse) return nb + 1;
    const Coarse &e = cell[c];
    return e.base + ((x & mask) >> e.shift);
  }
};

// fn(key) for the sampled entries of a unit (every step-th entry of its block, counted from the block's first), those not skipped
template <class Fn> void sample_unit(const Layout &L, const double *F, const GridData &gd, const Layout::Unit &u, size_t step, Fn fn){
  const size_t s = L.members[u.blk].size();
  const double *Fb = F + L.off[u.blk], *a = gd.a[u.blk].data(), *b = gd.b[u.blk].data();
  for(size_t q = (u.j0 * s + step - 1) / step * step; q < u.j1 * s; q += step){
    const size_t jj = q / s, ii = q % s;
    const double v = (a[ii] * b[jj]) / Fb[q];
    if(skipped(ii, jj, v)) continue;
    fn(sortable_key(v));
  }
}

// The coarse cells of the sampled key range [lo, hi]; none if nothing was sampled (every entry then goes to catch-all 1)
void set_range(Buckets &bk, uint64_t lo, uint64_t hi){
  bk.nb = 0;
  bk.ncoarse = 0;
  bk.lo = 0;
  bk.cshift = 0;
  bk.mask = 0;
  if(lo > hi) return;
  const int R = bit_width(hi - lo);
  bk.lo = lo;
  bk.cshift = R > COARSE_BITS ? R - COARSE_BITS : 0;
  bk.mask = bk.cshift ? (((uint64_t) 1 << bk.cshift) - 1) : 0;
  bk.ncoarse = (size_t) ((hi - lo) >> bk.cshift) + 1;
}

// Sub-buckets of each coarse cell from its sampled count h[c]: the smallest power of two 2^t with 2^t * total >= h[c] * target, in
// integers (every team member builds its own table, so the tables must not depend on a thread's floating-point mode), the counts
// scaled down to fit 64 bits (target <= 2^15, so t <= 15); halve the target until at most NB_CAP buckets (target 1 gives one bucket
// per coarse cell)
void set_cells(Buckets &bk, const uint64_t *h){
  uint64_t total = 0;
  for(size_t c = 0; c < bk.ncoarse; ++c) total += h[c];
  bk.cell.resize(bk.ncoarse);
  if(bk.ncoarse == 0) return;
  const int scale = std::max(0, bit_width(total) + 17 - 64);
  const uint64_t tot = total >> scale;
  for(uint64_t target = NB_TARGET; ; target /= 2){
    size_t next = 1;
    for(size_t c = 0; c < bk.ncoarse; ++c){
      const uint64_t want = (h[c] >> scale) * target;
      int t = 0;
      while(t < bk.cshift && (tot << t) < want) ++t;
      bk.cell[c].base = (uint32_t) next;
      bk.cell[c].shift = (uint32_t) (bk.cshift - t);
      next += (size_t) 1 << t;
    }
    bk.nb = next - 1;
    if(bk.nb <= NB_CAP || target <= 1) break;
  }
}

// First pass over a unit: entries per bucket and the id of each entry
template <class S> void count_unit(const Layout &L, const double *F, const GridData &gd, const Buckets &bk, const Layout::Unit &u,
                                   S *cnt, uint16_t *ids){
  const size_t s = L.members[u.blk].size(), off = L.off[u.blk];
  const double *Fb = F + off, *a = gd.a[u.blk].data(), *b = gd.b[u.blk].data();
  double q[CHUNK];
  for(size_t jj = u.j0; jj < u.j1; ++jj){
    const double *f = Fb + jj * s;
    uint16_t *id = ids + off + jj * s;
    for(size_t i0 = 0; i0 < s; i0 += CHUNK){
      const size_t m = std::min(CHUNK, s - i0);
      make_q(a + i0, b[jj], f + i0, m, q);
      const size_t diag = jj - i0;     // wraps above m when the diagonal is not in these rows
      for(size_t k = 0; k < m; ++k){
        uint64_t v;
        std::memcpy(&v, q + k, sizeof v);
        if(zero_or_nan(v) || k == diag){
          id[i0 + k] = SKIP;
          continue;
        }
        const size_t d = bk.locate(flip_key(v));
        ++cnt[d];
        id[i0 + k] = (uint16_t) d;
      }
    }
  }
}

// Second pass over a unit: the bounds of each entry's rank added to its row and column (written as receiver bounds, out by cell, lo
// then hi). Two integer adds per entry and side: the high sum gets S_{d+1}, the - 1 per entry is taken off once per counted entry at
// the end, and a third sum counts the skipped entries (rare: the diagonal), so the counts come from the block size. A row's three
// sums sit together, rows[3 i .. 3 i + 2] (low, high, skipped): one cache line per entry, where three arrays of n each start on a
// page boundary and alias
template <class S> void sum_unit(const Layout &L, const Layout::Unit &u, const S *st, const uint16_t *ids, double denom, uint64_t *rows,
                                 double *rlo, double *rhi){
  const std::vector<int> &mem = L.members[u.blk];
  const size_t s = mem.size(), off = L.off[u.blk];
  for(size_t jj = u.j0; jj < u.j1; ++jj){
    const uint16_t *id = ids + off + jj * s;
    uint64_t cl = 0, ch = 0, skip = 0;
    for(size_t ii = 0; ii < s; ++ii){
      const uint16_t d = id[ii];
      uint64_t *r = rows + 3 * (size_t) mem[ii];
      if(d == SKIP){
        ++r[2];
        ++skip;
        continue;
      }
      const uint64_t l = st[d], h = st[d + 1];
      r[0] += l;
      r[1] += h;
      cl += l;
      ch += h;
    }
    const uint64_t cn = s - skip;
    rlo[mem[jj]] = (double) cl / ((double) cn * denom);
    rhi[mem[jj]] = (double) (ch - cn) / ((double) cn * denom);
  }
}

struct Spin {
  std::atomic<size_t> count{0}, generation{0};
};

// All team threads wait here until every one has arrived; false once the run has failed (a thread could not start), which also
// ends the wait of members whose team will never be complete
bool team_wait(Spin &b, size_t nt, const std::atomic<bool> &failed){
  const size_t gen = b.generation.load(std::memory_order_acquire);
  if(b.count.fetch_add(1, std::memory_order_acq_rel) + 1 == nt){
    b.count.store(0, std::memory_order_relaxed);
    b.generation.fetch_add(1, std::memory_order_acq_rel);
    return !failed.load();
  }
  unsigned spins = 0;
  while(b.generation.load(std::memory_order_acquire) == gen){
    if(failed.load(std::memory_order_relaxed)) return false;
    if(++spins > SPIN) std::this_thread::yield();
  }
  return !failed.load();
}

struct Input {
  const double *F, *A, *B;
  const int *lig, *rec;
  size_t nA, nB, n, G;
};

struct Output {
  double *slo, *shi, *rlo, *rhi;
};

// A team of threads that screens one grid at a time: grids from a shared queue, units of columns shared out among the members,
// counts and row sums per member, merged. Every member works out the buckets itself from the members' samples (no wait for one
// member). Every buffer is made before the threads start (the leader gathers each grid's a and b into buffers already sized, and
// each member's bucket table has room for every coarse cell), so no team thread allocates.
template <class S> struct Team {
  size_t size = 0;
  std::vector<uint16_t> own;            // the ids when not kept
  uint16_t *ids = nullptr;
  std::vector<std::vector<S>> cnt;
  std::vector<S> start;
  std::vector<std::vector<uint64_t>> rows, hist, h;     // rows: each member's row sums, 3 per cell (sum_unit)
  std::vector<Buckets> bk;
  std::vector<uint64_t> lo, hi;
  std::vector<double> b;
  GridData gd;
  std::atomic<size_t> unit[4];
  Spin bar;
  size_t g = 0;
  double denom = 0;

  Team(const Layout &L, const Input &in, size_t members, size_t t, bool keep) : size(members), cnt(members, std::vector<S>(NB_CAP + 2, 0)),
    start(NB_CAP + 3, 0), rows(members, std::vector<uint64_t>(3 * in.n, 0)),
    hist(members, std::vector<uint64_t>((size_t) 1 << COARSE_BITS, 0)),
    h(members, std::vector<uint64_t>((size_t) 1 << COARSE_BITS, 0)), bk(members), lo(members), hi(members), b(in.n) {
    for(auto &x : bk) x.cell.reserve((size_t) 1 << COARSE_BITS);
    if(keep){
      if(kept_ids.size() <= t){
        kept_ids.resize(t + 1);
        kept_size.resize(t + 1, 0);
      }
      if(kept_size[t] < L.total){
        kept_ids[t].reset();
        kept_size[t] = 0;
        kept_ids[t].reset(new uint16_t[L.total]);
        kept_size[t] = L.total;
      }
      ids = kept_ids[t].get();
    } else {
      own.resize(L.total);
      ids = own.data();
    }
    for(size_t k = 0; k < in.n; ++k) b[k] = in.B[(size_t) in.rec[0] + k * in.nB];
    gather(L, in.A + (size_t) in.lig[0] * in.nA, b.data(), gd);
    for(auto &u : unit) u.store(0);
  }
};

// Member me of team tm: screen grids until the queue is empty
template <class S> void run_member(const Layout &L, const Input &in, const Output &out, size_t step, Team<S> &tm, size_t me,
                                   std::atomic<size_t> &next, const std::atomic<bool> &failed){
  const size_t n = in.n, U = L.units.size(), k = tm.size;
  Buckets &bk = tm.bk[me];
  std::vector<uint64_t> &h = tm.h[me];
  for(;;){
    if(me == 0){
      tm.g = next.fetch_add(1);
      if(tm.g < in.G){
        for(size_t c = 0; c < n; ++c) tm.b[c] = in.B[(size_t) in.rec[tm.g] + c * in.nB];
        gather(L, in.A + (size_t) in.lig[tm.g] * in.nA, tm.b.data(), tm.gd);
      }
      for(auto &u : tm.unit) u.store(0);
    }
    if(!team_wait(tm.bar, k, failed)) return;
    const size_t g = tm.g;
    if(g >= in.G) break;

    uint64_t lo = ~(uint64_t) 0, hi = 0;
    for(size_t x = tm.unit[0].fetch_add(1); x < U; x = tm.unit[0].fetch_add(1)){
      sample_unit(L, in.F, tm.gd, L.units[x], step, [&](uint64_t key){ lo = std::min(lo, key); hi = std::max(hi, key); });
    }
    tm.lo[me] = lo;
    tm.hi[me] = hi;
    if(!team_wait(tm.bar, k, failed)) return;
    set_range(bk, *std::min_element(tm.lo.begin(), tm.lo.end()), *std::max_element(tm.hi.begin(), tm.hi.end()));
    uint64_t *hm = tm.hist[me].data();
    std::fill(hm, hm + ((size_t) 1 << COARSE_BITS), 0);
    for(size_t x = tm.unit[1].fetch_add(1); x < U; x = tm.unit[1].fetch_add(1)){
      sample_unit(L, in.F, tm.gd, L.units[x], step, [&](uint64_t key){ ++hm[(key - bk.lo) >> bk.cshift]; });
    }
    if(!team_wait(tm.bar, k, failed)) return;
    for(size_t c = 0; c < bk.ncoarse; ++c){
      h[c] = 0;
      for(size_t t = 0; t < k; ++t) h[c] += tm.hist[t][c];
    }
    set_cells(bk, h.data());

    S *cm = tm.cnt[me].data();
    for(size_t x = tm.unit[2].fetch_add(1); x < U; x = tm.unit[2].fetch_add(1)){
      count_unit(L, in.F, tm.gd, bk, L.units[x], cm, tm.ids);
    }
    if(!team_wait(tm.bar, k, failed)) return;
    const size_t nd = bk.nb + 2;
    for(size_t d = nd * me / k; d < nd * (me + 1) / k; ++d){
      S c = 0;
      for(size_t t = 0; t < k; ++t){
        c += tm.cnt[t][d];
        tm.cnt[t][d] = 0;
      }
      tm.start[d + 1] = c;
    }
    if(!team_wait(tm.bar, k, failed)) return;
    if(me == 0){
      tm.start[0] = 0;
      for(size_t d = 1; d < nd + 1; ++d) tm.start[d] += tm.start[d - 1];
      tm.denom = (double) tm.start[nd] - 1.0;
    }
    if(!team_wait(tm.bar, k, failed)) return;
    for(size_t x = tm.unit[3].fetch_add(1); x < U; x = tm.unit[3].fetch_add(1)){
      sum_unit(L, L.units[x], tm.start.data(), tm.ids, tm.denom, tm.rows[me].data(), out.rlo + g * n, out.rhi + g * n);
    }
    if(!team_wait(tm.bar, k, failed)) return;
    for(size_t i = n * me / k; i < n * (me + 1) / k; ++i){
      uint64_t l = 0, hh = 0, c = 0;
      for(size_t t = 0; t < k; ++t){
        uint64_t *r = tm.rows[t].data() + 3 * i;
        l += r[0];
        hh += r[1];
        c += r[2];
        r[0] = r[1] = r[2] = 0;
      }
      c = L.size_of[i] - c;
      hh -= c;
      out.slo[g * n + i] = (double) l / ((double) c * tm.denom);
      out.shi[g * n + i] = (double) hh / ((double) c * tm.denom);
    }
    if(!team_wait(tm.bar, k, failed)) return;
  }
}

// Teams of about k threads (sizes differing by at most one), at most one per grid, all screening grids from one queue
template <class S> void screen_teams(const Layout &L, const Input &in, const Output &out, size_t nt, size_t k, size_t step, bool keep){
  const size_t nteams = std::max<size_t>(1, std::min(nt / k, in.G));
  std::vector<std::unique_ptr<Team<S>>> teams;
  for(size_t t = 0; t < nteams; ++t) teams.emplace_back(new Team<S>(L, in, nt / nteams + (t < nt % nteams), t, keep));
  std::atomic<size_t> next(0);
  std::atomic<bool> failed(false);
  const unsigned mxcsr = get_mxcsr();
  std::vector<std::thread> pool;
  try {
    pool.reserve(nt);
    for(size_t t = 0; t < nteams; ++t){
      for(size_t me = 0; me < teams[t]->size; ++me){
        if(t == 0 && me == 0) continue;
        pool.emplace_back([&, t, me]{
          set_mxcsr(mxcsr);
          run_member(L, in, out, step, *teams[t], me, next, failed);
        });
      }
    }
  } catch(...){
    failed.store(true);
  }
  if(!failed.load()) run_member(L, in, out, step, *teams[0], 0, next, failed);
  for(auto &th : pool) th.join();
  if(failed.load()) throw std::runtime_error("could not start the threads");
}

}  // namespace

// The within-sample blocks of a dense influence matrix in the kernels' layout (grid_layout.h), copied from F without a temporary:
// cell holds the 0-based node of each cell, sample its sample (from 1)
// [[Rcpp::export]]
Rcpp::NumericVector dense_blocks_cpp(Rcpp::NumericMatrix F, Rcpp::IntegerVector cell, Rcpp::IntegerVector sample){
  const size_t n = cell.size(), N = F.nrow();
  if((size_t) sample.size() != n || (size_t) F.ncol() != N) Rcpp::stop("dense blocks: inputs of different sizes");
  for(size_t k = 0; k < n; ++k){
    if(cell[k] < 0 || (size_t) cell[k] >= N) Rcpp::stop("dense blocks: cell index out of range");
    if(sample[k] < 1) Rcpp::stop("dense blocks: sample ids must be positive integers");
  }
  const std::vector<std::vector<int>> members = sample_blocks(INTEGER(sample), n);
  size_t total = 0;
  for(const auto &m : members) total += m.size() * m.size();
  Rcpp::NumericVector out = alloc_vector<REALSXP>((R_xlen_t) total, false);
  const double *f = REAL(F);
  double *o = REAL(out);
  for(const auto &m : members){
    const size_t s = m.size();
    for(size_t jj = 0; jj < s; ++jj){
      const double *col = f + (size_t) cell[m[jj]] * N;
      for(size_t ii = 0; ii < s; ++ii) *o++ = col[cell[m[ii]]];
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::List screen_blocks_cpp(Rcpp::NumericVector blocks, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample,
                             Rcpp::IntegerVector lig, Rcpp::IntegerVector rec, int threads, double id_budget, bool keep_ids=false){
  const size_t n = A.nrow(), G = lig.size();
  if((size_t) B.ncol() != n || (size_t) sample.size() != n || rec.size() != lig.size()){
    Rcpp::stop("grid screening: inputs of different sizes");
  }
  for(size_t k = 0; k < n; ++k){
    if(sample[k] < 1) Rcpp::stop("grid screening: sample ids must be positive integers");
  }
  for(size_t g = 0; g < G; ++g){
    if(lig[g] < 0 || lig[g] >= A.ncol() || rec[g] < 0 || rec[g] >= B.nrow()) Rcpp::stop("grid screening: gene index out of range");
  }
  const Layout L = make_layout(INTEGER(sample), n);
  if((size_t) blocks.size() != L.total) Rcpp::stop("grid screening: blocks do not match the samples");
  size_t biggest = 0;
  for(const auto &m : L.members) biggest = std::max(biggest, m.size());
  if((double) biggest * (double) L.total >= 9.2e18) Rcpp::stop("grid screening: a sample too big for 64-bit rank sums");
  Rcpp::NumericMatrix slo = alloc_numeric_matrix(n, G, false), shi = alloc_numeric_matrix(n, G, false),
    rlo = alloc_numeric_matrix(n, G, false), rhi = alloc_numeric_matrix(n, G, false);
  if(G == 0 || n == 0){
    return Rcpp::List::create(Rcpp::Named("sender_lo") = slo, Rcpp::Named("sender_hi") = shi, Rcpp::Named("receiver_lo") = rlo,
                              Rcpp::Named("receiver_hi") = rhi);
  }
  const Input in = {REAL(blocks), REAL(A), REAL(B), INTEGER(lig), INTEGER(rec), (size_t) A.nrow(), (size_t) B.nrow(), n, G};
  const Output out = {REAL(slo), REAL(shi), REAL(rlo), REAL(rhi)};
  const size_t step = std::max(std::max<size_t>(1, std::min<size_t>(SAMPLE, L.total / SAMPLE_MIN)), (L.total >> SAMPLE_BITS) | 1);
  const int nt = std::max(1, threads);
#ifdef HYPERCOM_SCREEN_WIDE
  const bool narrow = false;            // test builds: the uint64 counts that grids of 2^32 entries or more use
#else
  const bool narrow = L.total < ((uint64_t) 1 << 32);
#endif
  const double ids = (double) L.total * sizeof(uint16_t);
  const size_t k = id_budget > 0 ? (size_t) std::min<double>(nt, std::max(1.0, std::ceil(nt * ids / id_budget))) : (size_t) nt;
  try {
    if(narrow) screen_teams<uint32_t>(L, in, out, nt, k, step, keep_ids);
    else screen_teams<uint64_t>(L, in, out, nt, k, step, keep_ids);
  } catch(const std::bad_alloc &){
    Rcpp::stop("grid screening: not enough memory for the bucket ids (" + std::to_string((size_t) ids >> 20) + " MB per team of threads)");
  } catch(const std::exception &){
    Rcpp::stop("grid screening: could not start " + std::to_string(nt) + " threads");
  }

  return Rcpp::List::create(Rcpp::Named("sender_lo") = slo, Rcpp::Named("sender_hi") = shi, Rcpp::Named("receiver_lo") = rlo,
                            Rcpp::Named("receiver_hi") = rhi);
}

// Free the bucket ids that screen_blocks_cpp kept (keep_ids)
// [[Rcpp::export]]
void screen_release_cpp(){
  kept_ids.clear();
  kept_size.clear();
}
