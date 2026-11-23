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
// Buckets: from a sample of the entries (every SAMPLE_STEP-th entry of each block, at least SAMPLE_MIN entries in all, at most about
// 2^SCREEN_SAMPLE_BITS), the sampled key range is cut into 2^COARSE_BITS coarse cells, and each cell into the power of two of
// equal-width sub-buckets nearest above its sampled share of NB_TARGET buckets; at most IDS_MAX buckets, so their ids fit uint16,
// and keys outside the sampled range go to two catch-alls.
// Locating a key reads one 8-byte coarse entry. Two passes per grid: the first counts the entries of each bucket and keeps each
// entry's bucket id (2 B), the second adds the bounds of the ids in position order. Teams of threads screen grids from one queue, one
// grid per team at a time; the members of a team share its ids and split its columns, with per-member counts and row sums merged
// (integers: the bounds do not depend on the threads or teams). Teams are as small as the ids of all teams at once fit id_budget
// bytes: one thread each for small grids, all threads in one team for big ones. Memory: 2 B per within-sample entry per team, beside
// the 8 B of the blocks, and 24 B per cell per thread for the row sums. The ids can be kept from one call to the next (keep_ids):
// made once, never zero-filled (the count pass writes every id before the sum pass reads it), freed by screen_release_cpp: new
// memory is slow to first touch, so batches of a big data set screen faster with kept ids. The threads run as threads.h runs them (a
// thread that cannot start, an error or an interrupt ends the run; the barriers let the others go). Worker threads take the calling
// thread's SSE control word, as in score_grids.cpp.
//
// Memory mode (screen_lambda_cpp): the same buckets, passes and teams on the integer keys of f2.h, A[i] + B'[j] + Λ[i, j], read from
// Λ's tiles instead of F, so the bounds hold for the memory-mode scores as the exact ones do for the exact scores. The keys are
// integers made from integers, so they do not depend on the threads either. Only the diagonal is skipped (a and b are positive and
// finite, so no key stands for 0 or NaN). The count pass takes a unit's columns (one column of tiles) a chunk of rows at a time: their
// Λ, gathered once from whole tiles into a strip, serves all of them. Memory: Λ's ~3 B per entry instead of the blocks' 8, so R gives
// an id budget of 2 B per entry, one team for all threads, one id buffer.

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
#include "f2.h"
#include "grid_common.h"
#include "r_alloc.h"
#include "threads.h"

namespace {

// Tuning knobs: speed and tightness only, never a result
const size_t NB_TARGET = 32768;         // buckets aimed at (about 2x after rounding to powers of two, at most IDS_MAX)
const int COARSE_BITS = 12;             // coarse cells of the sampled key range
const size_t SAMPLE_MIN = (size_t) 1 << 16;  // bucket bounds from every SAMPLE_STEP-th entry of a block, or more often, so that a grid
                                        // gives at least this many,
const int SCREEN_SAMPLE_BITS = 20;      // or less often, every (entries >> 20)-th (odd), when that is sparser: big grids sample faster
const size_t CHUNK = 256;               // rows of Q made at once
const size_t COLS = 8;                  // columns per unit of work

// Bucket ids kept across calls: one buffer per team, grown when a call needs more, freed by screen_release_cpp
std::vector<std::unique_ptr<uint16_t[]>> kept_ids;
std::vector<size_t> kept_size;

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
  Layout lay;
  lay.members = sample_blocks(sample, n);
  lay.size_of.assign(n, 0);
  for(const auto &m : lay.members){
    for(int k : m) lay.size_of[k] = m.size();
  }
  for(size_t blk = 0; blk < lay.members.size(); ++blk){
    const size_t s = lay.members[blk].size();
    lay.off.push_back(lay.total);
    lay.total += s * s;
    for(size_t j0 = 0; j0 < s; j0 += COLS) lay.units.push_back({blk, j0, std::min(s, j0 + COLS)});
  }
  return lay;
}

// a and b of each block, contiguous; in memory mode also their keys (tile_keys)
struct GridData {
  std::vector<std::vector<double>> a, b;
  std::vector<std::vector<int64_t>> ka, kb;
};

void gather(const Layout &lay, const double *a, const double *b, GridData &gd){
  gd.a.resize(lay.members.size());
  gd.b.resize(lay.members.size());
  for(size_t blk = 0; blk < lay.members.size(); ++blk){
    const std::vector<int> &m = lay.members[blk];
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

struct ScreenBuckets {
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
template <class Fn> void sample_unit(const Layout &lay, const double *F, const GridData &gd, const Layout::Unit &u, size_t step, Fn fn){
  const size_t s = lay.members[u.blk].size();
  const double *Fb = F + lay.off[u.blk], *a = gd.a[u.blk].data(), *b = gd.b[u.blk].data();
  for(size_t q = (u.j0 * s + step - 1) / step * step; q < u.j1 * s; q += step){
    const size_t jj = q / s, ii = q % s;
    const double v = (a[ii] * b[jj]) / Fb[q];
    if(skipped(ii, jj, v)) continue;
    fn(sortable_key(v));
  }
}

// The coarse cells of the sampled key range [lo, hi]; none if nothing was sampled (every entry then goes to catch-all 1)
void set_range(ScreenBuckets &bk, uint64_t lo, uint64_t hi){
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
// scaled down to fit 64 bits (target <= 2^15, so t <= 15); halve the target until at most IDS_MAX buckets (target 1 gives one bucket
// per coarse cell)
void set_cells(ScreenBuckets &bk, const uint64_t *h){
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
    if(bk.nb <= IDS_MAX || target <= 1) break;
  }
}

// First pass over a unit: entries per bucket and the id of each entry
template <class Count> void count_unit(const Layout &lay, const double *F, const GridData &gd, const ScreenBuckets &bk, const Layout::Unit &u,
                                       Count *cnt, uint16_t *ids){
  const size_t s = lay.members[u.blk].size(), off = lay.off[u.blk];
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
template <class Count> void sum_unit(const Layout &lay, const Layout::Unit &u, const Count *st, const uint16_t *ids, double denom, uint64_t *rows,
                                     double *rlo, double *rhi){
  const std::vector<int> &mem = lay.members[u.blk];
  const size_t s = mem.size(), off = lay.off[u.blk];
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

struct Input {
  const double *F, *A, *B;
  const int *lig, *rec;
  size_t nA, nB, n, G;
  // Memory mode: Λ's tiles (f2.h) instead of F, which is then unused
  const char *lam = nullptr;
  std::vector<int64_t> base{};          // per block: Λ's base
  std::vector<size_t> toff{};           // per block: byte offset of its tiles
  std::vector<std::vector<int64_t>> lD{};  // per block: fixlog52(D) of its cells
  int sh = 0;                           // 52 - p: the keys' rounding shift
};

// Memory mode: the keys of each block's a and b, as score_grids.cpp makes them: ka = fix(ln a), kb = fix(ln b - ln D) + Λ's base
void tile_keys(const Layout &lay, const Input &in, GridData &gd){
  gd.ka.resize(lay.members.size());
  gd.kb.resize(lay.members.size());
  for(size_t blk = 0; blk < lay.members.size(); ++blk){
    const size_t s = lay.members[blk].size();
    gd.ka[blk].resize(s);
    gd.kb[blk].resize(s);
    for(size_t k = 0; k < s; ++k){
      gd.ka[blk][k] = f2::round_shift(f2::fixlog52(gd.a[blk][k]), in.sh);
      gd.kb[blk][k] = f2::round_shift(f2::fixlog52(gd.b[blk][k]) - in.lD[blk][k], in.sh) + in.base[blk];
    }
  }
}

// Memory mode: sample_unit on the keys
template <class Fn> void sample_tiles(const Layout &lay, const Input &in, const GridData &gd, const Layout::Unit &u, size_t step, Fn fn){
  const size_t s = lay.members[u.blk].size();
  const int64_t *ka = gd.ka[u.blk].data(), *kb = gd.kb[u.blk].data();
  const char *T = in.lam + in.toff[u.blk];
  for(size_t q = (u.j0 * s + step - 1) / step * step; q < u.j1 * s; q += step){
    const size_t jj = q / s, ii = q % s;
    if(ii == jj) continue;
    fn((uint64_t) (ka[ii] + kb[jj] + f2::tile_value(T, ii, jj)) ^ SIGN);
  }
}

// Memory mode: count_unit on the keys. Row chunks outer: a unit's columns are one column of tiles, so each chunk's strip of their
// Λ (column-major, from whole tiles; those below the diagonal tiles transposed) is gathered once for all of them. Ids by position
// and counts do not depend on the order.
template <class Count> void count_tiles(const Layout &lay, const Input &in, const GridData &gd, const ScreenBuckets &bk, const Layout::Unit &u,
                                        Count *cnt, uint16_t *ids){
  static_assert(COLS == f2::TILE, "a unit of work is one column of tiles");
  const size_t s = lay.members[u.blk].size(), off = lay.off[u.blk], J = u.j0 / f2::TILE, ncol = u.j1 - u.j0;
  const int64_t *ka = gd.ka[u.blk].data(), *kb = gd.kb[u.blk].data();
  const char *T = in.lam + in.toff[u.blk];
  int64_t strip[COLS * CHUNK];
  for(size_t i0 = 0; i0 < s; i0 += CHUNK){
    const size_t m = std::min(CHUNK, s - i0);
    f2::tile_strip(T, J, ncol, i0, i0 + m, strip);
    for(size_t jj = u.j0; jj < u.j1; ++jj){
      const int64_t *col = strip + (jj - u.j0) * m, bj = kb[jj];
      uint16_t *id = ids + off + jj * s + i0;
      const size_t diag = jj - i0;     // wraps above m when the diagonal is not in these rows
      for(size_t k = 0; k < m; ++k){
        if(k == diag){
          id[k] = SKIP;
          continue;
        }
        const size_t d = bk.locate((uint64_t) (ka[i0 + k] + bj + col[k]) ^ SIGN);
        ++cnt[d];
        id[k] = (uint16_t) d;
      }
    }
  }
}

struct Output {
  double *slo, *shi, *rlo, *rhi;
};

// A team of threads that screens one grid at a time: grids from a shared queue, units of columns shared out among the members,
// counts and row sums per member, merged. Every member works out the buckets itself from the members' samples (no wait for one
// member). Every buffer is made before the threads start (the leader gathers each grid's a and b into buffers already sized, and
// each member's bucket table has room for every coarse cell), so no team thread allocates.
template <class Count> struct ScreenTeam {
  size_t size = 0;
  std::vector<uint16_t> own;            // the ids when not kept
  uint16_t *ids = nullptr;
  std::vector<std::vector<Count>> cnt;
  std::vector<Count> start;
  std::vector<std::vector<uint64_t>> rows, hist, h;     // rows: each member's row sums, 3 per cell (sum_unit)
  std::vector<ScreenBuckets> bk;
  std::vector<uint64_t> lo, hi;
  std::vector<double> b;
  GridData gd;
  std::atomic<size_t> unit[4];
  hc::Barrier bar;
  size_t g = 0;
  double denom = 0;

  ScreenTeam(const Layout &lay, const Input &in, size_t members, size_t t, bool keep) : size(members),
    cnt(members, std::vector<Count>(IDS_MAX + 2, 0)), start(IDS_MAX + 3, 0), rows(members, std::vector<uint64_t>(3 * in.n, 0)),
    hist(members, std::vector<uint64_t>((size_t) 1 << COARSE_BITS, 0)),
    h(members, std::vector<uint64_t>((size_t) 1 << COARSE_BITS, 0)), bk(members), lo(members), hi(members), b(in.n), bar(members) {
    for(auto &x : bk) x.cell.reserve((size_t) 1 << COARSE_BITS);
    if(keep){
      if(kept_ids.size() <= t){
        kept_ids.resize(t + 1);
        kept_size.resize(t + 1, 0);
      }
      if(kept_size[t] < lay.total){
        kept_ids[t].reset();
        kept_size[t] = 0;
        kept_ids[t].reset(new uint16_t[lay.total]);
        kept_size[t] = lay.total;
      }
      ids = kept_ids[t].get();
    } else {
      own.resize(lay.total);
      ids = own.data();
    }
    for(size_t k = 0; k < in.n; ++k) b[k] = in.B[(size_t) in.rec[0] + k * in.nB];
    gather(lay, in.A + (size_t) in.lig[0] * in.nA, b.data(), gd);
    if(in.lam) tile_keys(lay, in, gd);
    for(auto &u : unit) u.store(0);
  }
};

// Member me of team tm: screen grids until the queue is empty
template <class Count> void run_member(const Layout &lay, const Input &in, const Output &out, size_t step, ScreenTeam<Count> &tm, size_t me,
                                       std::atomic<size_t> &next, const std::atomic<bool> &failed){
  const size_t n = in.n, U = lay.units.size(), k = tm.size;
  ScreenBuckets &bk = tm.bk[me];
  std::vector<uint64_t> &h = tm.h[me];
  for(;;){
    if(me == 0){
      tm.g = next.fetch_add(1);
      if(tm.g < in.G){
        for(size_t c = 0; c < n; ++c) tm.b[c] = in.B[(size_t) in.rec[tm.g] + c * in.nB];
        gather(lay, in.A + (size_t) in.lig[tm.g] * in.nA, tm.b.data(), tm.gd);
        if(in.lam) tile_keys(lay, in, tm.gd);
      }
      for(auto &u : tm.unit) u.store(0);
    }
    if(!tm.bar.wait(failed)) return;
    const size_t g = tm.g;
    if(g >= in.G) break;

    uint64_t lo = ~(uint64_t) 0, hi = 0;
    auto range = [&](uint64_t key){ lo = std::min(lo, key); hi = std::max(hi, key); };
    for(size_t x = tm.unit[0].fetch_add(1); x < U; x = tm.unit[0].fetch_add(1)){
      if(in.lam) sample_tiles(lay, in, tm.gd, lay.units[x], step, range);
      else sample_unit(lay, in.F, tm.gd, lay.units[x], step, range);
    }
    tm.lo[me] = lo;
    tm.hi[me] = hi;
    if(!tm.bar.wait(failed)) return;
    set_range(bk, *std::min_element(tm.lo.begin(), tm.lo.end()), *std::max_element(tm.hi.begin(), tm.hi.end()));
    uint64_t *hm = tm.hist[me].data();
    std::fill(hm, hm + ((size_t) 1 << COARSE_BITS), 0);
    auto tally = [&](uint64_t key){ ++hm[(key - bk.lo) >> bk.cshift]; };
    for(size_t x = tm.unit[1].fetch_add(1); x < U; x = tm.unit[1].fetch_add(1)){
      if(in.lam) sample_tiles(lay, in, tm.gd, lay.units[x], step, tally);
      else sample_unit(lay, in.F, tm.gd, lay.units[x], step, tally);
    }
    if(!tm.bar.wait(failed)) return;
    for(size_t c = 0; c < bk.ncoarse; ++c){
      h[c] = 0;
      for(size_t t = 0; t < k; ++t) h[c] += tm.hist[t][c];
    }
    set_cells(bk, h.data());

    Count *cm = tm.cnt[me].data();
    for(size_t x = tm.unit[2].fetch_add(1); x < U; x = tm.unit[2].fetch_add(1)){
      if(in.lam) count_tiles(lay, in, tm.gd, bk, lay.units[x], cm, tm.ids);
      else count_unit(lay, in.F, tm.gd, bk, lay.units[x], cm, tm.ids);
    }
    if(!tm.bar.wait(failed)) return;
    const size_t nd = bk.nb + 2;
    for(size_t d = nd * me / k; d < nd * (me + 1) / k; ++d){
      Count c = 0;
      for(size_t t = 0; t < k; ++t){
        c += tm.cnt[t][d];
        tm.cnt[t][d] = 0;
      }
      tm.start[d + 1] = c;
    }
    if(!tm.bar.wait(failed)) return;
    if(me == 0){
      tm.start[0] = 0;
      for(size_t d = 1; d < nd + 1; ++d) tm.start[d] += tm.start[d - 1];
      tm.denom = (double) tm.start[nd] - 1.0;
    }
    if(!tm.bar.wait(failed)) return;
    for(size_t x = tm.unit[3].fetch_add(1); x < U; x = tm.unit[3].fetch_add(1)){
      sum_unit(lay, lay.units[x], tm.start.data(), tm.ids, tm.denom, tm.rows[me].data(), out.rlo + g * n, out.rhi + g * n);
    }
    if(!tm.bar.wait(failed)) return;
    for(size_t i = n * me / k; i < n * (me + 1) / k; ++i){
      uint64_t l = 0, hh = 0, c = 0;
      for(size_t t = 0; t < k; ++t){
        uint64_t *r = tm.rows[t].data() + 3 * i;
        l += r[0];
        hh += r[1];
        c += r[2];
        r[0] = r[1] = r[2] = 0;
      }
      c = lay.size_of[i] - c;
      hh -= c;
      out.slo[g * n + i] = (double) l / ((double) c * tm.denom);
      out.shi[g * n + i] = (double) hh / ((double) c * tm.denom);
    }
    if(!tm.bar.wait(failed)) return;
  }
}

// Teams of about k threads (sizes differing by at most one), at most one per grid, all screening grids from one queue
template <class Count> void screen_teams(const Layout &lay, const Input &in, const Output &out, size_t nt, size_t k, size_t step, bool keep){
  const size_t nteams = std::max<size_t>(1, std::min(nt / k, in.G));
  std::vector<std::unique_ptr<ScreenTeam<Count>>> teams;
  std::vector<std::pair<size_t, size_t>> member;  // each thread's team and place in it
  try {
    for(size_t t = 0; t < nteams; ++t){
      teams.emplace_back(new ScreenTeam<Count>(lay, in, nt / nteams + (t < nt % nteams), t, keep));
      for(size_t me = 0; me < teams[t]->size; ++me) member.emplace_back(t, me);
    }
  } catch(const std::bad_alloc &){
    Rcpp::stop("grid screening: not enough memory for the bucket ids (" + std::to_string(lay.total * sizeof(uint16_t) >> 20) +
               " MB per team of threads)");
  }
  std::atomic<size_t> next(0);
  hc::Stop stop;
  const unsigned mxcsr = get_mxcsr();
  hc::run_threads(member.size(), stop, [&](size_t w){
    set_mxcsr(mxcsr);
    run_member(lay, in, out, step, *teams[member[w].first], member[w].second, next, stop.failed);
  });
  stop.raise("grid screening");
}

// The checks of the inputs both screening exports take
void check_screen(const Rcpp::NumericMatrix &A, const Rcpp::NumericMatrix &B, const Rcpp::IntegerVector &sample,
                  const Rcpp::IntegerVector &lig, const Rcpp::IntegerVector &rec){
  const size_t n = A.nrow(), G = lig.size();
  if((size_t) B.ncol() != n || (size_t) sample.size() != n || rec.size() != lig.size()){
    Rcpp::stop("grid screening: inputs of different sizes");
  }
  check_samples(sample, "grid screening");
  for(size_t g = 0; g < G; ++g){
    if(lig[g] < 0 || lig[g] >= A.ncol() || rec[g] < 0 || rec[g] >= B.nrow()) Rcpp::stop("grid screening: gene index out of range");
  }
}

// Screen every grid of in (exact or memory mode) and return the four bound matrices, sender_lo, sender_hi, receiver_lo, receiver_hi
Rcpp::List screen_run(const Layout &lay, const Input &in, int threads, double id_budget, bool keep_ids){
  const size_t n = in.n, G = in.G;
  size_t biggest = 0;
  for(const auto &m : lay.members) biggest = std::max(biggest, m.size());
  if((double) biggest * (double) lay.total >= RANK_SUM_MAX) Rcpp::stop("grid screening: a sample too big for 64-bit rank sums");
  Rcpp::NumericMatrix slo = alloc_numeric_matrix(n, G, false), shi = alloc_numeric_matrix(n, G, false),
    rlo = alloc_numeric_matrix(n, G, false), rhi = alloc_numeric_matrix(n, G, false);
  if(G == 0 || n == 0){
    return Rcpp::List::create(Rcpp::Named("sender_lo") = slo, Rcpp::Named("sender_hi") = shi, Rcpp::Named("receiver_lo") = rlo,
                              Rcpp::Named("receiver_hi") = rhi);
  }
  const Output out = {REAL(slo), REAL(shi), REAL(rlo), REAL(rhi)};
  const size_t step = std::max(std::max<size_t>(1, std::min<size_t>(SAMPLE_STEP, lay.total / SAMPLE_MIN)), (lay.total >> SCREEN_SAMPLE_BITS) | 1);
  const int nt = std::max(1, threads);
#ifdef HYPERCOM_SCREEN_WIDE
  const bool narrow = false;            // test builds: the uint64 counts that grids of 2^32 entries or more use
#else
  const bool narrow = lay.total < ((uint64_t) 1 << 32);
#endif
  const double ids = (double) lay.total * sizeof(uint16_t);
  const size_t k = id_budget > 0 ? (size_t) std::min<double>(nt, std::max(1.0, std::ceil(nt * ids / id_budget))) : (size_t) nt;
  if(narrow) screen_teams<uint32_t>(lay, in, out, nt, k, step, keep_ids);
  else screen_teams<uint64_t>(lay, in, out, nt, k, step, keep_ids);

  return Rcpp::List::create(Rcpp::Named("sender_lo") = slo, Rcpp::Named("sender_hi") = shi, Rcpp::Named("receiver_lo") = rlo,
                            Rcpp::Named("receiver_hi") = rhi);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List screen_blocks_cpp(Rcpp::NumericVector blocks, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample,
                             Rcpp::IntegerVector lig, Rcpp::IntegerVector rec, int threads, double id_budget, bool keep_ids=false){
  check_screen(A, B, sample, lig, rec);
  const size_t n = A.nrow(), G = lig.size();
  const Layout lay = make_layout(INTEGER(sample), n);
  if((size_t) blocks.size() != lay.total) Rcpp::stop("grid screening: blocks do not match the samples");
  const Input in = {REAL(blocks), REAL(A), REAL(B), INTEGER(lig), INTEGER(rec), (size_t) A.nrow(), (size_t) B.nrow(), n, G};

  return screen_run(lay, in, threads, id_budget, keep_ids);
}

// Memory mode: screen_blocks_cpp on the keys of f2.h, from Λ's tiles (bytes, base, toff as influence_lambda_cpp returns them, made at
// p fractional bits), D the adjacency's row sums of the cells in the order of sample; A and B positive and finite (else R screens
// the blocks). Bounds of the memory-mode scores; ids kept as screen_blocks_cpp keeps them, in the same buffers
// [[Rcpp::export]]
Rcpp::List screen_lambda_cpp(Rcpp::RawVector bytes, Rcpp::NumericVector base, Rcpp::NumericVector toff, Rcpp::NumericVector D,
                             Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample, Rcpp::IntegerVector lig,
                             Rcpp::IntegerVector rec, int threads, double id_budget, bool keep_ids, int p){
  check_screen(A, B, sample, lig, rec);
  const size_t n = A.nrow(), G = lig.size();
  const Layout lay = make_layout(INTEGER(sample), n);
  const size_t nblocks = lay.members.size();
  std::vector<size_t> sizes(nblocks);
  for(size_t blk = 0; blk < nblocks; ++blk) sizes[blk] = lay.members[blk].size();
  check_lambda(bytes, base, toff, D, A, B, sizes, n, p, "grid screening");
  Input in = {nullptr, REAL(A), REAL(B), INTEGER(lig), INTEGER(rec), (size_t) A.nrow(), (size_t) B.nrow(), n, G};
  in.lam = reinterpret_cast<const char *>(RAW(bytes));
  in.sh = 52 - p;
  in.base.resize(nblocks);
  if(nblocks > 0) std::memcpy(in.base.data(), REAL(base), nblocks * sizeof(int64_t));  // int64 bits in doubles
  in.toff.resize(nblocks);
  in.lD.resize(nblocks);
  for(size_t blk = 0; blk < nblocks; ++blk){
    const std::vector<int> &m = lay.members[blk];
    in.toff[blk] = (size_t) toff[blk];
    in.lD[blk].resize(m.size());
    for(size_t k = 0; k < m.size(); ++k) in.lD[blk][k] = f2::fixlog52(D[m[k]]);
  }

  return screen_run(lay, in, threads, id_budget, keep_ids);
}

// Free the bucket ids that screen_blocks_cpp and screen_lambda_cpp kept (keep_ids)
// [[Rcpp::export]]
void screen_release_cpp(){
  kept_ids.clear();
  kept_size.clear();
}
