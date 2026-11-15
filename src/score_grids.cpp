// score_grids.cpp: sender and receiver scores of cell x cell grids, bit-identical to master's score_sr
//
// For a ligand l and one receptor subunit r, master builds the cell x cell grid
//   Q[i, j] = (F[i, l] * F[r, j]) / F[i, j]      (F = influence matrix; Q[i, j] = F[i, j] for the network score)
// sets cross-sample entries, the diagonal, zeros and NaN to NA, percent-ranks the remaining K entries together
// ((min rank - 1) / (K - 1)) and averages the ranks over j (sender of i, colMeans) and over i (receiver of j, rowMeans).
//
// Here only the within-sample blocks are visited (every cross-sample entry is NA). Entries are ranked by an
// order-preserving integer image of their value: each entry falls in one of 2^11 (or more) value buckets in key order (bounds
// from a sample of the entries: equal width, or for a grid above the cache budget equal counts from a histogram of the sample,
// so every bucket sorts in L2; plus two catch-all buckets) and its key bits go to its bucket, then each bucket is sorted in
// cache as packed (key bits, index) words and its entries get their min ranks in place. Any partition in key order gives the
// same ranks. The means are accumulated in long double in master's order, so the results equal master's bit for bit. A team
// of threads scores one grid at a time; several teams run at once when their grids fit the cache together, else one team
// takes all threads and the row ranges of a big block sum as a pipeline. Memory per within-sample entry per team: its bucket id
// (2 B) and its key bits, 8 B for a grid that fits the cache, 6 B for a grid above the cache budget (a bucket then keeps at most 48
// key bits; the few entries of the catch-alls keep 8 B), plus 8 B per entry of samples whose cells are not consecutive in F while
// those copies fit the cache, plus each thread's sort buffers (16 B per entry of the biggest bucket it ranks). Do not compile with
// -ffast-math.
//
// score_grids_cpp reads F (a dense n x n influence matrix); score_blocks_cpp is given instead every within-sample block of F
// (made from the factored influence matrix, see influence.cpp) and F[cell, ligand] / F[receptor, cell] for the genes of the pairs,
// so no n x n matrix exists. Both rank and average the same values the same way.

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#define HYPERCOM_X86 1
#endif

namespace {

struct FpuState {
  unsigned short x87;
  unsigned int mxcsr;
};

FpuState get_fpu(){
  FpuState s = {0, 0};
#ifdef HYPERCOM_X86
  __asm__ volatile("fnstcw %0" : "=m"(s.x87));
  s.mxcsr = _mm_getcsr();
#endif
  return s;
}

void set_fpu(const FpuState &s){
#ifdef HYPERCOM_X86
  __asm__ volatile("fldcw %0" : : "m"(s.x87));
  _mm_setcsr(s.mxcsr);
#else
  (void) s;
#endif
}

typedef double v2d __attribute__((vector_size(16)));  // two doubles, one SSE2 instruction

const size_t SLAB = 64;                 // entries per pool slab on one thread: each bucket is a linked list of slabs
const uint16_t SKIP = 0xFFFF;           // bucket of a skipped entry
const int B1_MIN = 11, B1_MAX = 14;    // value buckets per grid: 2^b1 plus two catch-alls, b1 grows with the grid
const size_t CACHE_BUDGET = (size_t) 24 << 20;   // grid buffers that may be in flight at once; a bigger grid gets all threads
const size_t KEEP_SORT_BUFFER = (size_t) 1 << 20;  // per-thread sort buffers above this many entries are freed after each grid
const size_t PIPE = 8;                  // columns a pipeline stage hands on at once
const size_t SAMPLE = 61;               // bucket bounds from every 61st entry
const size_t RUN = 32;                  // a block whose nodes come in runs this long on average is read from F, others are copied
                                        // while all copies fit CACHE_BUDGET, else gathered from F
const int SLABS = 0, IDS = 1, SIX = 2;   // how a team stores entries: slabs, runs of 8-byte keys, runs of 6-byte keys
const size_t BUCKETS_IDS = 65533;        // bucket ids are uint16 below SKIP
const int NARROW_BITS = 48;             // key bits of a 6-byte pool slot: a bucket keeps at most this many
const int TABLE_BITS = 18;              // big grids: the sampled key range is cut into 2^18 cells whose counts set the buckets
const size_t BUCKET_TARGET = 65536;     // big grids: entries per bucket aimed at, so a bucket's sort buffers stay in L2
const size_t BUCKETS_MIN = 2048, BUCKETS_MAX = 32768;

inline uint64_t sortable_key(double v){
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | ((uint64_t) 1 << 63));
}

inline int bit_width(uint64_t x){
  return x ? 64 - __builtin_clzll(x) : 0;
}

inline bool skipped(size_t ii, size_t jj, double v){
  return ii == jj || v == 0 || std::isnan(v);
}

void insertion_sort(uint64_t *x, size_t m){
  for(size_t i = 1; i < m; ++i){
    const uint64_t v = x[i];
    size_t j = i;
    while(j > 0 && x[j - 1] > v){
      x[j] = x[j - 1];
      --j;
    }
    x[j] = v;
  }
}

// Sort packed words (key bits << P | index) into y, x is scratch. c: nb + 1 counts of the top digit (key bits >> shift),
// shifted by one, from the caller, or nullptr to take the key range from the data. One radix pass, deeper passes for
// crowded sub-buckets, then one insertion sort over the nearly sorted result.
void sort_words(uint64_t *x, uint64_t *y, size_t m, int P, int shift, size_t nb, uint32_t *c, std::vector<uint32_t> &cnt){
  uint64_t lo = 0;
  const size_t base = cnt.size();
  if(!c){
    uint64_t hi = lo = x[0] >> P;
    for(size_t p = 1; p < m; ++p){
      const uint64_t k = x[p] >> P;
      lo = std::min(lo, k);
      hi = std::max(hi, k);
    }
    if(lo == hi){
      std::memcpy(y, x, m * sizeof(uint64_t));
      return;
    }
    const int bits = bit_width(hi - lo), want = std::min(12, bit_width(m));
    shift = bits > want ? bits - want : 0;
    nb = (size_t) ((hi - lo) >> shift) + 1;
    cnt.resize(base + nb + 1);
    c = cnt.data() + base;
    std::fill(c, c + nb + 1, 0u);
    for(size_t p = 0; p < m; ++p) ++c[(((x[p] >> P) - lo) >> shift) + 1];
  }
  const size_t coff = c - cnt.data();
  uint32_t most = 0;
  for(size_t b = 0; b < nb; ++b){
    most = std::max(most, c[b + 1]);
    c[b + 1] += c[b];
  }
  for(size_t p = 0; p < m; ++p) y[c[((x[p] >> P) - lo) >> shift]++] = x[p];
  if(shift > 0){
    if(most > 32){
      for(size_t b = 0, s = 0; b < nb; ++b){
        const size_t e = cnt[coff + b];
        if(e - s > 32){
          std::memcpy(x + s, y + s, (e - s) * sizeof(uint64_t));
          sort_words(x + s, y + s, e - s, P, 0, 0, nullptr, cnt);
        }
        s = e;
      }
    }
    insertion_sort(y, m);
  }
  cnt.resize(base);
}

struct Cells {
  const double *Fb = nullptr;          // F[i, j] of the entries of the copied (or given) blocks, in position order
  std::vector<double> Fown;            // the copies, when made here
  const double *F = nullptr;           // the whole influence matrix (L x L), or null when every block is given
  const double *A = nullptr;           // without F: F[cell, ligand], cells x ligands, rows in cell order
  const double *B = nullptr;           // and F[receptor, cell], receptors x cells (nb rows)
  size_t nb = 0;
  size_t L = 0;
  std::vector<int> node;
  std::vector<int> out;
  std::vector<size_t> block_start;
  std::vector<size_t> block_off;
  std::vector<size_t> block_fb;        // offset of a copied block in Fb, SIZE_MAX for a block read from F
  std::vector<bool> block_gather;      // a block read from F one node at a time
  std::vector<size_t> run;             // per block read from F in runs: the first row of each run of consecutive nodes, then s
  std::vector<size_t> block_run;
  size_t cap;
  size_t max_block;
};

// Value buckets: 0 below lo, 1..nb in key order, nb + 1 above. Equal width: bucket h + 1 holds rel = k - lo with rel >> t = h.
// Equal count (big grids): rel >> t is a cell and map gives its bucket, cells in order, so the buckets hold about as many
// sampled entries each. low = the key bits a bucket keeps (rel minus its first key; the whole key in the catch-alls), bits(d)
// of them. Any partition in key order gives the same ranks.
struct Buckets {
  uint64_t lo = 0, mask = 0;
  int t = 0;
  size_t nb = 0, ncell = 0;
  const uint32_t *map = nullptr;
  const uint64_t *base = nullptr;
  const int *width = nullptr;

  template <bool TABLE> size_t locate(uint64_t k, uint64_t &low) const {
    low = k;
    if(k < lo) return 0;
    const uint64_t rel = k - lo, h = rel >> t;
    if(!TABLE){
      if(h >= nb) return nb + 1;
      low = rel & mask;
      return h + 1;
    }
    if(h >= ncell) return nb + 1;
    const size_t d = map[h];
    low = rel - base[d];
    return d;
  }
  int bits(size_t d) const {
    if(d == 0 || d > nb) return 64;
    return map == nullptr ? t : width[d];
  }
};

// Column jj of a block: Q[ii] = (a[ii] * bj) / f[ii], two at a time, or f itself for the network grid
const double *column_q(const double *f, const double *a, double bj, bool net, size_t s, double *q){
  if(net) return f;
  const v2d b2 = {bj, bj};
  size_t ii = 0;
  for(; ii + 2 <= s; ii += 2){
    v2d fv, av;
    std::memcpy(&fv, f + ii, sizeof fv);
    std::memcpy(&av, a + ii, sizeof av);
    const v2d qv = (av * b2) / fv;
    std::memcpy(q + ii, &qv, sizeof qv);
  }
  for(; ii < s; ++ii) q[ii] = (a[ii] * bj) / f[ii];
  return q;
}

// The same with f[ii] = col[node[ii]]
const double *column_q_gather(const double *col, const int *node, const double *a, double bj, bool net, size_t s, double *q){
  if(net){
    for(size_t ii = 0; ii < s; ++ii) q[ii] = col[node[ii]];
    return q;
  }
  const v2d b2 = {bj, bj};
  size_t ii = 0;
  for(; ii + 2 <= s; ii += 2){
    const v2d fv = {col[node[ii]], col[node[ii + 1]]};
    v2d av;
    std::memcpy(&av, a + ii, sizeof av);
    const v2d qv = (av * b2) / fv;
    std::memcpy(q + ii, &qv, sizeof qv);
  }
  for(; ii < s; ++ii) q[ii] = (a[ii] * bj) / col[node[ii]];
  return q;
}

// Q of rows [i0, i1) of column jj of block blk: from the block's copy in Fb, or from F one run of nodes at a time, or one
// node at a time
const double *block_column(const Cells &c, size_t blk, size_t jj, size_t i0, size_t i1, const double *a, double bj, bool net, double *q){
  const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
  if(c.block_fb[blk] != SIZE_MAX) return column_q(c.Fb + c.block_fb[blk] + jj * s + i0, a, bj, net, i1 - i0, q);
  const int *node = c.node.data() + first;
  const double *col = c.F + (size_t) node[jj] * c.L;
  if(c.block_gather[blk]) return column_q_gather(col, node + i0, a, bj, net, i1 - i0, q);
  const size_t *r = std::upper_bound(c.run.data() + c.block_run[blk], c.run.data() + c.block_run[blk + 1], i0);
  if(*r >= i1) return column_q(col + node[i0], a, bj, net, i1 - i0, q);
  for(size_t u = i0; u < i1; ++r){
    const size_t e = std::min(i1, *r);
    if(net) std::memcpy(q + (u - i0), col + node[u], (e - u) * sizeof(double));
    else column_q(col + node[u], a + (u - i0), bj, false, e - u, q + (u - i0));
    u = e;
  }
  return q;
}

class Barrier {
 public:
  explicit Barrier(size_t n) : n(n), count(0), phase(0) {}

  void wait(const std::atomic<bool> &failed){
    const size_t ph = phase.load();
    if(count.fetch_add(1) + 1 == n){
      count.store(0);
      phase.store(ph + 1);
      return;
    }
    for(unsigned spin = 0; phase.load() == ph; spin += spin <= 100){
      if(failed.load()) throw std::runtime_error("another thread failed");
      if(spin > 100) std::this_thread::yield();
    }
  }

 private:
  const size_t n;
  std::atomic<size_t> count, phase;
};

// Columns of a pipeline stage that are done, on its own cache line
struct alignas(64) Progress {
  std::atomic<size_t> done;
};

void wait_for(const std::atomic<size_t> &done, size_t need, const std::atomic<bool> &failed){
  for(unsigned spin = 0; done.load(std::memory_order_acquire) < need; spin += spin <= 100){
    if(failed.load()) throw std::runtime_error("another thread failed");
    if(spin > 100) std::this_thread::yield();
  }
}

// A team of threads scores one grid at a time. The blocks are cut into chunks of about cap / size entries: runs of whole
// blocks, or row ranges of a big block; a chunk's entries are in position order (column by column, the chunk's rows of
// each). On one thread, the team appends each entry to its bucket's list of slabs in one pass. Otherwise a first pass
// counts each chunk's entries per bucket and a second puts them in place, so a bucket is one run of the pool holding its
// chunks' entries in order, with no slack: smaller, and faster when several threads share the memory. Each chunk is a sum
// unit. The row ranges of a block are the stages of a pipeline: stage q continues the receiver sums of stage q - 1,
// column by column, so every sum keeps master's order.
//
// Every entry's bucket is stored (2 bytes). Grids that fit the cache, and every grid of a single thread, keep 8-byte keys. Big
// grids (above the cache budget) on several threads keep the key bits of a bucket in 6 bytes (at most 48 bits; the catch-alls, and
// any wider bucket, keep 8 in ext): 8 bytes per entry, and their buckets hold equal counts, so each sorts in L2 while the threads
// share L3.
struct Team {
  Team(const Cells &c, size_t size, bool alone) : size(size), bar(size), a(c.node.size()), b(c.node.size()),
                                                  rsum(c.node.size()), rcnt(c.node.size()), lo(size), hi(size) {
    const size_t nblocks = c.block_start.size() - 1, target = (c.cap + size - 1) / size;
    size_t open = 0, filled = 0;
    for(size_t blk = 0; blk < nblocks; ++blk){
      const size_t s = c.block_start[blk + 1] - c.block_start[blk];
      if(s * s > target){
        if(filled) add_chunk(c, open, blk, 0, SIZE_MAX, true, true);
        const size_t k = std::min(s, (s * s + target - 1) / target);
        for(size_t q = 0; q < k; ++q) add_chunk(c, blk, blk + 1, s * q / k, s * (q + 1) / k, q == 0, q + 1 == k);
        open = blk + 1;
        filled = 0;
      } else {
        filled += s * s;
        if(filled >= target){
          add_chunk(c, open, blk + 1, 0, SIZE_MAX, true, true);
          open = blk + 1;
          filled = 0;
        }
      }
    }
    if(filled) add_chunk(c, open, nblocks, 0, SIZE_MAX, true, true);
    b1 = std::max(B1_MIN, std::min(B1_MAX, bit_width(c.cap) - 8));
    nbmax = ((size_t) 1 << b1) + 2;
    stride = nbmax;
    const size_t nchunks = entries.size();
    slabbed = alone && nchunks == 1;
    big = !slabbed && c.cap * (sizeof(uint64_t) + sizeof(uint16_t)) > CACHE_BUDGET;
    if(slabbed){
      const size_t slabs = nbmax + entries[0] / SLAB + 1;
      pool.reset(new uint64_t[slabs * SLAB]);
      link.resize(slabs);
      cur.resize(nbmax);
      full.resize(nbmax);
      bucket.resize(c.cap);
    } else {
      cnt.resize(nchunks * nbmax);
      tot.resize(nbmax);
      if(big){
        pool6.reset(new char[c.cap * 6]);
        hist.resize(size);
      } else {
        pool.reset(new uint64_t[std::max<size_t>(1, c.cap)]);
      }
      bucket.resize(c.cap);
    }
    progress.reset(new Progress[nchunks]());
  }

  void add_chunk(const Cells &c, size_t b0, size_t b1, size_t i0, size_t i1, bool head, bool last){
    size_t e = 0;
    for(size_t blk = b0; blk < b1; ++blk){
      const size_t s = c.block_start[blk + 1] - c.block_start[blk];
      e += s * (std::min(i1, s) - i0);
    }
    chunk_b0.push_back(b0);
    chunk_b1.push_back(b1);
    chunk_i0.push_back(i0);
    chunk_i1.push_back(i1);
    chunk_head.push_back(head);
    chunk_last.push_back(last);
    entries.push_back(e);
  }

  const size_t size;
  Barrier bar;
  std::atomic<size_t> next_count, next_split, next_place, next_bucket, next_sum;
  size_t g, nbmax, stride;
  int l, r, b1;
  bool slabbed, big;                   // big: runs of 6-byte keys with equal-count buckets (a grid above the cache budget on threads)
  Buckets bk;                          // the grid's buckets, made by the leader
  std::unique_ptr<uint64_t[]> pool;    // slabs, or runs of a grid that fits the cache
  std::unique_ptr<char[]> pool6;       // big grids: 6 bytes per entry
  std::unique_ptr<Progress[]> progress;
  std::vector<size_t> link, cur, full;  // slabs: next slab of each slab, write position and full slabs of each bucket
  std::vector<uint32_t> cnt;           // runs: per chunk and bucket (stride apart), entries, then their offset in the bucket
  std::vector<size_t> tot, chunk_b0, chunk_b1, chunk_i0, chunk_i1, entries;
  std::vector<bool> chunk_head, chunk_last;
  std::vector<uint16_t> bucket;        // the bucket of each entry
  std::vector<uint64_t> ext;           // big grids: the keys of the buckets wider than 48 bits (the catch-alls at least)
  std::vector<size_t> sbase;           // big grids: where each bucket's keys start, in pool6 (6-byte slots) or in ext
  std::vector<char> swide;             // big grids: the bucket's keys are in ext
  std::vector<uint32_t> map;           // equal-count buckets: the bucket of each cell
  std::vector<uint64_t> base;          // and the first key (rel) and key bits of each bucket
  std::vector<int> width;
  std::vector<std::vector<uint32_t>> hist;  // big grids: each member's sampled entries per cell
  std::vector<double> a, b;
  std::vector<long double> rsum;
  std::vector<int64_t> rcnt;
  std::vector<uint64_t> lo, hi;
};

struct Work {
  std::vector<uint64_t> x, y;
  std::vector<uint32_t> z, cnt;
  std::vector<std::pair<size_t, size_t>> segs;
  std::vector<std::pair<uint64_t, uint64_t>> group;
  std::vector<size_t> start, at;
  std::vector<const char *> rp;
  std::vector<double> q;
  std::vector<long double> ssum;
  std::vector<int64_t> scnt, rk;

  Work(const Cells &c, const Team &tm) : start(tm.nbmax + 1), at(tm.nbmax), rp(tm.nbmax), q(c.max_block), ssum(c.max_block),
                                         scnt(c.max_block), rk(2 * c.max_block) {}
};

// A slab's entries hold key bits until ranked, then the first half of the slab holds their ranks as uint32
inline uint32_t get_rank(const uint64_t *pool, size_t p){
  uint32_t rank;
  std::memcpy(&rank, (const char *) (pool + p - p % SLAB) + p % SLAB * sizeof(uint32_t), sizeof rank);
  return rank;
}

// Keys of a bucket held in one run of the pool; their ranks go to the run's first half as uint32
struct RunKeys {
  uint64_t *keys;

  void pack(uint64_t *x, uint32_t *c, size_t m, int cut, int I, int shift) const {
    for(size_t i = 0; i < m; ++i){
      const uint64_t key = keys[i] >> cut;
      x[i] = (key << I) | i;
      ++c[(key >> shift) + 1];
    }
  }
  void copy(uint64_t *x, size_t m) const {
    std::memcpy(x, keys, m * sizeof(uint64_t));
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy((char *) keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(size_t){}
};

// The low 6 bytes of a key (little-endian), as the 6-byte pool slots hold it
inline uint64_t get48(const char *p){
  uint64_t v = 0;
  std::memcpy(&v, p, 6);
  return v;
}

// Keys of a bucket held in one run of 6-byte slots (at most NARROW_BITS key bits); their ranks go to the run's first 4m bytes
struct SixKeys {
  char *keys;

  void pack(uint64_t *x, uint32_t *c, size_t m, int cut, int I, int shift) const {
    for(size_t i = 0; i < m; ++i){
      const uint64_t key = get48(keys + 6 * i) >> cut;
      x[i] = (key << I) | i;
      ++c[(key >> shift) + 1];
    }
  }
  void copy(uint64_t *x, size_t m) const {
    for(size_t i = 0; i < m; ++i) x[i] = get48(keys + 6 * i);
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy(keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(size_t){}
};

// Keys of a bucket held in slab segments (offset, length); their ranks go to z, then to each slab's first half
struct SlabKeys {
  uint64_t *pool;
  const std::vector<std::pair<size_t, size_t>> &segs;
  uint32_t *z;

  void pack(uint64_t *x, uint32_t *c, size_t, int cut, int I, int shift) const {
    size_t i = 0;
    for(const auto &sg : segs){
      const uint64_t *src = pool + sg.first;
      for(size_t j = 0; j < sg.second; ++j, ++i){
        const uint64_t key = src[j] >> cut;
        x[i] = (key << I) | i;
        ++c[(key >> shift) + 1];
      }
    }
  }
  void copy(uint64_t *x, size_t) const {
    for(const auto &sg : segs){
      std::memcpy(x, pool + sg.first, sg.second * sizeof(uint64_t));
      x += sg.second;
    }
  }
  void put(size_t i, size_t rank){
    z[i] = (uint32_t) rank;
  }
  void finish(size_t){
    size_t i = 0;
    for(const auto &sg : segs){
      const size_t p = sg.first;
      std::memcpy((char *) (pool + p - p % SLAB) + p % SLAB * sizeof(uint32_t), z + i, sg.second * sizeof(uint32_t));
      i += sg.second;
    }
  }
};

// Ranks of one bucket of m entries, held in position order as their low t key bits: sort (key bits, index) words, then
// give each entry its min rank - 1
template <class K> void rank_bucket(K &keys, size_t m, size_t s0, int t, Work &w){
  const int I = bit_width(m - 1), cut = t + I > 64 ? t + I - 64 : 0, width = t - cut;
  const int want = std::min(12, bit_width(m)), shift = width > want ? width - want : 0;
  const uint64_t imask = I ? (~(uint64_t) 0 >> (64 - I)) : 0;
  const size_t nb = (size_t) 1 << (width - shift);
  if(w.x.size() < m){
    w.x.resize(m);
    w.y.resize(m);
  }
  uint64_t *x = w.x.data(), *y = w.y.data();
  w.cnt.assign(nb + 1, 0);
  uint32_t *c = w.cnt.data();
  keys.pack(x, c, m, cut, I, shift);
  if(m > 1) sort_words(x, y, m, I, shift, nb, c, w.cnt);
  else y[0] = x[0];
  if(cut == 0){
    uint64_t prev = y[0] >> I, first = s0;
    for(size_t v = 0; v < m; ++v){
      if((y[v] >> I) != prev){
        prev = y[v] >> I;
        first = s0 + v;
      }
      keys.put(y[v] & imask, first);
    }
  } else {
    uint64_t *key = x;
    keys.copy(key, m);
    for(size_t q = 0; q < m; ){
      size_t u = q + 1;
      while(u < m && (y[u] >> I) == (y[q] >> I)) ++u;
      w.group.clear();
      for(size_t v = q; v < u; ++v) w.group.emplace_back(key[y[v] & imask], y[v] & imask);
      std::sort(w.group.begin(), w.group.end());
      for(size_t g = 0, h = 0; g < w.group.size(); ++g){
        if(w.group[g].first != w.group[h].first) h = g;
        keys.put(w.group[g].second, s0 + q + h);
      }
      q = u;
    }
  }
  keys.finish(m);
}

// Member me of the team: a and b of its share of the cells
void prepare_ab(const Cells &c, Team &tm, size_t me){
  if(tm.l < 0) return;
  const size_t n = c.node.size(), k1 = n * (me + 1) / tm.size;
  if(c.F == nullptr){
    for(size_t k = n * me / tm.size; k < k1; ++k){
      tm.a[k] = c.A[(size_t) c.out[k] + (size_t) tm.l * n];
      tm.b[k] = c.B[(size_t) tm.r + (size_t) c.out[k] * c.nb];
    }
    return;
  }
  for(size_t k = n * me / tm.size; k < k1; ++k){
    tm.a[k] = c.F[(size_t) c.node[k] + (size_t) tm.l * c.L];
    tm.b[k] = c.F[(size_t) tm.r + (size_t) c.node[k] * c.L];
  }
}

// F[i, j] of entry (ii, jj) of block blk: from the block's copy, else from F
double block_value(const Cells &c, size_t blk, size_t ii, size_t jj){
  const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
  if(c.block_fb[blk] != SIZE_MAX) return c.Fb[c.block_fb[blk] + jj * s + ii];

  return c.F[(size_t) c.node[first + ii] + (size_t) c.node[first + jj] * c.L];
}

// fn(key) for member me's share of every 61st entry of each block that is not skipped
template <class Fn> void each_sample(const Cells &c, const Team &tm, size_t me, Fn fn){
  const size_t nblocks = c.block_start.size() - 1;
  for(size_t blk = 0; blk < nblocks; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
    const size_t ns = s * s > 3 ? (s * s - 4) / SAMPLE + 1 : 0, k1 = ns * (me + 1) / tm.size;
    const double *a = tm.a.data() + first, *b = tm.b.data() + first;
    for(size_t k = ns * me / tm.size; k < k1; ++k){
      const size_t q = 3 + k * SAMPLE, jj = q / s, ii = q % s;
      double v = block_value(c, blk, ii, jj);
      if(tm.l >= 0) v = (a[ii] * b[jj]) / v;
      if(skipped(ii, jj, v)) continue;
      fn(sortable_key(v));
    }
  }
}

// Member me: key range of its share of the sample
void sample_keys(const Cells &c, Team &tm, size_t me){
  uint64_t lo = ~(uint64_t) 0, hi = 0;
  each_sample(c, tm, me, [&](uint64_t k){
    lo = std::min(lo, k);
    hi = std::max(hi, k);
  });
  tm.lo[me] = lo;
  tm.hi[me] = hi;
}

// The sampled key range of the grid: lo, hi and its bit width R; false (lo = 0, R = 64) if nothing was sampled
bool key_range(const Team &tm, uint64_t &lo, uint64_t &hi, int &R){
  lo = ~(uint64_t) 0;
  hi = 0;
  for(size_t me = 0; me < tm.size; ++me){
    lo = std::min(lo, tm.lo[me]);
    hi = std::max(hi, tm.hi[me]);
  }
  if(lo > hi){
    lo = 0;
    R = 64;
    return false;
  }
  R = bit_width(hi - lo);
  return true;
}

// Member me: its share of the sample counted per cell of 2^tc keys from lo
void sample_cells(const Cells &c, Team &tm, size_t me, uint64_t lo, int tc, size_t ncell){
  std::vector<uint32_t> &h = tm.hist[me];
  h.assign(ncell, 0u);
  each_sample(c, tm, me, [&](uint64_t k){
    ++h[(k - lo) >> tc];
  });
}

// The leader: equal-count buckets from the members' cell counts. A bucket starts where the sampled entries below it reach the
// next of target equal shares, and spans at most 2^(NARROW_BITS - tc) cells so its keys fit 6 bytes. False if nothing was sampled.
bool count_buckets(Team &tm, size_t cap, uint64_t lo, int tc, size_t ncell){
  std::vector<uint32_t> &h = tm.hist[0];
  for(size_t me = 1; me < tm.size; ++me){
    const uint32_t *hm = tm.hist[me].data();
    for(size_t k = 0; k < ncell; ++k) h[k] += hm[k];
  }
  uint64_t total = 0;
  for(size_t k = 0; k < ncell; ++k) total += h[k];
  if(total == 0) return false;
  const uint64_t target = std::min<uint64_t>(BUCKETS_MAX, std::max<uint64_t>(BUCKETS_MIN, cap / BUCKET_TARGET));
  const size_t maxspan = (size_t) 1 << std::min(TABLE_BITS, NARROW_BITS - tc);
  tm.map.resize(ncell);
  tm.base.assign(1, 0);
  tm.width.assign(1, 64);
  size_t d = 0, c0 = 0;
  uint64_t cum = 0, prev = 0;
  for(size_t k = 0; k < ncell; ++k){
    const uint64_t share = cum * target / total;
    if(k == 0 || share != prev || k - c0 >= maxspan){
      if(k > 0) tm.width[d] = tc + bit_width(k - c0 - 1);
      ++d;
      c0 = k;
      tm.base.push_back((uint64_t) k << tc);
      tm.width.push_back(0);
    }
    tm.map[k] = (uint32_t) d;
    prev = share;
    cum += h[k];
  }
  tm.width[d] = tc + bit_width(ncell - c0 - 1);

  Buckets &bk = tm.bk;
  bk.lo = lo;
  bk.t = tc;
  bk.mask = 0;
  bk.nb = d;
  bk.ncell = ncell;
  bk.map = tm.map.data();
  bk.base = tm.base.data();
  bk.width = tm.width.data();
  return true;
}

// The leader: the grid's buckets for everyone in tm.bk, equal count for a big grid (at most BUCKETS_IDS), else equal width (2^b1
// buckets; for a big grid at least 2^(R - NARROW_BITS), at most 2^15, so the keys fit 6 bytes unless R is 64), and the per-bucket
// counts sized
void make_buckets(Team &tm, size_t cap, bool sampled, uint64_t lo, uint64_t hi, int R){
  const size_t nchunks = tm.entries.size();
  bool table = false;
  if(tm.big && sampled){
    const int tc = R > TABLE_BITS ? R - TABLE_BITS : 0;
    table = count_buckets(tm, cap, lo, tc, (size_t) ((hi - lo) >> tc) + 1) && tm.bk.nb <= BUCKETS_IDS;
  }
  if(!table){
    const int b1 = tm.big ? std::min(R, std::min(15, std::max(tm.b1, R - NARROW_BITS))) : std::min(R, tm.b1);
    Buckets &bk = tm.bk;
    bk.lo = lo;
    bk.t = R - b1;
    bk.mask = bk.t ? (~(uint64_t) 0 >> (64 - bk.t)) : 0;
    bk.nb = (size_t) 1 << b1;
    bk.ncell = 0;
    bk.map = nullptr;
    bk.base = nullptr;
    bk.width = nullptr;
  }
  if(!tm.slabbed){
    tm.stride = tm.bk.nb + 2;
    if(tm.cnt.size() < nchunks * tm.stride) tm.cnt.resize(nchunks * tm.stride);
    if(tm.tot.size() < tm.stride) tm.tot.resize(tm.stride);
  }
}

// One thread: Q of each entry in position order, appended to its bucket's list of slabs
void fill_slabs(const Cells &c, Team &tm, const Buckets &bk, Work &w){
  const size_t NB = bk.nb + 2;
  const bool net = tm.l < 0;
  uint64_t *pool = tm.pool.get();
  size_t *cur = tm.cur.data(), *full = tm.full.data(), *link = tm.link.data();
  size_t slabs = NB;
  for(size_t d = 0; d < NB; ++d){
    cur[d] = d * SLAB;
    full[d] = 0;
  }
  for(size_t blk = tm.chunk_b0[0]; blk < tm.chunk_b1[0]; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = tm.chunk_i0[0], m = std::min(tm.chunk_i1[0], s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const size_t pos = c.block_off[blk] + jj * s + i0;
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      uint16_t *bc = tm.bucket.data() + pos;
      for(size_t i = 0; i < m; ++i){
        const double v = q[i];
        if(skipped(i0 + i, jj, v)){
          bc[i] = SKIP;
          continue;
        }
        uint64_t low;
        const size_t d = bk.locate<false>(sortable_key(v), low);
        bc[i] = (uint16_t) d;
        size_t p = cur[d];
        pool[p++] = low;
        if(p % SLAB == 0){
          link[p / SLAB - 1] = slabs;
          p = slabs++ * SLAB;
          ++full[d];
        }
        cur[d] = p;
      }
    }
  }
}

// First pass over chunk ch: the bucket of each entry (SKIP if skipped) and the chunk's entries per bucket
template <bool TABLE> void count_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  std::fill(cnt, cnt + bk.nb + 2, 0u);
  for(size_t blk = tm.chunk_b0[ch]; blk < tm.chunk_b1[ch]; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = tm.chunk_i0[ch], m = std::min(tm.chunk_i1[ch], s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const double v = q[i];
        if(skipped(i0 + i, jj, v)){
          bc[i] = SKIP;
          continue;
        }
        uint64_t low;
        const size_t d = bk.locate<TABLE>(sortable_key(v), low);
        bc[i] = (uint16_t) d;
        ++cnt[d];
      }
    }
  }
  tm.progress[ch].done.store(0);
}

// Each bucket's chunk counts to the offsets of the chunks' entries in the bucket, a run of buckets at a time; tot = sizes
void split_buckets(Team &tm, size_t NB){
  const size_t nchunks = tm.entries.size(), step = 256;
  size_t run[step];
  for(size_t d0; (d0 = tm.next_split.fetch_add(step)) < NB; ){
    const size_t d1 = std::min(NB, d0 + step);
    std::fill(run, run + step, (size_t) 0);
    for(size_t ch = 0; ch < nchunks; ++ch){
      uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
      for(size_t d = d0; d < d1; ++d){
        const uint32_t e = cnt[d];
        cnt[d] = (uint32_t) run[d - d0];
        run[d - d0] += e;
      }
    }
    for(size_t d = d0; d < d1; ++d) tm.tot[d] = run[d - d0];
  }
}

// Second pass over chunk ch: the low key bits of each entry to its place in its bucket (a grid that fits the cache)
void place_chunk_ids(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  size_t *at = w.at.data();
  for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = w.start[d] + cnt[d];
  uint64_t *pool = tm.pool.get();
  for(size_t blk = tm.chunk_b0[ch]; blk < tm.chunk_b1[ch]; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = tm.chunk_i0[ch], m = std::min(tm.chunk_i1[ch], s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const size_t pos = c.block_off[blk] + jj * s + i0;
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      const uint16_t *bc = tm.bucket.data() + pos;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t k = sortable_key(q[i]);
        pool[at[d]++] = d == 0 || d > bk.nb ? k : (k - bk.lo) & bk.mask;
      }
    }
  }
}

// Second pass over chunk ch of a big grid: the kept key bits of each entry to its place in its bucket, 6 bytes in pool6 or 8 in
// ext (sbase / swide)
void place_chunk_six(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  const size_t NB = bk.nb + 2;
  size_t *at = w.at.data();
  for(size_t d = 0; d < NB; ++d) at[d] = tm.sbase[d] + cnt[d];
  const char *wide = tm.swide.data();
  char *pool = tm.pool6.get();
  uint64_t *ext = tm.ext.data();
  for(size_t blk = tm.chunk_b0[ch]; blk < tm.chunk_b1[ch]; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = tm.chunk_i0[ch], m = std::min(tm.chunk_i1[ch], s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t k = sortable_key(q[i]);
        uint64_t low = k;
        if(d > 0 && d <= bk.nb) low = bk.map == nullptr ? (k - bk.lo) & bk.mask : k - bk.lo - bk.base[d];
        if(wide[d]) ext[at[d]++] = low;
        else std::memcpy(pool + 6 * at[d]++, &low, 6);
      }
    }
  }
}

// The leader, big grid: where each bucket's keys go, buckets of at most NARROW_BITS key bits one after the other in pool6, the
// others (the catch-alls at least) in ext, sized now
void six_storage(Team &tm, const Buckets &bk){
  const size_t NB = bk.nb + 2;
  tm.sbase.resize(NB);
  tm.swide.resize(NB);
  size_t narrow = 0, wide = 0;
  for(size_t d = 0; d < NB; ++d){
    tm.swide[d] = bk.bits(d) > NARROW_BITS;
    tm.sbase[d] = tm.swide[d] ? wide : narrow;
    (tm.swide[d] ? wide : narrow) += tm.tot[d];
  }
  tm.ext.resize(wide);
}

// Rank the buckets, a run of about E / (16 x team size) entries at a time; start[d] = number of entries in the buckets below d
void rank_buckets(Team &tm, const Buckets &bk, Work &w){
  const size_t NB = bk.nb + 2;
  const size_t *start = w.start.data(), *full = tm.full.data(), *link = tm.link.data();
  const size_t unit = std::max<size_t>(1, start[NB] / (16 * tm.size)), nunits = (start[NB] + unit - 1) / unit;
  for(size_t u; (u = tm.next_bucket.fetch_add(1)) < nunits; ){
    const size_t d1 = std::lower_bound(start, start + NB, (u + 1) * unit) - start;
    for(size_t d = std::lower_bound(start, start + NB, u * unit) - start; d < d1; ++d){
      const size_t m = start[d + 1] - start[d];
      if(m == 0) continue;
      const int t = bk.bits(d);
      if(tm.slabbed){
        w.segs.clear();
        size_t slab = d;
        for(size_t f = 0; f < full[d]; ++f){
          w.segs.emplace_back(slab * SLAB, SLAB);
          slab = link[slab];
        }
        if(tm.cur[d] % SLAB) w.segs.emplace_back(slab * SLAB, tm.cur[d] % SLAB);
        if(w.z.size() < m) w.z.resize(m);
        SlabKeys keys = {tm.pool.get(), w.segs, w.z.data()};
        rank_bucket(keys, m, start[d], t, w);
      } else if(!tm.big){
        RunKeys keys = {tm.pool.get() + start[d]};
        rank_bucket(keys, m, start[d], t, w);
      } else if(tm.swide[d]){
        RunKeys keys = {tm.ext.data() + tm.sbase[d]};
        rank_bucket(keys, m, start[d], t, w);
      } else {
        SixKeys keys = {tm.pool6.get() + 6 * tm.sbase[d]};
        rank_bucket(keys, m, start[d], t, w);
      }
    }
  }
}

// Ranks of m entries in position order from their buckets' slab lists; at[d] = next slab entry of bucket d
inline void pull_slabs(const uint16_t *bc, size_t m, int64_t *rk, const uint64_t *pool, const size_t *link, size_t *at){
  for(size_t i = 0; i < m; ++i){
    const size_t d = bc[i];
    if(d == SKIP){
      rk[i] = -1;
      continue;
    }
    size_t p = at[d];
    rk[i] = get_rank(pool, p);
    if(++p % SLAB == 0) p = link[p / SLAB - 1] * SLAB;
    at[d] = p;
  }
}

// The same from their buckets' runs; at[d] = index of the next uint32 rank of bucket d (a grid that fits the cache)
inline void pull_runs_ids(const uint16_t *bc, size_t m, int64_t *rk, const char *ranks, size_t *at){
  for(size_t i = 0; i < m; ++i){
    const size_t d = bc[i];
    if(d == SKIP){
      rk[i] = -1;
      continue;
    }
    uint32_t rank;
    std::memcpy(&rank, ranks + at[d]++ * sizeof(uint32_t), sizeof rank);
    rk[i] = rank;
  }
}

// Ranks of m entries of a big grid in position order from their buckets' runs; rp[d] = the next uint32 rank of bucket d
inline void pull_runs_six(const uint16_t *bc, size_t m, int64_t *rk, const char **rp){
  for(size_t i = 0; i < m; ++i){
    const size_t d = bc[i];
    if(d == SKIP){
      rk[i] = -1;
      continue;
    }
    uint32_t rank;
    std::memcpy(&rank, rp[d], sizeof rank);
    rp[d] += sizeof rank;
    rk[i] = rank;
  }
}

// Sender and receiver sums of chunk ch, two columns at a time, each column's ranks taken from the chunk's entries of each
// bucket in position order. A stage after the first continues each column's receiver sum where the stage before it left it.
template <int MODE> void sum_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, double denom, Work &w, double *sender,
                                       double *receiver, const std::atomic<bool> &failed){
  const uint64_t *pool = tm.pool.get();
  const size_t *link = tm.link.data();
  size_t *at = w.at.data();
  if(MODE == SLABS){
    for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = d * SLAB;
  } else if(MODE == IDS){
    const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
    for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = 2 * w.start[d] + cnt[d];
  } else {
    const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
    const char *ext = (const char *) tm.ext.data(), *pool6 = tm.pool6.get();
    const char **rp = w.rp.data();
    for(size_t d = 0; d < bk.nb + 2; ++d){
      rp[d] = (tm.swide[d] ? ext + sizeof(uint64_t) * tm.sbase[d] : pool6 + 6 * tm.sbase[d]) + sizeof(uint32_t) * cnt[d];
    }
  }
  const bool head = tm.chunk_head[ch], last = tm.chunk_last[ch];
  std::atomic<size_t> &done = tm.progress[ch].done;
  const v2d den = {denom, denom};
  long double *ssum = w.ssum.data();
  int64_t *scnt = w.scnt.data(), *rk0 = w.rk.data();
  for(size_t blk = tm.chunk_b0[ch]; blk < tm.chunk_b1[ch]; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, off = c.block_off[blk];
    const size_t i0 = tm.chunk_i0[ch], m = std::min(tm.chunk_i1[ch], s) - i0;
    const int *out = c.out.data() + first;
    long double *rsum = tm.rsum.data() + first;
    int64_t *rcnt = tm.rcnt.data() + first, *rk1 = rk0 + m;
    for(size_t k = 0; k < m; ++k){
      ssum[k] = 0.0L;
      scnt[k] = 0;
    }
    for(size_t jj = 0; jj < s; jj += 2){
      const size_t ncol = std::min((size_t) 2, s - jj);
      for(size_t h = 0; h < ncol; ++h){
        if(MODE == SLABS) pull_slabs(tm.bucket.data() + off + (jj + h) * s + i0, m, h ? rk1 : rk0, pool, link, at);
        else if(MODE == IDS) pull_runs_ids(tm.bucket.data() + off + (jj + h) * s + i0, m, h ? rk1 : rk0, (const char *) pool, at);
        else pull_runs_six(tm.bucket.data() + off + (jj + h) * s + i0, m, h ? rk1 : rk0, w.rp.data());
      }
      if(ncol == 1) std::fill(rk1, rk1 + m, (int64_t) -1);
      long double r0 = 0.0L, r1 = 0.0L;
      int64_t c0 = 0, c1 = 0;
      if(!head){
        wait_for(tm.progress[ch - 1].done, jj + ncol, failed);
        r0 = rsum[jj];
        c0 = rcnt[jj];
        if(ncol == 2){
          r1 = rsum[jj + 1];
          c1 = rcnt[jj + 1];
        }
      }
      for(size_t i = 0; i < m; ++i){
        const int64_t k0 = rk0[i], k1 = rk1[i];
        const v2d num = {(double) k0, (double) k1}, pr = num / den;
        const double p0 = k0 >= 0 ? pr[0] : 0.0, p1 = k1 >= 0 ? pr[1] : 0.0;
        ssum[i] = (ssum[i] + p0) + p1;
        scnt[i] += (k0 >= 0) + (k1 >= 0);
        r0 += p0;
        r1 += p1;
        c0 += k0 >= 0;
        c1 += k1 >= 0;
      }
      if(last){
        receiver[out[jj]] = (double) (r0 / (long double) c0);
        if(ncol == 2) receiver[out[jj + 1]] = (double) (r1 / (long double) c1);
      } else {
        rsum[jj] = r0;
        rcnt[jj] = c0;
        if(ncol == 2){
          rsum[jj + 1] = r1;
          rcnt[jj + 1] = c1;
        }
        if((jj + ncol) % PIPE == 0 || jj + ncol == s) done.store(jj + ncol, std::memory_order_release);
      }
    }
    for(size_t k = 0; k < m; ++k) sender[out[i0 + k]] = (double) (ssum[k] / (long double) scnt[k]);
  }
}

// One team member: score grids until none are left; me = 0 leads
void run_team(const Cells &c, Team &tm, size_t me, const int *lig, const int *rec, size_t G, std::atomic<size_t> &next,
              const std::atomic<bool> &failed, double *S, double *R){
  const size_t n = c.node.size(), nchunks = tm.entries.size();
  Work w(c, tm);
  for(;;){
    tm.bar.wait(failed);
    if(me == 0){
      tm.g = failed.load() ? G : next.fetch_add(1);
      if(tm.g < G){
        tm.l = lig[tm.g];
        tm.r = rec[tm.g];
      }
      tm.next_count = 0;
      tm.next_split = 0;
      tm.next_place = 0;
      tm.next_bucket = 0;
      tm.next_sum = 0;
    }
    tm.bar.wait(failed);
    const size_t g = tm.g;
    if(g >= G) break;
    prepare_ab(c, tm, me);
    tm.bar.wait(failed);
    sample_keys(c, tm, me);
    tm.bar.wait(failed);
    uint64_t lo, hi;
    int range;
    const bool sampled = key_range(tm, lo, hi, range);
    if(tm.big && sampled){
      const int tc = range > TABLE_BITS ? range - TABLE_BITS : 0;
      sample_cells(c, tm, me, lo, tc, (size_t) ((hi - lo) >> tc) + 1);
      tm.bar.wait(failed);
    }
    if(me == 0) make_buckets(tm, c.cap, sampled, lo, hi, range);
    tm.bar.wait(failed);
    const Buckets &bk = tm.bk;
    const size_t NB = bk.nb + 2;
    if(w.start.size() < NB + 1){
      w.start.resize(NB + 1);
      w.at.resize(NB);
      w.rp.resize(NB);
    }
    w.start[0] = 0;
    if(tm.slabbed){
      if(tm.next_count.fetch_add(1) == 0) fill_slabs(c, tm, bk, w);
      tm.bar.wait(failed);
      for(size_t d = 0; d < NB; ++d) w.start[d + 1] = w.start[d] + tm.full[d] * SLAB + tm.cur[d] % SLAB;
    } else {
      for(size_t ch; (ch = tm.next_count.fetch_add(1)) < nchunks; ){
        if(bk.map) count_chunk<true>(c, tm, bk, ch, w);
        else count_chunk<false>(c, tm, bk, ch, w);
      }
      tm.bar.wait(failed);
      split_buckets(tm, NB);
      tm.bar.wait(failed);
      for(size_t d = 0; d < NB; ++d) w.start[d + 1] = w.start[d] + tm.tot[d];
      if(tm.big){
        if(me == 0) six_storage(tm, bk);
        tm.bar.wait(failed);
      }
      for(size_t ch; (ch = tm.next_place.fetch_add(1)) < nchunks; ){
        if(tm.big) place_chunk_six(c, tm, bk, ch, w);
        else place_chunk_ids(c, tm, bk, ch, w);
      }
      tm.bar.wait(failed);
    }
    rank_buckets(tm, bk, w);
    tm.bar.wait(failed);
    const double denom = (double) w.start[NB] - 1.0;
    for(size_t ch; (ch = tm.next_sum.fetch_add(1)) < nchunks; ){
      if(tm.slabbed) sum_chunk<SLABS>(c, tm, bk, ch, denom, w, S + g * n, R + g * n, failed);
      else if(tm.big) sum_chunk<SIX>(c, tm, bk, ch, denom, w, S + g * n, R + g * n, failed);
      else sum_chunk<IDS>(c, tm, bk, ch, denom, w, S + g * n, R + g * n, failed);
    }
    if(w.x.size() > KEEP_SORT_BUFFER){
      std::vector<uint64_t>().swap(w.x);
      std::vector<uint64_t>().swap(w.y);
      std::vector<uint32_t>().swap(w.z);
    }
    if(w.group.capacity() > KEEP_SORT_BUFFER) std::vector<std::pair<uint64_t, uint64_t>>().swap(w.group);
  }
}


// Cells grouped into within-sample blocks: samples ascending, the cells of a sample in the given order (node = their index in F)
void group_cells(Cells &c, const int *cell, const int *sample, size_t n){
  int n_samples = 0;
  for(size_t k = 0; k < n; ++k) n_samples = std::max(n_samples, sample[k]);
  std::vector<std::vector<int>> members(n_samples);
  for(size_t k = 0; k < n; ++k) members[sample[k] - 1].push_back((int) k);
  c.cap = 0;
  c.max_block = 0;
  c.block_start.push_back(0);
  for(int s = 0; s < n_samples; ++s){
    if(members[s].empty()) continue;
    for(int k : members[s]){
      c.node.push_back(cell[k]);
      c.out.push_back(k);
    }
    const size_t size = members[s].size();
    c.block_off.push_back(c.cap);
    c.block_start.push_back(c.node.size());
    c.cap += size * size;
    c.max_block = std::max(c.max_block, size);
  }
  if(c.cap >= (size_t) 4294967295u) Rcpp::stop("grid scoring: too many cells per sample (the sum of squared sample sizes must stay below 2^32)");
}

// Sender and receiver scores (cell x pair) of the grids of ligands lp and receptor subunits rp (negative: the network score)
Rcpp::List run_grids(const Cells &c, const int *lp, const int *rp, size_t G, int threads){
  const size_t n = c.node.size();
  Rcpp::NumericMatrix sender(n, G), receiver(n, G);
  double *S = REAL(sender), *R = REAL(receiver);

  const FpuState fpu = get_fpu();
  std::atomic<size_t> next(0);
  std::atomic<bool> failed(false);
  std::string error;

  const size_t cores = std::thread::hardware_concurrency();
  size_t nt = G == 0 || threads < 1 ? 1 : threads;
  if(cores > 0) nt = std::min(nt, cores);
  const size_t grid_bytes = std::max<size_t>(1, c.cap * (sizeof(uint64_t) + sizeof(uint16_t)));
  const size_t fit = std::max<size_t>(1, CACHE_BUDGET / grid_bytes);
  const size_t nteams = std::max<size_t>(1, std::min(std::min(nt, fit), G));
  std::vector<std::unique_ptr<Team>> teams;
  try {
    for(size_t k = 0; k < nteams; ++k) teams.emplace_back(new Team(c, nt / nteams + (k < nt % nteams ? 1 : 0), nt == 1));
  } catch(std::bad_alloc &){
    Rcpp::stop("grid scoring: not enough memory for the grid buffers (" + std::to_string(nteams * grid_bytes >> 20) + " MB)");
  }
  error.reserve(256);

  auto worker = [&](Team *tm, size_t me){
    try {
      set_fpu(fpu);
      run_team(c, *tm, me, lp, rp, G, next, failed, S, R);
    } catch(std::exception &e){
      if(!failed.exchange(true)){
        try {
          error = e.what();
        } catch(...){}
      }
    } catch(...){
      if(!failed.exchange(true)) error = "unknown error";
    }
  };

  if(nt == 1){
    worker(teams[0].get(), 0);
  } else {
    std::vector<std::thread> pool;
    bool started = true;
    try {
      pool.reserve(nt);
      for(auto &tm : teams){
        for(size_t me = 0; me < tm->size; ++me) pool.emplace_back(worker, tm.get(), me);
      }
    } catch(...){
      started = false;
      failed.store(true);
    }
    for(auto &th : pool) th.join();
    if(!started) error = "could not start threads";
  }
  if(failed.load()) Rcpp::stop("grid scoring: " + error);

  return Rcpp::List::create(Rcpp::Named("sender") = sender, Rcpp::Named("receiver") = receiver);
}


}  // namespace

// [[Rcpp::export]]
Rcpp::List score_grids_cpp(Rcpp::NumericMatrix F, Rcpp::IntegerVector cell, Rcpp::IntegerVector sample, Rcpp::IntegerVector lig,
                           Rcpp::IntegerVector rec, int threads){
  const size_t L = F.nrow(), n = cell.size(), G = lig.size();
  if((size_t) F.ncol() != L) Rcpp::stop("score_grids_cpp: the influence matrix must be square");
  if((size_t) sample.size() != n) Rcpp::stop("score_grids_cpp: one sample per cell");
  if((size_t) rec.size() != G) Rcpp::stop("score_grids_cpp: one receptor per ligand");
  for(size_t k = 0; k < n; ++k){
    if(cell[k] < 0 || (size_t) cell[k] >= L) Rcpp::stop("score_grids_cpp: cell index out of range");
    if(sample[k] < 1 || sample[k] == NA_INTEGER) Rcpp::stop("score_grids_cpp: sample ids must be positive integers");
  }
  for(size_t g = 0; g < G; ++g){
    if(lig[g] >= (int) L || rec[g] >= (int) L || (lig[g] < 0) != (rec[g] < 0)) Rcpp::stop("score_grids_cpp: bad ligand or receptor index");
  }

  Cells c;
  c.F = REAL(F);
  c.L = L;
  group_cells(c, INTEGER(cell), INTEGER(sample), n);
  size_t copied = 0;
  c.block_run.push_back(0);
  for(size_t blk = 0; blk + 1 < c.block_start.size(); ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
    const int *node = c.node.data() + first;
    size_t runs = 1;
    for(size_t k = 1; k < s; ++k) runs += node[k] != node[k - 1] + 1;
    if(runs * RUN <= s){
      c.block_fb.push_back(SIZE_MAX);
      for(size_t k = 1; k < s; ++k){
        if(node[k] != node[k - 1] + 1) c.run.push_back(k);
      }
      c.run.push_back(s);
    } else {
      c.block_fb.push_back(copied);
      copied += s * s;
    }
    c.block_run.push_back(c.run.size());
  }
  const bool gather = copied * sizeof(double) > CACHE_BUDGET;
  for(size_t &off : c.block_fb){
    c.block_gather.push_back(gather && off != SIZE_MAX);
    if(gather) off = SIZE_MAX;
  }
  if(gather) copied = 0;
  c.Fown.resize(copied);
  for(size_t blk = 0; blk + 1 < c.block_start.size(); ++blk){
    if(c.block_fb[blk] == SIZE_MAX) continue;
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, off = c.block_fb[blk];
    const int *node = c.node.data() + first;
    for(size_t jj = 0; jj < s; ++jj){
      const double *col = c.F + (size_t) node[jj] * L;
      for(size_t ii = 0; ii < s; ++ii) c.Fown[off + jj * s + ii] = col[node[ii]];
    }
  }
  c.Fb = c.Fown.data();

  return run_grids(c, INTEGER(lig), INTEGER(rec), G, threads);
}

// [[Rcpp::export]]
Rcpp::List score_blocks_cpp(Rcpp::NumericVector blocks, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample,
                            Rcpp::IntegerVector lig, Rcpp::IntegerVector rec, int threads){
  const size_t n = sample.size(), G = lig.size();
  if((size_t) A.nrow() != n || (size_t) B.ncol() != n) Rcpp::stop("score_blocks_cpp: A needs one row and B one column per cell");
  if((size_t) rec.size() != G) Rcpp::stop("score_blocks_cpp: one receptor per ligand");
  for(size_t k = 0; k < n; ++k){
    if(sample[k] < 1 || sample[k] == NA_INTEGER) Rcpp::stop("score_blocks_cpp: sample ids must be positive integers");
  }
  for(size_t g = 0; g < G; ++g){
    if(lig[g] >= A.ncol() || rec[g] >= B.nrow() || (lig[g] < 0) != (rec[g] < 0)) Rcpp::stop("score_blocks_cpp: bad ligand or receptor index");
  }

  Cells c;
  c.A = REAL(A);
  c.B = REAL(B);
  c.nb = B.nrow();
  std::vector<int> cell(n);
  for(size_t k = 0; k < n; ++k) cell[k] = (int) k;
  group_cells(c, cell.data(), INTEGER(sample), n);
  if((size_t) blocks.size() != c.cap) Rcpp::stop("score_blocks_cpp: the blocks must hold every within-sample entry");
  c.Fb = REAL(blocks);
  c.block_run.push_back(0);
  for(size_t blk = 0; blk + 1 < c.block_start.size(); ++blk){
    c.block_fb.push_back(c.block_off[blk]);
    c.block_gather.push_back(false);
    c.block_run.push_back(0);
  }

  return run_grids(c, INTEGER(lig), INTEGER(rec), G, threads);
}
