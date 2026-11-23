// score_grids.cpp: sender and receiver scores of cell x cell grids, within one ulp of those of HyperCom 0.0.0.9000
//
// For a ligand l and one receptor subunit r, HyperCom 0.0.0.9000 builds the cell x cell grid
//   Q[i, j] = (F[i, l] * F[r, j]) / F[i, j]      (F = influence matrix; Q[i, j] = F[i, j] for the network score)
// sets cross-sample entries, the diagonal, zeros and NaN to NA, percent-ranks the remaining K entries together
// ((min rank - 1) / (K - 1)) and averages the ranks over j (sender of i, colMeans) and over i (receiver of j, rowMeans).
//
// Here only the within-sample blocks are visited (every cross-sample entry is NA). Entries are ranked by an
// order-preserving integer image of their value: each entry falls in one of the grid's value buckets in key order, or in one
// of two catch-alls below and above the sampled key range. The bounds come from a sample of the entries: 2^b1 buckets of
// equal width (b1 = 11-15, growing with the grid; fewer when the sampled key range is narrow), or, for a grid above the
// cache budget scored by several threads, buckets of equal counts from a histogram of the sample, so every bucket sorts in
// L2. Each entry's key bits go to its bucket, then each bucket is sorted in cache as packed (key bits, index) words and its
// entries get their min ranks in place. Any partition in key order gives the same ranks. Each score is its exact integer rank
// sum divided once, sum / (count x (K - 1)), correctly rounded: the same on every computer, and within one ulp of the long double
// means of HyperCom 0.0.0.9000 (J. Moy, 2026-11-20). A team of threads scores one grid at a time; several teams run at once when
// their grids fit the cache together, else one team takes all threads and splits the row ranges of a big block among them.
// Memory per within-sample entry per team: its bucket id (2 B) and its key bits, 8 B, or 6 B
// for a grid above the cache budget on several threads (a bucket then keeps at most 48 key bits; the few entries of the
// catch-alls keep 8 B), plus 8 B per entry of samples whose cells are not consecutive in F while those copies fit the cache,
// plus each thread's sort buffers (16 B per entry of the biggest bucket it ranks). Do not compile with -ffast-math.
//
// score_grids_cpp reads F (a dense n x n influence matrix); score_blocks_cpp is given instead every within-sample block of F
// (made from the factored influence matrix, see influence.cpp) and F[cell, ligand] / F[receptor, cell] for the genes of the pairs,
// so no n x n matrix exists. Both rank and average the same values the same way.
//
// Memory mode (score_lambda_cpp): no F block either. Entries are ranked by integer keys ln(a b / F) = A[i] + B'[j] + Λ[i, j] (f2.h),
// Λ = -ln F of the within-sample blocks as 6-byte offsets in tiles of the upper triangle (influence_lambda_cpp, about 3 B per entry),
// A, B' the logs of F[cell, ligand] and F[receptor, cell] / D, all made from integers alone; the network grid is ranked by
// ln D[j] - Λ[i, j], the order of F. Only the diagonal is skipped. Scores are within ~1e-11 of the exact kernel's.
//
// Passes (score_blocks_cpp, score_lambda_cpp): a big grid on several threads can place, rank and sum its buckets in P ranges, one
// after the other, recomputing each entry's bucket at each pass instead of storing it: no bucket ids, and a key pool of ~6.75 / P B per
// entry instead of 6 B. Same integer sums, so the same scores at any P.

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
#include "f2.h"
#include "grid_layout.h"
#include "r_alloc.h"
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#include <emmintrin.h>
#define HYPERCOM_X86 1
#endif

namespace {

// The SSE control word (rounding, denormals) of the calling thread, given to the workers: Q's divisions round as R's would. No x87
// state: since 0.2.1 the sums are integers, not long double
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

typedef double v2d __attribute__((vector_size(16)));  // two doubles, one SSE2 instruction

// Hard limits
const uint16_t SKIP = 0xFFFF;           // bucket id of a skipped entry; the ids of a grid's buckets (0..nb + 1) stay below it
const int B1_LIMIT = 15;                // equal width: at most 2^15 buckets
const size_t BUCKETS_IDS = 65533;       // equal count: at most this many buckets
const int NARROW_BITS = 48;             // key bits of a 6-byte pool slot: a big grid's bucket keeps at most this many
static_assert(((size_t) 1 << B1_LIMIT) + 1 < SKIP, "equal-width bucket ids must stay below SKIP");
static_assert(BUCKETS_IDS + 1 < SKIP, "equal-count bucket ids must stay below SKIP");
static_assert(NARROW_BITS == 8 * 6, "a 6-byte pool slot holds 48 key bits");
// Ranks are uint32. A grid of LOCAL_FROM entries or more (run stores only) keeps each rank counted from its bucket's first and adds
// the bucket's start (64-bit) when reading it, so a bucket holds at most BUCKET_MAX entries. A chunk holds fewer than CHUNK_MAX
// entries, so its per-bucket counts (uint32) cannot wrap, and a slab-stored grid (one chunk) keeps global ranks. Test builds lower
// them (-DHYPERCOM_LOCAL_FROM=0 -DHYPERCOM_BUCKET_MAX=... -DHYPERCOM_CHUNK_MAX=...).
#ifndef HYPERCOM_LOCAL_FROM
#define HYPERCOM_LOCAL_FROM 0xFFFFFFFFu
#endif
#ifndef HYPERCOM_BUCKET_MAX
#define HYPERCOM_BUCKET_MAX 0xFFFFFFFFu
#endif
#ifndef HYPERCOM_CHUNK_MAX
#define HYPERCOM_CHUNK_MAX 0xFFFFFFFFu
#endif
const size_t LOCAL_FROM = HYPERCOM_LOCAL_FROM;
const size_t BUCKET_MAX = HYPERCOM_BUCKET_MAX;
const size_t CHUNK_MAX = HYPERCOM_CHUNK_MAX;
const size_t GRID_MAX = (size_t) 1 << 52;  // entries of a grid: R's longest vector; ranks and their count exact as doubles
const char *const BUCKET_ERROR = "more than 2^32 - 1 entries fall in one value bucket (a sample far bigger than memory allows, or "
                                 "very many equal values)";
static_assert(LOCAL_FROM <= UINT32_MAX && BUCKET_MAX <= UINT32_MAX && CHUNK_MAX <= UINT32_MAX && CHUNK_MAX >= 4,
              "ranks and bucket counts are uint32");

// Tuning knobs: speed and memory only, never a result
const size_t CACHE_BUDGET = (size_t) 24 << 20;  // team buffers that may be in flight at once (grid_bytes); a bigger grid gets all
                                                 // threads. R's batch_size() mirrors the 10 B per entry and the 24 MB
const int B1_MIN = 11, B1_MAX = 14;     // equal width: b1 grows with the grid between these
const int FILL_BITS = 8;                // equal width: about 2^8 entries per bucket (b1 = bit width of the entries - 8)
const size_t SAMPLE = 61;               // bucket bounds from every 61st entry of a block,
const int SAMPLE_BITS = 22;             // or every (K >> 22)-th (odd) when that is sparser, about 2^22 entries: big grids sample faster
const size_t SAMPLE_FIRST = 3;          // starting at its entry 3
const int TABLE_BITS = 18;              // big grids: the sampled key range is cut into 2^18 cells whose counts set the buckets
const size_t BUCKET_TARGET = 65536;     // big grids: entries per bucket aimed at, so a bucket's sort buffers stay in L2
const size_t BUCKETS_MIN = 2048, BUCKETS_MAX = 8192;  // binds above ~5.4e8 entries: fewer, bigger buckets place faster (one sample on
                                                      // 16 threads: 30k cells same, 40k -14%, 65k -15% vs 32768; 4096 -11% slower)
const size_t SLAB = 64;                 // entries per pool slab on one thread: each bucket is a linked list of slabs
const size_t RUN = 32;                  // a block whose nodes come in runs this long on average is read from F, others are copied
                                        // while all copies fit CACHE_BUDGET, else gathered from F
const int TOP_BITS = 16;                // sorting a bucket: at most 2^16 sub-buckets in its first radix pass (about one entry each, so
                                        // the insertion sort barely moves anything; wider spills the scatter out of cache)
const int RADIX_BITS = 12;              // and at most 2^12 in a deeper pass
const int PIECE_BITS = 10;              // a bucket of more than 2^12 entries is cut by its top key bits into pieces of ~2^10, each sorted
                                        // and ranked at once in L1 (bucket buffers bigger than a thread's share of L3 are not swept
                                        // again)
const size_t INSERTION_MAX = 32;        // sub-buckets of at most this many words are left to the final insertion sort
const size_t SPLIT_STEP = 256;          // buckets a thread takes at once when it turns the chunk counts into offsets
const size_t RANK_UNITS = 16;           // ranking: work units per team member
const size_t KEEP_SORT_BUFFER = (size_t) 1 << 20;  // per-thread sort buffers above this many entries are freed after each grid
const unsigned SPIN = 100;              // busy-wait rounds before a waiting thread yields

enum class Store { SLABS, RUNS };       // where the ranks are: slab lists (one thread), or one run of the pool per bucket

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

// Insertion sort with x[0] made the minimum first (searched in x[0, e0), where it is), so the inner loop needs no bounds test
void insertion_sort_sentinel(uint64_t *x, size_t m, size_t e0){
  size_t k = 0;
  for(size_t i = 1; i < e0; ++i){
    if(x[i] < x[k]) k = i;
  }
  std::swap(x[0], x[k]);
  for(size_t i = 1; i < m; ++i){
    const uint64_t v = x[i];
    if(x[i - 1] <= v) continue;
    size_t j = i;
    do {
      x[j] = x[j - 1];
      --j;
    } while(x[j - 1] > v);
    x[j] = v;
  }
}

// After a radix pass into y (cnt[coff + b] = end of sub-bucket b of nb, the biggest holds most): deeper passes for crowded
// sub-buckets, then one insertion sort over the nearly sorted result
void sort_rest(uint64_t *x, uint64_t *y, size_t m, int P, size_t nb, uint32_t most, size_t coff, std::vector<uint32_t> &cnt);

// Sort packed words (key bits << P | index) of x into y, x is scratch: one radix pass over the key range of the words, then sort_rest
void sort_words(uint64_t *x, uint64_t *y, size_t m, int P, std::vector<uint32_t> &cnt){
  uint64_t lo = x[0] >> P, hi = lo;
  for(size_t p = 1; p < m; ++p){
    const uint64_t k = x[p] >> P;
    lo = std::min(lo, k);
    hi = std::max(hi, k);
  }
  if(lo == hi){
    std::memcpy(y, x, m * sizeof(uint64_t));
    return;
  }
  const int bits = bit_width(hi - lo), want = std::min(RADIX_BITS, bit_width(m)), shift = bits > want ? bits - want : 0;
  const size_t nb = (size_t) ((hi - lo) >> shift) + 1, base = cnt.size();
  cnt.resize(base + nb + 1);
  uint32_t *c = cnt.data() + base;
  std::fill(c, c + nb + 1, 0u);
  for(size_t p = 0; p < m; ++p) ++c[(((x[p] >> P) - lo) >> shift) + 1];
  uint32_t most = 0;
  for(size_t b = 0; b < nb; ++b){
    most = std::max(most, c[b + 1]);
    c[b + 1] += c[b];
  }
  for(size_t p = 0; p < m; ++p) y[c[((x[p] >> P) - lo) >> shift]++] = x[p];
  if(shift > 0) sort_rest(x, y, m, P, nb, most, base, cnt);
  cnt.resize(base);
}

void sort_rest(uint64_t *x, uint64_t *y, size_t m, int P, size_t nb, uint32_t most, size_t coff, std::vector<uint32_t> &cnt){
  if(most > INSERTION_MAX){
    for(size_t b = 0, s = 0; b < nb; ++b){
      const size_t e = cnt[coff + b];
      if(e - s > INSERTION_MAX){
        std::memcpy(x + s, y + s, (e - s) * sizeof(uint64_t));
        sort_words(x + s, y + s, e - s, P, cnt);
      }
      s = e;
    }
  }
  size_t b = 0;
  while(cnt[coff + b] == 0) ++b;
  insertion_sort_sentinel(y, m, cnt[coff + b]);
}

struct Cells {
  const double *Fb = nullptr;          // F[i, j] of the entries of the copied (or given) blocks, in position order
  std::vector<double> Fown;            // the copies, when made here
  const double *F = nullptr;           // the whole influence matrix (L x L), or null when every block is given
  const double *A = nullptr;           // without F: F[cell, ligand], cells x ligands, rows in cell order
  const double *B = nullptr;           // and F[receptor, cell], receptors x cells (nrec rows)
  size_t nrec = 0;
  size_t L = 0;
  std::vector<int> node;
  std::vector<int> out;
  std::vector<size_t> block_start;
  std::vector<size_t> block_off;
  std::vector<size_t> block_fb;        // offset of a copied block in Fb, SIZE_MAX for a block read from F
  std::vector<bool> block_gather;      // a block read from F one node at a time
  std::vector<size_t> run;             // per block read from F in runs: the end (exclusive) of each run of consecutive nodes
  std::vector<size_t> block_run;
  size_t cap = 0;
  size_t max_block = 0;
  const char *LamT = nullptr;          // memory mode: Λ - lbase in tiles (f2.h), block blk's at toff[blk]; null in the exact kernel
  std::vector<size_t> toff;
  std::vector<int64_t> lbase;          // per block: its Λ base
  int p = 44;                          // fractional bits of the keys
  std::vector<int64_t> lD52;           // fixlog52(D) by position
  std::vector<int64_t> kD;             // ln D by position at 2^p (network grid)
  size_t call = 0;                     // this call's number, so a thread's tile strip of an earlier call is never reused
};

// Team buffers of one grid: a key (8 B) and a bucket id (2 B) per within-sample entry
inline size_t grid_bytes(const Cells &c){
  return c.cap * (sizeof(uint64_t) + sizeof(uint16_t));
}

// Value buckets: 0 below lo, 1..nb in key order, nb + 1 above. Equal width: bucket h + 1 holds rel = k - lo with rel >> t = h.
// Equal count (big grids): rel >> t is a cell and map gives its bucket, cells in order, so the buckets hold about as many
// sampled entries each. A bucket keeps low_bits of each key, bits(d) of them. Any partition in key order gives the same ranks.
struct Buckets {
  uint64_t lo = 0, mask = 0;
  int t = 0;
  size_t nb = 0, ncell = 0;
  const uint32_t *map = nullptr;
  const uint64_t *base = nullptr;
  const int *width = nullptr;

  // The bucket of key k, and in low the key bits it keeps (as low_bits)
  template <bool TABLE> size_t locate(uint64_t k, uint64_t &low) const {
    low = k;
    if(k < lo) return 0;
    const uint64_t h = (k - lo) >> t;
    if(h >= (TABLE ? ncell : nb)) return nb + 1;
    const size_t d = TABLE ? (size_t) map[h] : h + 1;
    low = rel_bits<TABLE>(k, d);
    return d;
  }
  template <bool TABLE> size_t locate(uint64_t k) const {
    uint64_t low;
    return locate<TABLE>(k, low);
  }
  // The key bits bucket d keeps of key k: the whole key in the catch-alls, else rel_bits
  template <bool TABLE> uint64_t low_bits(uint64_t k, size_t d) const {
    return d == 0 || d > nb ? k : rel_bits<TABLE>(k, d);
  }
  // Of bucket d in 1..nb: rel = k - lo minus the bucket's first key (equal width: the low t bits of rel)
  template <bool TABLE> uint64_t rel_bits(uint64_t k, size_t d) const {
    return TABLE ? k - lo - base[d] : (k - lo) & mask;
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

// Wait until done() holds: spin, then yield after SPIN rounds; throws when another thread has failed
template <class Pred> void spin_until(Pred done, const std::atomic<bool> &failed){
  for(unsigned spin = 0; !done(); spin += spin <= SPIN){
    if(failed.load()) throw std::runtime_error("another thread failed");
    if(spin > SPIN) std::this_thread::yield();
  }
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
    spin_until([&]{ return phase.load() != ph; }, failed);
  }

 private:
  const size_t n;
  std::atomic<size_t> count, phase;
};

// A sum unit: rows [i0, min(i1, s)) of each block in [b0, b1); head: the first row range of its blocks, last: the last one
struct Chunk {
  size_t b0, b1, i0, i1;
  bool head, last;
};

// A team of threads scores one grid at a time. The blocks are cut into chunks of about cap / size entries: runs of whole
// blocks, or row ranges of a big block; a chunk's entries are in position order (column by column, the chunk's rows of
// each). When the whole call has one thread, the team appends each entry to its bucket's list of slabs in one pass.
// Otherwise a first pass counts each chunk's entries per bucket and a second puts them in place, so a bucket is one run of
// the pool holding its chunks' entries in order, with no slack: smaller, and faster when several threads share the memory.
// Each chunk is a sum unit. Its rows' sender sums and, for whole blocks, its columns' receiver sums are complete in it; the row
// ranges of a block add their column sums to the team's (integers, so the order does not matter).
//
// Every entry's bucket is stored (2 bytes). Keys take 8 bytes, except in big grids (above the cache budget) on several
// threads: these keep the key bits of a bucket in 6 bytes (at most 48 bits; the catch-alls, and any wider bucket, keep 8 in
// ext), 8 bytes per entry in all, and their buckets hold equal counts, so each sorts in L2 while the threads share L3. A big grid
// in passes > 1 (multi) stores no bucket ids and keeps the keys of one range of buckets at a time: cap6 6-byte slots, about cap / passes
// plus an eighth for a range cut at a bucket's edge.
struct Team {
  Team(const Cells &c, size_t size, bool alone, size_t passes) : size(size), bar(size), step(std::max(SAMPLE, (c.cap >> SAMPLE_BITS) | 1)),
                                                                 a(c.node.size()), b(c.node.size()), ka(c.LamT ? c.node.size() : 0),
                                                                 kb(c.LamT ? c.node.size() : 0), lo(size), hi(size) {
    // a chunk holds at most 2 x target entries (a run of blocks) or target + s (rows of a block of s cells): only a grid of CHUNK_MAX
    // entries or more needs smaller chunks, so smaller grids are cut as before (one thread: one chunk, slabs)
    const size_t nblocks = c.block_start.size() - 1, even = (c.cap + size - 1) / size,
                 target = c.cap >= CHUNK_MAX ? std::min(even, CHUNK_MAX / 2) : even;
    size_t open = 0, filled = 0;
    for(size_t blk = 0; blk < nblocks; ++blk){
      const size_t s = c.block_start[blk + 1] - c.block_start[blk];
      if(s * s > target){
        if(filled) chunks.push_back({open, blk, 0, SIZE_MAX, true, true});
        const size_t k = std::min(s, (s * s + target - 1) / target);
        for(size_t q = 0; q < k; ++q) chunks.push_back({blk, blk + 1, s * q / k, s * (q + 1) / k, q == 0, q + 1 == k});
        open = blk + 1;
        filled = 0;
      } else {
        filled += s * s;
        if(filled >= target){
          chunks.push_back({open, blk + 1, 0, SIZE_MAX, true, true});
          open = blk + 1;
          filled = 0;
        }
      }
    }
    if(filled) chunks.push_back({open, nblocks, 0, SIZE_MAX, true, true});
    bucket_bits = std::max(B1_MIN, std::min(B1_MAX, bit_width(c.cap) - FILL_BITS));
    nbmax = ((size_t) 1 << bucket_bits) + 2;
    stride = nbmax;
    const size_t nchunks = chunks.size();
    slabbed = alone && nchunks == 1;
    big = !slabbed && grid_bytes(c) > CACHE_BUDGET;
    local = !slabbed && c.cap >= LOCAL_FROM;
    multi = big && passes > 1;
    if(slabbed){
      const size_t slabs = nbmax + c.cap / SLAB + 1;  // the one chunk holds every entry
      pool.reset(new uint64_t[slabs * SLAB]);
      link.resize(slabs);
      cur.resize(nbmax);
      full.resize(nbmax);
      bucket.resize(c.cap);
    } else {
      cnt.resize(nchunks * nbmax);
      tot.resize(nbmax);
      if(multi){
        cap6 = c.cap / passes + std::min<size_t>(c.cap / (8 * passes), (size_t) 1 << 26) + 1;
        pool6.reset(new char[cap6 * 6 + f2::SLACK]);
        hist.resize(size);
        msum.assign(2 * c.node.size(), 0);
      } else if(big){
        pool6.reset(new char[c.cap * 6 + f2::SLACK]);
        hist.resize(size);
      } else {
        pool.reset(new uint64_t[std::max<size_t>(1, c.cap)]);
      }
      if(!multi) bucket.resize(c.cap);
    }
    isum.reset(new std::atomic<uint64_t>[c.node.size()]());
    icnt.reset(new std::atomic<uint64_t>[c.node.size()]());
    split_cell.assign(c.node.size(), 0);
    for(const Chunk &ck : chunks){
      if(ck.head && ck.last) continue;
      has_split = true;
      for(size_t k = c.block_start[ck.b0]; k < c.block_start[ck.b1]; ++k) split_cell[k] = 1;
    }
  }

  const size_t size;
  Barrier bar;
  const size_t step;                   // sample every step-th entry of each block
  std::atomic<size_t> next_count{0}, next_split{0}, next_place{0}, next_bucket{0}, next_sum{0};
  size_t g = 0, nbmax, stride;
  int l = 0, r = 0, bucket_bits;        // bucket_bits: b1 of equal-width buckets before the key range is known
  bool slabbed, big;                   // big: runs of 6-byte keys with equal-count buckets (a grid above the cache budget on threads)
  bool local;                          // ranks counted from the bucket's first (a grid of LOCAL_FROM entries or more)
  bool multi = false;                  // a big grid in passes: ranges of buckets in turn, each entry's bucket recomputed (no ids)
  size_t cap6 = 0;                     // multi: entries pool6 holds
  std::vector<size_t> ranges;          // multi: bucket range r is [ranges[r], ranges[r + 1])
  std::vector<uint64_t> kmin, kend;    // multi: range r's keys are [kmin[r], kend[r]) (kend 0: no upper end)
  std::vector<uint64_t> msum;          // multi: sender rank sum and count of each position side by side, over the passes
  Buckets bk;                          // the grid's buckets, made by the leader
  std::vector<Chunk> chunks;
  std::unique_ptr<uint64_t[]> pool;    // slabs, or runs of a grid that fits the cache
  std::unique_ptr<char[]> pool6;       // big grids: 6 bytes per entry
  std::vector<size_t> link, cur, full;  // slabs: next slab of each slab, write position and full slabs of each bucket
  std::vector<uint32_t> cnt;           // runs: per chunk and bucket (stride apart), entries, then their offset in the bucket
  std::vector<size_t> tot;
  std::vector<uint16_t> bucket;        // the bucket of each entry (not in multi)
  std::vector<uint64_t> ext;           // big grids: the keys of the buckets wider than 48 bits (the catch-alls at least)
  std::vector<size_t> sbase;           // big grids: where each bucket's keys start, in pool6 (6-byte slots) or in ext
  std::vector<char> swide;             // big grids: the bucket's keys are in ext
  std::vector<uint32_t> map;           // equal-count buckets: the bucket of each cell
  std::vector<uint64_t> base;          // and the first key (rel) and key bits of each bucket
  std::vector<int> width;
  std::vector<std::vector<uint64_t>> hist;  // big grids: each member's sampled entries per cell
  std::vector<double> a, b;
  std::vector<int64_t> ka, kb;         // memory mode: ln a and ln b - ln D at 2^p by position
  std::vector<uint64_t> lo, hi;
  std::unique_ptr<std::atomic<uint64_t>[]> isum, icnt;  // receiver rank sums and counts of the columns of split blocks (multi: all)
  std::vector<char> split_cell;        // the cell's block is cut into row ranges
  bool has_split = false;
};

struct Work {
  std::vector<uint64_t> x, y;
  std::vector<uint32_t> z, cnt;
  std::vector<uint32_t> bnd, c2;       // a big bucket's top-digit offsets, a piece's second-digit offsets
  std::vector<uint64_t> pb;            // one piece's words
  std::vector<std::pair<size_t, size_t>> segs;
  std::vector<std::pair<uint64_t, uint64_t>> group;
  std::vector<size_t> start, at;
  std::vector<const char *> rp;
  std::vector<double> q;
  std::vector<uint64_t> kq;            // memory mode: a column's keys
  std::vector<int64_t> rk;
  std::vector<uint64_t> irc;           // integer row sums and counts side by side: one cache line per entry, where two arrays a page
                                       // multiple apart would alias
  std::vector<char> stage;             // big grids: each bucket's 64-byte line being filled
  std::vector<char *> cur;             // and the next byte of each bucket
  std::vector<uint16_t> bc;            // a column's buckets (multi: of its in-range entries, in ix order)
  std::vector<uint32_t> ix;            // multi: the rows of a column's in-range entries

  Work(const Cells &c, const Team &tm) : start(tm.nbmax + 1), at(tm.nbmax), rp(tm.nbmax), q(c.max_block), kq(c.LamT ? c.max_block : 0),
                                         rk(c.max_block), irc(2 * c.max_block), bc(tm.multi ? c.max_block : 0),
                                         ix(tm.multi ? c.max_block : 0) {}
};

// Byte offset in the pool of the rank of slab entry p: a slab's entries hold key bits until ranked, then the first half of the
// slab holds their ranks as uint32
inline size_t slab_rank_offset(size_t p){
  return (p - p % SLAB) * sizeof(uint64_t) + p % SLAB * sizeof(uint32_t);
}

inline uint32_t get_rank(const uint64_t *pool, size_t p){
  uint32_t rank;
  std::memcpy(&rank, (const char *) pool + slab_rank_offset(p), sizeof rank);
  return rank;
}

// Keys of a bucket held in one run of the pool; their ranks go to the run's first half as uint32
struct RunKeys {
  uint64_t *keys;

  void copy(uint64_t *x, size_t m) const {
    std::memcpy(x, keys, m * sizeof(uint64_t));
  }
  void hist(uint32_t *c, size_t m, int cut, int shift) const {
    for(size_t i = 0; i < m; ++i) ++c[(keys[i] >> cut >> shift) + 1];
  }
  void scatter(uint64_t *y, uint32_t *c, size_t m, int cut, int I, int shift) const {
    for(size_t i = 0; i < m; ++i){
      const uint64_t key = keys[i] >> cut;
      y[c[key >> shift]++] = (key << I) | i;
    }
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy((char *) keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(){}
};

// A 6-byte pool slot holds the low 48 bits of a key, least significant byte first on every platform (f2::put48; a 6-byte memcpy of a
// uint64_t would keep the high bytes on a big-endian machine), read by one 8-byte load (f2::get48: pool6 has f2::SLACK bytes of
// slack)
using f2::get48;
using f2::put48;

// Keys of a bucket held in one run of 6-byte slots (at most NARROW_BITS key bits); their ranks go to the run's first 4m bytes
struct SixKeys {
  char *keys;

  void copy(uint64_t *x, size_t m) const {
    for(size_t i = 0; i < m; ++i) x[i] = get48(keys + 6 * i);
  }
  void hist(uint32_t *c, size_t m, int cut, int shift) const {
    for(size_t i = 0; i < m; ++i) ++c[(get48(keys + 6 * i) >> cut >> shift) + 1];
  }
  void scatter(uint64_t *y, uint32_t *c, size_t m, int cut, int I, int shift) const {
    for(size_t i = 0; i < m; ++i){
      const uint64_t key = get48(keys + 6 * i) >> cut;
      y[c[key >> shift]++] = (key << I) | i;
    }
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy(keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(){}
};

// Keys of a bucket held in slab segments (offset, length); their ranks go to z, then to each slab's first half
struct SlabKeys {
  uint64_t *pool;
  const std::vector<std::pair<size_t, size_t>> &segs;
  uint32_t *z;

  void copy(uint64_t *x, size_t) const {
    for(const auto &sg : segs){
      std::memcpy(x, pool + sg.first, sg.second * sizeof(uint64_t));
      x += sg.second;
    }
  }
  void hist(uint32_t *c, size_t, int cut, int shift) const {
    for(const auto &sg : segs){
      const uint64_t *src = pool + sg.first;
      for(size_t j = 0; j < sg.second; ++j) ++c[(src[j] >> cut >> shift) + 1];
    }
  }
  void scatter(uint64_t *y, uint32_t *c, size_t, int cut, int I, int shift) const {
    size_t i = 0;
    for(const auto &sg : segs){
      const uint64_t *src = pool + sg.first;
      for(size_t j = 0; j < sg.second; ++j, ++i){
        const uint64_t key = src[j] >> cut;
        y[c[key >> shift]++] = (key << I) | i;
      }
    }
  }
  void put(size_t i, size_t rank){
    z[i] = (uint32_t) rank;
  }
  void finish(){
    size_t i = 0;
    for(const auto &sg : segs){
      std::memcpy((char *) pool + slab_rank_offset(sg.first), z + i, sg.second * sizeof(uint32_t));
      i += sg.second;
    }
  }
};

// Ranks of one bucket of m entries, held in position order as their low t key bits: sort (key bits, index) words, then
// give each entry its min rank - 1: s0 plus its rank in the bucket (s0 = the bucket's start, or 0 for local ranks). A bucket of
// more than 2^(PIECE_BITS + 2) entries is cut by a top digit into pieces (equal keys share one), then each piece is sorted by a
// second digit, deeper passes and an insertion sort, and ranked, in one reused buffer; a smaller one gets one radix pass of up
// to TOP_BITS straight from its keys, then sort_rest.
template <class K> void rank_bucket(K &keys, size_t m, size_t s0, int t, Work &w){
  const int I = bit_width(m - 1), cut = t + I > 64 ? t + I - 64 : 0, width = t - cut;
  const uint64_t imask = I ? (~(uint64_t) 0 >> (64 - I)) : 0;
  const bool pieces = cut == 0 && bit_width(m) > PIECE_BITS + 2;
  if(w.x.size() < m) w.x.resize(m);
  if(!pieces && w.y.size() < m) w.y.resize(m);
  uint64_t *x = w.x.data(), *y = w.y.data();
  if(pieces){
    const int b1 = std::min(width, bit_width(m) - PIECE_BITS), s1 = width - b1;
    const size_t nb1 = (size_t) 1 << b1;
    w.bnd.assign(nb1 + 1, 0);
    uint32_t *c1 = w.bnd.data();
    keys.hist(c1, m, 0, s1);
    for(size_t b = 0; b < nb1; ++b) c1[b + 1] += c1[b];
    keys.scatter(x, c1, m, 0, I, s1);
    w.cnt.clear();
    for(size_t b = 0, s = 0; b < nb1; ++b){
      const size_t e = c1[b], n = e - s;
      if(n == 0) continue;
      if(w.pb.size() < n) w.pb.resize(n);
      uint64_t *pc = w.pb.data(), *xs = x + s;
      size_t e0 = n;
      if(n == 1 || s1 == 0){
        std::memcpy(pc, xs, n * sizeof(uint64_t));
      } else {
        const int b2 = std::min(s1, bit_width(n)), s2 = s1 - b2, sh = I + s2;
        const uint64_t mask2 = ((uint64_t) 1 << b2) - 1;
        const size_t nb2 = (size_t) mask2 + 1;
        w.c2.assign(nb2 + 1, 0);
        uint32_t *c2 = w.c2.data();
        for(size_t q = 0; q < n; ++q) ++c2[((xs[q] >> sh) & mask2) + 1];
        uint32_t most = 0;
        for(size_t k = 0; k < nb2; ++k){
          most = std::max(most, c2[k + 1]);
          c2[k + 1] += c2[k];
        }
        for(size_t q = 0; q < n; ++q) pc[c2[(xs[q] >> sh) & mask2]++] = xs[q];
        for(size_t k = 0; k < nb2; ++k){
          if(c2[k] > 0){
            e0 = c2[k];
            break;
          }
        }
        if(most > INSERTION_MAX && s2 > 0){
          for(size_t k = 0, cs = 0; k < nb2; ++k){
            const size_t ce = c2[k];
            if(ce - cs > INSERTION_MAX){
              std::memcpy(xs + cs, pc + cs, (ce - cs) * sizeof(uint64_t));
              sort_words(xs + cs, pc + cs, ce - cs, I, w.cnt);
            }
            cs = ce;
          }
        }
      }
      if(n > 1 && s1 > 0) insertion_sort_sentinel(pc, n, e0);
      uint64_t prev = pc[0] >> I, first = s0 + s;
      for(size_t v = 0; v < n; ++v){
        if((pc[v] >> I) != prev){
          prev = pc[v] >> I;
          first = s0 + s + v;
        }
        keys.put(pc[v] & imask, first);
      }
      s = e;
    }
    keys.finish();
    return;
  }
  const int want = std::min(TOP_BITS, bit_width(m)), shift = width > want ? width - want : 0;
  const size_t nb = (size_t) 1 << (width - shift);
  w.cnt.assign(nb + 1, 0);
  uint32_t *c = w.cnt.data();
  keys.hist(c, m, cut, shift);
  uint32_t most = 0;
  for(size_t b = 0; b < nb; ++b){
    most = std::max(most, c[b + 1]);
    c[b + 1] += c[b];
  }
  keys.scatter(y, c, m, cut, I, shift);
  if(m > 1 && shift > 0) sort_rest(x, y, m, I, nb, most, 0, w.cnt);
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
  keys.finish();
}

// Memory mode: Λ[i, j] - lbase of block blk from its tiles (f2.h); a diagonal tile is stored whole, so only tiles below it swap
inline int64_t tile_value(const Cells &c, size_t blk, size_t i, size_t j){
  const size_t TL = f2::TILE;
  if(i / TL > j / TL) std::swap(i, j);
  const size_t I = i / TL, J = j / TL;
  return (int64_t) get48(c.LamT + c.toff[blk] + f2::TILE_BYTES * (J * (J + 1) / 2 + I) + 6 * (TL * (j % TL) + i % TL));
}

// Memory mode: a thread's strip of the TILE columns of tile column J of a block for rows [i0, i1), column-major, made from whole
// tiles (tile (I, J) at and above the diagonal, tile (J, I) transposed below it): each tile is read once for its TILE columns
struct TileStrip {
  size_t call = 0, blk = SIZE_MAX, J = SIZE_MAX, i0 = 0, i1 = 0;
  std::vector<int64_t> v;
};
thread_local TileStrip tile_strip;
std::atomic<size_t> tile_calls{0};      // memory-mode calls so far (Cells::call)

// Memory mode: Λ - lbase of rows [i0, i0 + m) of column jj of block blk, in the thread's strip (0 on the diagonal: the tiles hold 0)
const int64_t *tile_strip_column(const Cells &c, size_t blk, size_t jj, size_t i0, size_t m){
  const size_t TL = f2::TILE, s = c.block_start[blk + 1] - c.block_start[blk], i1 = i0 + m, J = jj / TL;
  TileStrip &st = tile_strip;
  if(st.call != c.call || st.blk != blk || st.J != J || st.i0 != i0 || st.i1 != i1){
    st.call = c.call;
    st.blk = blk;
    st.J = J;
    st.i0 = i0;
    st.i1 = i1;
    if(st.v.size() < TL * m) st.v.resize(TL * m);
    const char *T = c.LamT + c.toff[blk];
    const size_t ncol = std::min(TL, s - TL * J);
    for(size_t I = i0 / TL; TL * I < i1; ++I){
      const size_t r0 = std::max(i0, TL * I), r1 = std::min(i1, TL * I + TL);
      if(I <= J){
        const char *t = T + f2::TILE_BYTES * (J * (J + 1) / 2 + I);
        for(size_t b = 0; b < ncol; ++b){
          int64_t *col = st.v.data() + b * m;
          for(size_t r = r0; r < r1; ++r) col[r - i0] = (int64_t) get48(t + 6 * (TL * b + r % TL));
        }
      } else {
        const char *t = T + f2::TILE_BYTES * (I * (I + 1) / 2 + J);
        for(size_t b = 0; b < ncol; ++b){
          int64_t *col = st.v.data() + b * m;
          for(size_t r = r0; r < r1; ++r) col[r - i0] = (int64_t) get48(t + 6 * (TL * (r % TL) + b));
        }
      }
    }
  }
  return st.v.data() + (jj % TL) * m;
}

// Memory mode: the keys of rows [i0, i0 + m) of column jj of block blk into k, by integer adds: A[i] + B'[j] + Λ[i, j] (the network
// grid: ln D[j] - Λ[i, j], the order of F), signed, then flipped to unsigned order
const uint64_t *lambda_keys(const Cells &c, const Team &tm, size_t blk, size_t jj, size_t i0, size_t m, uint64_t *k){
  const size_t first = c.block_start[blk];
  const uint64_t SIGN = (uint64_t) 1 << 63;
  const int64_t *L = tile_strip_column(c, blk, jj, i0, m), base = c.lbase[blk];
  if(tm.l < 0){
    const int64_t dj = c.kD[first + jj] - base;
    for(size_t i = 0; i < m; ++i) k[i] = (uint64_t) (dj - L[i]) ^ SIGN;
    return k;
  }
  const int64_t *ka = tm.ka.data() + first + i0, bj = tm.kb[first + jj] + base;
  for(size_t i = 0; i < m; ++i) k[i] = (uint64_t) (ka[i] + bj + L[i]) ^ SIGN;
  return k;
}

// Memory mode: the key of entry (ii, jj) of block blk
inline uint64_t lambda_key(const Cells &c, const Team &tm, size_t blk, size_t ii, size_t jj){
  const size_t first = c.block_start[blk];
  const int64_t L = c.lbase[blk] + tile_value(c, blk, ii, jj);
  const int64_t v = tm.l < 0 ? c.kD[first + jj] - L : tm.ka[first + ii] + tm.kb[first + jj] + L;
  return (uint64_t) v ^ ((uint64_t) 1 << 63);
}

// Member me of the team: a and b of its share of the cells (memory mode: their keys ka, kb)
void prepare_ab(const Cells &c, Team &tm, size_t me){
  if(tm.l < 0) return;
  const size_t n = c.node.size(), k1 = n * (me + 1) / tm.size;
  if(c.F == nullptr){
    for(size_t k = n * me / tm.size; k < k1; ++k){
      tm.a[k] = c.A[(size_t) c.out[k] + (size_t) tm.l * n];
      tm.b[k] = c.B[(size_t) tm.r + (size_t) c.out[k] * c.nrec];
    }
    if(c.LamT){
      const int sh = 52 - c.p;
      for(size_t k = n * me / tm.size; k < k1; ++k){
        tm.ka[k] = f2::round_shift(f2::fixlog52(tm.a[k]), sh);
        tm.kb[k] = f2::round_shift(f2::fixlog52(tm.b[k]) - c.lD52[k], sh);
      }
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

// fn(key) for member me's share of every step-th entry of each block from entry SAMPLE_FIRST, those not skipped
template <class Fn> void each_sample(const Cells &c, const Team &tm, size_t me, Fn fn){
  const size_t nblocks = c.block_start.size() - 1;
  for(size_t blk = 0; blk < nblocks; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
    const size_t step = tm.step, ns = s * s > SAMPLE_FIRST ? (s * s - SAMPLE_FIRST - 1) / step + 1 : 0, k1 = ns * (me + 1) / tm.size;
    const double *a = tm.a.data() + first, *b = tm.b.data() + first;
    for(size_t k = ns * me / tm.size; k < k1; ++k){
      const size_t q = SAMPLE_FIRST + k * step, jj = q / s, ii = q % s;
      if(c.LamT){
        if(ii != jj) fn(lambda_key(c, tm, blk, ii, jj));
        continue;
      }
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
  std::vector<uint64_t> &h = tm.hist[me];
  h.assign(ncell, 0u);
  each_sample(c, tm, me, [&](uint64_t k){
    ++h[(k - lo) >> tc];
  });
}

// The leader: equal-count buckets from the members' cell counts. A bucket starts where the sampled entries below it reach the
// next of target equal shares, and spans at most 2^(NARROW_BITS - tc) cells so its keys fit 6 bytes. False if nothing was sampled.
bool count_buckets(Team &tm, size_t cap, uint64_t lo, int tc, size_t ncell){
  std::vector<uint64_t> &h = tm.hist[0];
  for(size_t me = 1; me < tm.size; ++me){
    const uint64_t *hm = tm.hist[me].data();
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

// The leader: the grid's buckets for everyone in tm.bk, equal count for a big grid (at most BUCKETS_IDS; histogram cells of
// 2^tc keys from lo), else equal width (2^b1 buckets; for a big grid at least 2^(R - NARROW_BITS), at most 2^B1_LIMIT, so the
// keys fit 6 bytes unless R is 64), and the per-bucket counts sized
void make_buckets(Team &tm, size_t cap, bool sampled, uint64_t lo, int R, int tc, size_t ncell){
  const size_t nchunks = tm.chunks.size();
  bool table = false;
  if(tm.big && sampled) table = count_buckets(tm, cap, lo, tc, ncell) && tm.bk.nb <= BUCKETS_IDS;
  if(!table){
    const int b1 = tm.big ? std::min(R, std::min(B1_LIMIT, std::max(tm.bucket_bits, R - NARROW_BITS))) : std::min(R, tm.bucket_bits);
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

// One thread: Q of each entry in position order, appended to its bucket's list of slabs. Kept out of line like pull_slabs and
// pull_runs: inlined into run_team, their loops lose registers to it (one thread ~4% slower, GCC 13 -O2).
__attribute__((noinline)) void fill_slabs(const Cells &c, Team &tm, const Buckets &bk, Work &w){
  const size_t NB = bk.nb + 2;
  const bool net = tm.l < 0;
  const Chunk &ck = tm.chunks[0];
  uint64_t *pool = tm.pool.get();
  size_t *cur = tm.cur.data(), *full = tm.full.data(), *link = tm.link.data();
  size_t slabs = NB;
  for(size_t d = 0; d < NB; ++d){
    cur[d] = d * SLAB;
    full[d] = 0;
  }
  // key k of an entry not skipped: its bucket to bc, its kept bits to the bucket's slab list
  auto append = [&](uint16_t &bc, uint64_t k){
    uint64_t low;
    const size_t d = bk.locate<false>(k, low);
    bc = (uint16_t) d;
    size_t p = cur[d];
    pool[p++] = low;
    if(p % SLAB == 0){
      link[p / SLAB - 1] = slabs;
      p = slabs++ * SLAB;
      ++full[d];
    }
    cur[d] = p;
  };
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const size_t pos = c.block_off[blk] + jj * s + i0;
      uint16_t *bc = tm.bucket.data() + pos;
      if(c.LamT){
        const uint64_t *kk = lambda_keys(c, tm, blk, jj, i0, m, w.kq.data());
        for(size_t i = 0; i < m; ++i){
          if(i0 + i == jj) bc[i] = SKIP;
          else append(bc[i], kk[i]);
        }
        continue;
      }
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      for(size_t i = 0; i < m; ++i){
        const double v = q[i];
        if(skipped(i0 + i, jj, v)) bc[i] = SKIP;
        else append(bc[i], sortable_key(v));
      }
    }
  }
}

// First pass over chunk ch: the bucket of each entry (SKIP if skipped; multi: to scratch only) and the chunk's entries per bucket
template <bool TABLE> void count_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  const Chunk &ck = tm.chunks[ch];
  uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  std::fill(cnt, cnt + bk.nb + 2, 0u);
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      uint16_t *bc = tm.multi ? w.bc.data() : tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      const size_t diag = jj - i0;       // wraps above m when the diagonal is not in the chunk
      if(c.LamT){
        const uint64_t *kk = lambda_keys(c, tm, blk, jj, i0, m, w.kq.data());
        for(size_t i = 0; i < m; ++i){
          if(i == diag){
            bc[i] = SKIP;
            continue;
          }
          const size_t d = bk.locate<TABLE>(kk[i]);
          bc[i] = (uint16_t) d;
          ++cnt[d];
        }
        continue;
      }
      const double *q = block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      for(size_t i = 0; i < m; ++i){
        uint64_t u;
        std::memcpy(&u, q + i, sizeof u);
        if(zero_or_nan(u) || i == diag){
          bc[i] = SKIP;
          continue;
        }
        const size_t d = bk.locate<TABLE>(flip_key(u));
        bc[i] = (uint16_t) d;
        ++cnt[d];
      }
    }
  }
}

// Each bucket's chunk counts to the offsets of the chunks' entries in the bucket, SPLIT_STEP buckets at a time; tot = sizes
void split_buckets(Team &tm, size_t NB){
  const size_t nchunks = tm.chunks.size();
  size_t run[SPLIT_STEP];
  for(size_t d0; (d0 = tm.next_split.fetch_add(SPLIT_STEP)) < NB; ){
    const size_t d1 = std::min(NB, d0 + SPLIT_STEP);
    std::fill(run, run + SPLIT_STEP, (size_t) 0);
    for(size_t ch = 0; ch < nchunks; ++ch){
      uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
      for(size_t d = d0; d < d1; ++d){
        const uint32_t e = cnt[d];
        cnt[d] = (uint32_t) run[d - d0];
        run[d - d0] += e;
      }
    }
    for(size_t d = d0; d < d1; ++d){
      if(run[d - d0] > BUCKET_MAX) throw std::runtime_error(BUCKET_ERROR);
      tm.tot[d] = run[d - d0];
    }
  }
}

// Second pass over chunk ch: the low key bits of each entry to its place in its bucket (a grid that fits the cache)
void place_chunk_ids(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  const Chunk &ck = tm.chunks[ch];
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  size_t *at = w.at.data();
  for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = w.start[d] + cnt[d];
  uint64_t *pool = tm.pool.get();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const size_t pos = c.block_off[blk] + jj * s + i0;
      const uint64_t *kk = c.LamT ? lambda_keys(c, tm, blk, jj, i0, m, w.kq.data()) : nullptr;
      const double *q = kk ? nullptr
                           : block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      const uint16_t *bc = tm.bucket.data() + pos;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        pool[at[d]++] = bk.low_bits<false>(kk ? kk[i] : sortable_key(q[i]), d);
      }
    }
  }
}

// Second pass over chunk ch of a big grid: the kept key bits of each entry to its place in its bucket, 6 bytes in pool6 or 8 in
// ext (sbase / swide)
template <bool TABLE> void place_chunk_six(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const bool net = tm.l < 0;
  const Chunk &ck = tm.chunks[ch];
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  const size_t NB = bk.nb + 2;
  size_t *at = w.at.data();
  for(size_t d = 0; d < NB; ++d) at[d] = tm.sbase[d] + cnt[d];
  const char *wide = tm.swide.data();
  char *pool = tm.pool6.get();
  uint64_t *ext = tm.ext.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const uint64_t *kk = c.LamT ? lambda_keys(c, tm, blk, jj, i0, m, w.kq.data()) : nullptr;
      const double *q = kk ? nullptr
                           : block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t low = bk.low_bits<TABLE>(kk ? kk[i] : sortable_key(q[i]), d);
        if(wide[d]) ext[at[d]++] = low;
        else put48(pool + 6 * at[d]++, low);
      }
    }
  }
}

// place_chunk_six with cache-line write-combining (x86). Each narrow bucket's buffer holds the bytes of the 64-byte line of pool6
// it is filling; a line that lies wholly inside this chunk's run of the bucket goes out with non-temporal stores (no read for
// ownership, not kept in cache), a line it shares with a neighbouring run (another chunk or bucket, maybe another thread) gets only
// its own bytes with ordinary stores. A bucket's keys keep their order.
template <bool TABLE> void place_chunk_six_nt(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
#ifdef HYPERCOM_X86
  const bool net = tm.l < 0;
  const Chunk &ck = tm.chunks[ch];
  const size_t NB = bk.nb + 2;
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  char *pool = tm.pool6.get();
  uint64_t *ext = tm.ext.data();
  const char *wide = tm.swide.data();
  size_t *at = w.at.data();
  if(w.stage.size() < NB * 64 + 64) w.stage.resize(NB * 64 + 64);
  char *stage = (char *) (((uintptr_t) w.stage.data() + 63) & ~(uintptr_t) 63);
  w.rp.resize(std::max(w.rp.size(), NB));
  const char **r0 = w.rp.data();           // first byte of each bucket's run in this chunk
  std::vector<char *> &cur = w.cur;      // next byte to write of each narrow bucket
  cur.resize(NB);
  for(size_t d = 0; d < NB; ++d){
    at[d] = tm.sbase[d] + cnt[d];
    cur[d] = pool + 6 * at[d];
    r0[d] = cur[d];
  }
  auto flush = [&](size_t d, char *line, const char *buf){
    if(line >= r0[d]){
      const __m128i *src = (const __m128i *) buf;
      __m128i *dst = (__m128i *) line;
      _mm_stream_si128(dst, _mm_load_si128(src));
      _mm_stream_si128(dst + 1, _mm_load_si128(src + 1));
      _mm_stream_si128(dst + 2, _mm_load_si128(src + 2));
      _mm_stream_si128(dst + 3, _mm_load_si128(src + 3));
    } else {
      std::memcpy((char *) r0[d], buf + (r0[d] - line), line + 64 - r0[d]);
    }
  };
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const uint64_t *kk = c.LamT ? lambda_keys(c, tm, blk, jj, i0, m, w.kq.data()) : nullptr;
      const double *q = kk ? nullptr
                           : block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t low = bk.low_bits<TABLE>(kk ? kk[i] : sortable_key(q[i]), d);
        if(wide[d]){
          ext[at[d]++] = low;
          continue;
        }
        char *o = cur[d], *buf = stage + 64 * d;
        const size_t off = (uintptr_t) o & 63;
        if(off + 6 <= 64){
          put48(buf + off, low);
          if(off + 6 == 64) flush(d, o - off, buf);
        } else {
          char b6[8];
          put48(b6, low);
          const size_t n1 = 64 - off;
          std::memcpy(buf + off, b6, n1);
          flush(d, o - off, buf);
          std::memcpy(buf, b6 + n1, 6 - n1);
        }
        cur[d] = o + 6;
      }
    }
  }
  for(size_t d = 0; d < NB; ++d){
    if(wide[d]) continue;
    char *o = cur[d], *line = (char *) ((uintptr_t) o & ~(uintptr_t) 63);
    if(o == line) continue;
    const char *from = std::max<const char *>(line, r0[d]);
    if(from < o) std::memcpy((char *) from, stage + 64 * d + (from - line), o - from);
  }
  _mm_sfence();
#else
  place_chunk_six<TABLE>(c, tm, bk, ch, w);
#endif
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

// Rank the buckets, a run of about E / (RANK_UNITS x team size) entries at a time; start[d] = number of entries in the buckets
// below d
void rank_buckets(Team &tm, const Buckets &bk, Work &w){
  const size_t NB = bk.nb + 2;
  const size_t *start = w.start.data(), *full = tm.full.data(), *link = tm.link.data();
  const size_t unit = std::max<size_t>(1, start[NB] / (RANK_UNITS * tm.size)), nunits = (start[NB] + unit - 1) / unit;
  for(size_t u; (u = tm.next_bucket.fetch_add(1)) < nunits; ){
    const size_t d1 = std::lower_bound(start, start + NB, (u + 1) * unit) - start;
    for(size_t d = std::lower_bound(start, start + NB, u * unit) - start; d < d1; ++d){
      const size_t m = start[d + 1] - start[d];
      if(m == 0) continue;
      if(m > BUCKET_MAX) throw std::runtime_error(BUCKET_ERROR);
      const size_t s0 = tm.local ? 0 : start[d];
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
        rank_bucket(keys, m, s0, t, w);
      } else if(tm.swide[d]){
        RunKeys keys = {tm.ext.data() + tm.sbase[d]};
        rank_bucket(keys, m, s0, t, w);
      } else {
        SixKeys keys = {tm.pool6.get() + 6 * tm.sbase[d]};
        rank_bucket(keys, m, s0, t, w);
      }
    }
  }
}

// Ranks of m entries in position order from their buckets' slab lists; at[d] = next slab entry of bucket d
__attribute__((noinline)) void pull_slabs(const uint16_t *bc, size_t m, int64_t *rk, const uint64_t *pool, const size_t *link, size_t *at){
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

// The same from their buckets' runs; rp[d] = the next uint32 rank of bucket d, plus its start base[d] when ranks are LOCAL
template <bool LOCAL> __attribute__((noinline)) void pull_runs(const uint16_t *bc, size_t m, int64_t *rk, const char **rp,
                                                         const size_t *base){
  for(size_t i = 0; i < m; ++i){
    const size_t d = bc[i];
    if(d == SKIP){
      rk[i] = -1;
      continue;
    }
    uint32_t rank;
    std::memcpy(&rank, rp[d], sizeof rank);
    rp[d] += sizeof rank;
    rk[i] = LOCAL ? (int64_t) (base[d] + rank) : (int64_t) rank;
  }
}

// A score as one correctly rounded division of its exact integer rank sum, (sum of min rank - 1) / (count x (K - 1)); exact operands
// while count x K < 2^53 (a sum fits uint64 while the cells of a sample, cubed, stay below 2^64: 2.6 million cells)
inline double rank_mean(uint64_t sum, uint64_t cnt, double denom){
  return (double) sum / ((double) cnt * denom);
}

// Sender and receiver sums of chunk ch as exact integers, each column's ranks taken from the chunk's entries of each bucket in
// position order: rows of the chunk complete here, columns too when the chunk holds whole blocks, else added to the team's column
// sums for run_team to finish
template <Store STORE> void sum_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, double denom, Work &w, double *sender,
                                      double *receiver){
  const uint64_t *pool = tm.pool.get();
  const size_t *link = tm.link.data();
  size_t *at = w.at.data();
  const char **rp = w.rp.data();
  if(STORE == Store::SLABS){
    for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = d * SLAB;
  } else {
    const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
    const char *ext = (const char *) tm.ext.data(), *pool6 = tm.pool6.get();
    for(size_t d = 0; d < bk.nb + 2; ++d){
      const char *keys = tm.big ? (tm.swide[d] ? ext + sizeof(uint64_t) * tm.sbase[d] : pool6 + 6 * tm.sbase[d])
                                : (const char *) (pool + w.start[d]);
      rp[d] = keys + sizeof(uint32_t) * cnt[d];
    }
  }
  const Chunk &ck = tm.chunks[ch];
  const bool whole = ck.head && ck.last;
  uint64_t *sc = w.irc.data();          // row sum and count of row k at 2 k, 2 k + 1
  int64_t *rk = w.rk.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first, off = c.block_off[blk];
    const size_t i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    const int *out = c.out.data() + first;
    std::fill(sc, sc + 2 * m, (uint64_t) 0);
    for(size_t jj = 0; jj < s; ++jj){
      const uint16_t *bc = tm.bucket.data() + off + jj * s + i0;
      if(STORE == Store::SLABS) pull_slabs(bc, m, rk, pool, link, at);
      else if(tm.local) pull_runs<true>(bc, m, rk, rp, w.start.data());
      else pull_runs<false>(bc, m, rk, rp, nullptr);
      uint64_t cs = 0, cn = 0;
      for(size_t i = 0; i < m; ++i){
        const int64_t k = rk[i];
        const uint64_t keep = ~(uint64_t) (k >> 63), v = (uint64_t) k & keep, one = keep & 1;
        sc[2 * i] += v;
        sc[2 * i + 1] += one;
        cs += v;
        cn += one;
      }
      if(whole){
        receiver[out[jj]] = rank_mean(cs, cn, denom);
      } else {
        tm.isum[first + jj].fetch_add(cs, std::memory_order_relaxed);
        tm.icnt[first + jj].fetch_add(cn, std::memory_order_relaxed);
      }
    }
    for(size_t k = 0; k < m; ++k) sender[out[i0 + k]] = rank_mean(sc[2 * k], sc[2 * k + 1], denom);
  }
}

// Multi: the in-range entries of rows [i0, i0 + m) of column jj of block blk, keys in [kmin, kend) (kend 0: no upper end) and not
// skipped, made as count_chunk makes them: their rows to w.ix and buckets to w.bc without a branch (each entry written, the count
// advanced by the in-range flag); returns how many. kk / q: the column's keys (memory mode) or quotients
template <bool TABLE> size_t column_buckets_in(const Cells &c, Team &tm, const Buckets &bk, size_t blk, size_t jj, size_t i0, size_t m,
                                               Work &w, uint64_t kmin, uint64_t kend, const uint64_t *&kk, const double *&q){
  const bool net = tm.l < 0;
  const size_t first = c.block_start[blk];
  kk = c.LamT ? lambda_keys(c, tm, blk, jj, i0, m, w.kq.data()) : nullptr;
  q = kk ? nullptr : block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data());
  uint16_t *bc = w.bc.data();
  uint32_t *ix = w.ix.data();
  const size_t diag = jj - i0;         // wraps above m when the diagonal is not in the chunk
  const uint64_t top = kend ? kend : ~(uint64_t) 0;  // kend 0: every key from kmin (the top key is never one: flipped keys)
  size_t n = 0;
  if(kk){
    for(size_t i = 0; i < m; ++i){
      const uint64_t k = kk[i];
      ix[n] = (uint32_t) i;
      bc[n] = (uint16_t) bk.locate<TABLE>(k);
      n += (size_t) ((k >= kmin) & (k < top) & (i != diag));
    }
    return n;
  }
  for(size_t i = 0; i < m; ++i){
    uint64_t u;
    std::memcpy(&u, q + i, sizeof u);
    const bool skip = zero_or_nan(u) || i == diag;
    const uint64_t k = flip_key(u);
    ix[n] = (uint32_t) i;
    bc[n] = (uint16_t) bk.locate<TABLE>(k);
    n += (size_t) ((k >= kmin) & (k < top) & !skip);
  }
  return n;
}

size_t column_buckets_any(const Cells &c, Team &tm, const Buckets &bk, size_t blk, size_t jj, size_t i0, size_t m, Work &w, uint64_t kmin,
                          uint64_t kend, const uint64_t *&kk, const double *&q){
  if(bk.map) return column_buckets_in<true>(c, tm, bk, blk, jj, i0, m, w, kmin, kend, kk, q);
  return column_buckets_in<false>(c, tm, bk, blk, jj, i0, m, w, kmin, kend, kk, q);
}

// Multi: the first key of bucket d (0 for the low catch-all; nb + 1, the high catch-all, starts after the last cell); the end of the
// keys of the buckets below d1 is range_first(d1), 0 for d1 = nb + 2 (no upper end)
uint64_t range_first(const Buckets &bk, size_t d){
  if(d == 0) return 0;
  if(d >= bk.nb + 2) return 0;
  const size_t ncells = bk.map ? bk.ncell : bk.nb;
  size_t h = d - 1;                    // equal width: bucket d is cell d - 1
  if(bk.map){
    h = 0;
    while(h < ncells && (size_t) bk.map[h] < d) ++h;
  }
  return bk.lo + ((uint64_t) std::min(h, ncells) << bk.t);
}

// Multi, the leader: the buckets cut into ranges of at most cap6 narrow entries and each range's keys. Wide buckets go to ext and do not
// count, so no cut falls at an empty or a wide bucket (the catch-alls among them)
void cut_ranges(Team &tm, const Buckets &bk){
  const size_t NB = bk.nb + 2;
  tm.ranges.assign(1, 0);
  size_t in = 0;
  for(size_t d = 0; d < NB; ++d){
    const size_t e = bk.bits(d) > NARROW_BITS ? 0 : tm.tot[d];
    if(e > tm.cap6) throw std::runtime_error("a value bucket holds more entries than the key pool of one pass; give the exact kernels' memory with the memory argument");
    if(in + e > tm.cap6){
      tm.ranges.push_back(d);
      in = 0;
    }
    in += e;
  }
  tm.ranges.push_back(NB);
  tm.kmin.clear();
  tm.kend.clear();
  for(size_t r = 0; r + 1 < tm.ranges.size(); ++r){
    tm.kmin.push_back(range_first(bk, tm.ranges[r]));
    tm.kend.push_back(range_first(bk, tm.ranges[r + 1]));
  }
}

// Multi, the leader: six_storage for buckets [d0, d1) only
void six_storage_range(Team &tm, const Buckets &bk, size_t d0, size_t d1){
  const size_t NB = bk.nb + 2;
  tm.sbase.assign(NB, 0);
  tm.swide.assign(NB, 0);
  size_t narrow = 0, wide = 0;
  for(size_t d = d0; d < d1; ++d){
    tm.swide[d] = bk.bits(d) > NARROW_BITS;
    tm.sbase[d] = tm.swide[d] ? wide : narrow;
    (tm.swide[d] ? wide : narrow) += tm.tot[d];
  }
  tm.ext.resize(wide);
}

// Multi: place_chunk_six for the entries of chunk ch in buckets [d0, d1) (keys [kmin, kend)), their buckets recomputed
void place_chunk_range(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w, size_t d0, size_t d1, uint64_t kmin, uint64_t kend){
  const Chunk &ck = tm.chunks[ch];
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  size_t *at = w.at.data();
  for(size_t d = d0; d < d1; ++d) at[d] = tm.sbase[d] + cnt[d];
  const char *wide = tm.swide.data();
  char *pool = tm.pool6.get();
  uint64_t *ext = tm.ext.data();
  const uint16_t *bc = w.bc.data();
  const uint32_t *ix = w.ix.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const uint64_t *kk;
      const double *q;
      const size_t n = column_buckets_any(c, tm, bk, blk, jj, i0, m, w, kmin, kend, kk, q);
      for(size_t t = 0; t < n; ++t){
        const size_t d = bc[t], i = ix[t];
        const uint64_t key = kk ? kk[i] : sortable_key(q[i]);
        const uint64_t low = bk.map ? bk.low_bits<true>(key, d) : bk.low_bits<false>(key, d);
        if(wide[d]) ext[at[d]++] = low;
        else put48(pool + 6 * at[d]++, low);
      }
    }
  }
}

// Multi: rank_buckets over buckets [d0, d1)
void rank_range(Team &tm, const Buckets &bk, Work &w, size_t d0, size_t d1){
  const size_t *start = w.start.data(), E = start[d1] - start[d0];
  const size_t unit = std::max<size_t>(1, E / (RANK_UNITS * tm.size)), nunits = (E + unit - 1) / unit;
  for(size_t u; (u = tm.next_bucket.fetch_add(1)) < nunits; ){
    const size_t e1 = std::lower_bound(start + d0, start + d1, start[d0] + (u + 1) * unit) - start;
    for(size_t d = std::lower_bound(start + d0, start + d1, start[d0] + u * unit) - start; d < e1; ++d){
      const size_t m = start[d + 1] - start[d];
      if(m == 0) continue;
      if(m > BUCKET_MAX) throw std::runtime_error(BUCKET_ERROR);
      const size_t s0 = tm.local ? 0 : start[d];
      const int t = bk.bits(d);
      if(tm.swide[d]){
        RunKeys keys = {tm.ext.data() + tm.sbase[d]};
        rank_bucket(keys, m, s0, t, w);
      } else {
        SixKeys keys = {tm.pool6.get() + 6 * tm.sbase[d]};
        rank_bucket(keys, m, s0, t, w);
      }
    }
  }
}

// Multi: the integer rank sums of chunk ch's entries in buckets [d0, d1), their buckets recomputed; rows to msum (a chunk's rows are
// its own), columns to the team's column sums (isum / icnt), all divided after the last pass
void sum_chunk_range(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w, size_t d0, size_t d1, uint64_t kmin, uint64_t kend){
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  const char *ext = (const char *) tm.ext.data(), *pool6 = tm.pool6.get();
  const char **rp = w.rp.data();
  for(size_t d = d0; d < d1; ++d){
    rp[d] = (tm.swide[d] ? ext + sizeof(uint64_t) * tm.sbase[d] : pool6 + 6 * tm.sbase[d]) + sizeof(uint32_t) * cnt[d];
  }
  const Chunk &ck = tm.chunks[ch];
  const size_t *base = w.start.data();
  const uint16_t *bc = w.bc.data();
  const uint32_t *ix = w.ix.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
    const size_t i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    uint64_t *sc = tm.msum.data() + 2 * (first + i0);  // row sum and count of row i at 2 i, 2 i + 1
    for(size_t jj = 0; jj < s; ++jj){
      const uint64_t *kk;
      const double *q;
      const size_t n = column_buckets_any(c, tm, bk, blk, jj, i0, m, w, kmin, kend, kk, q);
      uint64_t cs = 0;
      for(size_t t = 0; t < n; ++t){
        const size_t d = bc[t], i = ix[t];
        uint32_t rank;
        std::memcpy(&rank, rp[d], sizeof rank);
        rp[d] += sizeof rank;
        const uint64_t v = tm.local ? base[d] + rank : rank;
        sc[2 * i] += v;
        ++sc[2 * i + 1];
        cs += v;
      }
      tm.isum[first + jj].fetch_add(cs, std::memory_order_relaxed);
      tm.icnt[first + jj].fetch_add(n, std::memory_order_relaxed);
    }
  }
}

// Multi, one team member after the counts and their offsets: each range of buckets placed, ranked and summed in turn, then the
// member's share of the scores, one division each
void run_passes(const Cells &c, Team &tm, size_t me, const Buckets &bk, Work &w, const std::atomic<bool> &failed, double *snd, double *rcv){
  const size_t n = c.node.size(), nchunks = tm.chunks.size(), NB = bk.nb + 2;
  if(me == 0) cut_ranges(tm, bk);
  tm.bar.wait(failed);
  for(size_t r = 0; r + 1 < tm.ranges.size(); ++r){
    const size_t d0 = tm.ranges[r], d1 = tm.ranges[r + 1];
    if(me == 0){
      six_storage_range(tm, bk, d0, d1);
      tm.next_place = 0;
      tm.next_bucket = 0;
      tm.next_sum = 0;
    }
    tm.bar.wait(failed);
    for(size_t ch; (ch = tm.next_place.fetch_add(1)) < nchunks; ) place_chunk_range(c, tm, bk, ch, w, d0, d1, tm.kmin[r], tm.kend[r]);
    tm.bar.wait(failed);
    rank_range(tm, bk, w, d0, d1);
    tm.bar.wait(failed);
    for(size_t ch; (ch = tm.next_sum.fetch_add(1)) < nchunks; ) sum_chunk_range(c, tm, bk, ch, w, d0, d1, tm.kmin[r], tm.kend[r]);
    tm.bar.wait(failed);
  }
  const double denom = (double) w.start[NB] - 1.0;
  for(size_t k = n * me / tm.size, k1 = n * (me + 1) / tm.size; k < k1; ++k){
    snd[c.out[k]] = rank_mean(tm.msum[2 * k], tm.msum[2 * k + 1], denom);
    rcv[c.out[k]] = rank_mean(tm.isum[k].load(std::memory_order_relaxed), tm.icnt[k].load(std::memory_order_relaxed), denom);
    tm.msum[2 * k] = tm.msum[2 * k + 1] = 0;
    tm.isum[k].store(0, std::memory_order_relaxed);
    tm.icnt[k].store(0, std::memory_order_relaxed);
  }
}

// After each grid: frees per-thread sort buffers above KEEP_SORT_BUFFER entries
void trim_work(Work &w){
  if(w.x.size() > KEEP_SORT_BUFFER){
    std::vector<uint64_t>().swap(w.x);
    std::vector<uint64_t>().swap(w.y);
    std::vector<uint32_t>().swap(w.z);
  }
  if(w.pb.size() > KEEP_SORT_BUFFER) std::vector<uint64_t>().swap(w.pb);
  if(w.group.capacity() > KEEP_SORT_BUFFER) std::vector<std::pair<uint64_t, uint64_t>>().swap(w.group);
}

// One team member: score grids until none are left; me = 0 leads
void run_team(const Cells &c, Team &tm, size_t me, const int *lig, const int *rec, size_t G, std::atomic<size_t> &next,
              const std::atomic<bool> &failed, double *S, double *R){
  const size_t n = c.node.size(), nchunks = tm.chunks.size();
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
    const int tc = range > TABLE_BITS ? range - TABLE_BITS : 0;  // big grids: histogram cells of 2^tc keys from lo
    const size_t ncell = (size_t) ((hi - lo) >> tc) + 1;
    if(tm.big && sampled){
      sample_cells(c, tm, me, lo, tc, ncell);
      tm.bar.wait(failed);
    }
    if(me == 0) make_buckets(tm, c.cap, sampled, lo, range, tc, ncell);
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
      if(tm.multi){
        run_passes(c, tm, me, bk, w, failed, S + g * n, R + g * n);
        trim_work(w);
        continue;
      }
      if(tm.big){
        if(me == 0) six_storage(tm, bk);
        tm.bar.wait(failed);
      }
      for(size_t ch; (ch = tm.next_place.fetch_add(1)) < nchunks; ){
        if(!tm.big) place_chunk_ids(c, tm, bk, ch, w);
        else if(bk.map) place_chunk_six_nt<true>(c, tm, bk, ch, w);
        else place_chunk_six_nt<false>(c, tm, bk, ch, w);
      }
      tm.bar.wait(failed);
    }
    rank_buckets(tm, bk, w);
    tm.bar.wait(failed);
    const double denom = (double) w.start[NB] - 1.0;
    for(size_t ch; (ch = tm.next_sum.fetch_add(1)) < nchunks; ){
      if(tm.slabbed) sum_chunk<Store::SLABS>(c, tm, bk, ch, denom, w, S + g * n, R + g * n);
      else sum_chunk<Store::RUNS>(c, tm, bk, ch, denom, w, S + g * n, R + g * n);
    }
    if(tm.has_split){
      tm.bar.wait(failed);
      double *rcv = R + g * n;
      for(size_t k = n * me / tm.size, k1 = n * (me + 1) / tm.size; k < k1; ++k){
        if(!tm.split_cell[k]) continue;
        rcv[c.out[k]] = rank_mean(tm.isum[k].load(std::memory_order_relaxed), tm.icnt[k].load(std::memory_order_relaxed), denom);
        tm.isum[k].store(0, std::memory_order_relaxed);
        tm.icnt[k].store(0, std::memory_order_relaxed);
      }
    }
    trim_work(w);
  }
}


// Cells grouped into within-sample blocks as in grid_layout.h (node = their index in F, out = their position). Ranks and
// per-bucket counts are uint32, counted per bucket in grids of 2^32 - 1 entries or more (LOCAL_FROM).
void group_cells(Cells &c, const int *cell, const int *sample, size_t n){
  c.block_start.push_back(0);
  for(const std::vector<int> &members : sample_blocks(sample, n)){
    for(int k : members){
      c.node.push_back(cell[k]);
      c.out.push_back(k);
    }
    const size_t size = members.size();
    c.block_off.push_back(c.cap);
    c.block_start.push_back(c.node.size());
    c.cap += size * size;
    c.max_block = std::max(c.max_block, size);
  }
  if(c.cap >= GRID_MAX) Rcpp::stop("grid scoring: too many cells per sample (the sum of squared sample sizes must stay below 2^52)");
}

// Sender and receiver scores (cell x pair) of the grids of ligands lp and receptor subunits rp (negative: the network score); a big grid
// on several threads in `passes` ranges of buckets when passes > 1
Rcpp::List run_grids(const Cells &c, const int *lp, const int *rp, size_t G, int threads, size_t passes){
  const size_t n = c.node.size();
  Rcpp::NumericMatrix sender = alloc_numeric_matrix(n, G), receiver = alloc_numeric_matrix(n, G);
  double *S = REAL(sender), *R = REAL(receiver);

  const unsigned mxcsr = get_mxcsr();
  std::atomic<size_t> next(0);
  std::atomic<bool> failed(false);
  std::string error;

  const size_t cores = std::thread::hardware_concurrency();
  size_t nt = G == 0 || threads < 1 ? 1 : threads;
  if(cores > 0) nt = std::min(nt, cores);
  const size_t bytes = std::max<size_t>(1, grid_bytes(c));
  const size_t fit = std::max<size_t>(1, CACHE_BUDGET / bytes);
  const size_t nteams = std::max<size_t>(1, std::min(std::min(nt, fit), G));
  std::vector<std::unique_ptr<Team>> teams;
  try {
    for(size_t k = 0; k < nteams; ++k) teams.emplace_back(new Team(c, nt / nteams + (k < nt % nteams ? 1 : 0), nt == 1, passes));
  } catch(std::bad_alloc &){
    Rcpp::stop("grid scoring: not enough memory for the grid buffers (" + std::to_string(nteams * bytes >> 20) + " MB)");
  }
  error.reserve(256);

  auto worker = [&](Team *tm, size_t me){
    try {
      set_mxcsr(mxcsr);
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
    if(sample[k] < 1) Rcpp::stop("score_grids_cpp: sample ids must be positive integers");
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

  return run_grids(c, INTEGER(lig), INTEGER(rec), G, threads, 1);
}

// [[Rcpp::export]]
Rcpp::List score_blocks_cpp(Rcpp::NumericVector blocks, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample,
                            Rcpp::IntegerVector lig, Rcpp::IntegerVector rec, int threads, int passes = 1){
  const size_t n = sample.size(), G = lig.size();
  if((size_t) A.nrow() != n || (size_t) B.ncol() != n) Rcpp::stop("score_blocks_cpp: A needs one row and B one column per cell");
  if((size_t) rec.size() != G) Rcpp::stop("score_blocks_cpp: one receptor per ligand");
  if(passes < 1) Rcpp::stop("score_blocks_cpp: passes must be at least 1");
  for(size_t k = 0; k < n; ++k){
    if(sample[k] < 1) Rcpp::stop("score_blocks_cpp: sample ids must be positive integers");
  }
  for(size_t g = 0; g < G; ++g){
    if(lig[g] >= A.ncol() || rec[g] >= B.nrow() || (lig[g] < 0) != (rec[g] < 0)) Rcpp::stop("score_blocks_cpp: bad ligand or receptor index");
  }

  Cells c;
  c.A = REAL(A);
  c.B = REAL(B);
  c.nrec = B.nrow();
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

  return run_grids(c, INTEGER(lig), INTEGER(rec), G, threads, (size_t) passes);
}

// Memory mode: score_blocks_cpp ranking by the keys A + B' + Λ (f2.h). bytes, base, toff: Λ from influence_lambda_cpp for the same
// cells and samples at the same p; D: the adjacency's row sums of the cells, in the order of sample; A and B as score_blocks_cpp's,
// every entry positive and finite (their logs are taken)
// [[Rcpp::export]]
Rcpp::List score_lambda_cpp(Rcpp::RawVector bytes, Rcpp::NumericVector base, Rcpp::NumericVector toff, Rcpp::NumericVector D,
                            Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample, Rcpp::IntegerVector lig,
                            Rcpp::IntegerVector rec, int threads, int p, int passes){
  const size_t n = sample.size(), G = lig.size();
  if((size_t) A.nrow() != n || (size_t) B.ncol() != n || (size_t) D.size() != n){
    Rcpp::stop("score_lambda_cpp: one row of A, column of B and D per cell");
  }
  if((size_t) rec.size() != G) Rcpp::stop("score_lambda_cpp: one receptor per ligand");
  if(p < 0 || p > 52) Rcpp::stop("score_lambda_cpp: p must be in 0..52");
  if(passes < 1) Rcpp::stop("score_lambda_cpp: passes must be at least 1");
  for(size_t k = 0; k < n; ++k){
    if(sample[k] < 1) Rcpp::stop("score_lambda_cpp: sample ids must be positive integers");
    if(!(D[k] > 0) || !std::isfinite(D[k])) Rcpp::stop("score_lambda_cpp: D must be positive and finite");
  }
  for(size_t g = 0; g < G; ++g){
    if(lig[g] >= A.ncol() || rec[g] >= B.nrow() || (lig[g] < 0) != (rec[g] < 0)) Rcpp::stop("score_lambda_cpp: bad ligand or receptor index");
  }
  for(R_xlen_t k = 0; k < A.size(); ++k){
    if(!(A[k] > 0) || !std::isfinite(A[k])) Rcpp::stop("score_lambda_cpp: A must be positive and finite");
  }
  for(R_xlen_t k = 0; k < B.size(); ++k){
    if(!(B[k] > 0) || !std::isfinite(B[k])) Rcpp::stop("score_lambda_cpp: B must be positive and finite");
  }

  Cells c;
  c.A = REAL(A);
  c.B = REAL(B);
  c.nrec = B.nrow();
  std::vector<int> cell(n);
  for(size_t k = 0; k < n; ++k) cell[k] = (int) k;
  group_cells(c, cell.data(), INTEGER(sample), n);
  const size_t nblocks = c.block_start.size() - 1;
  if((size_t) base.size() != nblocks || (size_t) toff.size() != nblocks) Rcpp::stop("score_lambda_cpp: lambda does not match the samples");
  for(size_t blk = 0; blk < nblocks; ++blk){  // every block's tiles, and the slack get48 reads past the last, inside bytes
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], T = (s + f2::TILE - 1) / f2::TILE;
    if(!(toff[blk] >= 0) || toff[blk] + (double) (f2::TILE_BYTES * (T * (T + 1) / 2) + f2::SLACK) > (double) bytes.size()){
      Rcpp::stop("score_lambda_cpp: lambda does not match the samples");
    }
  }
  c.LamT = reinterpret_cast<const char *>(RAW(bytes));
  c.lbase.resize(nblocks);
  std::memcpy(c.lbase.data(), REAL(base), nblocks * sizeof(int64_t));
  c.toff.resize(nblocks);
  for(size_t blk = 0; blk < nblocks; ++blk) c.toff[blk] = (size_t) toff[blk];
  c.call = ++tile_calls;
  c.p = p;
  c.lD52.resize(n);
  c.kD.resize(n);
  for(size_t k = 0; k < n; ++k){
    c.lD52[k] = f2::fixlog52(D[c.out[k]]);
    c.kD[k] = f2::round_shift(c.lD52[k], 52 - p);
  }
  c.block_run.push_back(0);
  for(size_t blk = 0; blk < nblocks; ++blk){
    c.block_fb.push_back(c.block_off[blk]);
    c.block_gather.push_back(false);
    c.block_run.push_back(0);
  }

  return run_grids(c, INTEGER(lig), INTEGER(rec), G, threads, (size_t) passes);
}
