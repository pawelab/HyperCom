// grid_team.h: the grid kernel's limits and tuning knobs, the cells and value buckets of the grids, the buffers a team of threads
// shares (Team: the team scores one grid at a time, its members split each stage into chunks or buckets taken one at a time and
// meet at barriers; the leader, member 0, does what is done once per grid) and each thread's own (Work). Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_TEAM_H
#define HYPERCOM_GRID_TEAM_H

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>
#include "f2.h"
#include "grid_common.h"
#include "threads.h"

namespace {

// Hard limits
const int B1_LIMIT = 15;                // equal width: at most 2^15 buckets
const int NARROW_BITS = 48;             // key bits of a 6-byte pool slot: a big grid's bucket keeps at most this many
static_assert(((size_t) 1 << B1_LIMIT) + 1 < SKIP, "equal-width bucket ids must stay below SKIP");
static_assert(IDS_MAX + 1 < SKIP, "equal-count bucket ids must stay below SKIP");
static_assert(NARROW_BITS == 8 * 6, "a 6-byte pool slot holds 48 key bits");
// Ranks are uint32. A grid of LOCAL_FROM entries or more (run stores only) keeps each rank counted from its bucket's first and adds
// the bucket's start (64-bit) when reading it, so a bucket holds at most BUCKET_MAX entries. A chunk holds fewer than CHUNK_MAX
// entries, so its per-bucket counts (uint32) cannot wrap, and a slab-stored grid (one chunk) keeps global ranks. The three macros are
// test-only knobs: test builds lower them (-DHYPERCOM_LOCAL_FROM=0 -DHYPERCOM_BUCKET_MAX=... -DHYPERCOM_CHUNK_MAX=...) so that small
// grids take the paths of huge ones; a release build leaves them at these values.
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

// Tuning knobs: speed and memory only, never a result. HYPERCOM_CACHE_BUDGET is a test-only knob like the three above: test builds
// may lower it (-DHYPERCOM_CACHE_BUDGET=bytes) so that small grids take the paths of big ones.
#ifndef HYPERCOM_CACHE_BUDGET
#define HYPERCOM_CACHE_BUDGET ((size_t) 24 << 20)
#endif
const size_t CACHE_BUDGET = HYPERCOM_CACHE_BUDGET;  // team buffers that may be in flight at once (grid_bytes); a bigger grid gets all
                                                    // threads. R reads it and grid_bytes' 10 B per entry from kernel_constants_cpp
const int B1_MIN = 11, B1_MAX = 14;     // equal width: b1 grows with the grid between these
const int FILL_BITS = 8;                // equal width: about 2^8 entries per bucket (b1 = bit width of the entries - 8)
const int RANK_SAMPLE_BITS = 22;        // bucket bounds from every SAMPLE_STEP-th entry of a block, or every (K >> 22)-th (odd) when
                                        // that is sparser (about 2^22 entries): big grids sample faster
const size_t SAMPLE_FIRST = 3;          // starting at its entry 3
const int TABLE_BITS = 18;              // big grids: the sampled key range is cut into 2^18 cells whose counts set the buckets
const size_t BUCKET_TARGET = 65536;     // big grids: entries per bucket aimed at, so a bucket's sort buffers stay in L2
const size_t BUCKETS_MIN = 2048, BUCKETS_MAX = 8192;  // and between these many buckets: above ~5.4e8 entries fewer, bigger buckets
                                                      // place faster
const size_t PASS_SLACK_DIV = 8;        // passes > 1: a range of buckets is cut at a bucket's edge, so its key pool holds cap / passes
const size_t PASS_SLACK_MAX = (size_t) 1 << 26;  // entries plus an eighth of that, at most 2^26 (far more than one equal-count bucket)
const size_t SLAB = 64;                 // entries per pool slab on one thread: each bucket is a linked list of slabs
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

enum class Store { SLABS, RUNS };       // where the ranks are: slab lists (one thread), or one run of the pool per bucket

// The grids' cells and values. Cells are taken in position order: block by block (samples ascending), the cells of a block in the
// order given; out maps a position to the cell's row of the results
struct Cells {
  const double *Fb = nullptr;          // exact kernel: F[i, j] of every within-sample block, block blk's at block_off[blk], column-major
  const double *A = nullptr;           // F[cell, ligand], cells x ligands, rows in cell order
  const double *B = nullptr;           // and F[receptor, cell], receptors x cells (nrec rows)
  size_t nrec = 0;
  std::vector<int> out;                // the cell at each position
  std::vector<size_t> block_start;     // first position of each block, and the end
  std::vector<size_t> block_off;       // first entry of each block (sum of the squared sizes of the blocks before it)
  size_t cap = 0;                      // within-sample entries
  size_t max_block = 0;                // cells of the biggest block
  const char *LamT = nullptr;          // memory mode: Λ - lbase in tiles (f2.h), block blk's at toff[blk]; null in the exact kernel
  std::vector<size_t> toff;
  std::vector<int64_t> lbase;          // per block: its Λ base
  int p = f2::KEY_BITS;                // fractional bits of the keys
  std::vector<int64_t> lD52;           // fixlog52(D) by position
  std::vector<int64_t> kD;             // ln D by position at 2^p (network grid)
};

// Team buffers of one grid: a key (8 B) and a bucket id (2 B) per within-sample entry
const size_t GRID_ENTRY_BYTES = sizeof(uint64_t) + sizeof(uint16_t);
inline size_t grid_bytes(const Cells &c){
  return c.cap * GRID_ENTRY_BYTES;
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

// Q of rows [i0, i1) of column jj of block blk (into q), or F itself for the network grid
const double *block_column(const Cells &c, size_t blk, size_t jj, size_t i0, size_t i1, const double *a, double bj, bool net, double *q){
  const size_t s = c.block_start[blk + 1] - c.block_start[blk];
  const double *f = c.Fb + c.block_off[blk] + jj * s + i0;
  if(net) return f;
  make_q(a, bj, f, i1 - i0, q);
  return q;
}

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
  Team(const Cells &c, size_t size, bool alone, size_t passes) : size(size), bar(size), step(std::max(SAMPLE_STEP, (c.cap >> RANK_SAMPLE_BITS) | 1)),
                                                                 a(c.out.size()), b(c.out.size()), ka(c.LamT ? c.out.size() : 0),
                                                                 kb(c.LamT ? c.out.size() : 0), lo(size), hi(size) {
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
      slab_pos.resize(nbmax);
      full.resize(nbmax);
      bucket.resize(c.cap);
    } else {
      cnt.resize(nchunks * nbmax);
      tot.resize(nbmax);
      if(multi){
        cap6 = c.cap / passes + std::min(c.cap / (PASS_SLACK_DIV * passes), PASS_SLACK_MAX) + 1;
        pool6.reset(new char[cap6 * 6 + f2::SLACK]);
        hist.resize(size);
        msum.assign(2 * c.out.size(), 0);
      } else if(big){
        pool6.reset(new char[c.cap * 6 + f2::SLACK]);
        hist.resize(size);
      } else {
        pool.reset(new uint64_t[std::max<size_t>(1, c.cap)]);
      }
      if(!multi) bucket.resize(c.cap);
    }
    isum.reset(new std::atomic<uint64_t>[c.out.size()]());
    icnt.reset(new std::atomic<uint64_t>[c.out.size()]());
    split_cell.assign(c.out.size(), 0);
    for(const Chunk &ck : chunks){
      if(ck.head && ck.last) continue;
      has_split = true;
      for(size_t k = c.block_start[ck.b0]; k < c.block_start[ck.b1]; ++k) split_cell[k] = 1;
    }
  }

  const size_t size;
  hc::Barrier bar;
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
  std::vector<size_t> link, slab_pos, full;  // slabs: next slab of each slab, write position and full slabs of each bucket
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

// Memory mode: a thread's strip of the TILE columns of tile column J of a block for rows [i0, i1), column-major, made from whole
// tiles (tile (I, J) at and above the diagonal, tile (J, I) transposed below it): each tile is read once for its TILE columns. Λ does
// not change during a call, so a strip serves every grid of the call.
struct TileStrip {
  size_t blk = SIZE_MAX, J = SIZE_MAX, i0 = 0, i1 = 0;
  std::vector<int64_t> v;
};

// One thread's buffers, made in run_team for one call
struct Work {
  TileStrip strip;                     // memory mode
  std::vector<uint64_t> x, y;
  std::vector<uint32_t> z, cnt;
  std::vector<uint32_t> bnd;           // a big bucket's top-digit offsets
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
  std::vector<char *> next_byte;       // and the next byte of each bucket
  std::vector<uint16_t> bc;            // a column's buckets (multi: of its in-range entries, in ix order)
  std::vector<uint32_t> ix;            // multi: the rows of a column's in-range entries

  Work(const Cells &c, const Team &tm) : start(tm.nbmax + 1), at(tm.nbmax), rp(tm.nbmax), q(c.max_block), kq(c.LamT ? c.max_block : 0),
                                         rk(c.max_block), irc(2 * c.max_block), bc(tm.multi ? c.max_block : 0),
                                         ix(tm.multi ? c.max_block : 0) {}
};

}  // namespace

#endif
