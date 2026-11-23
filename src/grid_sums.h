// grid_sums.h: the scores from the ranks: each chunk's exact integer rank sums (one member per chunk), divided once per score.
// Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_SUMS_H
#define HYPERCOM_GRID_SUMS_H

#include "grid_sort.h"

namespace {

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

// The run stores: where the ranks of bucket d start (they overwrite the first 4 bytes per entry of the bucket's run of keys): in
// pool6 or ext for a big grid, else in pool
inline const char *bucket_ranks(const Team &tm, const Work &w, size_t d){
  if(!tm.big) return (const char *) (tm.pool.get() + w.start[d]);
  return tm.swide[d] ? (const char *) tm.ext.data() + sizeof(uint64_t) * tm.sbase[d] : tm.pool6.get() + 6 * tm.sbase[d];
}

// Sender and receiver sums of chunk ch as exact integers, each column's ranks taken from the chunk's entries of each bucket in
// position order: rows of the chunk complete here, columns too when the chunk holds whole blocks, else added to the team's column
// sums for run_team to finish
template <Store STORE> void sum_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, double denom, Work &w, double *snd,
                                      double *rcv){
  const uint64_t *pool = tm.pool.get();
  const size_t *link = tm.link.data();
  size_t *at = w.at.data();
  const char **rp = w.rp.data();
  if(STORE == Store::SLABS){
    for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = d * SLAB;
  } else {
    const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
    for(size_t d = 0; d < bk.nb + 2; ++d) rp[d] = bucket_ranks(tm, w, d) + sizeof(uint32_t) * cnt[d];
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
        rcv[out[jj]] = rank_mean(cs, cn, denom);
      } else {
        tm.isum[first + jj].fetch_add(cs, std::memory_order_relaxed);
        tm.icnt[first + jj].fetch_add(cn, std::memory_order_relaxed);
      }
    }
    for(size_t k = 0; k < m; ++k) snd[out[i0 + k]] = rank_mean(sc[2 * k], sc[2 * k + 1], denom);
  }
}

}  // namespace

#endif
