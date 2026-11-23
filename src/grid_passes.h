// grid_passes.h: passes > 1, a big grid's buckets placed, ranked and summed one range of buckets at a time, each entry's bucket
// recomputed at each pass, so that no bucket ids are stored. The members split each pass's chunks and buckets as for one pass,
// barriers between the stages. Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_PASSES_H
#define HYPERCOM_GRID_PASSES_H

#include "grid_place.h"
#include "grid_sums.h"

namespace {

// Multi: the in-range entries of rows [i0, i0 + m) of column jj of block blk, keys in [kmin, kend) (kend 0: no upper end) and not
// skipped, made as count_chunk makes them: their rows to w.ix and buckets to w.bc without a branch (each entry written, the count
// advanced by the in-range flag); returns how many. col: the column's keys (memory mode) or Q
template <bool TABLE> size_t column_buckets_in(const Cells &c, Team &tm, const Buckets &bk, size_t blk, size_t jj, size_t i0, size_t m,
                                               Work &w, uint64_t kmin, uint64_t kend, Column &col){
  col = column_values(c, tm, w, blk, jj, i0, m);
  const uint64_t *kk = col.kk;
  const double *q = col.q;
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
                          uint64_t kend, Column &col){
  if(bk.map) return column_buckets_in<true>(c, tm, bk, blk, jj, i0, m, w, kmin, kend, col);
  return column_buckets_in<false>(c, tm, bk, blk, jj, i0, m, w, kmin, kend, col);
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
    if(e > tm.cap6){
      throw std::runtime_error("a value bucket holds more entries than the key pool of one pass (very many equal values); use the "
                               "exact kernels (give more memory)");
    }
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
      Column col;
      const size_t n = column_buckets_any(c, tm, bk, blk, jj, i0, m, w, kmin, kend, col);
      for(size_t t = 0; t < n; ++t){
        const size_t d = bc[t], i = ix[t];
        const uint64_t key = col.kk ? col.kk[i] : sortable_key(col.q[i]);
        const uint64_t low = bk.map ? bk.low_bits<true>(key, d) : bk.low_bits<false>(key, d);
        if(wide[d]) ext[at[d]++] = low;
        else put48(pool + 6 * at[d]++, low);
      }
    }
  }
}

// Multi: the integer rank sums of chunk ch's entries in buckets [d0, d1), their buckets recomputed; rows to msum (a chunk's rows are
// its own), columns to the team's column sums (isum / icnt), all divided after the last pass
void sum_chunk_range(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w, size_t d0, size_t d1, uint64_t kmin, uint64_t kend){
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  const char **rp = w.rp.data();
  for(size_t d = d0; d < d1; ++d) rp[d] = bucket_ranks(tm, w, d) + sizeof(uint32_t) * cnt[d];
  const Chunk &ck = tm.chunks[ch];
  const size_t *base = w.start.data();
  const uint16_t *bc = w.bc.data();
  const uint32_t *ix = w.ix.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t first = c.block_start[blk], s = c.block_start[blk + 1] - first;
    const size_t i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    uint64_t *sc = tm.msum.data() + 2 * (first + i0);  // row sum and count of row i at 2 i, 2 i + 1
    for(size_t jj = 0; jj < s; ++jj){
      Column col;
      const size_t n = column_buckets_any(c, tm, bk, blk, jj, i0, m, w, kmin, kend, col);
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

// Member me's share of the cells: the receiver scores whose column sums were added up across chunks (isum / icnt: every cell in
// multi, else the cells of split blocks), one division each, and those sums reset for the next grid
void finish_columns(const Cells &c, Team &tm, size_t me, double denom, double *rcv){
  const size_t n = c.out.size();
  for(size_t k = n * me / tm.size, k1 = n * (me + 1) / tm.size; k < k1; ++k){
    if(!tm.multi && !tm.split_cell[k]) continue;
    rcv[c.out[k]] = rank_mean(tm.isum[k].load(std::memory_order_relaxed), tm.icnt[k].load(std::memory_order_relaxed), denom);
    tm.isum[k].store(0, std::memory_order_relaxed);
    tm.icnt[k].store(0, std::memory_order_relaxed);
  }
}

// Multi, one team member after the counts and their offsets: each range of buckets placed, ranked and summed in turn, then the
// member's share of the scores, one division each
void run_passes(const Cells &c, Team &tm, size_t me, const Buckets &bk, Work &w, const std::atomic<bool> &failed, double *snd, double *rcv){
  const size_t n = c.out.size(), nchunks = tm.chunks.size(), NB = bk.nb + 2;
  if(me == 0) cut_ranges(tm, bk);
  tm.bar.sync(failed);
  for(size_t r = 0; r + 1 < tm.ranges.size(); ++r){
    const size_t d0 = tm.ranges[r], d1 = tm.ranges[r + 1];
    if(me == 0){
      six_storage(tm, bk, d0, d1);
      tm.next_place = 0;
      tm.next_bucket = 0;
      tm.next_sum = 0;
    }
    tm.bar.sync(failed);
    for(size_t ch; (ch = tm.next_place.fetch_add(1)) < nchunks; ) place_chunk_range(c, tm, bk, ch, w, d0, d1, tm.kmin[r], tm.kend[r]);
    tm.bar.sync(failed);
    rank_buckets(tm, bk, w, d0, d1);
    tm.bar.sync(failed);
    for(size_t ch; (ch = tm.next_sum.fetch_add(1)) < nchunks; ) sum_chunk_range(c, tm, bk, ch, w, d0, d1, tm.kmin[r], tm.kend[r]);
    tm.bar.sync(failed);
  }
  const double denom = (double) w.start[NB] - 1.0;
  for(size_t k = n * me / tm.size, k1 = n * (me + 1) / tm.size; k < k1; ++k){
    snd[c.out[k]] = rank_mean(tm.msum[2 * k], tm.msum[2 * k + 1], denom);
    tm.msum[2 * k] = tm.msum[2 * k + 1] = 0;
  }
  finish_columns(c, tm, me, denom, rcv);
}

}  // namespace

#endif
