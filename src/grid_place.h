// grid_place.h: the entries of a grid put into their buckets and ranked. One thread appends them to slab lists; a team counts its
// chunks' entries per bucket, turns the counts into offsets and places each chunk's entries (chunks, buckets and ranking units each
// taken by one member at a time, barriers between the stages). Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_PLACE_H
#define HYPERCOM_GRID_PLACE_H

#include "grid_sort.h"
#include "grid_input.h"
#ifdef HYPERCOM_X86
#include <emmintrin.h>
#endif

namespace {

// One thread: Q of each entry in position order, appended to its bucket's list of slabs. Kept out of line like pull_slabs and
// pull_runs: inlined into run_team, their loops lose registers to it.
__attribute__((noinline)) void fill_slabs(const Cells &c, Team &tm, const Buckets &bk, Work &w){
  const size_t NB = bk.nb + 2;
  const Chunk &ck = tm.chunks[0];
  uint64_t *pool = tm.pool.get();
  size_t *cur = tm.slab_pos.data(), *full = tm.full.data(), *link = tm.link.data();
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
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      const Column col = column_values(c, tm, w, blk, jj, i0, m);
      if(col.kk){
        for(size_t i = 0; i < m; ++i){
          if(i0 + i == jj) bc[i] = SKIP;
          else append(bc[i], col.kk[i]);
        }
        continue;
      }
      for(size_t i = 0; i < m; ++i){
        const double v = col.q[i];
        if(skipped(i0 + i, jj, v)) bc[i] = SKIP;
        else append(bc[i], sortable_key(v));
      }
    }
  }
}

// First pass over chunk ch: the bucket of each entry (SKIP if skipped; multi: to scratch only) and the chunk's entries per bucket
template <bool TABLE> void count_chunk(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const Chunk &ck = tm.chunks[ch];
  uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  std::fill(cnt, cnt + bk.nb + 2, 0u);
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      uint16_t *bc = tm.multi ? w.bc.data() : tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      const size_t diag = jj - i0;       // wraps above m when the diagonal is not in the chunk
      const Column col = column_values(c, tm, w, blk, jj, i0, m);
      if(col.kk){
        for(size_t i = 0; i < m; ++i){
          if(i == diag){
            bc[i] = SKIP;
            continue;
          }
          const size_t d = bk.locate<TABLE>(col.kk[i]);
          bc[i] = (uint16_t) d;
          ++cnt[d];
        }
        continue;
      }
      for(size_t i = 0; i < m; ++i){
        uint64_t u;
        std::memcpy(&u, col.q + i, sizeof u);
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
  const Chunk &ck = tm.chunks[ch];
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  size_t *at = w.at.data();
  for(size_t d = 0; d < bk.nb + 2; ++d) at[d] = w.start[d] + cnt[d];
  uint64_t *pool = tm.pool.get();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const Column col = column_values(c, tm, w, blk, jj, i0, m);
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        pool[at[d]++] = bk.low_bits<false>(col.kk ? col.kk[i] : sortable_key(col.q[i]), d);
      }
    }
  }
}

// Second pass over chunk ch of a big grid: the kept key bits of each entry to its place in its bucket, 6 bytes in pool6 or 8 in
// ext (sbase / swide)
template <bool TABLE> void place_chunk_six(const Cells &c, Team &tm, const Buckets &bk, size_t ch, Work &w){
  const Chunk &ck = tm.chunks[ch];
  const uint32_t *cnt = tm.cnt.data() + ch * tm.stride;
  const size_t NB = bk.nb + 2;
  size_t *at = w.at.data();
  for(size_t d = 0; d < NB; ++d) at[d] = tm.sbase[d] + cnt[d];
  const char *wide = tm.swide.data();
  char *pool = tm.pool6.get();
  uint64_t *ext = tm.ext.data();
  for(size_t blk = ck.b0; blk < ck.b1; ++blk){
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const Column col = column_values(c, tm, w, blk, jj, i0, m);
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t low = bk.low_bits<TABLE>(col.kk ? col.kk[i] : sortable_key(col.q[i]), d);
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
  std::vector<char *> &cur = w.next_byte;  // next byte to write of each narrow bucket
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
    const size_t s = c.block_start[blk + 1] - c.block_start[blk], i0 = ck.i0, m = std::min(ck.i1, s) - i0;
    for(size_t jj = 0; jj < s; ++jj){
      const Column col = column_values(c, tm, w, blk, jj, i0, m);
      const uint16_t *bc = tm.bucket.data() + c.block_off[blk] + jj * s + i0;
      for(size_t i = 0; i < m; ++i){
        const size_t d = bc[i];
        if(d == SKIP) continue;
        const uint64_t low = bk.low_bits<TABLE>(col.kk ? col.kk[i] : sortable_key(col.q[i]), d);
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

// The leader, big grid: where the keys of buckets [d0, d1) go (all buckets, or a pass's range), buckets of at most NARROW_BITS key
// bits one after the other in pool6, the others (the catch-alls at least) in ext, sized now
void six_storage(Team &tm, const Buckets &bk, size_t d0, size_t d1){
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

// Rank buckets [d0, d1) (all buckets, or a pass's range), a run of about E / (RANK_UNITS x team size) of their E entries at a time;
// start[d] = number of entries in the buckets below d
void rank_buckets(Team &tm, const Buckets &bk, Work &w, size_t d0, size_t d1){
  const size_t *start = w.start.data(), *full = tm.full.data(), *link = tm.link.data(), E = start[d1] - start[d0];
  const size_t unit = std::max<size_t>(1, E / (RANK_UNITS * tm.size)), nunits = (E + unit - 1) / unit;
  for(size_t u; (u = tm.next_bucket.fetch_add(1)) < nunits; ){
    const size_t e1 = std::lower_bound(start + d0, start + d1, start[d0] + (u + 1) * unit) - start;
    for(size_t d = std::lower_bound(start + d0, start + d1, start[d0] + u * unit) - start; d < e1; ++d){
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
        if(tm.slab_pos[d] % SLAB) w.segs.emplace_back(slab * SLAB, tm.slab_pos[d] % SLAB);
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

}  // namespace

#endif
