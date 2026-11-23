// grid_buckets.h: the value buckets of a grid from a sample of its entries. Each member samples its share; the leader (member 0)
// sets the buckets for the whole team between two barriers. Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_BUCKETS_H
#define HYPERCOM_GRID_BUCKETS_H

#include "grid_input.h"

namespace {

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

// The sampled key range of the grid: lo, hi and its bit width; false (lo = 0, width 64) if nothing was sampled
bool key_range(const Team &tm, uint64_t &lo, uint64_t &hi, int &width){
  lo = ~(uint64_t) 0;
  hi = 0;
  for(size_t me = 0; me < tm.size; ++me){
    lo = std::min(lo, tm.lo[me]);
    hi = std::max(hi, tm.hi[me]);
  }
  if(lo > hi){
    lo = 0;
    width = 64;
    return false;
  }
  width = bit_width(hi - lo);
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

// The leader: the grid's buckets for everyone in tm.bk, equal count for a big grid (at most IDS_MAX; histogram cells of
// 2^tc keys from lo), else equal width (2^b1 buckets; for a big grid at least 2^(width - NARROW_BITS), at most 2^B1_LIMIT, so the
// keys fit 6 bytes unless the width is 64), and the per-bucket counts sized
void make_buckets(Team &tm, size_t cap, bool sampled, uint64_t lo, int width, int tc, size_t ncell){
  const size_t nchunks = tm.chunks.size();
  bool table = false;
  if(tm.big && sampled) table = count_buckets(tm, cap, lo, tc, ncell) && tm.bk.nb <= IDS_MAX;
  if(!table){
    const int b1 = tm.big ? std::min(width, std::min(B1_LIMIT, std::max(tm.bucket_bits, width - NARROW_BITS))) : std::min(width, tm.bucket_bits);
    Buckets &bk = tm.bk;
    bk.lo = lo;
    bk.t = width - b1;
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

}  // namespace

#endif
