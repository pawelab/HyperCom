// grid_input.h: the grid entries a thread reads, a column at a time: Q = (a b) / F from the blocks of F, or in memory mode integer
// keys from Λ's tiles (each thread gathers tiles into its own strip); a and b of the team's cells, made by each member for its
// share. Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_INPUT_H
#define HYPERCOM_GRID_INPUT_H

#include "grid_team.h"

namespace {

// Memory mode: Λ - lbase of rows [i0, i0 + m) of column jj of block blk, in the thread's strip (Work::strip; 0 on the diagonal: the
// tiles hold 0)
const int64_t *tile_strip_column(const Cells &c, TileStrip &st, size_t blk, size_t jj, size_t i0, size_t m){
  const size_t TL = f2::TILE, s = c.block_start[blk + 1] - c.block_start[blk], i1 = i0 + m, J = jj / TL;
  if(st.blk != blk || st.J != J || st.i0 != i0 || st.i1 != i1){
    st.blk = blk;
    st.J = J;
    st.i0 = i0;
    st.i1 = i1;
    if(st.v.size() < TL * m) st.v.resize(TL * m);
    f2::tile_strip(c.LamT + c.toff[blk], J, std::min(TL, s - TL * J), i0, i1, st.v.data());
  }
  return st.v.data() + (jj % TL) * m;
}

// Memory mode: the keys of rows [i0, i0 + m) of column jj of block blk into k, by integer adds: A[i] + B'[j] + Λ[i, j] (the network
// grid: ln D[j] - Λ[i, j], the order of F), signed, then flipped to unsigned order
const uint64_t *lambda_keys(const Cells &c, const Team &tm, TileStrip &st, size_t blk, size_t jj, size_t i0, size_t m, uint64_t *k){
  const size_t first = c.block_start[blk];
  const int64_t *lam = tile_strip_column(c, st, blk, jj, i0, m), base = c.lbase[blk];
  if(tm.l < 0){
    const int64_t dj = c.kD[first + jj] - base;
    for(size_t i = 0; i < m; ++i) k[i] = (uint64_t) (dj - lam[i]) ^ SIGN;
    return k;
  }
  const int64_t *ka = tm.ka.data() + first + i0, bj = tm.kb[first + jj] + base;
  for(size_t i = 0; i < m; ++i) k[i] = (uint64_t) (ka[i] + bj + lam[i]) ^ SIGN;
  return k;
}

// Memory mode: the key of entry (ii, jj) of block blk
inline uint64_t lambda_key(const Cells &c, const Team &tm, size_t blk, size_t ii, size_t jj){
  const size_t first = c.block_start[blk];
  const int64_t lam = c.lbase[blk] + f2::tile_value(c.LamT + c.toff[blk], ii, jj);
  const int64_t v = tm.l < 0 ? c.kD[first + jj] - lam : tm.ka[first + ii] + tm.kb[first + jj] + lam;
  return (uint64_t) v ^ SIGN;
}

// Rows [i0, i0 + m) of column jj of a block: in memory mode their keys (kk), else their Q (q); the other is null
struct Column {
  const uint64_t *kk;
  const double *q;
};

inline Column column_values(const Cells &c, const Team &tm, Work &w, size_t blk, size_t jj, size_t i0, size_t m){
  if(c.LamT) return {lambda_keys(c, tm, w.strip, blk, jj, i0, m, w.kq.data()), nullptr};
  const size_t first = c.block_start[blk];
  const bool net = tm.l < 0;
  return {nullptr, block_column(c, blk, jj, i0, i0 + m, tm.a.data() + first + i0, net ? 0 : tm.b[first + jj], net, w.q.data())};
}

// Member me of the team: a and b of its share of the cells (memory mode: their keys ka, kb)
void prepare_ab(const Cells &c, Team &tm, size_t me){
  if(tm.l < 0) return;
  const size_t n = c.out.size(), k1 = n * (me + 1) / tm.size;
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
}

// F[i, j] of entry (ii, jj) of block blk
inline double block_value(const Cells &c, size_t blk, size_t ii, size_t jj){
  const size_t s = c.block_start[blk + 1] - c.block_start[blk];
  return c.Fb[c.block_off[blk] + jj * s + ii];
}

}  // namespace

#endif
