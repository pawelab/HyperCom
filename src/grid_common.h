// grid_common.h: what the grid kernels (score_grids.cpp, screen_grids.cpp) and the blocks made for them (influence.cpp) share: the
// within-sample blocks of the cell x cell grids, the grid entries Q and their order-preserving keys, the SSE control word the worker
// threads take over, and the checks of the arguments
#ifndef HYPERCOM_GRID_COMMON_H
#define HYPERCOM_GRID_COMMON_H

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "HyperCom must not be compiled with -ffast-math or -ffinite-math-only: its NaN and Inf checks would be removed"
#endif

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "f2.h"
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#define HYPERCOM_X86 1
#endif

// A row's rank sum adds at most (cells of its sample) ranks below K (all within-sample entries): it fits uint64 while their product
// stays below this, 2^63 less a margin
const double RANK_SUM_MAX = 9.2e18;
const uint16_t SKIP = 0xFFFF;           // bucket id of a skipped entry (the diagonal, 0, NaN)
const size_t IDS_MAX = 65533;           // at most this many buckets, so the ids 0..nb + 1 (with the two catch-alls) stay below SKIP
const size_t SAMPLE_STEP = 61;          // bucket bounds from every 61st entry of a block (or a sparser step for big grids)
const uint64_t SIGN = (uint64_t) 1 << 63;  // memory mode: flips a key's int64 to an order-preserving uint64
// The bits of +Inf shifted left by one: the shifted bits of NaN lie above, those of 0 wrap to the top
const uint64_t INF_SHIFTED = (uint64_t) 0x7FF0000000000000ULL << 1;

// The SSE control word (rounding, denormals) of the calling thread, given to the workers so that Q's divisions round as R's would
inline unsigned get_mxcsr(){
#ifdef HYPERCOM_X86
  return _mm_getcsr();
#else
  return 0;
#endif
}

inline void set_mxcsr(unsigned m){
#ifdef HYPERCOM_X86
  _mm_setcsr(m);
#else
  (void) m;
#endif
}

typedef double v2d __attribute__((vector_size(16)));  // two doubles, one SSE2 instruction

// An order-preserving key of a double: unsigned order of the keys = numeric order of the values
inline uint64_t sortable_key(double v){
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | ((uint64_t) 1 << 63));
}

// sortable_key on the bits u of a value, the sign flip as a mask: ~u if negative, else u | sign bit
inline uint64_t flip_key(uint64_t u){
  return u ^ ((uint64_t) ((int64_t) u >> 63) | ((uint64_t) 1 << 63));
}

// The value with bits u is 0 or NaN (one compare: see INF_SHIFTED)
inline bool zero_or_nan(uint64_t u){
  return (u << 1) - 1 >= INF_SHIFTED;
}

inline int bit_width(uint64_t x){
  return x ? 64 - __builtin_clzll(x) : 0;
}

// An entry the grids leave out: the diagonal, 0 and NaN
inline bool skipped(size_t ii, size_t jj, double v){
  return ii == jj || v == 0 || std::isnan(v);
}

// Q[k] = (a[k] * bj) / f[k] for k < m, two at a time: the same doubles as R's (a * b) / F
inline void make_q(const double *a, double bj, const double *f, size_t m, double *q){
  const v2d b2 = {bj, bj};
  size_t k = 0;
  for(; k + 2 <= m; k += 2){
    v2d av, fv;
    std::memcpy(&av, a + k, sizeof av);
    std::memcpy(&fv, f + k, sizeof fv);
    const v2d qv = (av * b2) / fv;
    std::memcpy(q + k, &qv, sizeof qv);
  }
  for(; k < m; ++k) q[k] = (a[k] * bj) / f[k];
}

// The cells of each within-sample block, by their position 0..n-1: samples ascending (ids >= 1, checked by the caller), empty
// samples skipped, the cells of a sample in the given order. The blocks F[cells, cells] follow one another column-major in
// this order.
inline std::vector<std::vector<int>> sample_blocks(const int *sample, size_t n){
  int n_samples = 0;
  for(size_t k = 0; k < n; ++k) n_samples = std::max(n_samples, sample[k]);
  std::vector<std::vector<int>> members(n_samples);
  for(size_t k = 0; k < n; ++k) members[sample[k] - 1].push_back((int) k);
  members.erase(std::remove_if(members.begin(), members.end(), [](const std::vector<int> &m){ return m.empty(); }), members.end());

  return members;
}

// Stops (what: the error's prefix) unless every sample id is a positive integer
inline void check_samples(const Rcpp::IntegerVector &sample, const char *what){
  for(R_xlen_t k = 0; k < sample.size(); ++k){
    if(sample[k] < 1) Rcpp::stop(std::string(what) + ": sample ids must be positive integers");
  }
}

// Memory mode: stops (what: the error's prefix) unless Λ (bytes, base, toff from influence_lambda_cpp, made at p fractional bits) has
// every tile of blocks of the given sizes inside bytes, with the slack get48 reads past the last slot, and D (one per cell), A and B
// are positive and finite (their logs are taken)
inline void check_lambda(const Rcpp::RawVector &bytes, const Rcpp::NumericVector &base, const Rcpp::NumericVector &toff,
                         const Rcpp::NumericVector &D, const Rcpp::NumericMatrix &A, const Rcpp::NumericMatrix &B,
                         const std::vector<size_t> &sizes, size_t n, int p, const char *what){
  const std::string w(what);
  if(p < 0 || p > 52) Rcpp::stop(w + ": p must be in 0..52");
  if((size_t) D.size() != n) Rcpp::stop(w + ": one D per cell");
  for(size_t k = 0; k < n; ++k){
    if(!(D[k] > 0) || !std::isfinite(D[k])) Rcpp::stop(w + ": D must be positive and finite");
  }
  for(R_xlen_t k = 0; k < A.size(); ++k){
    if(!(A[k] > 0) || !std::isfinite(A[k])) Rcpp::stop(w + ": A must be positive and finite in memory mode");
  }
  for(R_xlen_t k = 0; k < B.size(); ++k){
    if(!(B[k] > 0) || !std::isfinite(B[k])) Rcpp::stop(w + ": B must be positive and finite in memory mode");
  }
  const size_t nblocks = sizes.size();
  if((size_t) base.size() != nblocks || (size_t) toff.size() != nblocks) Rcpp::stop(w + ": lambda does not match the samples");
  for(size_t blk = 0; blk < nblocks; ++blk){
    const double T = std::ceil((double) sizes[blk] / f2::TILE);
    if(!(toff[blk] >= 0) || toff[blk] != std::floor(toff[blk]) ||
       toff[blk] + (double) f2::TILE_BYTES * (T * (T + 1) / 2) + f2::SLACK > (double) bytes.size()){
      Rcpp::stop(w + ": lambda's tiles do not fit its bytes");
    }
  }
}

#endif
