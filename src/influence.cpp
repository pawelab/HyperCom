// influence.cpp: the influence matrix F = r (I - W)^-1 with W = (1 - r) A / rowSums(A), three ways
//
// lu (solve_influence_cpp): bit-identical to generate_infMat of HyperCom 0.0.0.9000. That version forms W from the sparse
// adjacency, densifies it, forms M = I - W and calls solve(M), which solves M X = I with LAPACK's dgesv on a copy of M and a
// dense identity: four or five n x n matrices at once. Here M is written straight from the adjacency (a base matrix, or the
// slots of a dgCMatrix or dsCMatrix) and dgesv factors it in place into an identity that becomes the result: two n x n
// matrices. Same arithmetic element by element, the same dgesv call and the same singularity checks as R's solve.
//
// chol (solve_influence_chol_cpp): the earlier Cholesky variant (method="chol" of HyperCom 0.1.0) with its arithmetic, step by
// step. With D = diag(rowSums(A)) and s = D^-1/2, M = I - W is similar to the symmetric positive definite
// S = I - (1 - r) diag(s) A diag(s) (eigenvalues in [r, 2 - r]), so F = r D^-1/2 S^-1 D^1/2. S is written into one n x n
// matrix, inverted in place by LAPACK dpotrf + dpotri (upper triangle), mirrored and scaled: one n x n matrix, about 3/8 of
// the flops of lu. Not bit-identical to lu (~1e-15 relative).
//
// woodbury (influence_factors_cpp, then influence_block_cpp / influence_blocks_cpp): the adjacency of a hypergraph is A = H H'
// with its diagonal set to 0, H the node x hyperedge incidence counts. So K = D - (1 - r) A = Delta - (1 - r) H H' with the diagonal
// Delta = D + (1 - r) diag(H H'), and by the Woodbury identity K^-1 = Delta^-1 + (1 - r) Y C^-1 Y' with Y = Delta^-1 H and the
// hyperedge x hyperedge matrix C = I - (1 - r) H' Y, symmetric positive definite. Since M = D^-1 K, F = r K^-1 D =
// r Delta^-1 D + r (1 - r) Y Z with Z = C^-1 Y' D (e x n). Every entry of Y, C^-1 and Z is >= 0, so each entry of F is a sum of
// nonnegative terms (no cancellation). F is never stored: the factors (H, Delta, Z) take O(e n) memory, and any block F[rows, cols]
// is made from them on demand at |rows| x |cols| x (hyperedges per row node), each entry the same at any `threads` and in any block.
// Cost O(e^3 + n nnz(H)) for the factors instead of O(n^3). Not bit-identical to lu (~1e-15 relative); the BLAS's own thread count
// can change the last bits of the e x e inverse, as it does for lu and chol. In memory mode (influence_lambda_cpp) the blocks are
// made a panel at a time and kept only as Λ = -ln F tiles of 6-byte integer keys (f2.h) for the grid kernel.
//
// The grid kernels read the within-sample blocks of F (grid_common.h): influence_blocks_cpp makes them from the factors,
// dense_blocks_cpp copies them from a dense F.
//
// chol and woodbury stop as solve() does when the system is computationally singular (LAPACK dpocon, reciprocal condition number
// below the machine epsilon, or NaN), e.g. for r near 0.
//
// Threads: fill_block, the hyperedge matrix and Z split their columns among worker threads (threads.h), each entry made by one
// thread in a fixed order, so the results do not depend on the thread count.
//
// No contraction of a * b + c into one fused multiply-add: an FMA rounds once instead of twice, so F would differ in the last bits
// between CPUs that have one (arm64) and those whose default build does not use it (x86-64). Clang honours the standard pragma;
// GCC ignores it in C++ and contracts on arm64 (its default -ffp-contract=fast), so there F can differ in the last bits.
#ifdef __clang__
#pragma STDC FP_CONTRACT OFF
#endif

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "f2.h"
#include "grid_common.h"
#include "r_alloc.h"
#include "threads.h"

namespace {

const size_t SOLVE_TILE = 64;           // chol: tile edge of the loops over both triangles, so the mirrored entries stay in cache
const size_t BLOCK_COLS = 32;           // fill_block: columns per task
const size_t Z_COLS = 64;               // influence_factors_cpp: columns of Z per task
const double ROWSUM_TOL = 1e-12;        // woodbury: relative tolerance of adj's row sums against the hyperedges
const size_t PANEL = 256;               // influence_lambda_cpp: columns per panel, a multiple of f2::TILE (a panel's rows end on a
                                        // tile's edge)
const int64_t MARGIN = (int64_t) 1 << 20;  // influence_lambda_cpp: Λ's base lies this far below the smallest libm estimate, far more
                                        // than libm's logs and fixlog52 can differ (a few units at 2^-44)
static_assert(PANEL % f2::TILE == 0, "a panel's rows end on a tile's edge");

// The calling thread's error for a reciprocal condition number below the machine epsilon, or NaN (a matrix with NaN), as solve()'s
void check_rcond(double rcond){
  if(!(rcond >= DBL_EPSILON)) Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);
}

// fn(i, j) for the entries i < j (i <= j with diagonal) of the upper triangle of an N x N matrix, tile by tile, so that the mirrored
// entries (j, i) of a tile stay in cache
template <class Fn> void upper_tiles(size_t N, bool diagonal, Fn fn){
  for(size_t jb = 0; jb < N; jb += SOLVE_TILE){
    for(size_t ib = 0; ib <= jb; ib += SOLVE_TILE){
      for(size_t j = jb; j < std::min(jb + SOLVE_TILE, N); ++j){
        for(size_t i = ib; i < std::min(ib + SOLVE_TILE, j + diagonal); ++i) fn(i, j);
      }
    }
  }
}

// order of a square adjacency given as a base matrix, dgCMatrix or dsCMatrix, with R's solve() messages
int adjacency_order(SEXP adj){
  const bool dense = Rf_isMatrix(adj);
  if(!dense && !Rf_inherits(adj, "dgCMatrix") && !Rf_inherits(adj, "dsCMatrix")){
    Rcpp::stop("adj must be a base matrix, a dgCMatrix or a dsCMatrix");
  }
  int n = 0;
  if(dense){
    n = Rf_nrows(adj);
    if(Rf_ncols(adj) != n) Rcpp::stop("'a' (%d x %d) must be square", n, Rf_ncols(adj));
  } else {
    const Rcpp::IntegerVector dim = Rf_getAttrib(adj, Rf_install("Dim"));  // the Matrix slot, not R's dim (R_DimSymbol)
    n = dim[0];
    if(dim[1] != n) Rcpp::stop("'a' (%d x %d) must be square", n, dim[1]);
  }
  if(n == 0) Rcpp::stop("'a' is 0-diml");

  return n;
}

// The n x n adjacency as a dense matrix: a base matrix copied, a dgCMatrix written from its slots, a dsCMatrix (either
// triangle stored) into the upper triangle, and into the lower one too when mirror
Rcpp::NumericMatrix dense_adjacency(SEXP adj, int n, bool mirror){
  const size_t N = n;
  const bool dense = Rf_isMatrix(adj);
  Rcpp::NumericMatrix M = alloc_numeric_matrix(n, n, !dense);
  double *m = REAL(M);
  if(dense){
    const Rcpp::NumericVector a(adj);
    std::copy(a.begin(), a.end(), m);
  } else {
    const Rcpp::IntegerVector i = Rf_getAttrib(adj, Rf_install("i")), p = Rf_getAttrib(adj, Rf_install("p"));
    const Rcpp::NumericVector x = Rf_getAttrib(adj, Rf_install("x"));
    const bool sym = Rf_inherits(adj, "dsCMatrix");
    for(size_t j = 0; j < N; ++j){
      for(int k = p[j]; k < p[j + 1]; ++k){
        const size_t row = i[k];
        if(!sym){
          m[row + j * N] = x[k];
        } else {
          const size_t lo = std::min(row, j), hi = std::max(row, j);
          m[lo + hi * N] = x[k];
          if(mirror) m[hi + lo * N] = x[k];
        }
      }
    }
  }

  return M;
}

// fn(worker, task) for tasks 0..ntasks-1, each once, on up to nthreads worker threads (worker 0..nthreads-1) while the calling thread
// polls for an interrupt (threads.h); with one worker the tasks run on the calling thread, which polls between them. An error or an
// interrupt stops the workers at their next task and is raised here.
template <class Fn>
void run_tasks(int nthreads, size_t ntasks, Fn fn){
  hc::Stop stop;
  const size_t nt = std::min<size_t>(std::max(1, nthreads), ntasks);
  if(nt <= 1){
    hc::Poll poll;
    for(size_t t = 0; t < ntasks; ++t){
      if(poll()){
        stop.interrupt();
        break;
      }
      fn(0, t);
    }
  } else {
    std::atomic<size_t> next(0);
    hc::run_threads(nt, stop, [&](size_t worker){
      for(size_t t; !stop.failed.load(std::memory_order_relaxed) && (t = next.fetch_add(1)) < ntasks; ) fn((int) worker, t);
    });
  }
  stop.raise("influence matrix");
}

}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_influence_cpp(SEXP adj, Rcpp::NumericVector rs, double r){
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("influence matrix: one row sum per row");
  const size_t N = n;
  const double *rsp = REAL(rs);

  Rcpp::NumericMatrix M = dense_adjacency(adj, n, true);
  double *m = REAL(M);
  const double c = 1 - r;
  for(size_t j = 0; j < N; ++j){
    for(size_t i = 0; i < N; ++i){
      double w = m[i + j * N] / rsp[i];
      w = c * w;
      m[i + j * N] = (i == j ? 1.0 : 0.0) - w;
    }
  }

  const char one = '1';
  const double anorm = F77_CALL(dlange)(&one, &n, &n, m, &n, nullptr FCONE);
  Rcpp::NumericMatrix X = alloc_numeric_matrix(n, n);
  double *x = REAL(X);
  for(size_t k = 0; k < N; ++k) x[k * N + k] = 1.0;
  std::vector<int> ipiv(N), iwork(N);
  int info = 0;
  F77_CALL(dgesv)(&n, &n, m, &n, ipiv.data(), x, &n, &info);
  if(info < 0) Rcpp::stop("argument %d of Lapack routine %s had invalid value", -info, "dgesv");
  if(info > 0) Rcpp::stop("Lapack routine %s: system is exactly singular: U[%d,%d] = 0", "dgesv", info, info);
  double rcond = 0;
  std::vector<double> work(4 * N);
  F77_CALL(dgecon)(&one, &n, m, &n, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
  check_rcond(rcond);

  for(size_t k = 0; k < N * N; ++k) x[k] = r * x[k];

  return X;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_influence_chol_cpp(SEXP adj, Rcpp::NumericVector rs, double r){
  const bool sym = Rf_inherits(adj, "dsCMatrix");
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("influence matrix: one row sum per row");
  const size_t N = n;

  Rcpp::NumericMatrix S = dense_adjacency(adj, n, false);
  double *m = REAL(S);
  if(!sym){
    upper_tiles(N, false, [&](size_t i, size_t j){
      if(!(m[i + j * N] == m[j + i * N])) Rcpp::stop("method \"chol\" needs a symmetric adjacency matrix (use method=\"lu\")");
    });
  }

  std::vector<double> s(N), row(N), col(N);
  for(size_t i = 0; i < N; ++i){
    s[i] = 1 / std::sqrt(rs[i]);
    row[i] = r * s[i];
    col[i] = 1 / s[i];
  }
  const double alpha = -(1 - r);
  for(size_t j = 0; j < N; ++j){
    for(size_t i = 0; i <= j; ++i) m[i + j * N] = ((m[i + j * N] * s[i]) * s[j]) * alpha;
    m[j + j * N] = m[j + j * N] + 1.0;
  }

  const char uplo = 'U', one = '1';
  int info = 0;
  std::vector<double> work(3 * N);
  std::vector<int> iwork(N);
  const double anorm = F77_CALL(dlansy)(&one, &uplo, &n, m, &n, work.data() FCONE FCONE);
  F77_CALL(dpotrf)(&uplo, &n, m, &n, &info FCONE);
  if(info < 0) Rcpp::stop("argument %d of Lapack routine %s had invalid value", -info, "dpotrf");
  if(info > 0) Rcpp::stop("the leading minor of order %d is not positive: the system is not positive definite (use method=\"lu\")", info);
  double rcond = 0;
  F77_CALL(dpocon)(&uplo, &n, m, &n, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
  check_rcond(rcond);
  F77_CALL(dpotri)(&uplo, &n, m, &n, &info FCONE);
  if(info != 0) Rcpp::stop("Lapack routine %s failed with info %d", "dpotri", info);

  upper_tiles(N, true, [&](size_t i, size_t j){
    const double u = m[i + j * N];
    m[i + j * N] = (u * row[i]) * col[j];
    if(i != j) m[j + i * N] = (u * row[j]) * col[i];
  });

  return S;
}

namespace {

// The incidence by node (hyperedges ascending, with their counts) from the slots of a node x hyperedge dgCMatrix of positive
// counts, and the size (summed counts) of each hyperedge
struct Incidence {
  size_t N, E;
  std::vector<size_t> nptr;
  std::vector<int> nedge;
  std::vector<double> nh, edge_size;
};

Incidence node_incidence(const Rcpp::IntegerVector &hp, const Rcpp::IntegerVector &hi, const Rcpp::NumericVector &hx, int n){
  if(n < 1) Rcpp::stop("the incidence needs at least one node");
  if(hp.size() < 2 || hp[0] != 0 || hp[hp.size() - 1] != hi.size() || hi.size() != hx.size()){
    Rcpp::stop("H must be the slots of a node x hyperedge dgCMatrix");
  }
  Incidence h;
  h.N = n;
  h.E = hp.size() - 1;
  const size_t NNZ = hi.size();
  for(size_t e = 0; e < h.E; ++e){
    if(hp[e + 1] < hp[e]) Rcpp::stop("H must be the slots of a node x hyperedge dgCMatrix");
  }
  for(size_t k = 0; k < NNZ; ++k){
    if(hi[k] < 0 || hi[k] >= n || !(hx[k] > 0) || !std::isfinite(hx[k])) Rcpp::stop("H must hold positive incidence counts of the nodes");
  }

  h.nptr.assign(h.N + 1, 0);
  for(size_t k = 0; k < NNZ; ++k) ++h.nptr[hi[k] + 1];
  for(size_t i = 0; i < h.N; ++i) h.nptr[i + 1] += h.nptr[i];
  h.nedge.resize(NNZ);
  h.nh.resize(NNZ);
  h.edge_size.assign(h.E, 0);
  std::vector<size_t> fill(h.nptr.begin(), h.nptr.end() - 1);
  for(size_t e = 0; e < h.E; ++e){
    for(int k = hp[e]; k < hp[e + 1]; ++k){
      h.nedge[fill[hi[k]]] = (int) e;
      h.nh[fill[hi[k]]++] = hx[k];
      h.edge_size[e] += hx[k];
    }
  }

  return h;
}

// Y = H / Delta by node, in the order of the node incidence
std::vector<double> scaled_incidence(const Incidence &h, const double *delta){
  std::vector<double> ny(h.nedge.size());
  for(size_t i = 0; i < h.N; ++i){
    for(size_t k = h.nptr[i]; k < h.nptr[i + 1]; ++k) ny[k] = h.nh[k] / delta[i];
  }

  return ny;
}

// The factors' sizes and the node indices of rows or columns
void check_factors(const Incidence &h, const Rcpp::NumericVector &rs, const Rcpp::NumericVector &delta, const Rcpp::NumericMatrix &Z){
  if((size_t) rs.size() != h.N || (size_t) delta.size() != h.N) Rcpp::stop("the factors need one row sum and one Delta per node");
  if((size_t) Z.nrow() != h.E || (size_t) Z.ncol() != h.N) Rcpp::stop("Z must have one row per hyperedge and one column per node");
}

void check_nodes(const Rcpp::IntegerVector &idx, size_t N){
  for(R_xlen_t k = 0; k < idx.size(); ++k){
    if(idx[k] < 0 || (size_t) idx[k] >= N) Rcpp::stop("node index out of range");
  }
}

// The within-sample blocks of the grid kernels' layout (grid_common.h), each the nodes of its cells: cell holds each cell's 0-based
// node, sample its sample (from 1); what: the error's prefix
std::vector<std::vector<int>> cell_blocks(const Rcpp::IntegerVector &cell, const Rcpp::IntegerVector &sample, size_t N, const char *what){
  check_nodes(cell, N);
  if(sample.size() != cell.size()) Rcpp::stop(std::string(what) + ": one sample per cell");
  check_samples(sample, what);
  std::vector<std::vector<int>> blocks = sample_blocks(INTEGER(sample), sample.size());
  for(std::vector<int> &m : blocks){
    for(int &k : m) k = cell[k];
  }

  return blocks;
}

// F[rows, cols] into out (column-major, leading dimension ld) from the factors: r (1 - r) times the sum over the row node's
// hyperedges, ascending, of Y * Z, plus r rs / Delta on the diagonal. Columns go in blocks of BLOCK_COLS, with Z of the block
// laid out hyperedge-major; each entry is made by one thread in a fixed order, so F is the same at any thread count.
void fill_block(const Incidence &h, const std::vector<double> &ny, const double *Z, const double *rs, const double *delta, double r,
                const int *rows, size_t nr, const int *cols, size_t nc, double *out, size_t ld, int threads){
  const size_t E = h.E, nblocks = (nc + BLOCK_COLS - 1) / BLOCK_COLS;
  if(nr == 0 || nblocks == 0) return;
  const int nt = std::max(1, std::min<int>(threads, (int) nblocks));
  std::vector<std::vector<double>> zt(nt, std::vector<double>(E * BLOCK_COLS));
  const double c = 1 - r, rc = r * c;
  run_tasks(nt, nblocks, [&](int w, size_t b){
    double *z = zt[w].data();
    const size_t c0 = b * BLOCK_COLS, nb = std::min(BLOCK_COLS, nc - c0);
    std::fill(z, z + E * BLOCK_COLS, 0.0);
    for(size_t jj = 0; jj < nb; ++jj){
      const double *zc = Z + (size_t) cols[c0 + jj] * E;
      for(size_t e = 0; e < E; ++e) z[e * BLOCK_COLS + jj] = zc[e];
    }
    for(size_t ii = 0; ii < nr; ++ii){
      const size_t i = rows[ii];
      double acc[BLOCK_COLS] = {0};
      for(size_t k = h.nptr[i]; k < h.nptr[i + 1]; ++k){
        const double y = ny[k];
        const double *ze = z + (size_t) h.nedge[k] * BLOCK_COLS;
        for(size_t jj = 0; jj < BLOCK_COLS; ++jj) acc[jj] += y * ze[jj];
      }
      for(size_t jj = 0; jj < nb; ++jj){
        const size_t j = cols[c0 + jj];
        double f = rc * acc[jj];
        if(i == j) f += r * rs[j] / delta[j];
        out[ii + (c0 + jj) * ld] = f;
      }
    }
  });
}

}

// [[Rcpp::export]]
Rcpp::List influence_factors_cpp(Rcpp::IntegerVector hp, Rcpp::IntegerVector hi, Rcpp::NumericVector hx, int n, Rcpp::NumericVector rs,
                                 double r, int threads){
  const Incidence h = node_incidence(hp, hi, hx, n);
  if(rs.size() != n) Rcpp::stop("influence matrix: one row sum per node");
  const size_t N = h.N, E = h.E;
  const double c = 1 - r;
  const double *rsp = REAL(rs), *hxp = REAL(hx);
  const int *hpp = INTEGER(hp), *hip = INTEGER(hi);

  // Delta = D + (1 - r) diag(H H') with D checked against the hyperedges; Y = H / Delta
  std::vector<double> delta(N);
  for(size_t i = 0; i < N; ++i){
    double d = 0, sq = 0;
    for(size_t k = h.nptr[i]; k < h.nptr[i + 1]; ++k){
      d += h.nh[k] * (h.edge_size[h.nedge[k]] - h.nh[k]);
      sq += h.nh[k] * h.nh[k];
    }
    if(!(rsp[i] > 0) || !std::isfinite(rsp[i])) Rcpp::stop("method \"woodbury\" needs a positive row sum for every node of adj");
    if(std::fabs(d - rsp[i]) > ROWSUM_TOL * rsp[i]){
      Rcpp::stop("hyperedges do not match adj: node %d has row sum %g in adj and %g from the hyperedges", (int) i + 1, rsp[i], d);
    }
    delta[i] = rsp[i] + c * sq;
  }
  const std::vector<double> ny = scaled_incidence(h, delta.data());

  // C = I - (1 - r) H' Y, lower triangle by column, then C^-1 in place (dpotrf + dpotri) and mirrored
  std::vector<double> C;
  try{
    C.assign(E * E, 0.0);
  } catch(std::bad_alloc &){
    Rcpp::stop("method \"woodbury\": not enough memory for the %d x %d hyperedge matrix", (int) E, (int) E);
  }
  double *cm = C.data();
  run_tasks(threads, E, [&](int, size_t f){
    double *col = cm + f * E;
    for(int k = hpp[f]; k < hpp[f + 1]; ++k){
      const size_t i = hip[k];
      const double yf = hxp[k] / delta[i];
      const int *first = h.nedge.data() + h.nptr[i], *last = h.nedge.data() + h.nptr[i + 1];
      for(const int *q = std::lower_bound(first, last, (int) f); q < last; ++q) col[*q] += h.nh[q - h.nedge.data()] * yf;
    }
    for(size_t e = f; e < E; ++e) col[e] = (e == f ? 1.0 : 0.0) - c * col[e];
  });
  const char uplo = 'L', one = '1';
  int ne = (int) E, info = 0;
  std::vector<double> work(3 * E);
  std::vector<int> iwork(E);
  const double anorm = F77_CALL(dlansy)(&one, &uplo, &ne, cm, &ne, work.data() FCONE FCONE);
  F77_CALL(dpotrf)(&uplo, &ne, cm, &ne, &info FCONE);
  if(info != 0) Rcpp::stop("method \"woodbury\": the hyperedge matrix is not positive definite (Lapack dpotrf info %d; use method=\"lu\")", info);
  double rcond = 0;
  F77_CALL(dpocon)(&uplo, &ne, cm, &ne, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
  check_rcond(rcond);
  F77_CALL(dpotri)(&uplo, &ne, cm, &ne, &info FCONE);
  if(info != 0) Rcpp::stop("Lapack routine %s failed with info %d", "dpotri", info);
  for(size_t f = 0; f < E; ++f){
    for(size_t e = f + 1; e < E; ++e) cm[f + e * E] = cm[e + f * E];
  }

  // Z = C^-1 Y' D, column j summed over the node's hyperedges in ascending order
  Rcpp::NumericMatrix Z = alloc_numeric_matrix(E, n);
  double *zm = REAL(Z);
  run_tasks(threads, (N + Z_COLS - 1) / Z_COLS, [&](int, size_t b){
    for(size_t j = b * Z_COLS; j < std::min(N, b * Z_COLS + Z_COLS); ++j){
      double *zc = zm + j * E;
      for(size_t k = h.nptr[j]; k < h.nptr[j + 1]; ++k){
        const double wj = ny[k] * rsp[j];
        const double *cf = cm + (size_t) h.nedge[k] * E;
        for(size_t e = 0; e < E; ++e) zc[e] += cf[e] * wj;
      }
    }
  });

  Rcpp::NumericVector delta_out = alloc_vector<REALSXP>(N, false);
  std::copy(delta.begin(), delta.end(), delta_out.begin());

  return Rcpp::List::create(Rcpp::Named("Z") = Z, Rcpp::Named("delta") = delta_out);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix influence_block_cpp(Rcpp::IntegerVector hp, Rcpp::IntegerVector hi, Rcpp::NumericVector hx, int n, Rcpp::NumericVector rs,
                                        Rcpp::NumericVector delta, Rcpp::NumericMatrix Z, double r, Rcpp::IntegerVector rows,
                                        Rcpp::IntegerVector cols, int threads){
  const Incidence h = node_incidence(hp, hi, hx, n);
  check_factors(h, rs, delta, Z);
  check_nodes(rows, h.N);
  check_nodes(cols, h.N);
  const std::vector<double> ny = scaled_incidence(h, REAL(delta));
  const size_t nr = rows.size(), nc = cols.size();

  Rcpp::NumericMatrix F = alloc_numeric_matrix(nr, nc, false);
  fill_block(h, ny, REAL(Z), REAL(rs), REAL(delta), r, INTEGER(rows), nr, INTEGER(cols), nc, REAL(F), nr, threads);

  return F;
}

// [[Rcpp::export]]
Rcpp::NumericVector influence_blocks_cpp(Rcpp::IntegerVector hp, Rcpp::IntegerVector hi, Rcpp::NumericVector hx, int n, Rcpp::NumericVector rs,
                                         Rcpp::NumericVector delta, Rcpp::NumericMatrix Z, double r, Rcpp::IntegerVector cell,
                                         Rcpp::IntegerVector sample, int threads){
  const Incidence h = node_incidence(hp, hi, hx, n);
  check_factors(h, rs, delta, Z);
  const std::vector<std::vector<int>> blocks = cell_blocks(cell, sample, h.N, "influence blocks");
  const std::vector<double> ny = scaled_incidence(h, REAL(delta));

  // each block column-major, one after the other
  size_t total = 0;
  for(const std::vector<int> &m : blocks) total += m.size() * m.size();
  Rcpp::NumericVector out = alloc_vector<REALSXP>((R_xlen_t) total, false);
  size_t off = 0;
  for(const std::vector<int> &m : blocks){
    fill_block(h, ny, REAL(Z), REAL(rs), REAL(delta), r, m.data(), m.size(), m.data(), m.size(), REAL(out) + off, m.size(), threads);
    off += m.size() * m.size();
  }

  return out;
}

// influence_blocks_cpp for a dense influence matrix F: the blocks copied from it, without a temporary
// [[Rcpp::export]]
Rcpp::NumericVector dense_blocks_cpp(Rcpp::NumericMatrix F, Rcpp::IntegerVector cell, Rcpp::IntegerVector sample){
  const size_t N = F.nrow();
  if((size_t) F.ncol() != N) Rcpp::stop("dense blocks: the influence matrix must be square");
  const std::vector<std::vector<int>> blocks = cell_blocks(cell, sample, N, "dense blocks");

  size_t total = 0;
  for(const std::vector<int> &m : blocks) total += m.size() * m.size();
  Rcpp::NumericVector out = alloc_vector<REALSXP>((R_xlen_t) total, false);
  const double *f = REAL(F);
  double *o = REAL(out);
  for(const std::vector<int> &m : blocks){
    for(int j : m){
      const double *col = f + (size_t) j * N;
      for(int i : m) *o++ = col[i];
    }
  }

  return out;
}

// Memory mode: Λ of the within-sample blocks straight from the factors, no F block held. Λ[i, j] = -round((ln(F[i, j] / D[j]) +
// ln(F[j, i] / D[i])) / 2 at 2^p), the logs by f2::fixlog52 (integers only, so the same bits on every computer), stored as 6-byte
// offsets from a per-block base in the tiles of f2.h: the upper triangle (F is symmetric up to rounding, so Λ is made symmetric by
// the average), diagonal tiles whole (their lower half too, read as stored), the diagonal and padding 0. toff: each block's byte
// offset; base: each block's int64 base, its bits in a double. A panel of PANEL columns at a time: F[cells, panel] and
// F[panel, cells] by fill_block. Each block is made twice: the offsets need the block's base before the first is written, and
// keeping F until then would take the 8 bytes per entry memory mode avoids. The first pass estimates each block's smallest Λ with
// libm logs and takes MARGIN below it as the base; the second makes the exact keys and writes them. Stops if a block's span does not
// fit 48 bits.
// [[Rcpp::export]]
Rcpp::List influence_lambda_cpp(Rcpp::IntegerVector hp, Rcpp::IntegerVector hi, Rcpp::NumericVector hx, int n, Rcpp::NumericVector rs,
                                Rcpp::NumericVector delta, Rcpp::NumericMatrix Z, double r, Rcpp::IntegerVector cell,
                                Rcpp::IntegerVector sample, int p, int threads){
  const Incidence h = node_incidence(hp, hi, hx, n);
  check_factors(h, rs, delta, Z);
  const std::vector<std::vector<int>> blocks = cell_blocks(cell, sample, h.N, "memory mode");
  if(p < 0 || p > 52) Rcpp::stop("memory mode: p must be in 0..52");
  const std::vector<double> ny = scaled_incidence(h, REAL(delta));
  const size_t nblocks = blocks.size(), TL = f2::TILE;
  size_t total = 0;                     // bytes
  Rcpp::NumericVector toff = alloc_vector<REALSXP>(nblocks), base = alloc_vector<REALSXP>(nblocks);
  for(size_t b = 0; b < nblocks; ++b){
    const size_t T = (blocks[b].size() + TL - 1) / TL;
    toff[b] = (double) total;
    total += f2::TILE_BYTES * (T * (T + 1) / 2);
  }
  Rcpp::RawVector bytes = alloc_vector<RAWSXP>((R_xlen_t) (total + f2::SLACK), true);  // zeroed: diagonal, padding, slack
  char *out = reinterpret_cast<char *>(RAW(bytes));
  const int nt = std::max(1, threads), sh = 52 - p;
  const double sc = std::ldexp(1.0, p);
  std::atomic<bool> wide(false), bad(false);
  for(size_t b = 0; b < nblocks; ++b){
    const std::vector<int> &m = blocks[b];
    const size_t s = m.size();
    if(s < 2) continue;
    std::vector<double> Dn(s);
    for(size_t k = 0; k < s; ++k) Dn[k] = rs[m[k]];
    char *o = out + (size_t) toff[b];
    std::vector<double> A, C;
    int64_t lo = INT64_MAX, hi = INT64_MIN;
    for(int pass = 0; pass < 2; ++pass){
      const int64_t bs = lo;
      for(size_t j0 = 0; j0 < s; j0 += PANEL){
        const size_t w = std::min(PANEL, s - j0), nrow = j0 + w;
        A.resize(nrow * w);
        C.resize(w * nrow);
        fill_block(h, ny, REAL(Z), REAL(rs), REAL(delta), r, m.data(), nrow, m.data() + j0, w, A.data(), nrow, nt);
        fill_block(h, ny, REAL(Z), REAL(rs), REAL(delta), r, m.data() + j0, w, m.data(), nrow, C.data(), w, nt);
        std::vector<int64_t> tlo(w, INT64_MAX), thi(w, INT64_MIN);
        run_tasks(nt, w, [&](int, size_t jl){
          const size_t j = j0 + jl, i1 = std::min(s, (j | (TL - 1)) + 1);  // rows to the end of j's diagonal tile
          const double *a = A.data() + jl * nrow;
          for(size_t i = 0; i < i1; ++i){
            if(i == j) continue;
            const double fij = a[i] / Dn[j], fji = C[jl + i * w] / Dn[i];
            if(pass == 0){
              if(!(fij > 0) || !(fji > 0) || !std::isfinite(fij) || !std::isfinite(fji)){
                bad.store(true);
                continue;
              }
              const int64_t v = (int64_t) std::llround(-0.5 * (std::log(fij) + std::log(fji)) * sc);
              tlo[jl] = std::min(tlo[jl], v);
              thi[jl] = std::max(thi[jl], v);
              continue;
            }
            const int64_t lam = -f2::round_shift(f2::fixlog52(fij) + f2::fixlog52(fji), sh + 1), d = lam - bs;
            if(d < 0 || (d >> 48)) wide.store(true);
            f2::put48(o + f2::tile_offset(i, j), (uint64_t) d);
          }
        });
        if(pass == 0){
          for(size_t jl = 0; jl < w; ++jl){
            lo = std::min(lo, tlo[jl]);
            hi = std::max(hi, thi[jl]);
          }
        }
      }
      if(pass == 0){
        if(bad.load()){
          Rcpp::stop("memory mode: cells of one sample are not connected in the hypergraph (influence 0); use the exact kernels (give "
                     "more memory)");
        }
        if(lo > hi) lo = hi = 0;
        lo -= MARGIN;
        if((uint64_t) (hi + MARGIN - lo) >> 48){
          Rcpp::stop("memory mode: the influence values of a sample span too wide a range for its 48-bit keys; use the exact kernels "
                     "(give more memory)");
        }
      }
    }
    std::memcpy(&base[b], &lo, sizeof lo);
  }
  if(wide.load()) Rcpp::stop("memory mode: an entry of lambda left its block's 48-bit range");

  return Rcpp::List::create(Rcpp::Named("bytes") = bytes, Rcpp::Named("base") = base, Rcpp::Named("toff") = toff);
}
