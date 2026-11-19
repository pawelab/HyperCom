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
// can change the last bits of the e x e inverse, as it does for lu and chol.
//
// chol and woodbury stop as solve() does when the system is computationally singular (LAPACK dpocon, reciprocal condition number
// below the machine epsilon), e.g. for r near 0.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <thread>
#include <vector>
#include "grid_layout.h"
#include "r_alloc.h"

namespace {

const size_t TILE = 64;                 // chol: tile edge of the loops over both triangles
const size_t BLOCK_COLS = 32;           // fill_block: columns per task
const size_t Z_COLS = 64;               // influence_factors_cpp: columns of Z per task
const double ROWSUM_TOL = 1e-12;        // woodbury: relative tolerance of adj's row sums against the hyperedges

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

// fn(worker, task) for tasks 0..ntasks-1 on up to nthreads workers (the calling thread is worker 0), each task once; fn must not
// throw. If a thread cannot start, the workers already running take its tasks.
template <class Fn>
void run_tasks(int nthreads, size_t ntasks, Fn fn){
  std::atomic<size_t> next(0);
  auto work = [&](int worker){
    for(size_t t = next.fetch_add(1); t < ntasks; t = next.fetch_add(1)) fn(worker, t);
  };
  std::vector<std::thread> pool;
  const int extra = (int) std::min<size_t>(nthreads > 1 ? nthreads - 1 : 0, ntasks > 1 ? ntasks - 1 : 0);
  try{
    for(int k = 1; k <= extra; ++k) pool.emplace_back(work, k);
  } catch(...){
  }
  work(0);
  for(auto &t : pool) t.join();
}

}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_influence_cpp(SEXP adj, Rcpp::NumericVector rs, double r){
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("solve_influence_cpp: one row sum per row");
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
  if(rcond < DBL_EPSILON) Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);

  for(size_t k = 0; k < N * N; ++k) x[k] = r * x[k];

  return X;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_influence_chol_cpp(SEXP adj, Rcpp::NumericVector rs, double r){
  const bool sym = Rf_inherits(adj, "dsCMatrix");
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("solve_influence_chol_cpp: one row sum per row");
  const size_t N = n;

  Rcpp::NumericMatrix S = dense_adjacency(adj, n, false);
  double *m = REAL(S);
  if(!sym){
    for(size_t jb = 0; jb < N; jb += TILE){
      for(size_t ib = 0; ib <= jb; ib += TILE){
        for(size_t j = jb; j < std::min(jb + TILE, N); ++j){
          for(size_t i = ib; i < std::min(ib + TILE, j); ++i){
            if(!(m[i + j * N] == m[j + i * N])) Rcpp::stop("method \"chol\" needs a symmetric adjacency matrix (use method=\"lu\")");
          }
        }
      }
    }
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
  if(rcond < DBL_EPSILON) Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);
  F77_CALL(dpotri)(&uplo, &n, m, &n, &info FCONE);
  if(info != 0) Rcpp::stop("Lapack routine %s failed with info %d", "dpotri", info);

  for(size_t jb = 0; jb < N; jb += TILE){
    for(size_t ib = 0; ib <= jb; ib += TILE){
      for(size_t j = jb; j < std::min(jb + TILE, N); ++j){
        for(size_t i = ib; i < std::min(ib + TILE, j + 1); ++i){
          const double u = m[i + j * N];
          m[i + j * N] = (u * row[i]) * col[j];
          if(i != j) m[j + i * N] = (u * row[j]) * col[i];
        }
      }
    }
  }

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
  if(rs.size() != n) Rcpp::stop("influence_factors_cpp: one row sum per node");
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
  if(rcond < DBL_EPSILON) Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);
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

  return Rcpp::List::create(Rcpp::Named("Z") = Z, Rcpp::Named("delta") = Rcpp::NumericVector(delta.begin(), delta.end()));
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
  check_nodes(cell, h.N);
  if(sample.size() != cell.size()) Rcpp::stop("influence_blocks_cpp: one sample per cell");
  for(R_xlen_t k = 0; k < sample.size(); ++k){
    if(sample[k] < 1) Rcpp::stop("influence_blocks_cpp: sample ids must be positive integers");
  }
  const std::vector<double> ny = scaled_incidence(h, REAL(delta));

  // the within-sample blocks in the kernel's layout (grid_layout.h), each column-major, by node
  std::vector<std::vector<int>> blocks = sample_blocks(INTEGER(sample), sample.size());
  size_t total = 0;
  for(std::vector<int> &m : blocks){
    for(int &k : m) k = cell[k];
    total += m.size() * m.size();
  }
  Rcpp::NumericVector out = alloc_vector<REALSXP>((R_xlen_t) total, false);
  size_t off = 0;
  for(const auto &m : blocks){
    fill_block(h, ny, REAL(Z), REAL(rs), REAL(delta), r, m.data(), m.size(), m.data(), m.size(), REAL(out) + off, m.size(), threads);
    off += m.size() * m.size();
  }

  return out;
}
