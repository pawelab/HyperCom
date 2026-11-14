// influence.cpp: the influence matrix F = r (I - W)^-1 with W = (1 - r) A / rowSums(A), three ways
//
// lu (solve_influence_cpp): bit-identical to master's generate_infMat. Master forms W from the sparse adjacency, densifies it,
// forms M = I - W and calls solve(M), which solves M X = I with LAPACK's dgesv on a copy of M and a dense identity: four or five
// n x n matrices at once. Here M is written straight from the adjacency (a base matrix, or the slots of a dgCMatrix or dsCMatrix)
// and dgesv factors it in place into an identity that becomes the result: two n x n matrices. Same arithmetic element by element,
// the same dgesv call and the same singularity checks as R's solve.
//
// chol (solve_influence_chol_cpp): v1-fast's method with its arithmetic, step by step. With D = diag(rowSums(A)) and s = D^-1/2,
// M = I - W is similar to the symmetric positive definite S = I - (1 - r) diag(s) A diag(s) (eigenvalues in [r, 2 - r]), so
// F = r D^-1/2 S^-1 D^1/2. S is written into one n x n matrix, inverted in place by LAPACK dpotrf + dpotri (upper triangle),
// mirrored and scaled: one n x n matrix, about 3/8 of the flops of lu. Not bit-identical to lu (~1e-15 relative).
//
// hypergraph (solve_influence_hypergraph_cpp): the adjacency of a hypergraph is A = H H' with its diagonal set to 0, H the node x
// hyperedge incidence counts. So K = D - (1 - r) A = Delta - (1 - r) H H' with the diagonal Delta = D + (1 - r) diag(H H'), and by
// the Woodbury identity K^-1 = Delta^-1 + (1 - r) Y C^-1 Y' with Y = Delta^-1 H and the hyperedge x hyperedge matrix
// C = I - (1 - r) H' Y, symmetric positive definite. Since M = D^-1 K, F = r K^-1 D = r Delta^-1 D + r (1 - r) Y Z with
// Z = C^-1 Y' D. Every entry of Y, C^-1 and Z is >= 0, so each entry of F is a sum of nonnegative terms (no cancellation).
// Cost O(e^3 + n nnz(H)) instead of O(n^3) and one n x n matrix (the result): fast when nodes far outnumber hyperedges (cells
// against genes and interactions). Not bit-identical to lu (~1e-15 relative); the same at any `threads` (the BLAS's own thread
// count can change the last bits of the e x e inverse, as it does for lu and chol).
//
// chol and hypergraph stop as solve() does when the system is computationally singular (LAPACK dpocon, reciprocal condition number
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
#ifndef FCONE
#define FCONE
#endif

namespace {

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
    const Rcpp::IntegerVector dim = Rf_getAttrib(adj, Rf_install("Dim"));
    n = dim[0];
    if(dim[1] != n) Rcpp::stop("'a' (%d x %d) must be square", n, dim[1]);
  }
  if(n == 0) Rcpp::stop("'a' is 0-diml");

  return n;
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
  const bool dense = Rf_isMatrix(adj);
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("solve_influence_cpp: one row sum per row");
  const size_t N = n;

  Rcpp::NumericMatrix M = dense ? Rcpp::NumericMatrix(Rcpp::no_init(n, n)) : Rcpp::NumericMatrix(n, n);
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
        m[row + j * N] = x[k];
        if(sym && row != j) m[j + row * N] = x[k];
      }
    }
  }
  const double c = 1 - r;
  for(size_t j = 0; j < N; ++j){
    for(size_t i = 0; i < N; ++i){
      double w = m[i + j * N] / rs[i];
      w = c * w;
      m[i + j * N] = (i == j ? 1.0 : 0.0) - w;
    }
  }

  const char one = '1';
  const double anorm = F77_CALL(dlange)(&one, &n, &n, m, &n, nullptr FCONE);
  Rcpp::NumericMatrix X(n, n);
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
  const bool dense = Rf_isMatrix(adj), sym = Rf_inherits(adj, "dsCMatrix");
  int n = adjacency_order(adj);
  if(rs.size() != n) Rcpp::stop("solve_influence_chol_cpp: one row sum per row");
  const size_t N = n, T = 64;

  Rcpp::NumericMatrix S = dense ? Rcpp::NumericMatrix(Rcpp::no_init(n, n)) : Rcpp::NumericMatrix(n, n);
  double *m = REAL(S);
  if(dense){
    const Rcpp::NumericVector a(adj);
    std::copy(a.begin(), a.end(), m);
  } else {
    const Rcpp::IntegerVector i = Rf_getAttrib(adj, Rf_install("i")), p = Rf_getAttrib(adj, Rf_install("p"));
    const Rcpp::NumericVector x = Rf_getAttrib(adj, Rf_install("x"));
    for(size_t j = 0; j < N; ++j){
      for(int k = p[j]; k < p[j + 1]; ++k){
        const size_t row = i[k];
        if(!sym || row <= j){
          m[row + j * N] = x[k];
        } else {
          m[j + row * N] = x[k];
        }
      }
    }
  }
  if(!sym){
    for(size_t jb = 0; jb < N; jb += T){
      for(size_t ib = 0; ib <= jb; ib += T){
        for(size_t j = jb; j < std::min(jb + T, N); ++j){
          for(size_t i = ib; i < std::min(ib + T, j); ++i){
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

  for(size_t jb = 0; jb < N; jb += T){
    for(size_t ib = 0; ib <= jb; ib += T){
      for(size_t j = jb; j < std::min(jb + T, N); ++j){
        for(size_t i = ib; i < std::min(ib + T, j + 1); ++i){
          const double u = m[i + j * N];
          m[i + j * N] = (u * row[i]) * col[j];
          if(i != j) m[j + i * N] = (u * row[j]) * col[i];
        }
      }
    }
  }

  return S;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_influence_hypergraph_cpp(Rcpp::IntegerVector hp, Rcpp::IntegerVector hi, Rcpp::NumericVector hx, int n,
                                                   Rcpp::NumericVector rs, double r, int threads){
  if(n < 1 || rs.size() != n) Rcpp::stop("solve_influence_hypergraph_cpp: one row sum per node");
  if(hp.size() < 2 || hp[0] != 0 || hp[hp.size() - 1] != hi.size() || hi.size() != hx.size()){
    Rcpp::stop("solve_influence_hypergraph_cpp: H must be the slots of a node x hyperedge dgCMatrix");
  }
  const size_t N = n, E = hp.size() - 1, NNZ = hi.size();
  for(size_t e = 0; e < E; ++e){
    if(hp[e + 1] < hp[e]) Rcpp::stop("solve_influence_hypergraph_cpp: H must be the slots of a node x hyperedge dgCMatrix");
  }
  for(size_t k = 0; k < NNZ; ++k){
    if(hi[k] < 0 || hi[k] >= n || !(hx[k] > 0) || !std::isfinite(hx[k])){
      Rcpp::stop("solve_influence_hypergraph_cpp: H must hold positive incidence counts of the nodes");
    }
  }
  const double c = 1 - r;

  // incidence by node, hyperedges ascending; hyperedge sizes; Delta; Y = H / Delta
  std::vector<size_t> nptr(N + 1, 0);
  for(size_t k = 0; k < NNZ; ++k) ++nptr[hi[k] + 1];
  for(size_t i = 0; i < N; ++i) nptr[i + 1] += nptr[i];
  std::vector<int> nedge(NNZ);
  std::vector<double> nh(NNZ), size(E, 0);
  std::vector<size_t> fill(nptr.begin(), nptr.end() - 1);
  for(size_t e = 0; e < E; ++e){
    for(int k = hp[e]; k < hp[e + 1]; ++k){
      nedge[fill[hi[k]]] = (int) e;
      nh[fill[hi[k]]++] = hx[k];
      size[e] += hx[k];
    }
  }
  std::vector<double> delta(N), ny(NNZ);
  for(size_t i = 0; i < N; ++i){
    double d = 0, sq = 0;
    for(size_t k = nptr[i]; k < nptr[i + 1]; ++k){
      d += nh[k] * (size[nedge[k]] - nh[k]);
      sq += nh[k] * nh[k];
    }
    if(!(rs[i] > 0) || !std::isfinite(rs[i])) Rcpp::stop("method \"hypergraph\" needs a positive row sum for every node of adj");
    if(std::fabs(d - rs[i]) > 1e-12 * rs[i]){
      Rcpp::stop("hyperedges do not match adj: node %d has row sum %g in adj and %g from the hyperedges", (int) i + 1, rs[i], d);
    }
    delta[i] = rs[i] + c * sq;
    for(size_t k = nptr[i]; k < nptr[i + 1]; ++k) ny[k] = nh[k] / delta[i];
  }

  // C = I - (1 - r) H' Y, lower triangle by column, then C^-1 in place (dpotrf + dpotri) and mirrored
  std::vector<double> C;
  try{
    C.assign(E * E, 0.0);
  } catch(std::bad_alloc &){
    Rcpp::stop("method \"hypergraph\": not enough memory for the %d x %d hyperedge matrix", (int) E, (int) E);
  }
  double *cm = C.data();
  run_tasks(threads, E, [&](int, size_t f){
    double *col = cm + f * E;
    for(int k = hp[f]; k < hp[f + 1]; ++k){
      const size_t i = hi[k];
      const double yf = hx[k] / delta[i];
      const int *first = nedge.data() + nptr[i], *last = nedge.data() + nptr[i + 1];
      for(const int *q = std::lower_bound(first, last, (int) f); q < last; ++q) col[*q] += nh[q - nedge.data()] * yf;
    }
    for(size_t e = f; e < E; ++e) col[e] = (e == f ? 1.0 : 0.0) - c * col[e];
  });
  const char uplo = 'L', one = '1';
  int ne = (int) E, info = 0;
  std::vector<double> work(3 * E);
  std::vector<int> iwork(E);
  const double anorm = F77_CALL(dlansy)(&one, &uplo, &ne, cm, &ne, work.data() FCONE FCONE);
  F77_CALL(dpotrf)(&uplo, &ne, cm, &ne, &info FCONE);
  if(info != 0) Rcpp::stop("method \"hypergraph\": the hyperedge matrix is not positive definite (Lapack dpotrf info %d; use method=\"lu\")", info);
  double rcond = 0;
  F77_CALL(dpocon)(&uplo, &ne, cm, &ne, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
  if(rcond < DBL_EPSILON) Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);
  F77_CALL(dpotri)(&uplo, &ne, cm, &ne, &info FCONE);
  if(info != 0) Rcpp::stop("Lapack routine %s failed with info %d", "dpotri", info);
  for(size_t f = 0; f < E; ++f){
    for(size_t e = f + 1; e < E; ++e) cm[f + e * E] = cm[e + f * E];
  }

  // F by blocks of B columns: Z = C^-1 Y' D for the block (hyperedge-major, B values per hyperedge), then every row of the block
  // as Y Z, the sum over the node's hyperedges in ascending order; the diagonal adds r D / Delta
  constexpr size_t B = 32;
  const size_t nblocks = (N + B - 1) / B;
  const int nt = std::max(1, std::min<int>(threads, (int) nblocks));
  std::vector<std::vector<double>> zt, zc;
  try{
    zt.assign(nt, std::vector<double>(E * B));
    zc.assign(nt, std::vector<double>(E));
  } catch(std::bad_alloc &){
    Rcpp::stop("method \"hypergraph\": not enough memory for the thread buffers");
  }
  Rcpp::NumericMatrix F(Rcpp::no_init(n, n));
  double *fm = REAL(F);
  const double rc = r * c;
  run_tasks(nt, nblocks, [&](int w, size_t b){
    double *z = zt[w].data(), *zcol = zc[w].data();
    const size_t j0 = b * B, nb = std::min(B, N - j0);
    std::fill(z, z + E * B, 0.0);
    for(size_t jj = 0; jj < nb; ++jj){
      const size_t j = j0 + jj;
      std::fill(zcol, zcol + E, 0.0);
      for(size_t k = nptr[j]; k < nptr[j + 1]; ++k){
        const double wj = ny[k] * rs[j];
        const double *cf = cm + (size_t) nedge[k] * E;
        for(size_t e = 0; e < E; ++e) zcol[e] += cf[e] * wj;
      }
      for(size_t e = 0; e < E; ++e) z[e * B + jj] = zcol[e];
    }
    for(size_t i = 0; i < N; ++i){
      double acc[B] = {0};
      for(size_t k = nptr[i]; k < nptr[i + 1]; ++k){
        const double y = ny[k];
        const double *ze = z + (size_t) nedge[k] * B;
        for(size_t jj = 0; jj < B; ++jj) acc[jj] += y * ze[jj];
      }
      for(size_t jj = 0; jj < nb; ++jj) fm[i + (j0 + jj) * N] = rc * acc[jj];
    }
    for(size_t jj = 0; jj < nb; ++jj){
      const size_t j = j0 + jj;
      fm[j + j * N] += r * rs[j] / delta[j];
    }
  });

  return F;
}
