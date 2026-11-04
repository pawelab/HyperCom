// infmat_lapack.cpp -- the Cholesky inverse of a symmetric positive definite matrix IN PLACE.
//
// R's chol() and chol2inv() each return a new matrix, so between them two dense L x L copies are
// alive (~165 GB at 100k cells). The LAPACK routines they call, dpotrf (factorise) and dpotri
// (invert from the factor), overwrite their input instead. Calling them here on the one buffer
// gives the same numbers -- the same routines on the same BLAS -- with a single copy.
//
// Steps, exactly as R does them in La_chol() and La_chol2inv():
//   dpotrf('U')  the upper triangle of S becomes U with S = U'U (the lower triangle is not read)
//   dpotri('U')  the upper triangle of U becomes the upper triangle of S^-1
//   mirror       the lower triangle is filled from the upper (S^-1 is symmetric)
// Threads: dpotrf / dpotri run on the BLAS's threads (set by generate_infMat); `nthreads` is for
// the mirror loop only.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
using namespace Rcpp;

// [[Rcpp::export]]
void chol_inverse_inplace_cpp(NumericMatrix S, int nthreads = 1) {
  int n = S.nrow();
  if (S.ncol() != n) stop("chol_inverse_inplace_cpp: S must be square");
  double *a = REAL(S);
  int info = 0;
  char uplo = 'U';
  F77_CALL(dpotrf)(&uplo, &n, a, &n, &info FCONE);
  if (info > 0) stop("the leading minor of order %d is not positive: the matrix is not positive definite", info);
  if (info < 0) stop("dpotrf: argument %d had an illegal value", -info);
  F77_CALL(dpotri)(&uplo, &n, a, &n, &info FCONE);
  if (info > 0) stop("dpotri: the (%d, %d) element of the factor is zero; the inverse could not be computed", info, info);
  if (info < 0) stop("dpotri: argument %d had an illegal value", -info);
  const R_xlen_t N = n;
  if (nthreads < 1) nthreads = 1;
  #pragma omp parallel for num_threads(nthreads) if(nthreads > 1) schedule(dynamic, 64)
  for (R_xlen_t j = 1; j < N; ++j)
    for (R_xlen_t i = 0; i < j; ++i) a[j + i * N] = a[i + j * N];   // lower <- upper
}

// fill_scaled_sparse_cpp -- write S = alpha * diag(s) A diag(s) + I into the dense, zero-filled S,
// straight from the column-compressed slots (p, i, x) of a general sparse A.
//
// generate_infMat used to build this as as.matrix(Diagonal(s) %*% A %*% Diagonal(s)), then scale and
// add the identity. The Matrix products are sparse and cheap, but as.matrix() allocates the dense
// result while the sparse product still exists, and R's sparse -> dense coercion briefly holds a
// second dense L x L copy (38 GB peak for an 18.6 GB matrix at 50k cells; ~160 GB at 100k). Here the
// caller allocates S once (matrix(0, L, L)) and every stored entry is written in place.
// Arithmetic, in the order the old code did it: Matrix scales rows first then columns,
// (x * s_i) * s_j, then scale_inplace_cpp multiplied by alpha, then 1 was added on the diagonal.
// [[Rcpp::export]]
void fill_scaled_sparse_cpp(NumericMatrix S, IntegerVector p, IntegerVector i, NumericVector x,
                            NumericVector s, double alpha, int nthreads = 1) {
  const int n = S.nrow();
  if (S.ncol() != n) stop("fill_scaled_sparse_cpp: S must be square");
  if (p.size() != n + 1 || s.size() != n) stop("fill_scaled_sparse_cpp: p must have length n + 1 and s length n");
  if (i.size() != x.size()) stop("fill_scaled_sparse_cpp: i and x differ in length");
  double *a = REAL(S);
  const int *pp = INTEGER(p), *ii = INTEGER(i);
  const double *xx = REAL(x), *ss = REAL(s);
  const R_xlen_t N = n;
  if (nthreads < 1) nthreads = 1;
  #pragma omp parallel for num_threads(nthreads) if(nthreads > 1) schedule(dynamic, 64)
  for (int j = 0; j < n; ++j) {
    double *col = a + (R_xlen_t) j * N;
    const double sj = ss[j];
    for (int k = pp[j]; k < pp[j + 1]; ++k) {
      const int r = ii[k];
      col[r] = ((xx[k] * ss[r]) * sj) * alpha;
    }
    col[j] = col[j] + 1.0;
  }
}
