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
