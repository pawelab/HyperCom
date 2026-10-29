// build_grid.cpp -- the cell x cell grid of one ligand-receptor (or network) score,
// built in a single pass over memory.
//
// For every pair of cells (i, j) it computes, in this order,
//
//   v = (a[i] * b[j]) / C[i, j]     ligand-receptor grid   (C[i, j] alone for the network score)
//   v = v * mask[i, j]              only if a same-sample mask is given (TRUE = 1, FALSE = 0)
//   v = NA  if v == 0 or i == j     zeros and the diagonal are left out of the ranking
//   K += 1  if v is a number        the count of values that will be ranked
//
// and writes v once. The R code this replaces did the same arithmetic in five separate passes
// over the n x n grid (Rfast::Tcrossprod, collapse::setop for the division and the mask,
// collapse::setv for zeros, collapse::fnobs for K). The floating-point operations and their order
// are unchanged, so every value is bit-identical; only the memory traffic is lower.
//
// Do not compile with -ffast-math: it lets the compiler reorder the arithmetic.
//
// Threads: columns are split across `nthreads` OpenMP threads, all writing into one grid, so the
// memory used does not grow with the thread count. nthreads = 1 (as used inside forked workers)
// runs without starting any threads.

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix build_grid_cpp(NumericMatrix C,
                             Nullable<NumericVector> a = R_NilValue,
                             Nullable<NumericVector> b = R_NilValue,
                             Nullable<LogicalMatrix> mask = R_NilValue,
                             int nthreads = 1) {
  const R_xlen_t n = C.nrow();
  if (C.ncol() != n) stop("build_grid_cpp: C must be square");

  const bool has_ab = a.isNotNull() && b.isNotNull();
  if (a.isNotNull() != b.isNotNull()) stop("build_grid_cpp: give both a and b, or neither");

  NumericVector av, bv;
  const double *pa = nullptr, *pb = nullptr;
  if (has_ab) {
    av = as<NumericVector>(a.get());
    bv = as<NumericVector>(b.get());
    if (av.size() != n || bv.size() != n) stop("build_grid_cpp: a and b must have length nrow(C)");
    pa = REAL(av);
    pb = REAL(bv);
  }

  LogicalMatrix mv;
  const int *pm = nullptr;
  if (mask.isNotNull()) {
    mv = as<LogicalMatrix>(mask.get());
    if (mv.nrow() != n || mv.ncol() != n) stop("build_grid_cpp: mask must have the dimensions of C");
    pm = LOGICAL(mv);
  }

  // no_init: skipping the zero fill, which would itself be a full pass
  NumericMatrix Q(no_init((int) n, (int) n));
  double *q = REAL(Q);
  const double *pc = REAL(C);
  double K = 0.0;
  if (nthreads < 1) nthreads = 1;

  #pragma omp parallel for num_threads(nthreads) if(nthreads > 1) reduction(+:K) schedule(static)
  for (R_xlen_t j = 0; j < n; ++j) {
    const double bj = has_ab ? pb[j] : 0.0;
    for (R_xlen_t i = 0; i < n; ++i) {
      const R_xlen_t k = i + j * n;
      double v = has_ab ? (pa[i] * bj) / pc[k] : pc[k];
      if (pm) v = v * (double) pm[k];
      if (v == 0 || i == j) {
        v = NA_REAL;
      } else if (!ISNAN(v)) {
        K += 1.0;
      }
      q[k] = v;
    }
  }

  Q.attr("K") = K;
  return Q;
}
