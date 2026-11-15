// r_alloc.h: R vectors and matrices for results that can be large (nodes or cells squared)
//
// When R cannot allocate, its error longjmps past the C++ frames: their destructors do not run, so std::vector buffers leak and
// every Rcpp argument stays preserved (pinned in memory) for the rest of the session. Under Rcpp::unwindProtect the error becomes
// a C++ exception (Rcpp::LongjumpException) that unwinds normally, and the export wrapper's END_RCPP then resumes R's error, so the
// caller still sees R's "cannot allocate vector of size ...". Use on the main thread only, never while worker threads run, and
// never inside a catch(...) that would swallow the exception.

#ifndef HYPERCOM_R_ALLOC_H
#define HYPERCOM_R_ALLOC_H

#include <Rcpp.h>
#include <algorithm>

// An nr x nc double matrix, zero-filled unless zero is false
inline Rcpp::NumericMatrix alloc_numeric_matrix(R_xlen_t nr, R_xlen_t nc, bool zero=true){
  Rcpp::NumericMatrix m(Rcpp::unwindProtect([&]{ return Rf_allocMatrix(REALSXP, (int) nr, (int) nc); }));
  if(zero) std::fill(m.begin(), m.end(), 0.0);

  return m;
}

// A vector of R type RTYPE (REALSXP, INTSXP) of length n, zero-filled unless zero is false
template <int RTYPE>
inline Rcpp::Vector<RTYPE> alloc_vector(R_xlen_t n, bool zero=true){
  Rcpp::Vector<RTYPE> v(Rcpp::unwindProtect([&]{ return Rf_allocVector(RTYPE, n); }));
  if(zero) std::fill(v.begin(), v.end(), typename Rcpp::Vector<RTYPE>::stored_type(0));

  return v;
}

#endif
