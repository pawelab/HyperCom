#' HyperCom: a hypergraph method to understand cell-cell communication
#'
#' Branch v1-fast: the HyperCom method (Moy & Przytycki) with the same outputs, computed faster
#' and in less memory -- the influence matrix by Cholesky on a multithreaded BLAS, and the per-LR
#' scoring by a C++ kernel that never materialises the cell x cell grid. See the README for the
#' benchmarks and the `FLAG:MODIFIED_20261025` blocks in `R/functions.R` for what changed where.
#'
#' @keywords internal
#' @useDynLib HyperCom, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
"_PACKAGE"
