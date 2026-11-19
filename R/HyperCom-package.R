# Cell grids are scored by compiled code (src/score_grids.cpp). Importing a Matrix class loads the Matrix namespace with HyperCom,
# so the bundled sparse counts can be subset right after library(HyperCom).

#' @useDynLib HyperCom, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importClassesFrom Matrix dgCMatrix dsCMatrix
#' @keywords internal
"_PACKAGE"
