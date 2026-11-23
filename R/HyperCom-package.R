# R only: cell grids are scored in R (R/grids.R), in parallel through foreach on a cluster (cl). Importing a Matrix class loads the
# Matrix namespace with HyperCom, so the bundled sparse counts can be subset right after library(HyperCom).

#' @importClassesFrom Matrix dgCMatrix dsCMatrix
#' @keywords internal
"_PACKAGE"
