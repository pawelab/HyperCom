# Pure R: cell grids are scored in R (R/grids.R), with no compiled code. Importing a Matrix class loads the Matrix namespace with
# HyperCom, so the bundled sparse counts can be subset right after library(HyperCom).

#' @importFrom dplyr .data
#' @importClassesFrom Matrix dgCMatrix dsCMatrix
#' @keywords internal
"_PACKAGE"
