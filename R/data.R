#' Example counts of 1000 mouse heart cells
#'
#' Log-normalised expression of 1000 cells from mouse hearts after sham surgery or transverse aortic constriction (TAC), a subset of
#' Alexanian et al. 2021, for trying HyperCom with `cellchat.mm.lrdb`.
#' @format A sparse gene x cell matrix (dgCMatrix) of 17083 genes (mouse gene symbols) and 1000 cells (barcodes with the sample as
#' suffix, e.g. "AACCGCGTCCAGAGGA-2"), in the order of the rows of `alexanian2021.metadata`
#' @name alexanian2021.counts
#' @docType data
#' @author Michael Alexanian \email{michael.alexanian@@gladstone.ucsf.edu}
#' @references \url{https://doi.org/10.1038/s41586-021-03674-1}
#' @keywords data
NULL

#' Example metadata of 1000 mouse heart cells
#'
#' Cell metadata of `alexanian2021.counts`: 500 cells of a sham heart (sample 2) and 500 of a TAC heart (sample 3).
#' @format A dataframe with row names being the cells of `alexanian2021.counts` and columns sample (integer), condition ("SHAM" or
#' "TAC"), celltype (10 cell types) and UMAP_1, UMAP_2 (UMAP coordinates)
#' @name alexanian2021.metadata
#' @docType data
#' @author Michael Alexanian \email{michael.alexanian@@gladstone.ucsf.edu}
#' @references \url{https://doi.org/10.1038/s41586-021-03674-1}
#' @keywords data
NULL

#' CellChat human ligand-receptor interactions
#'
#' The ligand-receptor interactions of the CellChat database for H. sapiens, in the format HyperCom reads as `lrdb`.
#' @format A dataframe of 2240 interactions with columns interaction (its name), ligand (a human gene symbol) and r1 to r4 (the
#' receptor genes, NA where a receptor complex has fewer subunits)
#' @name cellchat.hs.lrdb
#' @docType data
#' @author Suoqin Jin \email{sqjin@@whu.edu.cn}
#' @references \url{https://doi.org/10.1038/s41467-021-21246-9}
#' @keywords data
NULL

#' CellChat mouse ligand-receptor interactions
#'
#' The ligand-receptor interactions of the CellChat database for M. musculus, in the format HyperCom reads as `lrdb`.
#' @format A dataframe of 2377 interactions with columns interaction (its name), ligand (a mouse gene symbol) and r1 to r4 (the
#' receptor genes, NA where a receptor complex has fewer subunits)
#' @name cellchat.mm.lrdb
#' @docType data
#' @author Suoqin Jin \email{sqjin@@whu.edu.cn}
#' @references \url{https://doi.org/10.1038/s41467-021-21246-9}
#' @keywords data
NULL
