# Pipeline steps 0-1: counts of the ligand and receptor genes, and the interactions whose genes are all expressed

#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression. Missing (NA) counts are ignored.
#' A gene name used by several rows keeps its first row, with a warning. Raw and normalised counts give the same results (HyperCom uses
#' whether each gene is above 0 in each cell); counts with negative values (scaled data) and counts with genes as column names stop.
#' @param counts A matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions, e.g. `cellchat.mm.lrdb`
#' @param outfile A string naming the output file, or NULL to write none
#' @returns A matrix with row names being cells and column names being genes of only ligand and receptor genes with greater than 0 expression,
#' sparse if counts is sparse
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' dim(counts)
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  check_lrdb(lrdb)
  genes <- unique(unlist(lapply(lrdb[, -1, drop=FALSE], as.character), use.names=FALSE))
  features <- intersect(rownames(counts), genes)
  if(length(features) == 0 && any(colnames(counts) %in% genes)){
    stop("counts has genes as column names: give a gene x cell matrix, e.g. t(counts)")
  }
  dup.genes <- unique(rownames(counts)[duplicated(rownames(counts)) & rownames(counts) %in% features])
  if(length(dup.genes) > 0){
    warning("gene names used by several rows of counts, the first row is kept: ", paste(dup.genes, collapse=", "))
  }
  exp <- Matrix::t(counts[features, , drop=FALSE])
  if(length(exp) > 0 && min(exp, na.rm=TRUE) < 0){
    stop("counts has negative values (scaled data?): give raw or normalised counts; HyperCom uses whether each gene is above 0 in each cell")
  }
  genes.exp <- Matrix::colMeans(exp, na.rm=TRUE)
  genes.filter <- features[which(genes.exp > 0)]
  cc.counts <- exp[, genes.filter, drop=FALSE]

  if(!is.null(outfile)){
    saveRDS(cc.counts, outfile)
  }
  return(cc.counts)
}

#' Check a table of ligand-receptor interactions
#'
#' `check_lrdb` stops unless lrdb looks like a table of ligand-receptor interactions
#'
#' This function takes the lrdb argument and stops, naming the expected format, unless it is a dataframe with a column named ligand
#' and at least one receptor column, so that swapped arguments fail at once.
#' @param lrdb A dataframe of ligand-receptor interactions
#' @returns NULL, invisibly
#' @noRd
check_lrdb <- function(lrdb){
  if(!is.data.frame(lrdb) || !("ligand" %in% colnames(lrdb)) || ncol(lrdb) < 3){
    stop("lrdb must be a dataframe of ligand-receptor interactions (interaction, ligand and receptor columns) like cellchat.mm.lrdb")
  }

  return(invisible(NULL))
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns a dataframe of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptor interactions to produce
#' a table of the ligand-receptor interactions whose genes are all in the counts matrix.
#' @param counts A filtered matrix with row names being cells and column names being genes, from filter_counts
#' @inheritParams filter_counts
#' @returns A dataframe of ligand-receptor interactions in the format of lrdb with only interactions present in the counts matrix, and
#' without receptor columns that are empty for all of them
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' head(lrs)
#' @export
generate_lrs_table <- function(counts, lrdb, outfile="1.lrs.csv"){
  genes <- as.matrix(lrdb[, -1, drop=FALSE])
  found <- matrix(genes %in% colnames(counts), nrow=nrow(genes), ncol=ncol(genes)) | is.na(genes)
  keep <- rowSums(!found) == 0
  lrdb <- lrdb[keep, , drop=FALSE]
  if(nrow(lrdb) > 0){
    lrdb <- Filter(function(x)!all(is.na(x)), lrdb)
  }

  if(!is.null(outfile)){
    utils::write.csv(lrdb, outfile, row.names=FALSE)
  }
  return(lrdb)
}
