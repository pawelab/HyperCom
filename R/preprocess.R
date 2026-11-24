# Pipeline steps 0-1: counts of the ligand and receptor genes, and the interactions whose genes are all expressed

#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression. Missing (NA) counts are ignored.
#' A gene name used by several rows keeps its first row, with a warning. Raw and normalized counts give the same results (HyperCom uses
#' whether each gene is above 0 in each cell); counts with negative values (scaled data) and counts with genes as column names stop.
#' Counts held another way than a base or Matrix matrix (e.g. a DelayedArray of HDF5-backed data) are read for the genes of lrdb
#' only, into a sparse matrix.
#' @param counts A matrix (dense or sparse) with row names being genes and column names being cells: raw or normalized counts, not
#' scaled data
#' @param lrdb A data frame of ligand-receptor interactions where the first column names the interaction, the second column is named
#' ligand, and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of
#' ligand-receptor interactions, e.g. `cellchat.mm.lrdb`. Interactions without a ligand or a receptor are left out with a warning;
#' interaction names must be distinct and differ from the gene names
#' @param outfile A string naming the output file, or NULL to write none
#' @returns A matrix with row names being cells and column names being genes of only ligand and receptor genes with greater than 0 expression,
#' sparse if counts is sparse
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' dim(counts)
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  lrdb <- clean_lrdb(lrdb)
  genes <- lrdb_genes(lrdb)
  counts <- read_counts(counts, genes)
  features <- intersect(rownames(counts), genes)
  if(length(features) == 0 && any(colnames(counts) %in% genes)){
    stop("counts has genes as column names: give a gene x cell matrix, e.g. t(counts)")
  }
  dup.genes <- unique(rownames(counts)[duplicated(rownames(counts)) & rownames(counts) %in% features])
  if(length(dup.genes) > 0){
    warning("gene names used by several rows of counts, the first row is kept: ", paste(dup.genes, collapse=", "))
  }
  lr.exp <- Matrix::t(counts[features, , drop=FALSE])
  if(length(lr.exp) > 0 && min(lr.exp, na.rm=TRUE) < 0){
    stop("counts has negative values (scaled data?): give raw or normalized counts; HyperCom uses whether each gene is above 0 in each ",
         "cell")
  }
  genes.exp <- Matrix::colMeans(lr.exp, na.rm=TRUE)
  genes.filter <- features[which(genes.exp > 0)]
  cc.counts <- lr.exp[, genes.filter, drop=FALSE]

  if(!is.null(outfile)){
    saveRDS(cc.counts, outfile)
  }
  return(cc.counts)
}

#' Clean a table of ligand-receptor interactions
#'
#' `clean_lrdb` returns a table of ligand-receptor interactions with checked columns and clean gene names
#'
#' This function takes a table of ligand-receptor interactions (lrdb, or lrs from generate_lrs_table) and returns it as a data frame
#' of strings: factors become strings, names lose leading and trailing spaces, and empty names become NA. It stops unless the first
#' column names the interactions, the second is named ligand and the others are receptor columns named r1, r2, ...; unless the
#' interaction names are present and distinct and differ from every gene name (interactions and genes both name hyperedges, so a
#' shared name would merge them); and, unless empty.ok, when no interaction is left. Interactions without a ligand or without any
#' receptor are left out, with a warning.
#' @param lrdb A data frame of ligand-receptor interactions
#' @param arg A string naming the argument in messages, "lrdb" or "lrs"
#' @param empty.ok A boolean of whether a table without interactions is allowed
#' @returns A data frame of ligand-receptor interactions in the format of lrdb
#' @noRd
clean_lrdb <- function(lrdb, arg="lrdb", empty.ok=FALSE){
  if(!is.data.frame(lrdb)){
    stop(arg, " must be a data frame of ligand-receptor interactions (columns interaction, ligand, r1, r2, ...) like cellchat.mm.lrdb")
  }
  cols <- colnames(lrdb)
  unexpected <- c(if(length(cols) >= 2 && cols[2] != "ligand") cols[2], cols[-(1:2)][!grepl("^r[0-9]+$", cols[-(1:2)])])
  if(ncol(lrdb) < 3 || length(unexpected) > 0){
    hint <- ""
    if(any(c("receptor", "interaction_name", "pathway_name") %in% cols)){
      hint <- paste0(". CellChat's interaction table names receptor complexes instead of their genes: use cellchat.mm.lrdb or ",
                     "cellchat.hs.lrdb, or put the genes of each complex in r1, r2, ...")
    }
    other <- setdiff(unexpected, "ligand")
    other <- other[!grepl("^r[0-9]+$", other)]
    ligand.pos <- match("ligand", cols)
    stop(arg, " must have the columns interaction, ligand, r1, r2, ... in that order (any name for the first), like cellchat.mm.lrdb",
         if(is.na(ligand.pos)) "; it has no column named ligand",
         if(!is.na(ligand.pos) && ligand.pos != 2) paste0("; ligand is column ", ligand.pos, ", not 2"),
         if(length(other) > 0) paste0("; columns not expected: ", paste(utils::head(other, 10), collapse=", ")),
         ". Keep interaction, ligand, r1, r2, ... in that order, e.g. ", arg, "[, c(\"interaction\", \"ligand\", \"r1\", \"r2\")]", hint)
  }

  lrdb <- as.data.frame(lrdb)
  lrdb[] <- lapply(lrdb, function(x){
    x <- trimws(as.character(x))
    x[x %in% ""] <- NA
    return(x)
  })
  no.ligand <- is.na(lrdb$ligand)
  no.receptor <- rowSums(!is.na(lrdb[, -(1:2), drop=FALSE])) == 0
  if(any(no.ligand | no.receptor)){
    warning(sum(no.ligand | no.receptor), " interactions of ", arg, " have no ligand or no receptor and are left out: ",
            paste(utils::head(lrdb[[1]][no.ligand | no.receptor], 10), collapse=", "))
    lrdb <- lrdb[!(no.ligand | no.receptor), , drop=FALSE]
  }
  if(nrow(lrdb) == 0 && !empty.ok){
    stop(arg, " has no interactions")
  }

  names <- lrdb[[1]]
  if(anyNA(names)){
    stop("every interaction of ", arg, " needs a name in its first column")
  }
  dup.names <- unique(names[duplicated(names)])
  if(length(dup.names) > 0){
    stop("interaction names of ", arg, " must be distinct; used more than once: ", paste(utils::head(dup.names, 10), collapse=", "))
  }
  gene.names <- names[names %in% lrdb_genes(lrdb)]
  if(length(gene.names) > 0){
    stop("interaction names of ", arg, " must differ from its gene names (both name hyperedges); named like a gene: ",
         paste(utils::head(gene.names, 10), collapse=", "))
  }

  return(lrdb)
}

#' List the genes of a table of ligand-receptor interactions
#'
#' `lrdb_genes` returns the distinct ligand and receptor genes of a table of ligand-receptor interactions
#'
#' This function takes a table of ligand-receptor interactions and returns its gene names, without NA.
#' @param lrdb A data frame of ligand-receptor interactions from clean_lrdb
#' @returns A vector of strings of gene names
#' @noRd
lrdb_genes <- function(lrdb){
  genes <- unlist(lapply(lrdb[, -1, drop=FALSE], as.character), use.names=FALSE)

  return(unique(genes[!is.na(genes)]))
}

#' Read counts held another way than a matrix
#'
#' `read_counts` returns counts as a base matrix, a Matrix or a data frame
#'
#' This function takes a gene x cell counts object and the genes of lrdb. Base and Matrix matrices and data frames are returned as
#' they are; any other object (e.g. a DelayedArray of HDF5-backed data) is read for the given genes only, into a sparse matrix, so the
#' rest of it is never loaded.
#' @param counts A matrix-like object with row names being genes and column names being cells
#' @param genes A vector of strings of gene names
#' @returns A matrix with row names being genes and column names being cells
#' @noRd
read_counts <- function(counts, genes){
  if(is.matrix(counts) || methods::is(counts, "Matrix") || is.data.frame(counts)){
    return(counts)
  }
  if(is.null(rownames(counts)) || is.null(colnames(counts))){
    stop("counts needs row names (genes) and column names (cells)")
  }

  return(methods::as(counts[rownames(counts) %in% genes, , drop=FALSE], "dgCMatrix"))
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns a data frame of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptor interactions to produce
#' a table of the ligand-receptor interactions whose genes are all in the counts matrix.
#' @param counts A filtered matrix with row names being cells and column names being genes, from filter_counts
#' @inheritParams filter_counts
#' @returns A data frame of ligand-receptor interactions in the format of lrdb with only interactions present in the counts matrix, and
#' without receptor columns that are empty for all of them
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' head(lrs)
#' @export
generate_lrs_table <- function(counts, lrdb, outfile="1.lrs.csv"){
  lrdb <- clean_lrdb(lrdb)
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
