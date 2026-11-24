# Pipeline step 5: node metadata (cells with their sample and group; gene nodes labelled "gene")

#' Generate cleaned metadata table
#'
#' `generate_metadata` returns a minimal metadata data frame with one row per node of the influence matrix
#'
#' This function takes an influence matrix and a cell metadata table and returns a table of the node name (cell), sample and group of
#' every node. Without sample, every cell is in one sample, "all"; without group, every cell is in one group, "none". Nodes of the influence matrix that are not cells of the metadata
#' are genes: their sample, group and kept text and factor columns are "gene", and their kept numeric and logical columns NA (those
#' columns keep their type). Kept columns named cell, sample or group are renamed with make.unique (e.g. sample.1); kept dates and
#' other classed columns other than factors become text. Cells are compared only within their sample, so every sample needs at
#' least 2 cells in the hypergraph (cells that express a ligand or receptor gene); their sample and group must not be missing or
#' "gene".
#' @param infMat An influence matrix from generate_infMat (only its node names are read)
#' @param metadata A data frame containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample (library or replicate) in the provided metadata, or NULL (the default)
#' for one sample; cells are compared only within their sample
#' @param group A string that is the column name for group (e.g. condition) in the provided metadata, or NULL (the default) for one
#' group; groups feed the purity term of the priority (see `prioritize_lr`)
#' @param keep A vector of strings that are additional column names to keep from the provided metadata, or NULL (the default)
#' @inheritParams filter_counts
#' @returns A data frame with columns cell, sample, group and the kept columns, one row per node of the influence matrix, in its
#' order (gene rows marked "gene")
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' metadata <- generate_metadata(infMat, alexanian2021.metadata, sample="sample", group="condition",
#'                               outfile=NULL)
#' head(metadata)
#' table(metadata$sample)
#' @export
generate_metadata <- function(infMat, metadata, sample=NULL, group=NULL, keep=NULL, outfile="5.metadata.csv"){
  check_columns(metadata, sample, group, keep)
  index <- data.frame(cell=row.names(infMat), check.names=FALSE)
  cells <- rownames(metadata)

  sample.vals <- label_values(metadata, sample, "all")
  group.vals <- label_values(metadata, group, "none")

  in.graph <- cells %in% index$cell
  if(!any(in.graph)){
    stop("no node of infMat is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  check_labels(sample.vals[in.graph], group.vals[in.graph])
  sizes <- table(sample.vals[in.graph])
  if(any(sizes < 2)){
    stop("every sample needs at least 2 cells in the graph (cells that express a ligand or receptor gene), as cells are compared only ",
         "within their sample; remove the cells of: ", paste(names(sizes)[sizes < 2], collapse=", "))
  }

  clean.df <- data.frame(cell=cells, sample=sample.vals, group=group.vals, check.names=FALSE)
  clean.df$sample <- factor(clean.df$sample, levels=unique(c(clean.df$sample, "gene")))
  clean.df$group <- factor(clean.df$group, levels=unique(c(clean.df$group, "gene")))

  if(!is.null(keep)){
    additional <- metadata[, keep, drop=FALSE]
    colnames(additional) <- make.unique(c(colnames(clean.df), colnames(additional)))[-seq_len(ncol(clean.df))]
    for(col in colnames(additional)){
      if(is.factor(additional[[col]])){
        additional[[col]] <- factor(additional[[col]], levels=unique(c(levels(additional[[col]]), "gene")))
      } else if(is.object(additional[[col]])){
        additional[[col]] <- as.character(additional[[col]])
      }
    }
    clean.df <- cbind(clean.df, additional)
  }

  clean.df <- dplyr::left_join(index, clean.df, by=dplyr::join_by(cell))
  genes <- !(clean.df$cell %in% cells)
  for(col in colnames(clean.df)[-1]){
    if(is.factor(clean.df[[col]]) || is.character(clean.df[[col]])){
      clean.df[genes, col] <- "gene"
    } else {
      clean.df[genes, col] <- NA
    }
  }

  if(!is.null(outfile)){
    utils::write.csv(clean.df, outfile, row.names=FALSE)
  }

  return(clean.df)
}

#' Get the sample or group label of every cell
#'
#' `label_values` returns the labels of a metadata column as strings, or one label for every cell when no column is named
#' @param metadata A data frame containing cell metadata with row names being cells
#' @param column A string that is a column name of metadata, or NULL
#' @param single A string labelling every cell when column is NULL
#' @returns A vector of strings, one per row of metadata
#' @noRd
label_values <- function(metadata, column, single){
  if(is.null(column)){
    return(rep(single, nrow(metadata)))
  }

  return(as.character(metadata[[column]]))
}

#' Check the sample and group labels of the cells in the graph
#'
#' `check_labels` stops when a cell in the graph has a missing sample or group, or one named "gene"
#'
#' This function takes the samples and groups of the cells in the graph and stops, counting the cells, when one is missing, and when
#' one is "gene", the label of gene nodes.
#' @param sample.vals A vector of strings of the sample of each cell
#' @param group.vals A vector of strings of the group of each cell
#' @returns NULL, invisibly
#' @noRd
check_labels <- function(sample.vals, group.vals){
  missing.sample <- sum(is.na(sample.vals))
  missing.group <- sum(is.na(group.vals))
  if(missing.sample > 0 || missing.group > 0){
    stop("sample and group must not be missing for cells in the graph: ", missing.sample, " cells miss a sample and ", missing.group,
         " a group")
  }
  if(any(sample.vals == "gene") || any(group.vals == "gene")){
    stop("\"gene\" labels gene nodes and cannot name a sample or group")
  }

  return(invisible(NULL))
}

#' Check metadata made by generate_metadata
#'
#' `check_metadata` stops unless metadata has the columns of generate_metadata, and its gene rows when p-values need them
#'
#' This function takes the metadata argument of prioritize_lr or score_hypercom and stops unless it is a data frame with columns cell,
#' sample and group, and, with significance, gene rows (sample "gene"), from which the background pairs of the p-values are drawn. A
#' subset of the cells, in any order, is allowed.
#' @param metadata A data frame from generate_metadata
#' @param significance A boolean of whether p-values are calculated
#' @returns NULL, invisibly
#' @noRd
check_metadata <- function(metadata, significance=FALSE){
  msg <- "metadata must come from generate_metadata (columns cell, sample and group)"
  if(!is.data.frame(metadata) || !all(c("cell", "sample", "group") %in% colnames(metadata))){
    stop(msg)
  }
  if(significance && !any(metadata$sample == "gene", na.rm=TRUE)){
    stop(msg, ", with its gene rows: the background pairs of the p-values are drawn from them")
  }

  return(invisible(NULL))
}

#' Check the metadata columns named by the arguments
#'
#' `check_columns` stops unless sample, group and keep name columns of the metadata
#'
#' This function takes the metadata and the column arguments of generate_metadata and stops, naming the problem, unless metadata is a
#' data frame, sample and group are each one string or NULL and every named column exists.
#' @param metadata A data frame containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata, or NULL
#' @param group A string that is the column name for group in the provided metadata, or NULL
#' @param keep A vector of strings that are additional column names to keep in the provided metadata, or NULL
#' @returns NULL, invisibly
#' @noRd
check_columns <- function(metadata, sample=NULL, group=NULL, keep=NULL){
  if(!is.data.frame(metadata)){
    stop("metadata must be a data frame with row names being cells")
  }
  for(arg in list(sample, group)){
    if(!is.null(arg) && !(is.character(arg) && length(arg) == 1)){
      stop("sample and group must each be one column name of metadata, e.g. sample=\"sample\"")
    }
  }
  absent <- setdiff(c(sample, group, keep[is.character(keep)]), colnames(metadata))
  if(length(absent) > 0){
    stop("not columns of metadata: ", paste(absent, collapse=", "), " (its columns: ", paste(utils::head(colnames(metadata), 20), collapse=", "),
         ")")
  }

  return(invisible(NULL))
}
