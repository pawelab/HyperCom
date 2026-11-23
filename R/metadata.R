# Pipeline step 5: node metadata (cells with their sample and group; gene nodes labelled "gene")

#' Generate cleaned metadata table
#'
#' `generate_metadata` returns a minimal metadata dataframe
#'
#' This function takes an influence matrix and metadata table and creates a metadata data table containing the
#' node name, sample, and group. If group is unspecified, it is assumed all nodes are part of the same group.
#' Nodes of the influence matrix that are not cells of the metadata are genes: their sample, group and kept columns are "gene".
#' Kept columns named cell, sample or group are renamed with make.unique (e.g. sample.1); kept dates and other classed columns
#' other than factors become strings. Cells are compared only within their sample, so every sample needs at least 2 cells in the
#' graph (cells that express a ligand or receptor gene).
#' @param infMat An influence matrix from generate_infMat (only its node names are read)
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @inheritParams filter_counts
#' @returns A dataframe with columns cell, sample, group and the kept columns, one row per node of the influence matrix
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
generate_metadata <- function(infMat, metadata, sample, group=NULL, keep=NULL, outfile="5.metadata.csv"){
  check_columns(metadata, sample, group, keep)
  index <- data.frame(cell=row.names(infMat), check.names=FALSE)
  cells <- rownames(metadata)

  if(is.null(group)){
    group.vals <- "none"
  } else {
    group.vals <- as.character(metadata[[group]])
  }
  sample.vals <- as.character(metadata[[sample]])

  in.graph <- cells %in% index$cell
  if(!any(in.graph)){
    stop("no node of infMat is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  missing.sample <- sum(is.na(sample.vals[in.graph]))
  missing.group <- sum(is.na(rep(group.vals, length.out=nrow(metadata))[in.graph]))
  if(missing.sample > 0 || missing.group > 0){
    stop("sample and group must not be missing for cells in the graph: ", missing.sample, " cells miss a sample and ", missing.group,
         " a group")
  }
  if(any(sample.vals == "gene", na.rm=TRUE) || any(group.vals == "gene", na.rm=TRUE)){
    stop("\"gene\" labels gene nodes and cannot name a sample or group")
  }
  sizes <- table(sample.vals[in.graph])
  if(any(sizes < 2)){
    stop("every sample needs at least 2 cells in the graph (cells that express a ligand or receptor gene), as cells are compared only ",
         "within their sample; remove the cells of: ", paste(names(sizes)[sizes < 2], collapse=", "))
  }

  clean.df <- data.frame(cell=cells, sample=sample.vals, group=group.vals, check.names=FALSE)
  clean.df$sample <- factor(clean.df$sample, levels=c(unique(clean.df$sample), "gene"))
  clean.df$group <- factor(clean.df$group, levels=c(unique(clean.df$group), "gene"))

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
    clean.df[genes, col] <- "gene"
  }

  if(!is.null(outfile)){
    utils::write.csv(clean.df, outfile, row.names=FALSE)
  }

  return(clean.df)
}

#' Check the metadata columns named by the arguments
#'
#' `check_columns` stops unless sample, group and keep name columns of the metadata
#'
#' This function takes the metadata and the column arguments of generate_metadata and stops, naming the problem, unless metadata is a
#' dataframe, sample and group are each one string and every named column exists.
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata, or NULL
#' @param keep A vector of strings that are additional column names to keep in the provided metadata, or NULL
#' @returns NULL, invisibly
#' @noRd
check_columns <- function(metadata, sample, group=NULL, keep=NULL){
  if(!is.data.frame(metadata)){
    stop("metadata must be a dataframe with row names being cells")
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
