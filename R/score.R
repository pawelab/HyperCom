# Cell scores: sender and receiver scores of cell grids (src/score_grids.cpp) and the HyperCom score of a pair

#' Check that genes are nodes of the influence matrix
#'
#' `check_genes` stops if a gene is not a node of the influence matrix
#'
#' This function takes gene names and the node names of an influence matrix and stops naming the genes that are not among them.
#' @param genes A vector of strings of gene names, NA allowed
#' @param nodes A vector of strings of the node names of the influence matrix
#' @returns NULL, invisibly
#' @noRd
check_genes <- function(genes, nodes){
  missing <- setdiff(stats::na.omit(genes), nodes)
  if(length(missing) > 0){
    stop("not nodes of infMat: ", paste(missing, collapse=", "))
  }

  return(invisible(NULL))
}

#' Find the cells of the metadata in an influence matrix
#'
#' `grid_cells` returns the positions of the cells in the influence matrix, ascending, with their samples
#'
#' This function takes the node names of an influence matrix and the metadata and returns the 0-based positions of the metadata's
#' cells among the nodes in increasing order, the sample number of each and the order that sorted them.
#' @param nodes A vector of strings of the node names
#' @param metadata A dataframe containing columns named cell, sample and group
#' @returns A list with idx (positions), sample (integers from 1) and o (the order of the metadata's cells)
#' @noRd
grid_cells <- function(nodes, metadata){
  if(anyNA(metadata$sample)){
    stop("metadata has missing samples: make it with generate_metadata")
  }
  cells <- metadata$sample != "gene"
  if(!any(cells)){
    stop("metadata has no cells: make it with generate_metadata from metadata whose row names are the cell names")
  }
  cell.idx <- match(metadata$cell[cells], nodes) - 1L
  if(anyNA(cell.idx) || anyDuplicated(cell.idx) > 0){
    stop("metadata cells must be distinct nodes of infMat: make it with generate_metadata")
  }
  samples <- as.integer(factor(metadata$sample[cells]))
  o <- order(cell.idx)

  return(list(idx=cell.idx[o], sample=samples[o], o=o))
}

#' Calculate sender and receiver scores of cell grids
#'
#' `score_grids` returns the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, and vectors of ligands and receptor subunits.
#' For each pair it percent-ranks all entries of the cell x cell transfer grid between cells of the same sample together and averages
#' the ranks of each cell as sender and as receiver. A pair with ligand and receptor NA gives the network score.
#' Cells and genes are found in the influence matrix by name, so the metadata rows may be in any order. A factored influence
#' matrix is read through its within-sample blocks (made here unless given) and the genes' rows and columns.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligands A vector of strings of ligand names, NA for the network score
#' @param receptors A vector of strings of receptor subunit names, NA for the network score
#' @param threads An integer number of threads
#' @param blocks A numeric vector from influence_blocks for a hypercom_influence, NULL to make it here
#' @returns A list of two matrices, sender and receiver, with one row per cell and one column per pair
#' @noRd
score_grids <- function(infMat, metadata, ligands, receptors, threads=1, blocks=NULL){
  ligands <- as.character(ligands)
  receptors <- as.character(receptors)
  nodes <- rownames(infMat)
  check_genes(c(ligands, receptors), nodes)
  if(any(is.na(ligands) != is.na(receptors))){
    stop("each pair needs both a ligand and a receptor, or neither for the network score")
  }
  g <- grid_cells(nodes, metadata)
  lig <- match(ligands, nodes) - 1L
  rec <- match(receptors, nodes) - 1L
  lig[is.na(lig)] <- -1L
  rec[is.na(rec)] <- -1L

  if(inherits(infMat, "hypercom_influence")){
    if(is.null(blocks)){
      blocks <- influence_blocks(infMat, g, threads)
    }
    lig.genes <- unique(lig[lig >= 0])
    rec.genes <- unique(rec[rec >= 0])
    A <- influence_slice(infMat, g$idx + 1L, lig.genes + 1L, threads)
    B <- influence_slice(infMat, rec.genes + 1L, g$idx + 1L, threads)
    sr <- score_blocks_cpp(blocks, A, B, g$sample, match(lig, lig.genes, nomatch=0L) - 1L, match(rec, rec.genes, nomatch=0L) - 1L, threads)
  } else {
    sr <- score_grids_cpp(infMat, g$idx, g$sample, lig, rec, threads)
  }
  if(is.unsorted(g$o)){
    sr <- list(sender=sr$sender[order(g$o), , drop=FALSE], receiver=sr$receiver[order(g$o), , drop=FALSE])
  }

  return(sr)
}

#' Average the cell scores of the subunits of a receptor complex
#'
#' `average_scores` returns the sender and receiver scores of a ligand and receptor complex
#'
#' This function takes the scores of grids and the columns of a complex's subunit grids and returns the sum of the columns,
#' added in order, divided by their number (one column is returned as is).
#' @param sr A list of two matrices, sender and receiver, with one row per cell and one column per grid
#' @param cols An integer vector of the columns of the complex's subunits
#' @returns A list of two vectors, sender and receiver
#' @noRd
average_scores <- function(sr, cols){
  ss <- sr$sender[, cols[1]]
  rr <- sr$receiver[, cols[1]]
  if(length(cols) > 1){
    for(k in cols[-1]){
      ss <- ss + sr$sender[, k]
      rr <- rr + sr$receiver[, k]
    }
    ss <- ss / length(cols)
    rr <- rr / length(cols)
  }

  return(list(sender=ss, receiver=rr))
}

#' Calculate the HyperCom score for each cell
#'
#' `score_hypercom` returns a dataframe containing the HyperCom score for each cell
#'
#' This function takes a named influence matrix, the corresponding metadata, the ligand, and receptor
#' to give each cell a HyperCom score. With a ligand and receptors this is the ligand-receptor pair score; with neither, the network
#' score. For a receptor complex the scores of its subunits are averaged. Scores are the same at any number of threads.
#' @inheritParams prioritize_lr
#' @param ligand A string of the ligand name, NULL for the network score
#' @param receptors A vector containing string elements of the receptor names, NULL for the network score
#' @param threads An integer number of threads of at least 1, overriding parallel; more than the available cores warns and uses
#' them all
#' @returns A dataframe with one row per cell, sorted by hypercom (highest first), with columns cell, sender and receiver (the cell's
#' sender and receiver scores), dif (receiver - sender), avg (their mean), max (the larger of the two) and hypercom (the HyperCom score:
#' the sign of dif times the percent rank of its size times avg)
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' metadata <- generate_metadata(infMat, alexanian2021.metadata, sample="sample", group="condition",
#'                               outfile=NULL)
#' scores <- score_hypercom(infMat, metadata, ligand="Thbs1", receptors="Sdc4")
#' head(scores)
#' # a receptor complex: the scores of its subunits are averaged
#' complex <- score_hypercom(infMat, metadata, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2"))
#' # the network score: no ligand and no receptors
#' network <- score_hypercom(infMat, metadata)
#' @export
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL, parallel=FALSE, threads=NULL){
  threads <- get_threads(parallel, threads=threads)
  complex.size <- length(receptors)
  if(is.null(ligand) != is.null(receptors)){
    stop("give both ligand and receptors for a pair score, or neither for the network score")
  }
  if(!is.null(ligand) && length(ligand) != 1){
    stop("ligand must be one gene")
  }
  if(!is.null(receptors) && complex.size == 0){
    stop("receptors must name at least one gene")
  }
  if(is.null(ligand)){
    sr <- score_grids(infMat, metadata, NA, NA, threads)
    cols <- 1
  } else {
    sr <- score_grids(infMat, metadata, rep(ligand, complex.size), receptors, threads)
    cols <- seq_len(complex.size)
  }

  sr <- average_scores(sr, cols)
  sr <- data.frame(cell=metadata$cell[metadata$sample != "gene"], sender=sr$sender, receiver=sr$receiver)
  if(!is.null(receptors) && all(ligand == receptors)){
    sr$dif <- 0
  } else {
    sr$dif <- sr$receiver - sr$sender
  }
  sr$avg <- (sr$receiver + sr$sender) / 2
  sr$max <- pmax(sr$sender, sr$receiver)
  sr$hypercom <- dplyr::percent_rank(abs(sr$dif))
  sr$hypercom <- sign(sr$dif) * sr$hypercom * sr$avg
  sr <- dplyr::arrange(sr, dplyr::desc(hypercom))

  return(sr)
}
