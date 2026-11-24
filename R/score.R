# Cell scores: sender and receiver scores of cell grids (score_blocks in R/grids.R) and the HyperCom score of a pair

#' Check that genes are nodes of the influence matrix
#'
#' `check_genes` stops if a gene is not a node of the influence matrix
#'
#' This function takes gene names and the node names of an influence matrix and stops naming the first 10 genes that are not among
#' them and counting the rest.
#' @param genes A vector of strings of gene names, NA allowed
#' @param nodes A vector of strings of the node names of the influence matrix
#' @returns NULL, invisibly
#' @noRd
check_genes <- function(genes, nodes){
  absent <- setdiff(genes[!is.na(genes)], nodes)
  if(length(absent) > 0){
    stop("not nodes of infMat (not expressed, not in lrdb, not connected to the rest of the hypergraph, or another species or case): ",
         paste(utils::head(absent, 10), collapse=", "),
         if(length(absent) > 10) paste0(" and ", length(absent) - 10, " more"))
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
#' @param metadata A data frame containing columns named cell, sample and group
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

#' Gather what score_blocks reads
#'
#' `grid_inputs` returns the cells, genes and influence values of cell grids, as score_grids passes them to score_blocks
#'
#' This function takes a named influence matrix, the metadata and vectors of ligands and receptor subunits (NA in both for the network
#' grid), and returns the cells (grid_cells), the 0-based node index of each ligand and receptor (-1 for the network grid), the values of
#' the cells to the distinct ligands (A) and of the distinct receptors to the cells (B) with the column of A and row of B of each grid
#' (-1 for the network grid).
#' @inheritParams score_grids
#' @returns A list with g, lig, rec, A, B, lc and rc
#' @noRd
grid_inputs <- function(infMat, metadata, ligands, receptors){
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
  gi <- list(g=g, lig=lig, rec=rec)
  lig.genes <- unique(lig[lig >= 0])
  rec.genes <- unique(rec[rec >= 0])
  if(inherits(infMat, "hypercom_influence")){
    gi$A <- influence_slice(infMat, g$idx + 1L, lig.genes + 1L)
    gi$B <- influence_slice(infMat, rec.genes + 1L, g$idx + 1L)
  } else {
    gi$A <- infMat[g$idx + 1L, lig.genes + 1L, drop=FALSE]
    gi$B <- infMat[rec.genes + 1L, g$idx + 1L, drop=FALSE]
  }
  gi$lc <- match(lig, lig.genes, nomatch=0L) - 1L
  gi$rc <- match(rec, rec.genes, nomatch=0L) - 1L

  return(gi)
}

#' Put cell rows in the order of the metadata
#'
#' `metadata_order` returns matrices with one row per cell, rows reordered from the order of grid_cells to that of the metadata
#' @param x A list of matrices with one row per cell, in the order of g$idx
#' @param g A list from grid_cells
#' @returns The list of matrices, rows in the order of the metadata's cells
#' @noRd
metadata_order <- function(x, g){
  if(!is.unsorted(g$o)){
    return(x)
  }
  back <- order(g$o)

  return(lapply(x, function(m) m[back, , drop=FALSE]))
}

#' Calculate sender and receiver scores of cell grids
#'
#' `score_grids` returns the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, and vectors of ligands and receptor subunits.
#' For each pair it percent-ranks all entries of the cell x cell transfer grid between cells of the same sample together and averages
#' the ranks of each cell as sender and as receiver. A pair with ligand and receptor NA gives the network score.
#' Cells and genes are found in the influence matrix by name, so the metadata rows may be in any order. The influence matrix is read
#' through its within-sample blocks (made here unless given) and the genes' rows and columns.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A data frame containing columns named cell, sample and group
#' @param ligands A vector of strings of ligand names, NA for the network score
#' @param receptors A vector of strings of receptor subunit names, NA for the network score
#' @param blocks The within-sample blocks from grid_blocks, or NULL to make them here
#' @returns A list of two matrices, sender and receiver, with one row per cell and one column per pair
#' @noRd
score_grids <- function(infMat, metadata, ligands, receptors, blocks=NULL){
  gi <- grid_inputs(infMat, metadata, ligands, receptors)
  if(is.null(blocks)){
    blocks <- grid_blocks(infMat, gi$g)
  }
  sr <- score_blocks(blocks, gi$A, gi$B, gi$g$sample, gi$lc, gi$rc)

  return(metadata_order(sr, gi$g))
}

#' Calculate the HyperCom score for each cell
#'
#' `score_hypercom` returns a data frame containing the HyperCom score for each cell
#'
#' This function takes a named influence matrix, the corresponding metadata (e.g. `res$infMat` and `res$metadata` of `HyperCom`), a
#' ligand and its receptor genes, and scores every cell of the hypergraph as sender and receiver of the pair. For each pair of cells
#' i and j of the same sample, a score says how strongly cell i reaches the ligand and the receptor reaches cell j, relative to how
#' strongly i reaches j: `F[i, ligand] * F[receptor, j] / F[i, j]` in the influence matrix F. These scores (of every sample) are
#' percent-ranked together; a cell's sender score is its mean rank as cell i, and its receiver score its mean rank as cell j.
#'
#' For a receptor complex the scores of its subunits are averaged. With neither ligand nor receptors, the network score ranks `F[i, j]`
#' itself: how strongly each cell reaches, and is reached by, the other cells of its sample. The scores are computed in this R process;
#' it stops at the start, giving the estimate, when they do not fit in memory.gb (see the Memory section of `prioritize_lr`).
#' @inheritParams prioritize_lr
#' @param metadata A data frame from generate_metadata (e.g. `res$metadata` of `HyperCom`), with columns cell, sample and group and
#' one row per node of infMat
#' @param ligand A string of the ligand gene, or NULL (the default) for the network score
#' @param receptors A vector of strings of the receptor genes (every subunit of a receptor complex), or NULL (the default) for the
#' network score
#' @param verbose A boolean of whether to report the memory the scoring needs as a message; FALSE by default
#' @returns A data frame with one row per cell of the hypergraph, sorted by hypercom (highest first), with columns:
#'
#' - cell: the cell name;
#' - sender, receiver: the cell's sender and receiver scores, from 0 to 1;
#' - dif: receiver - sender;
#' - avg: the mean of sender and receiver;
#' - max: the larger of the two;
#' - hypercom: one score per cell, negative for senders and positive for receivers: the sign of dif times the percent rank of its
#'   size among the cells, times avg.
#'
#' A cell without a score (no path of hyperedges to the pair's genes) gets NaN.
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
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL, memory.gb=NULL, verbose=FALSE){
  check_metadata(metadata)
  check_memory(memory.gb)
  check_flag(verbose, "verbose")
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
  plan_workers(metadata, 1, memory.gb, verbose)
  if(is.null(ligand)){
    sr <- score_grids(infMat, metadata, NA, NA)
    cols <- 1
  } else {
    sr <- score_grids(infMat, metadata, rep(ligand, complex.size), receptors)
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
  sr <- dplyr::arrange(sr, dplyr::desc(.data$hypercom))

  return(sr)
}
