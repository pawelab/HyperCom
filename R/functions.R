#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression. Missing counts are left out of the means.
#' A gene name used by several rows keeps its first row, with a warning.
#' @param counts A matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A matrix with row names being cells and column names being genes of only ligand and receptor genes with greater than 0 expression,
#' sparse if counts is sparse
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  genes <- unique(unlist(lapply(lrdb[, -1, drop=FALSE], as.character), use.names=FALSE))
  features <- intersect(rownames(counts), genes)
  dup.genes <- unique(rownames(counts)[duplicated(rownames(counts)) & rownames(counts) %in% features])
  if(length(dup.genes) > 0){
    warning("gene names used by several rows of counts, the first row is kept: ", paste(dup.genes, collapse=", "))
  }
  exp <- Matrix::t(counts[features, , drop=FALSE])
  genes.exp <- Matrix::colMeans(exp, na.rm=TRUE)
  genes.filter <- features[which(genes.exp > 0)]
  cc.counts <- exp[, genes.filter, drop=FALSE]

  if(!is.null(outfile)){
    saveRDS(cc.counts, outfile)
  }
  return(cc.counts)
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns dataframe of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptors interactions to produce
#' a table of ligand-receptors interactions in the counts matrix.
#' @param counts A filtered matrix with row names being cells and column names being genes
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' with only interactions present in the counts matrix
#' @export
generate_lrs_table <- function(counts, lrdb, outfile="1.lrs.csv"){
  genes <- as.matrix(lrdb[, 2:dim(lrdb)[2], drop=FALSE])
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

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns dataframe of hyperedges
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptors interactions to produce
#' a table of hyperedges: one per ligand-receptor interaction, one per gene, and one per expressed gene holding the cells that express it.
#' Cells are nodes named by the row names of counts, so they must be unique (prefix barcodes with their sample when merging samples).
#' @param counts A filtered matrix with row names being cells and column names being genes
#' @param lrs A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param outfile A string naming the output file
#' @returns A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @export
generate_hyperedges <- function(counts, lrs, outfile="2.hyperedges.csv"){
  dup.cells <- unique(rownames(counts)[duplicated(rownames(counts))])
  if(length(dup.cells) > 0){
    stop("cell names used by several cells of counts (prefix them with their sample): ", paste(utils::head(dup.cells, 10), collapse=", "))
  }

  lr.edges <- unique(tidyr::pivot_longer(lrs, -1, values_to="node", names_to=NULL, values_drop_na=TRUE))
  lr.edges$weight <- 1
  colnames(lr.edges) <- c("edge", "node", "weight")

  gene.edges <- data.frame(edge=colnames(counts), node=colnames(counts), weight=rep(1, ncol(counts)), check.names=FALSE)

  X <- methods::as(methods::as(methods::as(Matrix::t(counts), "dMatrix"), "generalMatrix"), "CsparseMatrix")
  cell.idx <- rep(seq_len(ncol(X)), diff(X@p))
  expressed <- X@x != 0 & !is.na(X@x)
  exp.edges <- data.frame(edge=rownames(X)[X@i[expressed] + 1], node=colnames(X)[cell.idx[expressed]], weight=rep(1, sum(expressed)),
                          check.names=FALSE)

  edges <- rbind(lr.edges, gene.edges, exp.edges)

  if(!is.null(outfile)){
    utils::write.csv(edges, outfile, row.names=FALSE)
  }
  return(edges)
}

#' Index the nodes and hyperedges of a hyperedge table
#'
#' `index_hyperedges` returns the node names and the node and hyperedge index of each row
#'
#' This function takes a dataframe of hyperedges and numbers its nodes by their names in the C locale and its hyperedges in order
#' of appearance, leaving out rows with a missing node or hyperedge.
#' @param hyperedges A dataframe of hyperedges with columns named edge and node
#' @returns A list with nodes (sorted node names), n.edges (number of hyperedges), node and edge (1-based indices of the rows kept)
#' @noRd
index_hyperedges <- function(hyperedges){
  if(is.null(hyperedges[["node"]]) | is.null(hyperedges[["edge"]])){
    stop("hyperedges must have columns named node and edge")
  }
  node <- hyperedges[["node"]]
  if(is.factor(node)){
    node <- as.character(node)
  }
  nodes <- as.character(sort(unique(node), method="radix"))
  edge <- as.character(hyperedges[["edge"]])
  edges <- unique(edge[!is.na(edge)])
  node.idx <- match(as.character(node), nodes)
  edge.idx <- match(edge, edges)
  found <- !is.na(node.idx) & !is.na(edge.idx)

  return(list(nodes=nodes, n.edges=length(edges), node=node.idx[found], edge=edge.idx[found]))
}

#' Generate adjacency matrix of hypergraph
#'
#' `generate_adjacency_matrix` returns the adjacency matrix of the hypergraph
#'
#' This function takes a dataframe of hyperedges and returns its adjacency matrix. Nodes are ordered by their names in the C locale,
#' so the order does not depend on the session's collation. Nodes that share no hyperedge with another node are dropped.
#' Only the upper triangle is computed and stored.
#' @param hyperedges A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @param outfile A string naming the output file
#' @returns An adjacency matrix where each entry is the number of shared hyperedges between nodes
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds"){
  h <- index_hyperedges(hyperedges)

  A <- adjacency_cpp(h$node - 1L, h$edge - 1L, length(h$nodes), h$n.edges)
  nodes <- h$nodes[A$keep]
  adj <- methods::new("dsCMatrix", i=A$i, p=A$p, x=A$x, Dim=rep(length(nodes), 2), Dimnames=list(node=nodes, node=nodes), uplo="U")

  if(!is.null(outfile)){
    saveRDS(adj, outfile)
  }
  return(adj)
}

#' Generate the incidence matrix behind an adjacency matrix
#'
#' `incidence_matrix` returns the node x hyperedge incidence counts of the hyperedges an adjacency matrix was made from
#'
#' This function takes a dataframe of hyperedges and its adjacency matrix and returns how often each node is in each hyperedge, for
#' the hyperedges with at least two distinct nodes (the only ones that add to the adjacency), nodes in the order of the adjacency.
#' It stops unless the incidence gives the adjacency (the incidence times its transpose with the diagonal set to zero): the same
#' nodes, and adj times two fixed integer vectors exactly equal to the products through the incidence (integers small enough for
#' every sum to be exact in double precision; no random numbers are drawn).
#' @param hyperedges A dataframe of hyperedges with columns named edge and node
#' @param adj An adjacency matrix with node names
#' @returns A sparse node x hyperedge matrix of counts
#' @noRd
incidence_matrix <- function(hyperedges, adj){
  h <- index_hyperedges(hyperedges)
  H <- Matrix::sparseMatrix(i=h$node, j=h$edge, x=1, dims=c(length(h$nodes), h$n.edges))
  H <- H[, diff(H@p) >= 2, drop=FALSE]
  kept <- which(Matrix::rowSums(H) > 0)
  idx <- match(rownames(adj), h$nodes[kept])
  if(length(idx) != length(kept) | anyNA(idx) | anyDuplicated(idx) > 0){
    stop("hyperedges do not match adj (other nodes): make adj with generate_adjacency_matrix(hyperedges)")
  }
  H <- H[kept[idx], , drop=FALSE]

  bound <- max(as.vector(H %*% Matrix::colSums(H)))
  top <- 2^max(1, min(30, floor(52 - log2(bound))))
  V <- floor(outer(seq_len(nrow(H)), c(0.6180339887498949, 0.7548776662466927)) %% 1 * top) + 1
  av <- as.matrix(adj %*% V)
  hv <- as.matrix(H %*% Matrix::crossprod(H, V)) - Matrix::rowSums(H^2) * V
  if(!identical(unname(av), unname(hv))){
    stop("hyperedges do not match adj (other entries): make adj with generate_adjacency_matrix(hyperedges)")
  }

  return(H)
}

#' Calculate influence matrix of hypergraph
#'
#' `generate_infMat` returns the influence matrix of the hypergraph
#'
#' This function takes an adjacency matrix and returns its influence matrix using heat diffusion, r (I - (1 - r) adj / D)^-1 with D
#' the row sums of adj. Three methods give the same matrix. "lu" solves the linear system with LAPACK in place as master's solve(),
#' bit for bit: two dense node x node matrices at once. "chol" is v1-fast's symmetric Cholesky inverse: one dense matrix and about
#' 3/8 of the arithmetic. "hypergraph" inverts through the hyperedges adj was made from (Woodbury identity): one dense matrix, the
#' result, and work growing with nodes x incidences instead of nodes cubed, so it is far faster when cells outnumber genes and
#' interactions (with at least as many hyperedges as nodes it gains nothing and computes "chol"). "chol" and "hypergraph" differ from
#' "lu" by about 1e-15 relative, which can break exact ties in the rankings, and stop like "lu" when the system is computationally
#' singular (r near 0).
#' @param adj A named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @param outfile A string naming the output file
#' @param method A string naming the method: "lu", "chol" or "hypergraph"
#' @param hyperedges A dataframe of the hyperedges adj was made from, needed by method "hypergraph"
#' @param threads An integer number of threads for method "hypergraph", NULL for all available cores ("lu" and "chol" run on the
#' threads of the BLAS)
#' @returns An influence matrix
#' @export
generate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds", method=c("lu", "chol", "hypergraph"), hyperedges=NULL, threads=NULL){
  method <- match.arg(method)
  if(is.data.frame(adj)){
    adj <- as.matrix(adj)
  }
  if(nrow(adj) == 0){
    stop("the adjacency matrix has no nodes: no ligand or receptor gene is expressed in any cell")
  }
  rs <- Matrix::rowSums(adj)
  if(!is.matrix(adj) & !methods::is(adj, "dsCMatrix")){
    adj <- methods::as(methods::as(methods::as(adj, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  }
  if(method != "lu"){
    if(!(is.numeric(r) && length(r) == 1 && isTRUE(r > 0 & r <= 1))){
      stop("r must be one number in (0, 1] for method \"", method, "\"")
    }
    if(!all(is.finite(rs))){
      stop("method \"", method, "\" needs finite adjacency weights")
    }
    if(!all(rs > 0)){
      stop("method \"", method, "\" needs every node of adj to have neighbours (positive row sums)")
    }
  }
  if(method == "hypergraph" & is.null(hyperedges)){
    stop("method \"hypergraph\" needs the hyperedges adj was made from")
  }

  if(method == "lu"){
    infMat <- solve_influence_cpp(adj, rs, r)
  } else if(method == "chol"){
    infMat <- solve_influence_chol_cpp(adj, rs, r)
  } else {
    H <- incidence_matrix(hyperedges, adj)
    if(ncol(H) < nrow(H)){
      infMat <- solve_influence_hypergraph_cpp(H@p, H@i, H@x, nrow(H), rs, r, get_threads(parallel=TRUE, threads=threads))
    } else {
      infMat <- solve_influence_chol_cpp(adj, rs, r)
    }
  }
  if(!is.null(rownames(adj))){
    dimnames(infMat) <- list(rownames(adj), rownames(adj))
  }

  if(!is.null(outfile)){
    saveRDS(infMat, outfile, compress=FALSE)
  }
  return(infMat)
}

#' Generate cleaned metadata table
#'
#' `generate_metadata` returns a minimal metadata dataframe
#'
#' This function takes an influence matrix and metadata table and creates a metadata data table containing the
#' node name, sample, and group. If group is unspecified, it is assumed all nodes are part of the same group.
#' Nodes of the influence matrix that are not cells of the metadata are genes: their sample, group and kept columns are "gene".
#' Kept columns named cell, sample or group are renamed with make.unique (e.g. sample.1); kept dates and other classed columns
#' other than factors become strings.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param outfile A string naming the output file
#' @returns A dataframe with columns cell, sample, group and the kept columns, one row per node of the influence matrix
#' @export
generate_metadata <- function(infMat, metadata, sample, group=NULL, keep=NULL, outfile="5.metadata.csv"){
  index <- data.frame(cell=row.names(infMat), check.names=FALSE)
  cells <- rownames(metadata)

  if(is.null(group)){
    group.vals <- "none"
  } else {
    group.vals <- as.character(metadata[, group])
  }
  sample.vals <- as.character(metadata[, sample])

  in.graph <- cells %in% index$cell
  if(!any(in.graph)){
    stop("no node of infMat is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  missing.sample <- sum(is.na(sample.vals[in.graph]))
  missing.group <- sum(is.na(rep(group.vals, length.out=nrow(metadata))[in.graph]))
  if(missing.sample > 0 | missing.group > 0){
    stop("sample and group must not be missing for cells in the graph: ", missing.sample, " cells miss a sample and ", missing.group,
         " a group")
  }
  if(any(sample.vals == "gene", na.rm=TRUE) | any(group.vals == "gene", na.rm=TRUE)){
    stop("\"gene\" labels gene nodes and cannot name a sample or group")
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

#' Choose the number of threads
#'
#' `get_threads` returns the number of threads used to score cell grids
#'
#' This function takes the parallel switch, an optional cluster and an optional thread count. An explicit thread count wins;
#' otherwise one thread without parallel, and with parallel the available cores (no more than the size of the cluster if one is given,
#' or than cl itself if it is a number). Never more than the available cores.
#' @param parallel A boolean of whether to run in parallel
#' @param cl A cluster object, a number of workers, or NULL
#' @param threads An integer number of threads or NULL
#' @returns An integer number of threads
#' @noRd
get_threads <- function(parallel=FALSE, cl=NULL, threads=NULL){
  cores <- parallelly::availableCores()
  if(!is.null(threads)){
    threads <- suppressWarnings(as.numeric(threads))
    if(length(threads) != 1 || is.na(threads)){
      stop("threads must be one number")
    }
    return(as.integer(max(1, min(floor(threads), cores))))
  }
  if(!parallel){
    return(1)
  }
  if(is.numeric(cl)){
    cores <- min(cores, max(1, floor(cl[1])))
  } else if(!is.null(cl)){
    cores <- min(cores, length(cl))
  }

  return(as.integer(cores))
}

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

#' Calculate sender and receiver scores of cell grids
#'
#' `score_grids` returns the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, and vectors of ligands and receptor subunits.
#' For each pair it percent-ranks all entries of the cell x cell transfer grid between cells of the same sample together and averages
#' the ranks of each cell as sender and as receiver. A pair with ligand and receptor NA gives the network score.
#' Cells and genes are found in the influence matrix by name, so the metadata rows may be in any order.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligands A vector of strings of ligand names, NA for the network score
#' @param receptors A vector of strings of receptor subunit names, NA for the network score
#' @param threads An integer number of threads
#' @returns A list of two matrices, sender and receiver, with one row per cell and one column per pair
#' @noRd
score_grids <- function(infMat, metadata, ligands, receptors, threads=1){
  ligands <- as.character(ligands)
  receptors <- as.character(receptors)
  nodes <- rownames(infMat)
  check_genes(c(ligands, receptors), nodes)
  if(any(is.na(ligands) != is.na(receptors))){
    stop("each pair needs both a ligand and a receptor, or neither for the network score")
  }
  if(anyNA(metadata$sample)){
    stop("metadata has missing samples: make it with generate_metadata")
  }
  cells <- metadata$sample != "gene"
  if(!any(cells)){
    stop("metadata has no cells: make it with generate_metadata from metadata whose row names are the cell names")
  }
  cell.idx <- match(metadata$cell[cells], nodes) - 1
  if(anyNA(cell.idx) | anyDuplicated(cell.idx) > 0){
    stop("metadata cells must be distinct nodes of infMat: make it with generate_metadata")
  }

  samples <- as.integer(factor(metadata$sample[cells]))
  lig <- match(ligands, nodes) - 1
  rec <- match(receptors, nodes) - 1
  lig[is.na(lig)] <- -1
  rec[is.na(rec)] <- -1
  o <- order(cell.idx)

  sr <- score_grids_cpp(infMat, cell.idx[o], samples[o], lig, rec, threads)
  if(is.unsorted(cell.idx)){
    sr <- list(sender=sr$sender[order(o), , drop=FALSE], receiver=sr$receiver[order(o), , drop=FALSE])
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
#' to give each cell a HyperCom score. If both the ligand and receptor are specified this is the LR Pair Score
#' otherwise the score is the Network Score. If a receptor complex is specified, the average of the score vectors is used.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligand A string of the ligand name
#' @param receptors A vector containing string elements of the receptor names
#' @returns A dataframe containing the sender score, receiver score, avg, dif, and HyperCom score
#' @export
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL){
  complex.size <- length(receptors)
  if(!is.null(ligand) & !is.null(receptors) & complex.size == 0){
    stop("receptors must name at least one gene")
  }
  if(is.null(ligand) | is.null(receptors)){
    sr <- score_grids(infMat, metadata, NA, NA)
    cols <- rep(1, max(complex.size, 1))
  } else {
    sr <- score_grids(infMat, metadata, rep(ligand, complex.size), receptors)
    cols <- seq_len(complex.size)
  }

  sr <- average_scores(sr, cols)
  sr <- data.frame(cell=metadata$cell[metadata$sample != "gene"], sender=sr$sender, receiver=sr$receiver)
  if(all(ligand == receptors) & !is.null(receptors)){
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

#' Generate a background set of random gene interactions
#'
#' `generate_background_lrs` returns a dataframe of random gene interactions.
#'
#' This generates a table of random gene interactions where the number of genes in an interaction is equal to the complex size + 1.
#' @param lrs A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param count An integer of the number of interactions to generate
#' @return A dataframe of random gene interactions
#' @noRd
generate_background_lrs <- function(lrs, metadata, count=1000){
  background <- data.frame(interaction_name=rep("background", count), check.names=FALSE)
  nodes <- metadata$cell[metadata$sample == "gene"]
  for (i in 2:dim(lrs)[2]) {
    col <- sample(nodes, count, replace=TRUE)
    background <- cbind(background, col)
  }

  colnames(background) <- colnames(lrs)

  return(background)
}

#' Perform a permutation test
#'
#' `permutation_test` returns a numeric p-value
#'
#' This function takes a numeric and list of numeric vector and returns how many values in the vector are less than the numeric.
#' @param x A numeric
#' @param background A numeric vector
#' @return A numeric p-value
#' @noRd
permutation_test <- function(x, background){
  sum(abs(background) >= abs(x)) / length(background)
}

#' Calculate the Gini impurity of a set of labels
#'
#' `gini_impurity` returns the Gini impurity of a vector of labels
#'
#' This function takes a vector of labels and returns the probability that a label drawn at random is mislabelled
#' when labelled according to the distribution of labels (as mltools::gini_impurity does for labels, not counts).
#' @param vals A vector of labels
#' @returns A numeric Gini impurity
#' @noRd
gini_impurity <- function(vals){
  if(is.factor(vals)){
    counts <- tabulate(vals, nbins=nlevels(vals))
  } else if(is.logical(vals)){
    counts <- c(sum(!vals, na.rm=TRUE), sum(vals, na.rm=TRUE))
    counts <- counts[counts > 0]
  } else {
    counts <- table(vals)
  }
  total <- sum(counts)
  return(sum((counts/total)*(1-counts/total)))
}

#' Group ligand-receptor pairs into batches of whole ligands
#'
#' `batch_pairs` returns the rows of each batch of ligand-receptor pairs
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and splits the pairs into batches of whole ligands
#' with about size distinct ligand and receptor subunit grids each. A grid belongs to one ligand, so no grid is needed by two batches.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @param size An integer number of grids per batch
#' @returns A list of integer vectors of rows, one per batch, rows in their original order
#' @noRd
batch_pairs <- function(ligands, receptors, size){
  grids <- unique(data.frame(ligand=rep(ligands, ncol(receptors)), receptor=as.vector(receptors)))
  grids <- grids[!is.na(grids$receptor), ]
  ligand.order <- unique(ligands)
  per.ligand <- tabulate(match(grids$ligand, ligand.order), nbins=length(ligand.order))

  batch <- integer(length(ligand.order))
  b <- 1
  filled <- 0
  for(k in seq_along(ligand.order)){
    if(filled > 0 & filled + per.ligand[k] > size){
      b <- b + 1
      filled <- 0
    }
    batch[k] <- b
    filled <- filled + per.ligand[k]
  }

  return(split(seq_along(ligands), batch[match(ligands, ligand.order)]))
}

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns the terms of the priority score of each ligand-receptor pair
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to score ligand-receptor pairs. Pairs are scored in batches of whole ligands: each ligand and receptor subunit grid is scored once,
#' whatever the number of pairs using it, and only one batch of cell scores is in memory at a time. The gene x cell counts of a
#' symmetric sparse adjacency are read from its stored triangle, without expanding the whole matrix.
#' @param lrs A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param threads An integer number of threads
#' @return A list with terms, a matrix with one row per pair and columns frac, pct.min, max, balance and purity, and title, the pair names
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, threads=1){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells.idx <- metadata$sample != "gene"
  cells <- metadata$cell[cells.idx]
  groups <- metadata$group[cells.idx]

  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, 3:dim(lrs)[2], drop=FALSE])
  genes <- unique(c(ligands, receptors[!is.na(receptors)]))
  check_genes(genes, rownames(infMat))

  if(methods::is(adj, "dsCMatrix") & !any(genes %in% cells)){
    U <- methods::new("dgCMatrix", i=adj@i, p=adj@p, x=adj@x, Dim=adj@Dim, Dimnames=dimnames(adj))
    gene.counts <- as.matrix(Matrix::t(U[genes, cells, drop=FALSE]) + U[cells, genes, drop=FALSE])
    rm(U)
  } else {
    gene.counts <- as.matrix(Matrix::t(adj[genes, cells, drop=FALSE]))
  }
  gene.frac <- colSums(gene.counts) / length(cells)
  expressed <- gene.counts > 0
  rm(gene.counts)

  terms <- matrix(NA_real_, nrow(lrs), 5, dimnames=list(NULL, c("frac", "pct.min", "max", "balance", "purity")))
  title <- character(nrow(lrs))
  for(rows in batch_pairs(ligands, receptors, max(256, 32 * threads))){
    pairs <- unique(data.frame(ligand=rep(ligands[rows], ncol(receptors)), receptor=as.vector(receptors[rows, , drop=FALSE])))
    pairs <- pairs[!is.na(pairs$receptor), ]
    keys <- paste(pairs$ligand, pairs$receptor, sep="\t")
    sr <- score_grids(infMat, metadata, pairs$ligand, pairs$receptor, threads)

    for(i in rows){
      l <- ligands[i]
      rs <- stats::na.omit(receptors[i, ])
      cols <- match(paste(l, rs, sep="\t"), keys)

      pair.sr <- average_scores(sr, cols)
      ss <- pair.sr$sender
      rr <- pair.sr$receiver
      if(all(l == rs)){
        type <- rep(FALSE, length(ss))
      } else {
        type <- (rr - ss) > 0
      }

      pair.genes <- c(l, rs)
      present <- expressed[, pair.genes[1]]
      for(g in pair.genes[-1]){
        present <- present | expressed[, g]
      }
      pct.min <- min(gene.frac[pair.genes])
      frac <- sum(present) / length(cells)
      max <- mean(pmax(ss, rr)[present])

      if(gini.max > 0){
        purity <- gini.max - gini_impurity(groups[present])
      }else{purity <- 1}

      balance <- gini_impurity(type[present]) / 0.5

      terms[i, ] <- c(frac, pct.min, max, balance, purity)
      title[i] <- paste(c(l, rs), collapse="_")
    }
    rm(sr)
  }

  return(list(terms=terms, title=title))
}

#' Combine priority terms into priority scores
#'
#' `score_priority` returns the priority score of each ligand-receptor pair for one weight
#'
#' This function takes the priority terms of calculate_priority and returns their weighted geometric mean.
#' @param priority.terms A list from calculate_priority
#' @param weight a value between 0 to 1 for the importance of number of cells involved in an interaction
#' @param verbose A boolean of whether to print the score of each pair
#' @return A numeric vector of priority scores
#' @noRd
score_priority <- function(priority.terms, weight=1, verbose=TRUE){
  w <- c(weight, rep(1, 4))
  terms <- priority.terms$terms
  priority <- numeric(nrow(terms))
  for(i in seq_len(nrow(terms))){
    priority[i] <- exp(stats::weighted.mean(log(terms[i, ]), w))
  }

  if(verbose){
    for(i in seq_along(priority)){
      print(paste0(priority.terms$title[i], ",", priority[i]))
    }
  }
  return(priority)
}

#' Order ligand receptor pairs based on signalling strength and specificity
#'
#' `prioritize_lr` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to prioritize ligand-receptor pairs. Cell grids are scored in C++ threads, once for all weights.
#' @param lrs A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix.
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel: threads on all available cores (`parallelly::availableCores()`)
#' @param cl A cluster object, kept for compatibility: it is not used or stopped, only its size caps the number of threads
#' @param outdir A string naming the output directory
#' @param threads An integer number of threads, overriding parallel and cl (at most the available cores)
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state
#' @param verbose A boolean of whether to print each weight and the score of each pair
#' @returns A dataframe of ligand-receptor pairs ordered by priority score for each weight
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="6.priority",
                          threads=NULL, seed=NULL, verbose=TRUE){
  if(anyDuplicated(as.character(weights)) > 0){
    stop("weights must differ: each names a result and an output file")
  }
  if(significance){
    if(!is.null(seed)){
      set.seed(seed)
    }
    background.lrs <- generate_background_lrs(lrs, metadata)
  }
  threads <- get_threads(parallel, cl, threads)

  if(significance){
    all.terms <- calculate_priority(rbind(lrs, background.lrs), adj, infMat, metadata, threads)
    fg <- seq_len(nrow(lrs))
    priority.terms <- list(terms=all.terms$terms[fg, , drop=FALSE], title=all.terms$title[fg])
    background.terms <- list(terms=all.terms$terms[-fg, , drop=FALSE], title=all.terms$title[-fg])
  } else {
    priority.terms <- calculate_priority(lrs, adj, infMat, metadata, threads)
  }

  lrs.list <- list()

  for(i in seq(1, length(weights))){
    weight <- weights[i]
    if(verbose){
      print(weight)
    }
    priority <- score_priority(priority.terms, weight, verbose)

    current.lrs <- lrs

    if(significance){
      background.priority <- score_priority(background.terms, weight, verbose)
    }

    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      current.lrs$p <- as.numeric(lapply(current.lrs$priority, FUN=permutation_test, background=background.priority))
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method="BH")
    }

    if(!is.null(outdir)){
      dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)
      outfile <- paste0(outdir, "/priority.", weight, ".csv")
      utils::write.csv(current.lrs, outfile, row.names=FALSE)
    }

    lrs.list[[i]] <- current.lrs
  }

  names(lrs.list) <- weights

  return(lrs.list)
}

#' Run the HyperCom pipeline
#'
#' `HyperCom` Runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`,
#' `generate_metadata`, and `prioritize_lr`, writing every step to outdir.
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param counts A matrix with row names being genes and column names being cells
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel: threads on all available cores
#' @param cl A cluster object, kept for compatibility: only its size caps the number of threads
#' @param outdir output directory of pipeline run
#' @param threads An integer number of threads, overriding parallel and cl (at most the available cores)
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state
#' @param verbose A boolean of whether to print each weight and the score of each pair
#' @param infmat.method A string naming the method of `generate_infMat`: "lu" (master's numbers bit for bit), "chol" or "hypergraph"
#' (faster, one dense matrix instead of two, about 1e-15 apart)
#' @returns 0 if run successfully
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL,
                     outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="lu"){
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  no.metadata <- setdiff(colnames(counts), rownames(metadata))
  if(length(no.metadata) > 0){
    warning(length(no.metadata), " cells of counts are not in metadata; those that are nodes of the graph count as genes")
  }

  counts <- filter_counts(counts, lrdb, outfile=file.path(outdir, "0.counts.rds"))
  lrs <- generate_lrs_table(counts, lrdb, outfile=file.path(outdir, "1.lrs.csv"))
  hyperedges <- generate_hyperedges(counts, lrs, outfile=file.path(outdir, "2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=file.path(outdir, "3.adj.rds"))
  infMat <- generate_infMat(adj, outfile=file.path(outdir, "4.infMat.rds"), method=infmat.method, hyperedges=hyperedges, threads=threads)
  metadata <- generate_metadata(infMat, metadata, sample, group, keep, outfile=file.path(outdir, "5.metadata.csv"))
  priority <- prioritize_lr(lrs, adj, infMat, metadata, weights, significance, parallel, cl, outdir=file.path(outdir, "6.priority"),
                            threads=threads, seed=seed, verbose=verbose)

  return(0)
}

#' Get the counts matrix of a Seurat assay layer
#'
#' `get_seurat_layer` returns the matrix with dimnames for given assay and layer
#'
#' This function takes a Seurat object and returns one layer of one assay, joining split layers of v5 assays.
#' @param seurat A seurat object
#' @param assay A string naming the assay to get data from
#' @param layer A string naming the layer to get data from
#' @returns A matrix with row names being genes and column names being cells
#' @noRd
get_seurat_layer <- function(seurat, assay="RNA", layer="counts"){
  if(!requireNamespace("SeuratObject", quietly=TRUE)){
    stop("reading a Seurat object needs the SeuratObject package")
  }
  a <- seurat[[assay]]

  if(!inherits(a, "Assay5")){
    m <- SeuratObject::GetAssayData(a, layer=layer)
    if(!length(m)){
      stop("layer '", layer, "' is empty or absent in assay '", assay, "'")
    }
    return(m)
  }

  lyrs <- SeuratObject::Layers(a, search=layer)
  if(!length(lyrs)){
    stop("no layer matching '", layer, "' in assay '", assay, "'; available: ", paste(SeuratObject::Layers(a), collapse=", "))
  }
  if(length(lyrs) == 1){
    return(SeuratObject::LayerData(a, layer=lyrs))
  }

  return(SeuratObject::LayerData(SeuratObject::JoinLayers(a), layer=layer))
}

#' Run the HyperCom pipeline on a seurat object
#'
#' `HyperCom_seurat` Runs the HyperCom pipeline on a seurat object
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`,
#' `generate_metadata`, and `prioritize_lr` on one assay layer and the cell metadata of a Seurat object. Needs the SeuratObject package.
#' @param seurat A seurat object
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param sample A string that is the column name for sample in the provided metadata
#' @param assay assay to get data from
#' @param layer layer to get data from
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel: threads on all available cores
#' @param cl A cluster object, kept for compatibility: only its size caps the number of threads
#' @param outdir output directory of pipeline run
#' @param threads An integer number of threads, overriding parallel and cl (at most the available cores)
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state
#' @param verbose A boolean of whether to print each weight and the score of each pair
#' @param infmat.method A string naming the method of `generate_infMat`: "lu" (master's numbers bit for bit), "chol" or "hypergraph"
#' (faster, one dense matrix instead of two, about 1e-15 apart)
#' @returns 0 if run successfully
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, weights=c(1), significance=FALSE,
                            parallel=FALSE, cl=NULL, outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="lu"){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat[[]]

  status <- HyperCom(lrdb, counts, metadata, sample, group, keep, weights, significance, parallel, cl, outdir, threads=threads, seed=seed,
                     verbose=verbose, infmat.method=infmat.method)

  return(status)
}
