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

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns a dataframe of hyperedges
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptor interactions to produce
#' a table of hyperedges: one per ligand-receptor interaction (its genes) and one per gene (the gene and the cells that express it).
#' Cells are nodes named by the row names of counts, so they must be unique (prefix barcodes with their sample when merging samples).
#' @inheritParams generate_lrs_table
#' @param lrs A dataframe of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @returns A dataframe with one row per node of each hyperedge: columns edge (the hyperedge, an interaction or a gene), node (a gene or
#' cell in it) and weight (1)
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
#' @returns A list with nodes (sorted node names), edges (hyperedge names in order of appearance), node and edge (1-based indices
#' of the rows kept)
#' @noRd
index_hyperedges <- function(hyperedges){
  if(is.null(hyperedges[["node"]]) || is.null(hyperedges[["edge"]])){
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

  return(list(nodes=nodes, edges=edges, node=node.idx[found], edge=edge.idx[found]))
}

#' Count how often each node is in each hyperedge
#'
#' `hyperedge_incidence` returns the node x hyperedge incidence counts of the nodes that share a hyperedge with another node
#'
#' This function takes indexed hyperedges and returns how often each node is in each hyperedge, for the hyperedges with at least two
#' distinct nodes (the only ones that add to the adjacency) and the nodes in one of them (the nodes of the adjacency), nodes in the
#' order of index_hyperedges.
#' @param h A list from index_hyperedges
#' @returns A sparse node x hyperedge matrix of counts with node and hyperedge names
#' @noRd
hyperedge_incidence <- function(h){
  H <- Matrix::sparseMatrix(i=h$node, j=h$edge, x=1, dims=c(length(h$nodes), length(h$edges)), dimnames=list(node=h$nodes, edge=h$edges))
  H <- H[, diff(H@p) >= 2, drop=FALSE]
  H <- H[Matrix::rowSums(H) > 0, , drop=FALSE]

  return(H)
}

#' Generate adjacency matrix of hypergraph
#'
#' `generate_adjacency_matrix` returns the adjacency matrix of the hypergraph
#'
#' This function takes a dataframe of hyperedges and returns its adjacency matrix: each entry is the number of hyperedges two nodes
#' share, and the diagonal is zero. Nodes are ordered by their names in the C locale, so the order does not depend on the session's
#' collation. Nodes that share no hyperedge with another node are dropped. The full node x node matrix is not stored: the result
#' records which nodes are in which hyperedge, so it takes memory in proportion to the hyperedge table instead of the number of nodes
#' squared. It works like a matrix for dim, dimnames and rownames; `adj[i, j]` gives those entries and as.matrix the whole matrix.
#' @param hyperedges A dataframe of hyperedges from generate_hyperedges, with columns edge (the hyperedge) and node (a node in it)
#' @inheritParams filter_counts
#' @returns A hypercom_adjacency: a list with H, the sparse node x hyperedge matrix of how often each node is in each hyperedge
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds"){
  adj <- list(H=hyperedge_incidence(index_hyperedges(hyperedges)))
  class(adj) <- "hypercom_adjacency"

  if(!is.null(outfile)){
    saveRDS(adj, outfile)
  }
  return(adj)
}

#' Expand an adjacency into the sparse matrix of HyperCom 0.0.0.9000
#'
#' `adjacency_sparse` returns the adjacency matrix of a hypercom_adjacency as the dsCMatrix of HyperCom 0.0.0.9000
#'
#' This function takes an adjacency from generate_adjacency_matrix and returns the upper triangle of its symmetric matrix with the
#' explicit zero diagonal, the object generate_adjacency_matrix made in HyperCom 0.0.0.9000. It takes memory in proportion to the node pairs that
#' share a hyperedge.
#' @param adj A hypercom_adjacency
#' @returns A dsCMatrix with node names
#' @noRd
adjacency_sparse <- function(adj){
  H <- adj$H
  node <- rep(H@i, H@x)
  edge <- rep(rep(seq_len(ncol(H)) - 1L, diff(H@p)), H@x)
  A <- adjacency_cpp(node, edge, nrow(H), ncol(H))
  nodes <- rownames(H)[A$keep]

  return(methods::new("dsCMatrix", i=A$i, p=A$p, x=A$x, Dim=rep(length(nodes), 2L), Dimnames=list(node=nodes, node=nodes), uplo="U"))
}

#' Find nodes by name, position or logical
#'
#' `node_index` returns the positions of indexed nodes
#'
#' This function takes an index as for a matrix (names, positive or negative positions, logical; a factor is read by name) and the
#' node names and returns the positions, stopping on a name or position that is not a node, on NA and on a logical index longer than
#' the nodes.
#' @param i A vector indexing nodes
#' @param nodes A vector of strings of the node names
#' @returns An integer vector of positions
#' @noRd
node_index <- function(i, nodes){
  if(anyNA(i)){
    stop("indices must not be NA")
  }
  if(is.logical(i) && length(i) > length(nodes)){
    stop("(subscript) logical subscript too long")
  }
  if(is.factor(i)){
    i <- as.character(i)
  }
  if(is.character(i)){
    idx <- match(i, nodes)
  } else {
    idx <- seq_along(nodes)[i]
  }
  if(anyNA(idx)){
    stop("subscript out of bounds")
  }

  return(idx)
}

#' Adjacency matrix of a hypergraph, stored as its incidence
#'
#' `hypercom_adjacency` objects are made by generate_adjacency_matrix
#'
#' A hypercom_adjacency keeps the node x hyperedge incidence H of a hypergraph; its adjacency matrix is H H' with the diagonal set
#' to zero. dim and dimnames give those of the matrix, `adj[i, j]` makes the indexed entries as a dense matrix (indices as for a
#' matrix: names, positions or logical; factors are read by name, not by their codes) and as.matrix the whole dense matrix (nodes
#' squared doubles).
#' @param x A hypercom_adjacency
#' @param i,j Indices of rows and columns
#' @param drop A boolean of whether to drop dimensions of extent one
#' @param ... Ignored
#' @returns dim: the number of nodes twice; dimnames: the node names twice; `[` and as.matrix: a dense matrix; print: x, invisibly
#' @name hypercom_adjacency
NULL

#' @rdname hypercom_adjacency
#' @export
dim.hypercom_adjacency <- function(x){
  return(rep(nrow(x$H), 2L))
}

#' @rdname hypercom_adjacency
#' @export
dimnames.hypercom_adjacency <- function(x){
  return(list(node=rownames(x$H), node=rownames(x$H)))
}

#' @rdname hypercom_adjacency
#' @export
`[.hypercom_adjacency` <- function(x, i, j, drop=TRUE){
  narg <- nargs() - !missing(drop)
  if(narg < 3){
    stop("index an adjacency as adj[i, j]")
  }
  nodes <- rownames(x$H)
  rows <- if(missing(i)) seq_along(nodes) else node_index(i, nodes)
  cols <- if(missing(j)) seq_along(nodes) else node_index(j, nodes)
  A <- as.matrix(Matrix::tcrossprod(x$H[rows, , drop=FALSE], x$H[cols, , drop=FALSE]))
  if(anyDuplicated(cols) == 0){
    same <- match(rows, cols)
    A[cbind(which(!is.na(same)), same[!is.na(same)])] <- 0
  } else {
    for(a in which(rows %in% cols)){
      A[a, cols == rows[a]] <- 0
    }
  }
  dimnames(A) <- list(node=nodes[rows], node=nodes[cols])
  if(drop && any(dim(A) == 1)){
    A <- A[, , drop=TRUE]
  }

  return(A)
}

#' @rdname hypercom_adjacency
#' @export
as.matrix.hypercom_adjacency <- function(x, ...){
  return(as.matrix(adjacency_sparse(x)))
}

#' @rdname hypercom_adjacency
#' @export
print.hypercom_adjacency <- function(x, ...){
  cat("HyperCom adjacency matrix of ", nrow(x$H), " nodes, stored as their incidence in ", ncol(x$H), " hyperedges\n", sep="")

  return(invisible(x))
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
#' @returns A sparse node x hyperedge matrix of counts with node and hyperedge names
#' @noRd
incidence_matrix <- function(hyperedges, adj){
  H <- hyperedge_incidence(index_hyperedges(hyperedges))
  idx <- match(rownames(adj), rownames(H))
  if(length(idx) != nrow(H) || anyNA(idx) || anyDuplicated(idx) > 0){
    stop("hyperedges do not match adj (other nodes): make adj with generate_adjacency_matrix(hyperedges)")
  }
  H <- H[idx, , drop=FALSE]

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
#' This function takes an adjacency matrix and returns its influence matrix by heat diffusion (a random walk with restart),
#' r (I - (1 - r) adj / D)^-1 with D the row sums of adj. Three methods give the same matrix up to rounding (about 1e-15 relative,
#' which can reorder exactly tied scores):
#'
#' - "woodbury" (the default) works through the hyperedges, so no node x node matrix is made: memory grows with the number of nodes
#'   instead of its square. The result works like a matrix for dim, dimnames and rownames; `infMat[i, j]` gives those entries and
#'   as.matrix the whole matrix. An adjacency from generate_adjacency_matrix holds its hyperedges; an adjacency matrix needs
#'   `hyperedges`, and without them (and without a method given) gets "lu". Small data with at least as many hyperedges as nodes get
#'   the dense matrix of "chol", which is then no bigger.
#' - "lu" makes the dense matrix with exactly the numbers of HyperCom 0.0.0.9000; it holds two node x node matrices at once.
#' - "chol" makes the dense matrix faster and with one node x node matrix.
#'
#' All methods stop when the matrix cannot be inverted (r near 0).
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @param method A string naming the method: "woodbury" (the default; "lu" for an adjacency matrix given without hyperedges), "lu"
#' or "chol"
#' @param hyperedges A dataframe of the hyperedges an adjacency matrix was made from, needed by method "woodbury" with an adjacency
#' matrix; ignored with an adjacency from generate_adjacency_matrix, which holds its own
#' @param threads An integer number of threads for method "woodbury", NULL for one per physical core; not used by "lu" and "chol"
#' @inheritParams filter_counts
#' @returns An influence matrix with node names: a hypercom_influence (method "woodbury", unless there are at least as many hyperedges
#' as nodes) or a dense matrix
#' @export
generate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds", method=c("woodbury", "lu", "chol"), hyperedges=NULL, threads=NULL){
  # woodbury: adj = H H' with a zero diagonal (H the node x hyperedge incidence), so the Woodbury identity inverts it through a
  # hyperedge x hyperedge matrix and the result stays factored (influence_factors). lu: LAPACK dgesv in place, bit for bit the
  # solve() of HyperCom 0.0.0.9000. chol: LAPACK dpotrf + dpotri, as HyperCom 0.1.0's method="chol".
  graph <- inherits(adj, "hypercom_adjacency")
  if(missing(method) && !graph && is.null(hyperedges)){
    method <- "lu"
  }
  method <- match_method(method)
  if(is.data.frame(adj)){
    adj <- as.matrix(adj)
  }
  if(nrow(adj) == 0){
    stop("the adjacency matrix has no nodes: no ligand or receptor gene is expressed in any cell")
  }
  if(graph && method != "woodbury"){
    adj <- adjacency_sparse(adj)
    graph <- FALSE
  }
  if(graph){
    H <- adj$H
    rs <- as.vector(H %*% Matrix::colSums(H) - Matrix::rowSums(H^2))
  } else {
    rs <- Matrix::rowSums(adj)
    if(!is.matrix(adj) && !methods::is(adj, "dsCMatrix")){
      adj <- methods::as(methods::as(methods::as(adj, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    }
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
  if(method == "woodbury" && !graph && is.null(hyperedges)){
    stop("method \"woodbury\" needs the hyperedges adj was made from")
  }

  if(method == "lu"){
    infMat <- solve_influence_cpp(adj, rs, r)
  } else if(method == "chol"){
    infMat <- solve_influence_chol_cpp(adj, rs, r)
  } else {
    if(!graph){
      H <- incidence_matrix(hyperedges, adj)
    }
    if(ncol(H) < nrow(H)){
      infMat <- influence_factors(H, rs, r, threads)
    } else {
      if(graph){
        adj <- adjacency_sparse(adj)
      }
      infMat <- solve_influence_chol_cpp(adj, rs, r)
    }
  }
  if(is.matrix(infMat) && !is.null(rownames(adj))){
    dimnames(infMat) <- list(rownames(adj), rownames(adj))
  }

  if(!is.null(outfile)){
    saveRDS(infMat, outfile, compress=FALSE)
  }
  return(infMat)
}

#' Match the name of an influence matrix method
#'
#' `match_method` returns the full name of an influence matrix method
#'
#' This function takes a method name as generate_infMat's method argument and returns it matched to "woodbury", "lu" or "chol"
#' (a unique abbreviation is completed; any other name stops).
#' @param method A string naming the method, or the vector of all three for the first
#' @returns A string naming the method
#' @noRd
match_method <- function(method){
  return(match.arg(method, c("woodbury", "lu", "chol")))
}

#' Factor the influence matrix through the hyperedges
#'
#' `influence_factors` returns the influence matrix of a hypergraph as a hypercom_influence
#'
#' This function takes the node x hyperedge incidence, the row sums of the adjacency and the restart probability and returns the
#' factors of the influence matrix F = r Delta^-1 D + r (1 - r) Y Z (Y = H / Delta, Z = C^-1 Y' D): H, the row sums D, Delta and the
#' hyperedge x node matrix Z.
#' @param H A sparse node x hyperedge matrix of counts with node and hyperedge names
#' @param rs A numeric vector of the row sums of the adjacency
#' @param r A double between 0 and 1 representing the restart probability
#' @param threads An integer number of threads, NULL for one per physical core
#' @returns A hypercom_influence: a list with H, rs, delta, Z and r
#' @noRd
influence_factors <- function(H, rs, r, threads=NULL){
  rs <- as.numeric(rs)
  f <- influence_factors_cpp(H@p, H@i, H@x, nrow(H), rs, r, get_threads(parallel=TRUE, threads=threads))
  Z <- f$Z
  dimnames(Z) <- list(edge=colnames(H), node=rownames(H))
  infMat <- list(H=H, rs=rs, delta=f$delta, Z=Z, r=r)
  class(infMat) <- "hypercom_influence"

  return(infMat)
}

#' Make a block of a factored influence matrix
#'
#' `influence_slice` returns the entries of a hypercom_influence at given rows and columns
#'
#' This function takes a hypercom_influence and positions of rows and columns and returns those entries as a dense matrix, each
#' entry the same as in the whole matrix at any thread count.
#' @param infMat A hypercom_influence
#' @param rows An integer vector of row positions
#' @param cols An integer vector of column positions
#' @param threads An integer number of threads
#' @returns A dense matrix without names
#' @noRd
influence_slice <- function(infMat, rows, cols, threads=1){
  H <- infMat$H
  block <- influence_block_cpp(H@p, H@i, H@x, nrow(H), infMat$rs, infMat$delta, infMat$Z, infMat$r, as.integer(rows) - 1L,
                               as.integer(cols) - 1L, get_threads(threads=threads))

  return(block)
}

#' Make the within-sample blocks of a factored influence matrix
#'
#' `influence_blocks` returns every within-sample cell x cell block of a hypercom_influence, as the grid kernel reads them
#'
#' This function takes a hypercom_influence and the cells from grid_cells and returns the blocks one after the other (samples
#' ascending, cells of a sample in the order given, each block by column): the sum of squared sample sizes of doubles.
#' @param infMat A hypercom_influence
#' @param g A list from grid_cells
#' @param threads An integer number of threads
#' @returns A numeric vector
#' @noRd
influence_blocks <- function(infMat, g, threads=1){
  if(sum(as.numeric(tabulate(g$sample))^2) >= 4294967295){
    stop("too many cells per sample: the sum of squared sample sizes must stay below 2^32")
  }
  H <- infMat$H

  return(influence_blocks_cpp(H@p, H@i, H@x, nrow(H), infMat$rs, infMat$delta, infMat$Z, infMat$r, g$idx, g$sample, threads))
}

#' Influence matrix of a hypergraph, stored as its factors
#'
#' `hypercom_influence` objects are made by generate_infMat with method "woodbury"
#'
#' A hypercom_influence keeps the factors of the influence matrix F = r Delta^-1 D + r (1 - r) Y Z of a hypergraph (H, the node x
#' hyperedge incidence; rs, the row sums D of the adjacency; delta, Delta = D + (1 - r) diag(H H'); Z, hyperedge x node; Y = H /
#' Delta) instead of F. dim and dimnames give those of F, `infMat[i, j]` makes the indexed entries as a dense matrix (indices as for a
#' matrix: names, positions or logical; factors are read by name, not by their codes) on one thread and as.matrix the whole dense
#' matrix (nodes squared doubles).
#' @param x A hypercom_influence
#' @param i,j Indices of rows and columns
#' @param drop A boolean of whether to drop dimensions of extent one
#' @param threads An integer number of threads for as.matrix
#' @param ... Ignored
#' @returns dim: the number of nodes twice; dimnames: the node names twice; `[` and as.matrix: a dense matrix; print: x, invisibly
#' @name hypercom_influence
NULL

#' @rdname hypercom_influence
#' @export
dim.hypercom_influence <- function(x){
  return(rep(nrow(x$H), 2L))
}

#' @rdname hypercom_influence
#' @export
dimnames.hypercom_influence <- function(x){
  return(list(rownames(x$H), rownames(x$H)))
}

#' @rdname hypercom_influence
#' @export
`[.hypercom_influence` <- function(x, i, j, drop=TRUE){
  narg <- nargs() - !missing(drop)
  if(narg < 3){
    stop("index an influence matrix as infMat[i, j]")
  }
  nodes <- rownames(x$H)
  rows <- if(missing(i)) seq_along(nodes) else node_index(i, nodes)
  cols <- if(missing(j)) seq_along(nodes) else node_index(j, nodes)
  block <- influence_slice(x, rows, cols)
  dimnames(block) <- list(nodes[rows], nodes[cols])
  if(drop && any(dim(block) == 1)){
    block <- block[, , drop=TRUE]
  }

  return(block)
}

#' @rdname hypercom_influence
#' @export
as.matrix.hypercom_influence <- function(x, threads=1, ...){
  nodes <- rownames(x$H)
  block <- influence_slice(x, seq_along(nodes), seq_along(nodes), threads)
  dimnames(block) <- list(nodes, nodes)

  return(block)
}

#' @rdname hypercom_influence
#' @export
print.hypercom_influence <- function(x, ...){
  cat("HyperCom influence matrix of ", nrow(x$H), " nodes (r = ", x$r, "), stored as its factors through ", ncol(x$H), " hyperedges\n",
      sep="")

  return(invisible(x))
}

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

#' Choose the number of threads
#'
#' `get_threads` returns the number of threads used to score cell grids
#'
#' This function takes the parallel switch and an optional thread count. An explicit thread count wins; otherwise one thread
#' without parallel, and with parallel the available physical cores: hyperthreads and big grids compete for the same caches. Never
#' more than the available (logical) cores: a larger thread count warns. A thread count below 1 stops.
#' @param parallel A boolean of whether to run in parallel
#' @param threads An integer number of threads or NULL
#' @returns An integer number of threads
#' @noRd
get_threads <- function(parallel=FALSE, threads=NULL){
  cores <- parallelly::availableCores()
  if(!is.null(threads)){
    threads <- suppressWarnings(as.numeric(threads))
    if(length(threads) != 1 || is.na(threads)){
      stop("threads must be one number")
    }
    if(threads < 1){
      stop("threads must be at least 1")
    }
    if(floor(threads) > cores){
      warning("threads=", threads, " is more than the ", cores, " available cores: using ", cores)
    }
    return(as.integer(min(floor(threads), cores)))
  }
  if(!parallel){
    return(1L)
  }
  cores <- min(cores, parallelly::availableCores(logical=FALSE), na.rm=TRUE)

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

#' Generate a background set of random gene interactions
#'
#' `generate_background_lrs` returns a dataframe of random gene interactions.
#'
#' This generates a table of random gene interactions where the number of genes in an interaction is equal to the complex size + 1.
#' @inheritParams prioritize_lr
#' @param count An integer of the number of interactions to generate
#' @returns A dataframe of random gene interactions
#' @noRd
generate_background_lrs <- function(lrs, metadata, count=1000){
  background <- data.frame(interaction_name=rep("background", count), check.names=FALSE)
  nodes <- metadata$cell[metadata$sample == "gene"]
  for(i in seq_len(ncol(lrs))[-1]){
    background[[i]] <- sample(nodes, count, replace=TRUE)
  }
  colnames(background) <- colnames(lrs)

  return(background)
}

#' Perform a permutation test
#'
#' `permutation_test` returns a numeric p-value
#'
#' This function takes a priority and the priorities of the background pairs and returns the fraction of background priorities at
#' least as high in absolute value (NA if any of them is NA or x is).
#' @param x A numeric priority
#' @param background A numeric vector of background priorities
#' @returns A numeric p-value
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
  grids <- pair_grids(ligands, receptors)
  ligand.order <- unique(ligands)
  per.ligand <- tabulate(match(grids$ligand, ligand.order), nbins=length(ligand.order))

  batch <- integer(length(ligand.order))
  b <- 1
  filled <- 0
  for(k in seq_along(ligand.order)){
    if(filled > 0 && filled + per.ligand[k] > size){
      b <- b + 1
      filled <- 0
    }
    batch[k] <- b
    filled <- filled + per.ligand[k]
  }

  return(split(seq_along(ligands), batch[match(ligands, ligand.order)]))
}

#' List the distinct grids of ligand-receptor pairs
#'
#' `pair_grids` returns the distinct ligand and receptor subunit pairs, one per cell grid
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and returns each distinct (ligand, receptor
#' subunit) pair once, in order of first appearance by receptor column: a complex is scored through one grid per subunit.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @returns A dataframe with columns ligand and receptor, one row per grid
#' @noRd
pair_grids <- function(ligands, receptors){
  grids <- unique(data.frame(ligand=rep(ligands, ncol(receptors)), receptor=as.vector(receptors)))

  return(grids[!is.na(grids$receptor), ])
}

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns the terms of the priority score of each ligand-receptor pair
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata
#' and returns the five terms of each pair's priority score. Pairs are scored in batches of whole ligands: each ligand and receptor
#' subunit grid is scored once, whatever the number of pairs using it, and only one batch of cell scores (about 64 MB at most, at least
#' one grid per thread) is in memory at a time. The gene x cell counts of a symmetric sparse adjacency are read from its stored
#' triangle, without expanding the whole matrix.
#' @inheritParams prioritize_lr
#' @param threads An integer number of threads
#' @param verbose A boolean of whether to report the pairs scored, about every 5 s, with an estimate of the time left (pairs are
#' scored in about 20 batches or more when grids are bigger than the kernel's 24 MB cache budget; see batch_size)
#' @returns A list with terms, a matrix with one row per pair and columns frac, pct.min, max, balance and purity, and title, the pair names
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, threads=1, verbose=FALSE){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells.idx <- metadata$sample != "gene"
  cells <- metadata$cell[cells.idx]
  groups <- metadata$group[cells.idx]

  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, -(1:2), drop=FALSE])
  genes <- unique(c(ligands, receptors[!is.na(receptors)]))
  check_genes(genes, rownames(infMat))

  if(inherits(adj, "hypercom_adjacency")){
    gene.counts <- adj[cells, genes, drop=FALSE]
  } else if(methods::is(adj, "dsCMatrix") && !any(genes %in% cells)){
    U <- methods::new("dgCMatrix", i=adj@i, p=adj@p, x=adj@x, Dim=adj@Dim, Dimnames=dimnames(adj))
    gene.counts <- as.matrix(Matrix::t(U[genes, cells, drop=FALSE]) + U[cells, genes, drop=FALSE])
    rm(U)
  } else {
    gene.counts <- as.matrix(Matrix::t(adj[genes, cells, drop=FALSE]))
  }
  gene.frac <- colSums(gene.counts) / length(cells)
  expressed <- lapply(seq_len(ncol(gene.counts)), function(k) which(gene.counts[, k] > 0))
  names(expressed) <- colnames(gene.counts)
  rm(gene.counts)

  blocks <- NULL
  if(inherits(infMat, "hypercom_influence")){
    blocks <- influence_blocks(infMat, grid_cells(rownames(infMat), metadata), threads)
  }

  terms <- matrix(NA_real_, nrow(lrs), 5, dimnames=list(NULL, c("frac", "pct.min", "max", "balance", "purity")))
  title <- character(nrow(lrs))
  n.grids <- nrow(pair_grids(ligands, receptors))
  size <- batch_size(threads, length(cells), tabulate(as.integer(factor(metadata$sample[cells.idx]))), n.grids)
  done <- 0
  grids.done <- 0
  t0 <- Sys.time()
  shown <- t0
  for(rows in batch_pairs(ligands, receptors, size)){
    pairs <- pair_grids(ligands[rows], receptors[rows, , drop=FALSE])
    keys <-paste(pairs$ligand, pairs$receptor, sep="\t")
    sr <- score_grids(infMat, metadata, pairs$ligand, pairs$receptor, threads, blocks)

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
      present <- logical(length(cells))
      for(g in pair.genes){
        present[expressed[[g]]] <- TRUE
      }
      pct.min <- min(gene.frac[pair.genes])
      frac <- sum(present) / length(cells)
      max.score <- mean(pmax(ss, rr)[present])

      purity <- group_purity(groups[present], gini.max)
      balance <- gini_impurity(type[present]) / 0.5

      terms[i, ] <- c(frac, pct.min, max.score, balance, purity)
      title[i] <- paste(c(l, rs), collapse="_")
    }
    rm(sr)
    done <- done + length(rows)
    grids.done <- grids.done + nrow(pairs)
    if(verbose && (as.numeric(difftime(Sys.time(), shown, units="secs")) >= 5 | done == nrow(lrs))){
      left <- as.numeric(difftime(Sys.time(), t0, units="secs")) / grids.done * (n.grids - grids.done)
      message("  ", done, " of ", nrow(lrs), " pairs scored", if(done < nrow(lrs)) paste0(", about ", format_duration(left), " left"))
      shown <- Sys.time()
    }
  }

  return(list(terms=terms, title=title))
}

#' Choose the number of grids per batch
#'
#' `batch_size` returns how many grids calculate_priority scores per kernel call
#'
#' This function takes the thread count, the number of cells, the cells of each sample and the number of grids and returns the batch
#' size: at least one grid per thread and at most about 64 MB of cell scores; when a grid is bigger than the kernel's 24 MB cache
#' budget (10 bytes per within-sample cell pair; CACHE_BUDGET in src/score_grids.cpp), all threads score one grid at a time, so smaller batches keep every thread busy and
#' there are about 20 batches or more, for progress messages. Each batch sets up the kernel's buffers again: on 20,000 cells, 104
#' pairs in 20 batches instead of 2 took 5% longer (90 instead of 86 s) and peaked 0.16 GB higher; the full CellChat database already
#' gives about 20 batches there.
#' @param threads An integer number of threads
#' @param cells An integer number of cells
#' @param sample.sizes An integer vector of the number of cells of each sample
#' @param n.grids An integer number of grids to score
#' @returns An integer number of grids per batch
#' @noRd
batch_size <- function(threads, cells, sample.sizes, n.grids){
  size <- max(threads, min(max(256, 32 * threads), floor(2^26 / (16 * cells))))
  if(10 * sum(as.numeric(sample.sizes)^2) > 24 * 2^20){
    size <- max(threads, min(size, ceiling(n.grids / 20)))
  }

  return(size)
}

#' Format a duration for a progress message
#'
#' `format_duration` returns a duration in seconds as a short string
#'
#' This function takes a number of seconds and returns it rounded in seconds below 90 s, in minutes below 90 min and in hours above.
#' @param seconds A double number of seconds
#' @returns A string, e.g. "40 s", "12 min" or "2.5 h"
#' @noRd
format_duration <- function(seconds){
  if(seconds < 90){
    return(paste0(max(1, round(seconds)), " s"))
  }
  if(seconds < 5400){
    return(paste0(round(seconds / 60), " min"))
  }

  return(paste0(round(seconds / 3600, 1), " h"))
}

#' Calculate the group purity of the cells of a pair
#'
#' `group_purity` returns how far the groups of a pair's cells are from evenly mixed
#'
#' This function takes the groups of the cells expressing a pair's genes and the largest possible Gini impurity and returns their
#' difference, at least 0: cells perfectly balanced over 5, 10 or 11 groups round the difference to -1.1e-16, whose log would make the
#' priority NaN. One group (largest impurity 0) gives 1.
#' @param groups A factor of the groups of the cells
#' @param gini.max A double of the largest possible Gini impurity, 1 - 1 / the number of groups
#' @returns A numeric purity
#' @noRd
group_purity <- function(groups, gini.max){
  if(gini.max <= 0){
    return(1)
  }

  return(max(0, gini.max - gini_impurity(groups)))
}

#' Combine priority terms into priority scores
#'
#' `score_priority` returns the priority score of each ligand-receptor pair for one weight
#'
#' This function takes the priority terms of calculate_priority and returns their weighted geometric mean.
#' @param priority.terms A list from calculate_priority
#' @param weight A weight of at least 0 for the fraction of cells expressing a pair's genes (frac), the first of the five terms
#' @returns A numeric vector of priority scores
#' @noRd
score_priority <- function(priority.terms, weight=1){
  w <- c(weight, rep(1, 4))
  terms <- priority.terms$terms
  priority <- numeric(nrow(terms))
  for(i in seq_len(nrow(terms))){
    priority[i] <- exp(stats::weighted.mean(log(terms[i, ]), w))
  }

  return(priority)
}

#' Check the weights of the priority scores
#'
#' `check_weights` stops unless the weights are distinct numbers of at least 0
#'
#' This function takes the weights argument of prioritize_lr and stops unless every weight is a finite number of at least 0 and no two
#' weights print alike (each names a result and an output file).
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_weights <- function(weights){
  if(!is.numeric(weights) || length(weights) == 0 || !all(is.finite(weights)) || any(weights < 0)){
    stop("weights must be numbers of at least 0, e.g. c(1, 0.5)")
  }
  if(anyDuplicated(as.character(weights)) > 0){
    stop("weights must differ: each names a result and an output file")
  }

  return(invisible(NULL))
}

#' Name the priority file of a weight
#'
#' `priority_filename` returns the name of the file prioritize_lr writes the priorities of a weight to
#'
#' This function takes weights and returns the file name of each, "priority.<weight>.csv".
#' @param weights A vector of weights
#' @returns A vector of strings of file names
#' @noRd
priority_filename <- function(weights){
  return(paste0("priority.", weights, ".csv"))
}

#' Order ligand-receptor pairs by signaling strength and specificity
#'
#' `prioritize_lr` returns an ordered list of ligand-receptor pairs
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata
#' and ranks the pairs by a priority score: the weighted geometric mean of five terms of the cells that express a pair's genes (their
#' fraction of all cells, the smallest fraction of cells expressing one of its genes, the mean of each cell's larger score as sender
#' or receiver, how evenly they split into senders and receivers, and how much they concentrate in some groups). Cells are scored once
#' for all weights; parallel or threads use several CPU cores.
#' With significance, the p-value of a pair is the fraction of 1000 background pairs of random genes whose priority is at least as
#' high, so it moves in steps of 0.001 and is 0 when no background pair scored as high; p.adj is its Benjamini-Hochberg adjustment.
#' @param lrs A dataframe of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param infMat An influence matrix from generate_infMat (a hypercom_influence or a named dense matrix)
#' @param metadata A dataframe from generate_metadata, containing columns named cell, sample and group
#' @param weights A vector of weights of at least 0, one result each, for the fraction of cells expressing a pair's genes in the priority
#' score: 1 weighs it like the other four terms, less than 1 favors pairs of fewer cells, 0 leaves it out
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel: one thread per physical core (one per logical core on Linux, where R does
#' not tell them apart)
#' @param outdir A string naming the output directory, or NULL to write no files
#' @param threads An integer number of threads of at least 1, overriding parallel; more than the available cores warns and uses
#' them all
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state
#' @param verbose A boolean of whether to report progress as messages: the pairs scored about every 5 s with
#' the time left, and the top pair of each weight
#' @returns A list of dataframes named by weight, each the ligand-receptor pairs ordered by their priority column, with columns p and
#' p.adj when significance is TRUE
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, outdir="6.priority",
                          threads=NULL, seed=NULL, verbose=TRUE){
  check_weights(weights)
  if(significance){
    if(!is.null(seed)){
      set.seed(seed)
    }
    background.lrs <- generate_background_lrs(lrs, metadata)
  }
  threads <- get_threads(parallel, threads)
  if(verbose){
    message("Scoring ", nrow(lrs), " ligand-receptor pairs", if(significance) paste0(" and ", nrow(background.lrs), " background pairs"),
            " on ", threads, " thread", if(threads > 1) "s")
  }

  if(significance){
    all.terms <- calculate_priority(rbind(lrs, background.lrs), adj, infMat, metadata, threads, verbose)
    fg <- seq_len(nrow(lrs))
    priority.terms <- list(terms=all.terms$terms[fg, , drop=FALSE], title=all.terms$title[fg])
    background.terms <- list(terms=all.terms$terms[-fg, , drop=FALSE], title=all.terms$title[-fg])
  } else {
    priority.terms <- calculate_priority(lrs, adj, infMat, metadata, threads, verbose)
  }

  if(!is.null(outdir)){
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  }
  lrs.list <- list()

  for(i in seq_along(weights)){
    weight <- weights[i]
    priority <- score_priority(priority.terms, weight)
    current.lrs <- lrs
    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      background.priority <- score_priority(background.terms, weight)
      current.lrs$p <- vapply(current.lrs$priority, permutation_test, numeric(1), background=background.priority)
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method="BH")
    }

    outfile <- NULL
    if(!is.null(outdir)){
      outfile <- file.path(outdir, priority_filename(weight))
      utils::write.csv(current.lrs, outfile, row.names=FALSE)
    }
    if(verbose){
      message("Weight ", weight, ": top pair ", current.lrs[[1]][1], if(!is.null(outfile)) paste0(", written to ", outfile))
    }

    lrs.list[[i]] <- current.lrs
  }

  names(lrs.list) <- weights

  return(lrs.list)
}

#' Run the HyperCom pipeline
#'
#' `HyperCom` runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`,
#' `generate_metadata`, and `prioritize_lr`, writing every step to outdir. By default no matrix with a row and a column for every
#' cell and gene is made, so memory grows with the number of cells rather than its square (see `generate_infMat`).
#' The arguments are checked before any step runs, and the run stops when no interaction of lrdb has all its genes expressed (e.g.
#' lrdb of another species). Priority files of other weights already in outdir are kept, with a warning. The results are returned as
#' well as written, so they need not be read back: `res <- HyperCom(...)`, then `res$priority[["1"]]` or
#' `score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")`.
#' @inheritParams filter_counts
#' @inheritParams generate_metadata
#' @inheritParams prioritize_lr
#' @param outdir A string naming the output directory of the pipeline run
#' @param verbose A boolean of whether to report progress as messages: what each step found, the pairs scored and the top pair of
#' each weight
#' @param infmat.method A string naming the method of `generate_infMat`: "woodbury" (the default; fastest and least memory), "chol" or
#' "lu" (a full matrix of all cells and genes; "lu" gives exactly the numbers of HyperCom 0.0.0.9000). The methods agree to about 1e-15
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency) and lrs (the interactions with all genes expressed)
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE,
                     outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="woodbury"){
  check_lrdb(lrdb)
  check_columns(metadata, sample, group, keep)
  check_weights(weights)
  infmat.method <- match_method(infmat.method)
  threads <- get_threads(parallel, threads)
  stale <- setdiff(list.files(file.path(outdir, "6.priority"), pattern="^priority\\..*\\.csv$"), priority_filename(weights))
  if(length(stale) > 0){
    warning("outdir already holds priorities of other weights, which are kept: ", paste(stale, collapse=", "))
  }
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  no.metadata <- setdiff(colnames(counts), rownames(metadata))
  if(length(no.metadata) > 0){
    warning(length(no.metadata), " cells of counts are not in metadata; those that are nodes of the graph count as genes")
  }

  counts <- filter_counts(counts, lrdb, outfile=file.path(outdir, "0.counts.rds"))
  lrs <- generate_lrs_table(counts, lrdb, outfile=file.path(outdir, "1.lrs.csv"))
  if(nrow(lrs) == 0){
    stop("no interaction of lrdb has all its genes expressed in counts (", ncol(counts), " of its genes are): check that lrdb is of the ",
         "species of counts (cellchat.mm.lrdb mouse, cellchat.hs.lrdb human) and that counts is gene x cell")
  }
  if(verbose){
    message(ncol(counts), " ligand and receptor genes expressed in ", nrow(counts), " cells; ", nrow(lrs), " interactions of lrdb have all ",
            "their genes expressed")
  }
  hyperedges <- generate_hyperedges(counts, lrs, outfile=file.path(outdir, "2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=file.path(outdir, "3.adj.rds"))
  rm(counts, hyperedges)
  infMat <- generate_infMat(adj, outfile=file.path(outdir, "4.infMat.rds"), method=infmat.method, threads=threads)
  if(verbose){
    message("Hypergraph of ", nrow(adj), " cells and genes; influence matrix made (method \"", infmat.method, "\")")
  }
  metadata <- generate_metadata(infMat, metadata, sample, group, keep, outfile=file.path(outdir, "5.metadata.csv"))
  priority <- prioritize_lr(lrs, adj, infMat, metadata, weights=weights, significance=significance, outdir=file.path(outdir, "6.priority"),
                            threads=threads, seed=seed, verbose=verbose)
  if(verbose){
    message("Done: results in ", outdir)
  }

  return(invisible(list(priority=priority, metadata=metadata, infMat=infMat, adj=adj, lrs=lrs)))
}

#' Get the counts matrix of a Seurat assay layer
#'
#' `get_seurat_layer` returns the matrix with dimnames for given assay and layer
#'
#' This function takes a Seurat object and returns one layer of one assay, joining split layers of v5 assays. Assays of other classes
#' (v3 assays, SeuratObject 4 and 5) are read from the slot of that name, as GetAssayData does in SeuratObject 5 (in SeuratObject 4
#' GetAssayData has no layer argument and returns the data slot).
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
    if(!(layer %in% c("counts", "data", "scale.data"))){
      stop("layer '", layer, "' is not a layer of assay '", assay, "': use counts, data or scale.data")
    }
    m <- methods::slot(a, layer)
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
#' `HyperCom_seurat` runs the HyperCom pipeline on a Seurat object
#'
#' This function runs `HyperCom` on one assay layer and the cell metadata of a Seurat object (v5 or older assays). Needs the
#' SeuratObject package.
#' @param seurat A Seurat object
#' @param sample A string that is the column name for sample in the cell metadata of the Seurat object
#' @param assay A string naming the assay to read counts from
#' @param layer A string naming the layer of the assay to read counts from
#' @param group A string that is the column name for group in the cell metadata of the Seurat object
#' @param keep A vector of strings that are additional column names to keep from the cell metadata of the Seurat object
#' @inheritParams HyperCom
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency) and lrs (the interactions with all genes expressed)
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, weights=c(1), significance=FALSE,
                            parallel=FALSE, outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="woodbury"){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat[[]]

  res <- HyperCom(lrdb, counts, metadata, sample=sample, group=group, keep=keep, weights=weights, significance=significance,
                  parallel=parallel, outdir=outdir, threads=threads, seed=seed, verbose=verbose, infmat.method=infmat.method)

  return(invisible(res))
}
