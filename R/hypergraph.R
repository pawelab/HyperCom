# Pipeline steps 2-3: hyperedges, and the adjacency matrix of the hypergraph kept as its incidence (class hypercom_adjacency)

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
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' head(hyperedges)
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
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' adj
#' dim(adj)
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds"){
  adj <- list(H=hyperedge_incidence(index_hyperedges(hyperedges)))
  class(adj) <- "hypercom_adjacency"

  if(!is.null(outfile)){
    saveRDS(adj, outfile)
  }
  return(adj)
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
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' adj
#' dim(adj)
#' adj[c("Thbs1", "Sdc4"), 1:3]
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
