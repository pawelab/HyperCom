# Pipeline steps 2-3: hyperedges, and the adjacency matrix of the hypergraph kept as its incidence (class hypercom_adjacency)

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns a data frame of hyperedges
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptor interactions to produce
#' a table of hyperedges (a hyperedge links any number of nodes): one per ligand-receptor interaction (its genes) and one per gene
#' (the gene and the cells that express it). Cells are nodes named by the row names of counts, so they must be unique (prefix
#' barcodes with their sample when merging samples).
#'
#' Only the largest connected part of the hypergraph is kept. Cells and genes that no chain of shared hyperedges joins to it (e.g. the
#' genes of an interaction that shares no gene with any other, and cells that express only those genes) would have no influence on
#' the other cells. They are
#' left out, with a warning that counts them ("... cells and ... genes are not connected to the rest of the hypergraph and are left
#' out"); a few such cells are normal and need no action.
#' @inheritParams generate_lrs_table
#' @param lrs A data frame of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @returns A tibble (a data frame) with one row per node of each hyperedge: columns edge (the hyperedge, an interaction or a gene), node (a gene or
#' cell in it) and weight (1)
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' head(hyperedges)
#' @export
generate_hyperedges <- function(counts, lrs, outfile="2.hyperedges.csv"){
  lrs <- clean_lrdb(lrs, "lrs", empty.ok=TRUE)
  dup.cells <- unique(rownames(counts)[duplicated(rownames(counts))])
  if(length(dup.cells) > 0){
    stop("cell names used by several cells of counts (prefix them with their sample): ", paste(utils::head(dup.cells, 10), collapse=", "))
  }

  lr.edges <- unique(tidyr::pivot_longer(lrs, -1, values_to="node", names_to=NULL, values_drop_na=TRUE))
  lr.edges$weight <- 1
  colnames(lr.edges) <- c("edge", "node", "weight")

  X <- methods::as(methods::as(methods::as(Matrix::t(counts), "dMatrix"), "generalMatrix"), "CsparseMatrix")
  cell.idx <- rep(seq_len(ncol(X)), diff(X@p))
  expressed <- X@x != 0 & !is.na(X@x)

  genes <- unique(c(colnames(counts), lr.edges$node))
  gene.idx <- match(colnames(counts), genes)
  edge <- c(match(lr.edges$edge, lrs[[1]]), nrow(lrs) + gene.idx, nrow(lrs) + gene.idx[X@i[expressed] + 1])
  node <- c(match(lr.edges$node, genes), gene.idx, length(genes) + cell.idx[expressed])
  left <- outside_nodes(edge, node, nrow(lrs) + length(genes), length(genes) + ncol(X))
  genes.left <- genes[left[seq_along(genes)]]
  cells.left <- left[length(genes) + seq_len(ncol(X))]
  if(any(left)){
    warning(sum(cells.left), " cells and ", length(genes.left), " genes are not connected to the rest of the hypergraph and are left out")
    lr.edges <- lr.edges[!(lr.edges$node %in% genes.left), , drop=FALSE]
    expressed <- expressed & !cells.left[cell.idx]
  }

  gene.kept <- !(colnames(counts) %in% genes.left)
  gene.edges <- data.frame(edge=colnames(counts)[gene.kept], node=colnames(counts)[gene.kept], weight=rep(1, sum(gene.kept)),
                           check.names=FALSE)
  exp.edges <- data.frame(edge=rownames(X)[X@i[expressed] + 1], node=colnames(X)[cell.idx[expressed]], weight=rep(1, sum(expressed)),
                          check.names=FALSE)

  edges <- rbind(lr.edges, gene.edges, exp.edges)

  if(!is.null(outfile)){
    utils::write.csv(edges, outfile, row.names=FALSE)
  }
  return(edges)
}

#' Find the nodes outside the largest connected part of a hypergraph
#'
#' `outside_nodes` returns which nodes of a hypergraph are not in its largest connected component
#'
#' This function takes the hyperedge and node index of each row of a hypergraph's incidence and returns, for each node, whether it is
#' a node of the hypergraph (in a hyperedge of at least two distinct nodes) outside its largest connected component (the first of
#' equally large ones). Nodes that are in no such hyperedge are not nodes of the hypergraph and are never outside.
#' @param edge An integer vector of hyperedge indices, one per row
#' @param node An integer vector of node indices, one per row
#' @param n.edge An integer number of hyperedges
#' @param n.node An integer number of nodes
#' @returns A logical vector, one per node
#' @noRd
outside_nodes <- function(edge, node, n.edge, n.node){
  distinct <- !duplicated(as.numeric(edge) * (n.node + 1) + node)
  linked <- tabulate(edge[distinct], nbins=n.edge)[edge] >= 2
  in.graph <- logical(n.node)
  in.graph[node[linked]] <- TRUE
  label <- incidence_components(node[linked], edge[linked], n.node, n.edge)
  size <- tabulate(label[in.graph], nbins=n.node)

  return(in.graph & label != which.max(size))
}

#' Label the connected components of a hypergraph
#'
#' `incidence_components` returns the connected component of each node of a hypergraph, from its incidence
#'
#' This function takes the node and hyperedge index of each row of an incidence and labels each node with the smallest node index of
#' its component: labels pass from nodes to their hyperedges and back, each step keeping the smallest, and each label then jumps to the
#' label of the node it names, until no label changes. Nodes in no row keep their own index.
#' @param node An integer vector of node indices, one per row
#' @param edge An integer vector of hyperedge indices, one per row
#' @param n.node An integer number of nodes
#' @param n.edge An integer number of hyperedges
#' @returns An integer vector of component labels, one per node
#' @noRd
incidence_components <- function(node, edge, n.node, n.edge){
  label <- seq_len(n.node)
  repeat{
    edge.label <- group_min(label[node], edge, n.edge)
    new <- pmin(label, group_min(edge.label[edge], node, n.node), na.rm=TRUE)
    new <- new[new]
    if(identical(new, label)){
      break
    }
    label <- new
  }

  return(label)
}

#' Take the smallest value of each group
#'
#' `group_min` returns the smallest value of each group, NA for groups without a value
#'
#' This function takes values and their group indices and returns the smallest value per group: the values are written in decreasing
#' order, so the last write to each group is its smallest.
#' @param x An integer vector of values
#' @param group An integer vector of group indices, one per value
#' @param n An integer number of groups
#' @returns An integer vector, one per group
#' @noRd
group_min <- function(x, group, n){
  out <- rep(NA_integer_, n)
  o <- order(x, decreasing=TRUE, method="radix")
  out[group[o]] <- x[o]

  return(out)
}

#' Index the nodes and hyperedges of a hyperedge table
#'
#' `index_hyperedges` returns the node names and the node and hyperedge index of each row
#'
#' This function takes a data frame of hyperedges and numbers its nodes by their names in the C locale and its hyperedges in order
#' of appearance, leaving out rows with a missing node or hyperedge.
#' @param hyperedges A data frame of hyperedges with columns named edge and node
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
#' This function takes a data frame of hyperedges and returns its adjacency matrix: each entry is the number of hyperedges two nodes
#' share, and the diagonal is zero. The nodes are the cells, ligand genes and receptor genes, ordered by name (in the C locale, so the
#' order does not depend on the session's language settings). Nodes that share no hyperedge with another node are dropped. The full
#' node x node matrix is not stored: the result records which nodes are in which hyperedge, so it takes memory in proportion to the
#' hyperedge table instead of the number of nodes squared. It works like a matrix: dim, dimnames, `adj[i, j]` and as.matrix (the whole
#' matrix; see `hypercom_adjacency`).
#' @param hyperedges A data frame of hyperedges from generate_hyperedges, with columns edge (the hyperedge) and node (a node in it)
#' @inheritParams filter_counts
#' @returns A hypercom_adjacency (a list with H, the sparse node x hyperedge matrix of how often each node is in each hyperedge),
#' which works like the adjacency matrix
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

#' Expand an adjacency into a sparse matrix
#'
#' `adjacency_sparse` returns the adjacency matrix of a hypercom_adjacency as a symmetric sparse matrix (dsCMatrix)
#'
#' This function takes an adjacency from generate_adjacency_matrix and returns the upper triangle of its symmetric matrix with the
#' explicit zero diagonal, as the methods "lu" and "chol" read it. It takes memory in proportion to the node pairs that share a
#' hyperedge.
#' @param adj A hypercom_adjacency
#' @returns A dsCMatrix with node names
#' @noRd
adjacency_sparse <- function(adj){
  H <- adj$H
  A <- adjacency(H)
  nodes <- rownames(H)[A$keep]

  return(methods::new("dsCMatrix", i=A$i, p=A$p, x=A$x, Dim=rep(length(nodes), 2L), Dimnames=list(node=nodes, node=nodes), uplo="U"))
}

#' Adjacency matrix of a hypergraph, stored as its incidence
#'
#' `hypercom_adjacency` objects are made by generate_adjacency_matrix
#'
#' A hypercom_adjacency keeps which nodes are in which hyperedge (the node x hyperedge incidence H of the hypergraph) instead of the
#' adjacency matrix, which is H H' with the diagonal set to zero. It works like a matrix:
#'
#' - dim and dimnames give the number and the names of the nodes (cells and genes);
#' - `adj[i, j]` computes the indexed entries as a dense matrix (indices as for a matrix: names, positions or logical; factors are
#'   read by name, not by their codes);
#' - as.matrix computes the whole dense matrix, 8 bytes times the number of nodes squared: slice instead on large data.
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
#' @exportS3Method dim hypercom_adjacency
dim.hypercom_adjacency <- function(x){
  return(node_dim(x))
}

#' @rdname hypercom_adjacency
#' @exportS3Method dimnames hypercom_adjacency
dimnames.hypercom_adjacency <- function(x){
  return(node_dimnames(x))
}

#' @rdname hypercom_adjacency
#' @export
`[.hypercom_adjacency` <- function(x, i, j, drop=TRUE){
  idx <- block_index(rownames(x$H), i, j, nargs() - !missing(drop), "an adjacency as adj[i, j]")
  rows <- idx$rows
  cols <- idx$cols
  A <- as.matrix(Matrix::tcrossprod(x$H[rows, , drop=FALSE], x$H[cols, , drop=FALSE]))
  if(anyDuplicated(cols) == 0){
    same <- match(rows, cols)
    A[cbind(which(!is.na(same)), same[!is.na(same)])] <- 0
  } else {
    for(a in which(rows %in% cols)){
      A[a, cols == rows[a]] <- 0
    }
  }

  return(name_block(A, idx, drop))
}

#' @rdname hypercom_adjacency
#' @method as.matrix hypercom_adjacency
#' @export
as.matrix.hypercom_adjacency <- function(x, ...){
  return(as.matrix(adjacency_sparse(x)))
}

#' @rdname hypercom_adjacency
#' @method print hypercom_adjacency
#' @export
print.hypercom_adjacency <- function(x, ...){
  cat("HyperCom adjacency matrix of ", nrow(x$H), " nodes, stored as their incidence in ", ncol(x$H), " hyperedges\n", sep="")

  return(invisible(x))
}

#' Generate the incidence matrix behind an adjacency matrix
#'
#' `incidence_matrix` returns the node x hyperedge incidence counts of the hyperedges an adjacency matrix was made from
#'
#' This function takes a data frame of hyperedges and its adjacency matrix and returns how often each node is in each hyperedge, for
#' the hyperedges with at least two distinct nodes (the only ones that add to the adjacency), nodes in the order of the adjacency.
#' It stops unless the incidence gives the adjacency (the incidence times its transpose with the diagonal set to zero): the same
#' nodes, and adj times two fixed integer vectors exactly equal to the products through the incidence (integers small enough for
#' every sum to be exact in double precision; no random numbers are drawn).
#' @param hyperedges A data frame of hyperedges with columns named edge and node
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
