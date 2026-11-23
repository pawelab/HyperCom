# Pipeline step 4: the influence matrix (methods woodbury, lu, chol; woodbury kept as its factors, class hypercom_influence)

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
#' - "lu" makes the dense matrix with exactly the numbers of HyperCom 0.0.0.9000; it holds four node x node matrices at once.
#' - "chol" makes the dense matrix faster, with three node x node matrices at once.
#'
#' All methods stop when the matrix cannot be inverted (r near 0).
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @param method A string naming the method: "woodbury" (the default; "lu" for an adjacency matrix given without hyperedges), "lu"
#' or "chol"
#' @param hyperedges A dataframe of the hyperedges an adjacency matrix was made from, needed by method "woodbury" with an adjacency
#' matrix; ignored with an adjacency from generate_adjacency_matrix, which holds its own
#' @inheritParams filter_counts
#' @returns An influence matrix with node names: a hypercom_influence (method "woodbury", unless there are at least as many hyperedges
#' as nodes) or a dense matrix
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' infMat
#' infMat[c("Thbs1", "Sdc4"), 1:3]
#' @export
generate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds", method=c("woodbury", "lu", "chol"), hyperedges=NULL){
  # woodbury: adj = H H' with a zero diagonal (H the node x hyperedge incidence), so the Woodbury identity inverts it through a
  # hyperedge x hyperedge matrix and the result stays factored (influence_factors). lu: the arithmetic and solve() of HyperCom
  # 0.0.0.9000, bit for bit. chol: chol() and chol2inv() of a symmetric matrix similar to I - W.
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
    infMat <- solve_influence_lu(adj, rs, r)
  } else if(method == "chol"){
    infMat <- solve_influence_chol(adj, rs, r)
  } else {
    if(!graph){
      H <- incidence_matrix(hyperedges, adj)
    }
    if(ncol(H) < nrow(H)){
      infMat <- influence_factors(H, rs, r)
    } else {
      if(graph){
        adj <- adjacency_sparse(adj)
      }
      infMat <- solve_influence_chol(adj, rs, r)
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

#' Influence matrix of a hypergraph, stored as its factors
#'
#' `hypercom_influence` objects are made by generate_infMat with method "woodbury"
#'
#' A hypercom_influence keeps the factors of the influence matrix F = r Delta^-1 D + r (1 - r) Y Z of a hypergraph (H, the node x
#' hyperedge incidence; rs, the row sums D of the adjacency; delta, Delta = D + (1 - r) diag(H H'); Z, hyperedge x node; Y = H /
#' Delta) instead of F. dim and dimnames give those of F, `infMat[i, j]` makes the indexed entries as a dense matrix (indices as for a
#' matrix: names, positions or logical; factors are read by name, not by their codes) and as.matrix the whole dense
#' matrix (nodes squared doubles).
#' @param x A hypercom_influence
#' @param i,j Indices of rows and columns
#' @param drop A boolean of whether to drop dimensions of extent one
#' @param ... Ignored
#' @returns dim: the number of nodes twice; dimnames: the node names twice; `[` and as.matrix: a dense matrix; print: x, invisibly
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' infMat
#' dim(infMat)
#' infMat[c("Thbs1", "Sdc4"), 1:3]
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
as.matrix.hypercom_influence <- function(x, ...){
  nodes <- rownames(x$H)
  block <- influence_slice(x, seq_along(nodes), seq_along(nodes))
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
