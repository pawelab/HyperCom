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
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' infMat
#' infMat[c("Thbs1", "Sdc4"), 1:3]
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
