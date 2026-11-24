# Pipeline step 4: the influence matrix (methods woodbury, lu, chol, in R/grids.R; woodbury kept as its factors, class
# hypercom_influence)

#' Calculate influence matrix of hypergraph
#'
#' `generate_infMat` returns the influence matrix of the hypergraph
#'
#' This function takes an adjacency and returns its influence matrix by a random walk with restart: entry (i, j) is the share of a
#' walk started at node i that ends at node j, so each row sums to 1. At every step the walk returns to its start with probability r.
#' In matrix form it is r (I - (1 - r) adj / D)^-1, with D the row sums of adj. Three methods give the same matrix up to rounding
#' (about 1e-15 relative, which can reorder exactly tied scores):
#'
#' - "woodbury" (the default) works through the hyperedges and keeps the matrix as its factors (a hypercom_influence), so its memory
#'   grows with the number of nodes instead of its square. Data with at least as many hyperedges (of two or more nodes) as nodes get
#'   the dense matrix of "chol" instead, which is then no bigger: that is, data with no more cells than interactions of two or more
#'   distinct genes, all expressed (the bundled 1,000 cells stay factored: 1,661 hyperedges, 1,663 nodes). An adjacency from
#'   generate_adjacency_matrix holds its hyperedges; a plain adjacency matrix needs `hyperedges`.
#' - "lu" makes the dense matrix with exactly the numbers of the original release, HyperCom 0.0.0.9000; it holds two node x node
#'   matrices at once.
#' - "chol" makes the dense matrix faster, with one node x node matrix.
#'
#' A hypercom_influence works like a matrix: dim, dimnames, `infMat[i, j]` and as.matrix (the whole matrix; see
#' `hypercom_influence`). All methods stop when r is not in (0, 1] or the matrix cannot be inverted.
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param r A number greater than 0 and at most 1, the restart probability of the random walk (0.5 by default): about r of each
#' node's influence stays on itself and the rest spreads over the nodes it shares hyperedges with
#' @param method A string naming the method, "woodbury", "lu" or "chol", or NULL (the default) for "woodbury" when the hyperedges are
#' known (adj from generate_adjacency_matrix, or hyperedges given) and "lu" otherwise
#' @param hyperedges A data frame of the hyperedges a plain adjacency matrix was made from, needed by method "woodbury" with such a
#' matrix; ignored with an adjacency from generate_adjacency_matrix, which holds its own
#' @inheritParams filter_counts
#' @returns An influence matrix with node names: a hypercom_influence (method "woodbury", except on small data, see Details) or a
#' dense matrix
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' infMat
#' infMat[c("Thbs1", "Sdc4"), 1:3]
#' @export
generate_infMat <- function(adj, r=0.5, method=NULL, hyperedges=NULL, outfile="4.infMat.rds"){
  graph <- inherits(adj, "hypercom_adjacency")
  if(is.null(method)){
    method <- if(graph || !is.null(hyperedges)) "woodbury" else "lu"
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
  check_restart(r)
  if(method != "lu"){
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
    infMat <- solve_influence(adj, rs, r)
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

#' Check the restart probability
#'
#' `check_restart` stops unless r is one number in (0, 1]
#' @param r The r argument of generate_infMat or HyperCom
#' @returns NULL, invisibly
#' @noRd
check_restart <- function(r){
  if(!(is.numeric(r) && length(r) == 1 && isTRUE(r > 0 && r <= 1))){
    stop("r must be one number in (0, 1]")
  }

  return(invisible(NULL))
}

#' Match the name of an influence matrix method
#'
#' `match_method` returns the full name of an influence matrix method
#'
#' This function takes a method name as generate_infMat's method argument and returns it matched to "woodbury", "lu" or "chol"
#' (a unique abbreviation is completed; any other name stops).
#' @param method A string naming the method
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
#' @returns A hypercom_influence: a list with H, rs, delta, Z and r
#' @noRd
influence_factors <- function(H, rs, r){
  rs <- as.numeric(rs)
  f <- influence_factors_r(H, rs, r)
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
#' This function takes a hypercom_influence and positions of rows and columns and returns those entries as a dense matrix:
#' r (1 - r) Y[rows, ] Z[, cols], plus r rs / Delta where a row and a column are the same node.
#' @param infMat A hypercom_influence
#' @param rows An integer vector of row positions
#' @param cols An integer vector of column positions
#' @returns A dense matrix without names
#' @noRd
influence_slice <- function(infMat, rows, cols){
  rows <- as.integer(rows)
  cols <- as.integer(cols)
  r <- infMat$r
  Y <- infMat$H[rows, , drop=FALSE] / infMat$delta[rows]
  block <- Y %*% infMat$Z[, cols, drop=FALSE]
  block <- (r * (1 - r)) * methods::as(block, "matrix")
  dimnames(block) <- NULL

  if(anyDuplicated(cols) == 0){
    same <- match(rows, cols)
    hit <- which(!is.na(same))
    at <- cbind(hit, same[hit])
    block[at] <- block[at] + r * infMat$rs[rows[hit]] / infMat$delta[rows[hit]]
  } else {
    for(a in which(rows %in% cols)){
      b <- which(cols == rows[a])
      block[a, b] <- block[a, b] + r * infMat$rs[rows[a]] / infMat$delta[rows[a]]
    }
  }

  return(block)
}

#' Influence matrix of a hypergraph, stored as its factors
#'
#' `hypercom_influence` objects are made by generate_infMat with method "woodbury"
#'
#' A hypercom_influence is an influence matrix kept as smaller pieces from which any entry can be computed, so the whole node x node
#' matrix is never stored. It works like a matrix:
#'
#' - dim and dimnames give the number and the names of the nodes (cells and genes);
#' - `infMat[i, j]` computes the indexed entries as a dense matrix (indices as for a matrix: names, positions or logical; factors are
#'   read by name, not by their codes);
#' - as.matrix computes the whole dense matrix, 8 bytes times the number of nodes squared: slice instead on large data.
#'
#' For the record, the pieces are the factors of F = r Delta^-1 D + r (1 - r) Y Z: H (the node x hyperedge incidence), rs (the row sums
#' D of the adjacency), delta (Delta = D + (1 - r) diag(H H')), Z (hyperedge x node) and Y = H / Delta.
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
#' @exportS3Method dim hypercom_influence
dim.hypercom_influence <- function(x){
  return(node_dim(x))
}

#' @rdname hypercom_influence
#' @exportS3Method dimnames hypercom_influence
dimnames.hypercom_influence <- function(x){
  return(node_dimnames(x))
}

#' @rdname hypercom_influence
#' @export
`[.hypercom_influence` <- function(x, i, j, drop=TRUE){
  idx <- block_index(rownames(x$H), i, j, nargs() - !missing(drop), "an influence matrix as infMat[i, j]")

  return(name_block(influence_slice(x, idx$rows, idx$cols), idx, drop))
}

#' @rdname hypercom_influence
#' @method as.matrix hypercom_influence
#' @export
as.matrix.hypercom_influence <- function(x, ...){
  nodes <- rownames(x$H)
  block <- influence_slice(x, seq_along(nodes), seq_along(nodes))
  dimnames(block) <- node_dimnames(x)

  return(block)
}

#' @rdname hypercom_influence
#' @method print hypercom_influence
#' @export
print.hypercom_influence <- function(x, ...){
  cat("HyperCom influence matrix of ", nrow(x$H), " nodes (r = ", x$r, "), stored as its factors through ", ncol(x$H), " hyperedges\n",
      sep="")

  return(invisible(x))
}
