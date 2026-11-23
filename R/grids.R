# R versions of the compiled parts of HyperCom 0.2: the adjacency from the incidence, the influence matrix (lu, chol, woodbury
# factors and blocks of them), and the scoring and screening of cell grids, which see only the within-sample cell x cell blocks.

#' Expand an adjacency into the sparse matrix of HyperCom 0.0.0.9000
#'
#' `adjacency_sparse` returns the adjacency matrix of a hypercom_adjacency as the dsCMatrix of HyperCom 0.0.0.9000
#'
#' This function takes an adjacency from generate_adjacency_matrix and returns the upper triangle of its symmetric matrix H H' with the
#' explicit zero diagonal, the object generate_adjacency_matrix made in HyperCom 0.0.0.9000. It takes memory in proportion to the node
#' pairs that share a hyperedge.
#' @param adj A hypercom_adjacency
#' @returns A dsCMatrix with node names
#' @noRd
adjacency_sparse <- function(adj){
  H <- adj$H
  A <- methods::as(Matrix::tcrossprod(H), "symmetricMatrix")
  A <- methods::as(methods::as(A, "CsparseMatrix"), "dMatrix")
  if(A@uplo != "U"){
    A <- Matrix::t(A)
  }
  col <- rep(seq_len(ncol(A)) - 1L, diff(A@p))
  A@x[A@i == col] <- 0
  A@Dimnames <- list(node=rownames(H), node=rownames(H))

  return(A)
}

#' Solve for the influence matrix as HyperCom 0.0.0.9000 did
#'
#' `solve_influence_lu` returns the dense influence matrix r (I - (1 - r) adj / rs)^-1 with the numbers of HyperCom 0.0.0.9000
#'
#' This function takes an adjacency matrix, its row sums and the restart probability and returns the influence matrix with the same
#' arithmetic as HyperCom 0.0.0.9000 (R's solve), writing I - W without a dense identity: four node x node matrices at once instead
#' of six.
#' @param adj An adjacency matrix: a base matrix or a sparse Matrix
#' @param rs A numeric vector of the row sums of adj
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A dense matrix
#' @noRd
solve_influence_lu <- function(adj, rs, r){
  W <- adj/rs
  W <- (1 - r) * W
  W <- as.matrix(W)
  M <- 0 - W
  diag(M) <- 1 - diag(W)
  rm(W)

  return(r * Matrix::solve(M))
}

#' Solve for the influence matrix by Cholesky
#'
#' `solve_influence_chol` returns the dense influence matrix r (I - (1 - r) adj / rs)^-1 through a symmetric matrix
#'
#' This function takes a symmetric adjacency matrix, its row sums and the restart probability and returns the influence matrix as
#' r D^-1/2 S^-1 D^1/2 with S = I - (1 - r) D^-1/2 adj D^-1/2 (D the row sums), S inverted by Cholesky (chol, chol2inv). It stops when
#' S is computationally singular.
#' @inheritParams solve_influence_lu
#' @returns A dense matrix
#' @noRd
solve_influence_chol <- function(adj, rs, r){
  s <- 1 / sqrt(rs)
  S <- as.matrix(adj)
  S <- -(1 - r) * (s * S)
  S <- t(s * t(S))
  diag(S) <- diag(S) + 1
  R <- factor_spd(S)
  rm(S)
  infMat <- chol2inv(R)
  rm(R)
  infMat <- (r * s) * infMat
  infMat <- t(t(infMat) / s)

  return(infMat)
}

#' Factor a symmetric positive definite matrix
#'
#' `factor_spd` returns the upper Cholesky factor of a matrix, stopping when it is computationally singular
#'
#' This function takes a symmetric positive definite matrix and returns chol's upper triangular factor; it stops, as solve does, when
#' the matrix is not positive definite or its reciprocal condition number (estimated from the factor) is below the machine epsilon,
#' e.g. for a restart probability near 0.
#' @param S A symmetric matrix
#' @returns An upper triangular matrix
#' @noRd
factor_spd <- function(S){
  R <- tryCatch(chol(S), error=function(e) NULL)
  if(is.null(R) || rcond(R, triangular=TRUE)^2 < .Machine$double.eps){
    stop("the influence matrix is computationally singular: the restart probability r is too close to 0")
  }

  return(R)
}

#' Factor the influence matrix through the hyperedges
#'
#' `influence_factors` returns the influence matrix of a hypergraph as a hypercom_influence
#'
#' This function takes the node x hyperedge incidence, the row sums of the adjacency and the restart probability and returns the
#' factors of the influence matrix. The adjacency is H H' with a zero diagonal, so K = D - (1 - r) adj = Delta - (1 - r) H H' with the
#' diagonal Delta = D + (1 - r) diag(H H'), and by the Woodbury identity F = r K^-1 D = r Delta^-1 D + r (1 - r) Y Z with Y = H / Delta,
#' C = I - (1 - r) H' Y (hyperedge x hyperedge, symmetric positive definite) and Z = C^-1 Y' D. Every entry is a sum of terms of at
#' least 0. Cost: hyperedges cubed plus hyperedges times the nonzeros of H, instead of nodes cubed.
#' @param H A sparse node x hyperedge matrix of counts with node and hyperedge names
#' @param rs A numeric vector of the row sums of the adjacency
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A hypercom_influence: a list with H, rs, delta, Z and r
#' @noRd
influence_factors <- function(H, rs, r){
  rs <- as.numeric(rs)
  delta <- rs + (1 - r) * Matrix::rowSums(H^2)
  Y <- H / delta
  C <- as.matrix(Matrix::crossprod(H, Y))
  C <- diag(nrow(C)) - (1 - r) * C
  R <- factor_spd(C)
  rm(C)
  Z <- as.matrix(chol2inv(R) %*% Matrix::t(Y * rs))
  dimnames(Z) <- list(edge=colnames(H), node=rownames(H))
  infMat <- list(H=H, rs=rs, delta=delta, Z=Z, r=r)
  class(infMat) <- "hypercom_influence"

  return(infMat)
}

#' Make a block of a factored influence matrix
#'
#' `influence_slice` returns the entries of a hypercom_influence at given rows and columns
#'
#' This function takes a hypercom_influence and positions of rows and columns and returns those entries as a dense matrix:
#' r (1 - r) Y[rows, ] Z[, cols], plus r rs / Delta where a row and a column are the same node. Each entry is summed over the row
#' node's hyperedges in ascending order, so it is the same in any block.
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

#' Make the within-sample blocks of an influence matrix
#'
#' `grid_blocks` returns every within-sample cell x cell block of an influence matrix
#'
#' This function takes an influence matrix and the cells from grid_cells and returns one block per sample (samples ascending, cells of a
#' sample in the order given): made from the factors of a hypercom_influence, copied from a dense matrix. They hold the sum of
#' squared sample sizes of doubles.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param g A list from grid_cells
#' @returns A list of dense matrices, one per sample
#' @noRd
grid_blocks <- function(infMat, g){
  if(sum(as.numeric(tabulate(g$sample))^2) >= 2^52){
    stop("too many cells per sample: the sum of squared sample sizes must stay below 2^52")
  }
  blocks <- lapply(split(g$idx + 1L, g$sample), function(pos){
    if(inherits(infMat, "hypercom_influence")){
      return(influence_slice(infMat, pos, pos))
    }
    block <- infMat[pos, pos, drop=FALSE]
    dimnames(block) <- NULL
    return(block)
  })

  return(unname(blocks))
}

#' Make the transfer grid of a ligand and receptor in each sample
#'
#' `grid_values` returns the entries of a pair's cell x cell grid within each sample, as HyperCom 0.0.0.9000 ranks them
#'
#' This function takes the within-sample blocks of the influence matrix, the influence of the cells on the ligand (one value per cell)
#' and of the receptor on the cells, and the cells of each sample, and returns all within-sample entries A[i] B[j] / F[i, j] one after
#' the other (each block by column) with the diagonal, zeros and NaN as NA: the entries of HyperCom 0.0.0.9000's n x n grid that it did
#' not set to NA, as it computed them. With a and b NULL the grid is the influence matrix itself (the network score).
#' @param blocks A list of blocks from grid_blocks
#' @param a A numeric vector of the influence of each cell on the ligand, or NULL
#' @param b A numeric vector of the influence of the receptor on each cell, or NULL
#' @param samples A list of the positions of each sample's cells
#' @returns A numeric vector
#' @noRd
grid_values <- function(blocks, a, b, samples){
  vals <- vector("list", length(samples))
  for(s in seq_along(samples)){
    idx <- samples[[s]]
    if(is.null(a)){
      trans <- blocks[[s]]
    } else {
      trans <- outer(a[idx], b[idx]) / blocks[[s]]
    }
    diag(trans) <- 0
    trans[trans == 0] <- NA
    vals[[s]] <- trans
  }

  return(unlist(vals, use.names=FALSE))
}

#' Read the influence of cells on ligands and of receptors on cells
#'
#' `grid_genes` returns the columns of the ligands and the rows of the receptors of an influence matrix at the cells
#'
#' This function takes an influence matrix, the cells from grid_cells and the distinct ligand and receptor genes and returns
#' F[cells, ligands] and F[receptors, cells], made from the factors of a hypercom_influence.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param g A list from grid_cells
#' @param lig.genes A vector of strings of ligand names
#' @param rec.genes A vector of strings of receptor names
#' @returns A list of two dense matrices, A (cells x ligands) and B (receptors x cells)
#' @noRd
grid_genes <- function(infMat, g, lig.genes, rec.genes){
  nodes <- rownames(infMat)
  if(inherits(infMat, "hypercom_influence")){
    A <- influence_slice(infMat, g$idx + 1L, match(lig.genes, nodes))
    B <- influence_slice(infMat, match(rec.genes, nodes), g$idx + 1L)
  } else {
    A <- unname(infMat[g$idx + 1L, match(lig.genes, nodes), drop=FALSE])
    B <- unname(infMat[match(rec.genes, nodes), g$idx + 1L, drop=FALSE])
  }

  return(list(A=A, B=B))
}

#' Bound the ranks of grid entries without ranking them
#'
#' `bound_ranks` returns, for each entry, bounds of the number of entries below it
#'
#' This function takes grid entries (NA left out) and places each in one of 2^20 value buckets of equal width on the 64th root of
#' the values, nearly equal width on their logarithm (infinite values in one more bucket at the top). Each step (square roots,
#' each rounded exactly; subtraction; multiplication by a positive number; truncation) never decreases, so every entry of a lower
#' bucket is smaller: an entry of bucket d has between S_d and S_d + n_d - 1 smaller entries (S_d entries below the bucket, n_d in
#' it), the bounds of its min rank - 1.
#' @param v A numeric vector of grid entries of at least 0, NA where an entry is not ranked
#' @returns A list of two numeric vectors, lo and hi, NA where v is NA, and the number k of entries ranked
#' @noRd
bound_ranks <- function(v){
  buckets <- 2^20
  y <- sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(v))))))
  finite <- is.finite(y)
  lo <- if(any(finite)) min(y[finite]) else 0
  hi <- if(any(finite)) max(y[finite]) else 0
  scale <- if(hi > lo) (buckets - 1) / (hi - lo) else 0
  d <- y
  d[finite] <- (y[finite] - lo) * scale
  rm(y)
  d <- trunc(d) + 1
  d[is.infinite(v)] <- buckets + 1
  d <- as.integer(d)
  counts <- tabulate(d, nbins=buckets + 1L)
  S <- c(0, cumsum(as.numeric(counts)))

  return(list(lo=S[d], hi=S[d + 1L] - 1, k=sum(counts)))
}

#' Calculate sender and receiver scores of cell grids
#'
#' `score_grids` returns the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, and vectors of ligands and receptor subunits.
#' For each pair it percent-ranks all entries of the cell x cell transfer grid between cells of the same sample together and averages
#' the ranks of each cell as sender and as receiver, as HyperCom 0.0.0.9000 did on the whole grid, whose entries between samples it
#' left out. A pair with ligand and receptor NA gives the network score. Cells and genes are found in the influence matrix by name, so
#' the metadata rows may be in any order.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligands A vector of strings of ligand names, NA for the network score
#' @param receptors A vector of strings of receptor subunit names, NA for the network score
#' @param blocks A list from grid_blocks, NULL to make it here
#' @returns A list of two matrices, sender and receiver, with one row per cell and one column per pair
#' @noRd
score_grids <- function(infMat, metadata, ligands, receptors, blocks=NULL){
  ligands <- as.character(ligands)
  receptors <- as.character(receptors)
  nodes <- rownames(infMat)
  check_genes(c(ligands, receptors), nodes)
  if(any(is.na(ligands) != is.na(receptors))){
    stop("each pair needs both a ligand and a receptor, or neither for the network score")
  }
  g <- grid_cells(nodes, metadata)
  if(is.null(blocks)){
    blocks <- grid_blocks(infMat, g)
  }
  samples <- unname(split(seq_along(g$idx), g$sample))
  lig.genes <- unique(ligands[!is.na(ligands)])
  rec.genes <- unique(receptors[!is.na(receptors)])
  gg <- grid_genes(infMat, g, lig.genes, rec.genes)

  sender <- receiver <- matrix(NA_real_, length(g$idx), length(ligands))
  for(k in seq_along(ligands)){
    if(is.na(ligands[k])){
      v <- grid_values(blocks, NULL, NULL, samples)
    } else {
      v <- grid_values(blocks, gg$A[, match(ligands[k], lig.genes)], gg$B[match(receptors[k], rec.genes), ], samples)
    }
    v <- dplyr::percent_rank(v)
    start <- 0
    for(idx in samples){
      m <- length(idx)
      R <- v[start + seq_len(m * m)]
      dim(R) <- c(m, m)
      sender[idx, k] <- rowMeans(R, na.rm=TRUE)
      receiver[idx, k] <- colMeans(R, na.rm=TRUE)
      start <- start + m * m
    }
  }
  if(is.unsorted(g$o)){
    back <- order(g$o)
    sender <- sender[back, , drop=FALSE]
    receiver <- receiver[back, , drop=FALSE]
  }

  return(list(sender=sender, receiver=receiver))
}

#' Bound the sender and receiver scores of cell grids
#'
#' `screen_grids` returns lower and upper bounds of the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, and vectors of ligands and receptor subunits, and returns
#' for each grid bounds of the scores score_grids would return, from one pass without ranking: each entry is only placed in a value
#' bucket, which bounds its rank (bound_ranks). Each bound is a cell's exact integer sum of rank bounds divided once by the entries
#' ranked minus 1 and once by the cell's entries.
#' @inheritParams score_grids
#' @returns A list of four matrices, sender.lo, sender.hi, receiver.lo and receiver.hi, with one row per cell and one column per grid
#' @noRd
screen_grids <- function(infMat, metadata, ligands, receptors, blocks=NULL){
  ligands <- as.character(ligands)
  receptors <- as.character(receptors)
  nodes <- rownames(infMat)
  check_genes(c(ligands, receptors), nodes)
  g <- grid_cells(nodes, metadata)
  if(is.null(blocks)){
    blocks <- grid_blocks(infMat, g)
  }
  samples <- unname(split(seq_along(g$idx), g$sample))
  lig.genes <- unique(ligands)
  rec.genes <- unique(receptors)
  gg <- grid_genes(infMat, g, lig.genes, rec.genes)

  sender.lo <- sender.hi <- receiver.lo <- receiver.hi <- matrix(NA_real_, length(g$idx), length(ligands))
  for(k in seq_along(ligands)){
    v <- grid_values(blocks, gg$A[, match(ligands[k], lig.genes)], gg$B[match(receptors[k], rec.genes), ], samples)
    b <- bound_ranks(v)
    rm(v)
    start <- 0
    for(idx in samples){
      m <- length(idx)
      pos <- start + seq_len(m * m)
      lo <- b$lo[pos]
      hi <- b$hi[pos]
      dim(lo) <- dim(hi) <- c(m, m)
      row.n <- rowSums(!is.na(lo))
      col.n <- colSums(!is.na(lo))
      sender.lo[idx, k] <- rowSums(lo, na.rm=TRUE) / (b$k - 1) / row.n
      sender.hi[idx, k] <- rowSums(hi, na.rm=TRUE) / (b$k - 1) / row.n
      receiver.lo[idx, k] <- colSums(lo, na.rm=TRUE) / (b$k - 1) / col.n
      receiver.hi[idx, k] <- colSums(hi, na.rm=TRUE) / (b$k - 1) / col.n
      start <- start + m * m
    }
  }
  back <- order(g$o)

  return(list(sender.lo=sender.lo[back, , drop=FALSE], sender.hi=sender.hi[back, , drop=FALSE],
              receiver.lo=receiver.lo[back, , drop=FALSE], receiver.hi=receiver.hi[back, , drop=FALSE]))
}
