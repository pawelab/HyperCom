# R versions of the compiled parts of HyperCom (src/ on branch master): the adjacency from the incidence, the influence matrix (lu,
# chol, woodbury factors and blocks of them) and the scoring of cell grids, which see only the within-sample cell x cell blocks.
# Each stands in for master's function of the same name without _cpp, so the other files stay close to master's.

#' Make the upper triangle of an adjacency from its incidence
#'
#' `adjacency` returns the slots of the upper triangle of H H' with an explicit zero diagonal, the dsCMatrix of HyperCom 0.0.0.9000
#'
#' This function takes a node x hyperedge incidence and returns the upper triangle of H H' with the diagonal set to explicit zeros, and
#' which nodes share a hyperedge with another node (the others are left out, as HyperCom 0.0.0.9000 left them out). Entries are sums
#' of products of integer counts, exact in any order.
#' @param H A sparse node x hyperedge matrix of counts
#' @returns A list with i, p and x (the slots of the kept nodes' upper triangle) and keep (a logical vector, one per node of H)
#' @noRd
adjacency <- function(H){
  A <- methods::as(Matrix::tcrossprod(H), "symmetricMatrix")
  A <- methods::as(methods::as(A, "CsparseMatrix"), "dMatrix")
  if(A@uplo != "U"){
    A <- Matrix::t(A)
  }
  col <- rep(seq_len(ncol(A)) - 1L, diff(A@p))
  diagonal <- A@i == col
  keep <- tabulate(c(A@i[!diagonal & A@x != 0], col[!diagonal & A@x != 0]) + 1L, nbins=nrow(A)) > 0
  A@x[diagonal] <- 0
  if(!all(keep)){
    A <- A[keep, keep]
  }

  return(list(i=A@i, p=A@p, x=A@x, keep=keep))
}

#' Solve for the influence matrix as HyperCom 0.0.0.9000 did
#'
#' `solve_influence` returns the dense influence matrix r (I - (1 - r) adj / rs)^-1 with the numbers of HyperCom 0.0.0.9000
#'
#' This function takes an adjacency matrix, its row sums and the restart probability and returns the influence matrix with the same
#' arithmetic as HyperCom 0.0.0.9000 (R's solve), writing I - W without a dense identity: four node x node matrices at once instead of
#' six.
#' @param adj An adjacency matrix: a base matrix or a sparse Matrix
#' @param rs A numeric vector of the row sums of adj
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A dense matrix
#' @noRd
solve_influence <- function(adj, rs, r){
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
#' adj is not symmetric or S is not positive definite or computationally singular.
#' @inheritParams solve_influence
#' @returns A dense matrix
#' @noRd
solve_influence_chol <- function(adj, rs, r){
  S <- as.matrix(adj)
  if(!isSymmetric(unname(S), tol=0)){
    stop("method \"chol\" needs a symmetric adjacency matrix (use method=\"lu\")")
  }
  s <- 1 / sqrt(rs)
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
#' `factor_spd` returns the upper Cholesky factor of a matrix, stopping when it is not positive definite or computationally singular
#'
#' This function takes a symmetric matrix and returns chol's upper triangular factor; it stops when the matrix is not positive definite
#' or its reciprocal condition number (estimated from the factor) is below the machine epsilon, e.g. for a restart probability near 0.
#' @param S A symmetric matrix
#' @returns An upper triangular matrix
#' @noRd
factor_spd <- function(S){
  R <- tryCatch(chol(S), error=function(e) NULL)
  if(is.null(R)){
    stop("the system is not positive definite (use method=\"lu\")")
  }
  rc <- rcond(R, triangular=TRUE)^2
  if(!(rc >= .Machine$double.eps)){
    stop("system is computationally singular: reciprocal condition number = ", signif(rc, 6))
  }

  return(R)
}

#' Factor the influence matrix through the hyperedges
#'
#' `influence_factors_r` returns Delta and Z of the factored influence matrix
#'
#' This function takes the node x hyperedge incidence, the row sums of the adjacency and the restart probability and returns the
#' factors of F = r Delta^-1 D + r (1 - r) Y Z. The adjacency is H H' with a zero diagonal, so K = D - (1 - r) adj = Delta - (1 - r) H H'
#' with the diagonal Delta = D + (1 - r) diag(H H'), and by the Woodbury identity F = r K^-1 D with Y = H / Delta, C = I - (1 - r) H' Y
#' (hyperedge x hyperedge, symmetric positive definite) and Z = C^-1 Y' D. Cost: hyperedges cubed plus hyperedges times the nonzeros
#' of H, instead of nodes cubed.
#' @param H A sparse node x hyperedge matrix of counts
#' @param rs A numeric vector of the row sums of the adjacency
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A list with delta (a numeric vector, one per node) and Z (a dense hyperedge x node matrix)
#' @noRd
influence_factors_r <- function(H, rs, r){
  delta <- rs + (1 - r) * Matrix::rowSums(H^2)
  Y <- H / delta
  C <- as.matrix(Matrix::crossprod(H, Y))
  C <- diag(nrow(C)) - (1 - r) * C
  R <- factor_spd(C)
  rm(C)
  Z <- as.matrix(chol2inv(R) %*% Matrix::t(Y * rs))
  dimnames(Z) <- NULL

  return(list(delta=as.numeric(delta), Z=Z))
}

#' Make the within-sample blocks of an influence matrix
#'
#' `grid_blocks` returns every within-sample cell x cell block of an influence matrix
#'
#' This function takes a named influence matrix and the cells from grid_cells and returns one block per sample (samples ascending,
#' cells of a sample in the order given): made from the factors of a hypercom_influence, copied from a dense matrix. They hold the sum
#' of squared sample sizes of doubles.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param g A list from grid_cells
#' @returns A list of dense matrices, one per sample
#' @noRd
grid_blocks <- function(infMat, g){
  check_cell_pairs(tabulate(g$sample))
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

#' Count the entries of within-sample blocks
#'
#' `block_entries` returns the within-sample cell pairs of the blocks from grid_blocks
#' @param blocks A list of dense matrices, or NULL
#' @returns A number
#' @noRd
block_entries <- function(blocks){
  return(sum(as.numeric(lengths(blocks))))
}

#' Make the transfer grid of a ligand and receptor in each sample
#'
#' `grid_values` returns the entries of a pair's cell x cell grid within each sample, as HyperCom 0.0.0.9000 ranks them
#'
#' This function takes the within-sample blocks of the influence matrix, the influence of the cells on the ligand (one value per cell)
#' and of the receptor on the cells, and the cells of each sample, and returns all within-sample entries A[i] B[j] / F[i, j] one after
#' the other (each block by column) with the diagonal, zeros and NaN as NA: the entries of HyperCom 0.0.0.9000's n x n grid that it did
#' not set to NA, as it computed them. With a and b NULL the grid is the influence matrix itself (the network score). One sample's
#' block is returned without a copy; several are written into one vector.
#' @param blocks A list of blocks from grid_blocks
#' @param a A numeric vector of the influence of each cell on the ligand, or NULL
#' @param b A numeric vector of the influence of the receptor on each cell, or NULL
#' @param samples A list of the positions of each sample's cells
#' @returns A numeric vector
#' @noRd
grid_values <- function(blocks, a, b, samples){
  vals <- if(length(samples) > 1) numeric(block_entries(blocks)) else NULL
  start <- 0
  for(s in seq_along(samples)){
    idx <- samples[[s]]
    if(is.null(a)){
      trans <- blocks[[s]]
    } else {
      trans <- outer(a[idx], b[idx]) / blocks[[s]]
    }
    diag(trans) <- 0
    trans[trans == 0] <- NA
    if(is.null(vals)){
      dim(trans) <- NULL
      return(trans)
    }
    vals[(start + 1):(start + length(trans))] <- trans
    start <- start + length(trans)
  }

  return(vals)
}

#' Calculate sender and receiver scores of cell grids from the within-sample blocks
#'
#' `score_blocks` returns the sender and receiver scores of the cells for each grid, cells in the order of grid_cells
#'
#' This function takes the within-sample blocks, the values of the cells to the distinct ligands (A, cells x ligands) and of the
#' distinct receptors to the cells (B, receptors x cells), the sample of each cell, and the column of A and row of B of each grid
#' (0-based, -1 for the network grid). For each grid it percent-ranks all entries of the transfer grid between cells of the same sample
#' together (dplyr::percent_rank, min ranks over the entries that are not NA) and averages the ranks of each cell as sender (its row)
#' and as receiver (its column), as HyperCom 0.0.0.9000 did on the whole grid, whose entries between samples it left out. With 2^24
#' within-sample cell pairs or more (a grid takes seconds), garbage is collected after each grid, as R would otherwise let several
#' grids' temporaries pile up.
#' @param blocks A list from grid_blocks
#' @param A A numeric matrix, cells x distinct ligands
#' @param B A numeric matrix, distinct receptors x cells
#' @param sample An integer vector of the sample of each cell (from 1)
#' @param lc An integer vector of the 0-based column of A of each grid, -1 for the network grid
#' @param rc An integer vector of the 0-based row of B of each grid, -1 for the network grid
#' @returns A list of two matrices, sender and receiver, with one row per cell and one column per grid
#' @noRd
score_blocks <- function(blocks, A, B, sample, lc, rc){
  samples <- unname(split(seq_along(sample), sample))
  big <- block_entries(blocks) >= 2^24
  sender <- receiver <- matrix(NA_real_, length(sample), length(lc))
  for(k in seq_along(lc)){
    if(lc[k] < 0){
      v <- grid_values(blocks, NULL, NULL, samples)
    } else {
      v <- grid_values(blocks, A[, lc[k] + 1L], B[rc[k] + 1L, ], samples)
    }
    v <- dplyr::percent_rank(v)
    if(length(samples) == 1){
      dim(v) <- rep(length(sample), 2)
      sender[, k] <- rowMeans(v, na.rm=TRUE)
      receiver[, k] <- colMeans(v, na.rm=TRUE)
    } else {
      start <- 0
      for(idx in samples){
        m <- length(idx)
        R <- v[(start + 1):(start + m * m)]
        dim(R) <- c(m, m)
        sender[idx, k] <- rowMeans(R, na.rm=TRUE)
        receiver[idx, k] <- colMeans(R, na.rm=TRUE)
        start <- start + m * m
      }
    }
    rm(v)
    if(big){
      gc()
    }
  }

  return(list(sender=sender, receiver=receiver))
}
