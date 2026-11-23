# Memory mode: when the exact grid kernels' 16 bytes per within-sample cell pair do not fit, prioritize_lr scores from log-domain
# integer keys of the factored influence matrix kept as 8 x 8 tiles of each block's upper triangle (src/f2.h), in passes over the
# key range. Scores are within ~1e-11 of the exact kernels' (approved by J. Moy 2026-11-20); runs that fit are unchanged.

#' Find the memory a run may use
#'
#' `memory_budget` returns the bytes prioritize_lr may use
#'
#' This function takes the memory argument of prioritize_lr and returns it in bytes, or 90% of the computer's physical memory when it is
#' NULL. A cluster job should give its request: the scheduler may not stop a job that uses more, and the node may be shared.
#' @param memory A number of gigabytes (10^9 bytes), Inf for no limit, or NULL
#' @returns A number of bytes
#' @noRd
memory_budget <- function(memory=NULL){
  if(is.null(memory)){
    total <- physical_memory_cpp()
    return(if(is.finite(total) && total > 0) 0.9 * total else Inf)
  }

  return(memory * 1e9)
}

#' Estimate the memory of a run's grid kernels
#'
#' `kernel_memory` returns the bytes the grid kernels need, exact and in memory mode
#'
#' This function takes the number of within-sample cell pairs and returns the peak of each kernel: exact, 16 bytes per pair (blocks 8,
#' keys 6, bucket ids 2); memory mode, the tiles' 3 bytes per pair plus scoring's keys for a third of the pairs (3 passes) plus an
#' eighth as slack, at most 2^26 keys more, 6 bytes each. Each is raised by 10% and 1 GiB for the run's other objects (measured:
#' 30,000 cells in one sample peak at 14.8 GiB against 15.6 estimated; 88,900 cells in memory mode at 42.2 GiB against 41.5).
#' @param K A number of within-sample cell pairs (the sum of squared sample sizes)
#' @returns A named numeric vector of bytes: exact and tiles
#' @noRd
kernel_memory <- function(K){
  per.pair <- c(16, 3 + 6 * (K / 3 + min(K / 24, 2^26)) / K)

  return(stats::setNames(1.1 * per.pair * K + 2^30, c("exact", "tiles")))
}

#' Choose the grid kernels of a run
#'
#' `choose_kernel` returns how prioritize_lr scores grids: exact, or memory mode
#'
#' This function takes a named influence matrix, the metadata, the memory argument and verbose and returns the exact kernels when their
#' estimate (kernel_memory) fits the budget (memory_budget), else memory mode (a hypercom_influence only: its tiles are made from the
#' factors). Memory mode always scores in 3 passes: at 100,000 cells 2 passes were no faster and took 8.8 GiB more, and the one crash
#' seen (2026-11-23, not reproduced) was at 2. It stops, giving the estimate, when nothing fits. The tiles are made once, by the
#' first grids that need them, and kept in the returned environment for the rest of the run.
#' @inheritParams prioritize_lr
#' @returns An environment with type ("exact" or "tiles"), passes, need and budget (bytes) and lambda (NULL until made)
#' @noRd
choose_kernel <- function(infMat, metadata, memory=NULL, verbose=FALSE){
  cells.idx <- metadata$sample != "gene"
  K <- sum(as.numeric(table(metadata$sample[cells.idx]))^2)
  need <- kernel_memory(K)
  budget <- memory_budget(memory)
  kernel <- new.env(parent=emptyenv())
  kernel$budget <- budget
  kernel$lambda <- NULL
  if(need[["exact"]] <= budget){
    kernel$type <- "exact"
    kernel$passes <- 0L
    kernel$need <- need[["exact"]]
  } else if(inherits(infMat, "hypercom_influence") && need[["tiles"]] <= budget){
    kernel$type <- "tiles"
    kernel$passes <- 3L
    kernel$need <- need[["tiles"]]
  } else {
    stop("scoring needs about ", format_gb(need[["exact"]]), " of memory", if(inherits(infMat, "hypercom_influence"))
           paste0(" (", format_gb(need[["tiles"]]), " in memory mode)"), " but ", format_gb(budget), " is available",
         if(!inherits(infMat, "hypercom_influence")) ": memory mode needs an influence matrix made with method \"woodbury\"",
         "; give more with the memory argument, or split the cells into more samples")
  }
  if(verbose){
    message("Grid kernels need about ", format_gb(kernel$need), " of memory",
            if(kernel$type == "tiles") paste0(" (memory mode, ", kernel$passes, " passes; exact kernels would need ",
                                              format_gb(need[["exact"]]), ")"))
  }

  return(kernel)
}

#' Format bytes as gigabytes
#'
#' `format_gb` returns a number of bytes as a string of gigabytes
#' @param bytes A number
#' @returns A string
#' @noRd
format_gb <- function(bytes){
  return(paste0(format(signif(bytes / 1e9, 3), big.mark=","), " GB"))
}

#' Make or reuse the grids' influence values of a run
#'
#' `grid_blocks` returns what the grid kernels read: the within-sample blocks of F (exact) or the tiles of Λ (memory mode)
#'
#' This function takes a named influence matrix, the cells from grid_cells, the thread count and a kernel from choose_kernel (NULL for
#' the exact kernels) and returns the blocks (influence_blocks, or dense_blocks_cpp for a dense matrix) or the tiles (lambda_blocks),
#' made once per run and kept in the kernel.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param g A list from grid_cells
#' @param threads An integer number of threads
#' @param kernel An environment from choose_kernel, or NULL
#' @returns A numeric vector or a hypercom_lambda
#' @noRd
grid_blocks <- function(infMat, g, threads=1, kernel=NULL){
  if(!is.null(kernel) && kernel$type == "tiles"){
    if(is.null(kernel$lambda)){
      kernel$lambda <- lambda_blocks(infMat, g, threads, kernel$passes)
    }
    return(kernel$lambda)
  }
  if(inherits(infMat, "hypercom_influence")){
    return(influence_blocks(infMat, g, threads))
  }

  return(dense_blocks_cpp(infMat, g$idx, g$sample))
}

#' Make the log-domain tiles of a factored influence matrix
#'
#' `lambda_blocks` returns Λ = -ln F of every within-sample block, rounded to integers, in 8 x 8 tiles of each block's upper triangle
#'
#' This function takes a hypercom_influence, the cells from grid_cells, the thread count and the passes of the scoring kernel and
#' returns the tiles as 6-byte offsets from a base per block (influence_lambda_cpp: 3 bytes per within-sample cell pair), with what
#' the kernels need besides: the row sums D of the cells and the key precision p (2^-44 in ln).
#' @param infMat A hypercom_influence
#' @param g A list from grid_cells
#' @param threads An integer number of threads
#' @param passes An integer number of passes of the scoring kernel
#' @returns A hypercom_lambda: a list with bytes, base, toff, D, p, passes and entries
#' @noRd
lambda_blocks <- function(infMat, g, threads=1, passes=3L){
  if(sum(as.numeric(tabulate(g$sample))^2) >= 2^52){
    stop("too many cells per sample: the sum of squared sample sizes must stay below 2^52")
  }
  H <- infMat$H
  lb <- influence_lambda_cpp(H@p, H@i, H@x, nrow(H), infMat$rs, infMat$delta, infMat$Z, infMat$r, g$idx, g$sample, 44L, threads)
  lb$D <- infMat$rs[g$idx + 1]
  lb$p <- 44L
  lb$passes <- as.integer(passes)
  lb$entries <- sum(as.numeric(tabulate(g$sample))^2)

  return(structure(lb, class="hypercom_lambda"))
}

#' Number of within-sample entries of blocks or tiles
#'
#' `block_entries` returns the entries the grid kernels visit, from the F blocks or a hypercom_lambda
#' @param blocks A numeric vector or a hypercom_lambda
#' @returns A number
#' @noRd
block_entries <- function(blocks){
  if(inherits(blocks, "hypercom_lambda")){
    return(blocks$entries)
  }

  return(length(blocks))
}

#' Check that memory mode can score some grids
#'
#' `check_tiles` stops unless the influence values of the grids' genes are all positive and finite, as their logarithms need
#'
#' This function takes the values of the cells to the ligands (A) and of the receptors to the cells (B). A zero (a gene and a cell
#' not connected in the hypergraph) cannot be scored in memory mode, and memory mode runs only when the exact kernels do not fit, so the
#' run stops with their estimate.
#' @param A A numeric matrix
#' @param B A numeric matrix
#' @param entries A number of within-sample cell pairs
#' @returns NULL, invisibly
#' @noRd
check_tiles <- function(A, B, entries){
  if(!all(A > 0 & is.finite(A)) || !all(B > 0 & is.finite(B))){
    stop("some cells have an influence of 0 (not connected in the hypergraph) to a gene of these pairs, which memory mode cannot ",
         "score; the exact kernels need about ", format_gb(kernel_memory(entries)[["exact"]]), ": give that with the memory argument")
  }

  return(invisible(NULL))
}
