# Grid kernels of prioritize_lr and score_hypercom: the memory a run may use, and the choice between the exact kernels (the
# within-sample blocks of the influence matrix) and memory mode (integer keys of -ln F kept as 8 x 8 tiles of each block's upper
# triangle, src/f2.h, ranked in passes over the key range; scores within ~1e-11 of the exact kernels')

#' Number of passes of memory mode
#'
#' `memory.passes` is the number of passes over the key range in which memory mode ranks each grid
#' @noRd
memory.passes <- 3L

#' Find the memory a run may use
#'
#' `memory_budget` returns the bytes prioritize_lr may use
#'
#' This function takes the memory.gb argument of prioritize_lr and returns it in bytes, or, when it is NULL, 90% of the computer's
#' physical memory or of the memory limit of the process's control group on Linux (a container or a cluster job), whichever is
#' smaller.
#' @param memory.gb A number of gigabytes (10^9 bytes), Inf for no limit, or NULL
#' @returns A number of bytes
#' @noRd
memory_budget <- function(memory.gb=NULL){
  if(!is.null(memory.gb)){
    return(memory.gb * 1e9)
  }
  total <- c(physical_memory_cpp(), cgroup_memory())
  total <- total[is.finite(total) & total > 0]
  if(length(total) == 0){
    return(Inf)
  }

  return(0.9 * min(total))
}

#' Read the memory limit of the process's control group
#'
#' `cgroup_memory` returns the smallest memory limit of the control groups of the process, Inf when there is none
#'
#' This function reads the process's control groups (cgroup v2 and the memory controller of v1) and returns the smallest limit set
#' on its group or any group above it ("max", or v1's huge number, means no limit). Where there are no control groups (Windows,
#' macOS) it returns Inf.
#' @param proc A string naming the file that lists the process's control groups
#' @param root A string naming the folder the control groups are mounted on
#' @returns A number of bytes, or Inf
#' @noRd
cgroup_memory <- function(proc="/proc/self/cgroup", root="/sys/fs/cgroup"){
  lines <- character(0)
  if(file.exists(proc)){
    lines <- tryCatch(readLines(proc, warn=FALSE), error=function(e) character(0))
  }
  v1.pattern <- "^[0-9]+:([^:]*,)?memory(,[^:]*)?:"
  v2 <- sub("^0::", "", lines[startsWith(lines, "0::")])
  v1 <- sub(v1.pattern, "", grep(v1.pattern, lines, value=TRUE))
  files <- c(file.path(root, cgroup_paths(v2), "memory.max"), file.path(root, "memory", cgroup_paths(v1), "memory.limit_in_bytes"),
             file.path(root, "memory.max"), file.path(root, "memory", "memory.limit_in_bytes"))
  limits <- vapply(unique(files), read_limit, numeric(1))

  return(min(c(limits, Inf)))
}

#' List a control group and the groups above it
#'
#' `cgroup_paths` returns the paths of control groups and of every group above them
#'
#' This function takes paths of control groups, e.g. "/a/b/c", and returns them with their parents below the root, e.g. "a", "a/b" and
#' "a/b/c".
#' @param paths A vector of strings of control group paths
#' @returns A vector of strings
#' @noRd
cgroup_paths <- function(paths){
  out <- character(0)
  for(p in paths){
    parts <- strsplit(p, "/", fixed=TRUE)[[1]]
    parts <- parts[parts != ""]
    out <- c(out, vapply(seq_along(parts), function(k) paste(parts[seq_len(k)], collapse="/"), character(1)))
  }

  return(out)
}

#' Read a memory limit file
#'
#' `read_limit` returns the number of bytes in a control group's memory limit file, Inf for no limit
#'
#' This function reads the first line of the file; a missing or unreadable file, "max" and numbers of 2^60 or more (cgroup v1 writes
#' a huge number for no limit) give Inf.
#' @param f A string naming the file
#' @returns A number of bytes, or Inf
#' @noRd
read_limit <- function(f){
  if(!file.exists(f)){
    return(Inf)
  }
  value <- tryCatch(suppressWarnings(as.numeric(readLines(f, n=1, warn=FALSE))), error=function(e) NA_real_)
  if(length(value) != 1 || is.na(value) || value <= 0 || value >= 2^60){
    return(Inf)
  }

  return(value)
}

#' Count the within-sample cell pairs
#'
#' `cell_pairs` returns the number of within-sample cell pairs, the sum of squared sample sizes
#' @param sizes A numeric vector of the number of cells of each sample
#' @returns A number
#' @noRd
cell_pairs <- function(sizes){
  return(sum(as.numeric(sizes)^2))
}

#' Check that the within-sample cell pairs can be indexed
#'
#' `check_cell_pairs` stops when the sum of squared sample sizes reaches 2^52
#'
#' This function takes the sample sizes and stops before any block is made when their within-sample cell pairs cannot all be indexed
#' exactly by doubles.
#' @param sizes A numeric vector of the number of cells of each sample
#' @returns NULL, invisibly
#' @noRd
check_cell_pairs <- function(sizes){
  if(cell_pairs(sizes) >= 2^52){
    stop("too many cells per sample: the sum of squared sample sizes must stay below 2^52")
  }

  return(invisible(NULL))
}

#' Estimate the memory of a run's grid kernels
#'
#' `kernel_memory` returns the bytes the grid kernels need, exact and in memory mode
#'
#' This function takes the number of within-sample cell pairs and returns the peak of each kernel: exact, 16 bytes per pair (blocks 8,
#' keys 6, bucket ids 2); memory mode, the tiles' 3 bytes per pair plus scoring's keys for a third of the pairs (3 passes) plus an
#' eighth as slack, at most 2^26 keys more, 6 bytes each. Each is raised by 10% and 1 GiB for the run's other objects, so no run needs
#' less than about 1.08 GB; the estimates lie above the peaks measured from 1,000 to 100,000 cells.
#' @param K A number of within-sample cell pairs (the sum of squared sample sizes)
#' @returns A named numeric vector of bytes: exact and memory
#' @noRd
kernel_memory <- function(K){
  kc <- kernel_constants_cpp()
  per.pair <- c(16, 3 + 6 * (K / memory.passes + min(K / (kc$pass.slack.div * memory.passes), kc$pass.slack.max)) / K)

  return(stats::setNames(1.1 * per.pair * K + 2^30, c("exact", "memory")))
}

#' Choose the kind of grid kernel that fits in memory
#'
#' `kernel_type` returns "exact" when the exact kernels fit in memory, else "memory" (memory mode) when it fits, else stops
#'
#' This function takes the within-sample cell pairs, the memory.gb argument and whether the influence matrix is factored (memory mode
#' reads the factors) and returns the kernel type, or stops with the estimates and what to do.
#' @param K A number of within-sample cell pairs
#' @param memory.gb A number of gigabytes, Inf for no limit, or NULL for the computer's memory
#' @param factored A boolean of whether the influence matrix is a hypercom_influence
#' @returns A string, "exact" or "memory"
#' @noRd
kernel_type <- function(K, memory.gb=NULL, factored=TRUE){
  need <- kernel_memory(K)
  budget <- memory_budget(memory.gb)
  if(need[["exact"]] <= budget){
    return("exact")
  }
  if(factored && need[["memory"]] <= budget){
    return("memory")
  }

  digits <- 3
  while(digits < 8 && format_gb(need[["exact"]], digits) == format_gb(budget, digits)){
    digits <- digits + 1
  }
  stop("scoring needs about ", format_gb(need[["exact"]], digits), " of memory", if(factored) paste0(" (", format_gb(need[["memory"]], digits),
       " in memory mode)"), " but ", format_gb(budget, digits), " is available",
       if(!factored) paste0(": memory mode needs an influence matrix factored by method \"woodbury\" (data with no more cells than ",
                            "interactions of two or more genes get a full matrix)"),
       "; give more with the memory.gb argument, or split the cells into more samples (this changes the results: cells are ",
       "compared only within their sample)")
}

#' Choose the grid kernels of a run
#'
#' `choose_kernel` returns how prioritize_lr and score_hypercom score grids: exact, or memory mode
#'
#' This function takes a named influence matrix, the metadata, the memory.gb argument and verbose and returns the exact kernels when their
#' estimate (kernel_memory) fits the budget (memory_budget), else memory mode (a hypercom_influence only: its tiles are made from the
#' factors), which needs a connected hypergraph (check_connected). It stops, giving the estimate, when nothing fits. The tiles are made
#' once, by the first grids that need them, and kept in the returned environment for the rest of the run.
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A data frame from generate_metadata
#' @param memory.gb A number of gigabytes, Inf for no limit, or NULL for the computer's memory
#' @param verbose A boolean of whether to report the estimate as a message
#' @returns An environment with type ("exact" or "memory"), need (bytes) and lambda (the tiles, NULL until made)
#' @noRd
choose_kernel <- function(infMat, metadata, memory.gb=NULL, verbose=FALSE){
  K <- cell_pairs(table(metadata$sample[metadata$sample != "gene"]))
  need <- kernel_memory(K)
  kernel <- new.env(parent=emptyenv())
  kernel$type <- kernel_type(K, memory.gb, inherits(infMat, "hypercom_influence"))
  kernel$need <- need[[kernel$type]]
  kernel$lambda <- NULL
  if(kernel$type == "memory"){
    check_connected(infMat, need[["exact"]])
  }
  if(verbose){
    message("Scoring needs at most about ", format_gb(kernel$need), " of memory",
            if(kernel$type == "memory") paste0(" (memory mode; exact scoring would need ", format_gb(need[["exact"]]), ")"))
  }

  return(kernel)
}

#' Check that a factored influence matrix has no zero
#'
#' `check_connected` stops unless the hypergraph of a hypercom_influence is connected
#'
#' This function takes a hypercom_influence and the memory of the exact kernels and stops when some nodes have no path of hyperedges
#' to the others: their influence on each other is 0, whose logarithm memory mode cannot take. generate_hyperedges keeps only the
#' largest connected part, so this happens only with hyperedges made another way.
#' @param infMat A hypercom_influence
#' @param exact A number of bytes the exact kernels need
#' @returns NULL, invisibly
#' @noRd
check_connected <- function(infMat, exact){
  H <- infMat$H
  label <- incidence_components(H@i + 1L, rep(seq_len(ncol(H)), diff(H@p)), nrow(H), ncol(H))
  if(any(label != 1L)){
    stop("memory mode needs a connected hypergraph, but ", sum(label != 1L), " of its ", nrow(H), " nodes have no path to the others ",
         "(influence 0), which memory mode cannot score: make the hyperedges with generate_hyperedges, which keeps the largest connected ",
         "part, or give the exact kernels' ", format_gb(exact), " with the memory.gb argument")
  }

  return(invisible(NULL))
}

#' Format bytes as gigabytes
#'
#' `format_gb` returns a number of bytes as a string of gigabytes
#' @param bytes A number
#' @param digits An integer number of significant digits
#' @returns A string
#' @noRd
format_gb <- function(bytes, digits=3){
  return(paste0(format(signif(bytes / 1e9, digits), big.mark=","), " GB"))
}

#' Make or reuse the grids' influence values of a run
#'
#' `grid_blocks` returns what the grid kernels read: the within-sample blocks of F (exact) or the tiles of Λ (memory mode)
#'
#' This function takes a named influence matrix, the cells from grid_cells, the thread count and a kernel from choose_kernel (NULL for
#' the exact kernels) and returns the blocks (influence_blocks, or dense_blocks_cpp for a dense matrix) or the tiles (lambda_blocks,
#' made once per run and kept in the kernel).
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param g A list from grid_cells
#' @param threads An integer number of threads
#' @param kernel An environment from choose_kernel, or NULL
#' @returns A numeric vector or a hypercom_lambda
#' @noRd
grid_blocks <- function(infMat, g, threads=1, kernel=NULL){
  if(!is.null(kernel) && kernel$type == "memory"){
    if(is.null(kernel$lambda)){
      kernel$lambda <- lambda_blocks(infMat, g, threads)
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
lambda_blocks <- function(infMat, g, threads=1, passes=memory.passes){
  key.bits <- kernel_constants_cpp()$key.bits
  sizes <- tabulate(g$sample)
  check_cell_pairs(sizes)
  H <- infMat$H
  lb <- influence_lambda_cpp(H@p, H@i, H@x, nrow(H), infMat$rs, infMat$delta, infMat$Z, infMat$r, g$idx, g$sample, key.bits, threads)
  lb$D <- infMat$rs[g$idx + 1]
  lb$p <- key.bits
  lb$passes <- as.integer(passes)
  lb$entries <- cell_pairs(sizes)

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
#' that no path of hyperedges joins) cannot be scored in memory mode, and memory mode runs only when the exact kernels do not fit, so
#' the run stops with their estimate. check_connected makes this stop unreachable for runs on a connected hypergraph.
#' @param A A numeric matrix
#' @param B A numeric matrix
#' @param entries A number of within-sample cell pairs
#' @returns NULL, invisibly
#' @noRd
check_tiles <- function(A, B, entries){
  if(!all(A > 0 & is.finite(A)) || !all(B > 0 & is.finite(B))){
    stop("some cells have an influence of 0 on a gene of these pairs (no path of hyperedges joins them), which memory mode cannot ",
         "score; the exact kernels need about ", format_gb(kernel_memory(entries)[["exact"]]), ": give that with the memory.gb argument")
  }

  return(invisible(NULL))
}
