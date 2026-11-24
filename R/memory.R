# Memory of prioritize_lr and score_hypercom: the memory a run may use, the estimate of a run's memory per R process, and the number
# of worker processes that fit (each holds its own copy of the within-sample blocks of the influence matrix)

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
  total <- c(physical_memory(), cgroup_memory())
  total <- total[is.finite(total) & total > 0]
  if(length(total) == 0){
    return(Inf)
  }

  return(0.9 * min(total))
}


#' Find the physical memory of the computer
#'
#' `physical_memory` returns the bytes of physical memory, NA when they cannot be read
#'
#' This function reads /proc/meminfo on Linux, sysctl hw.memsize on macOS and the Win32_ComputerSystem class (through PowerShell, about
#' 0.7 s) on Windows, once per session.
#' @returns A number of bytes, or NA
#' @noRd
physical_memory <- function(){
  if(is.null(session$physical.memory)){
    session$physical.memory <- tryCatch(read_physical_memory(), error=function(e) NA_real_, warning=function(w) NA_real_)
  }

  return(session$physical.memory)
}

#' Read the physical memory of the computer
#'
#' `read_physical_memory` returns the bytes of physical memory as the operating system reports them
#' @returns A number of bytes, or NA
#' @noRd
read_physical_memory <- function(){
  os <- Sys.info()[["sysname"]]
  if(os == "Linux"){
    line <- grep("^MemTotal:", readLines("/proc/meminfo", warn=FALSE), value=TRUE)
    return(as.numeric(gsub("[^0-9]", "", line[1])) * 1024)
  }
  if(os == "Darwin"){
    return(as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout=TRUE)[1]))
  }
  if(os == "Windows"){
    out <- system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
                   stdout=TRUE, stderr=FALSE)
    return(as.numeric(trimws(out[length(out)])))
  }

  return(NA_real_)
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

#' Bytes per within-sample cell pair of one R process that scores grids
#'
#' `pair.bytes` is the memory per within-sample cell pair of the process that scores grids: the blocks of the influence matrix (8) and
#' the peak of ranking one grid (grid_values and dplyr::percent_rank). Measured 2026-11-24 (one sample of 2,000 to 8,000 cells, peak
#' working set of worker processes): 73 bytes per pair above a worker's 0.4 GB
#' @noRd
pair.bytes <- 72

#' Memory of a worker process besides its blocks and grids
#'
#' `worker.bytes` is the memory of a worker process without its share of the cell pairs: R with HyperCom and its imports loaded, a
#' batch's inputs and cell scores (at most about 64 MB)
#' @noRd
worker.bytes <- 4e8

#' Estimate the memory of a run
#'
#' `run_memory` returns the bytes a run needs in the main R process and in each worker process
#'
#' This function takes the number of within-sample cell pairs and returns the estimates of scoring in the main process (pair.bytes per
#' pair, raised by 10%, plus 1 GiB for the run's other objects, so no run needs less than about 1.08 GB) and of each worker process
#' (pair.bytes per pair raised by 10%, plus worker.bytes). With workers the main process holds the blocks (8 bytes per pair) while it
#' hands them out.
#' @param K A number of within-sample cell pairs (the sum of squared sample sizes)
#' @returns A named numeric vector of bytes: main (scoring in the main process), handout (the main process while it hands the blocks to
#' the workers) and worker (each worker process)
#' @noRd
run_memory <- function(K){
  return(c(main=1.1 * pair.bytes * K + 2^30, handout=1.1 * 8 * K + 2^30, worker=1.1 * pair.bytes * K + worker.bytes))
}

#' Choose the number of R processes that score grids
#'
#' `fit_workers` returns how many worker processes fit in memory, 1 for scoring in the main process, or stops when that does not fit
#'
#' This function takes the within-sample cell pairs, the thread count and the memory.gb argument and returns the number of processes:
#' the thread count, lowered so that the main process and the workers fit the budget (memory_budget). Fewer than two workers that fit
#' means scoring in the main process alone. It stops, giving the estimate, when that does not fit either.
#' @param K A number of within-sample cell pairs
#' @param threads An integer number of processes asked for
#' @param memory.gb A number of gigabytes, Inf for no limit, or NULL for the computer's memory
#' @returns An integer number of processes
#' @noRd
fit_workers <- function(K, threads=1, memory.gb=NULL){
  need <- run_memory(K)
  budget <- memory_budget(memory.gb)
  if(need[["main"]] > budget){
    digits <- 3
    while(digits < 8 && format_gb(need[["main"]], digits) == format_gb(budget, digits)){
      digits <- digits + 1
    }
    stop("scoring needs about ", format_gb(need[["main"]], digits), " of memory but ", format_gb(budget, digits), " is available; give ",
         "more with the memory.gb argument, split the cells into more samples (this changes the results: cells are compared only ",
         "within their sample), or use HyperCom from branch master, whose compiled scoring needs about 4 times less")
  }
  fit <- floor((budget - need[["handout"]]) / need[["worker"]])
  workers <- min(threads, fit)
  if(workers < 2){
    return(1L)
  }

  return(as.integer(workers))
}

#' Plan the processes of a run
#'
#' `plan_workers` returns the number of R processes that score the grids of a run, reporting the memory they need
#'
#' This function takes the metadata, the thread count and the memory.gb argument and returns fit_workers for the run's within-sample
#' cell pairs; with verbose it reports the estimate, and the thread count asked for when fewer processes fit.
#' @param metadata A data frame from generate_metadata
#' @param threads An integer number of processes asked for
#' @param memory.gb A number of gigabytes, Inf for no limit, or NULL for the computer's memory
#' @param verbose A boolean of whether to report the estimate as a message
#' @returns An integer number of processes
#' @noRd
plan_workers <- function(metadata, threads=1, memory.gb=NULL, verbose=FALSE){
  K <- cell_pairs(table(metadata$sample[metadata$sample != "gene"]))
  workers <- fit_workers(K, threads, memory.gb)
  if(verbose){
    need <- run_memory(K)
    total <- if(workers == 1) need[["main"]] else need[["handout"]] + workers * need[["worker"]]
    message("Scoring needs at most about ", format_gb(total), " of memory",
            if(workers < threads) paste0(" (", workers, " of the ", threads, " processes asked for fit in ", format_gb(memory_budget(memory.gb)),
                                         ")"))
  }

  return(workers)
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
