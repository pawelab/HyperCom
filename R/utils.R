# Helpers shared by several steps: thread counts, node indexing, the output directory and the random number state

#' Find nodes by name, position or logical
#'
#' `node_index` returns the positions of indexed nodes
#'
#' This function takes an index as for a matrix (names, positive or negative positions, logical; a factor is read by name) and the
#' node names and returns the positions, stopping on a name or position that is not a node, on NA and on a logical index longer than
#' the nodes.
#' @param i A vector indexing nodes
#' @param nodes A vector of strings of the node names
#' @returns An integer vector of positions
#' @noRd
node_index <- function(i, nodes){
  if(anyNA(i)){
    stop("indices must not be NA")
  }
  if(is.logical(i) && length(i) > length(nodes)){
    stop("(subscript) logical subscript too long")
  }
  if(is.factor(i)){
    i <- as.character(i)
  }
  if(is.character(i)){
    idx <- match(i, nodes)
  } else {
    idx <- seq_along(nodes)[i]
  }
  if(anyNA(idx)){
    stop("subscript out of bounds")
  }

  return(idx)
}

#' Choose the number of threads
#'
#' `get_threads` returns the number of threads used to score cell grids
#'
#' This function takes the parallel switch and an optional thread count. An explicit thread count wins; otherwise one thread
#' without parallel, and with parallel the available physical cores: hyperthreads and big grids compete for the same caches. Never
#' more than the available (logical) cores: a larger thread count warns. A thread count below 1 stops.
#' @param parallel A boolean of whether to run in parallel
#' @param threads An integer number of threads or NULL
#' @returns An integer number of threads
#' @noRd
get_threads <- function(parallel=FALSE, threads=NULL){
  cores <- parallelly::availableCores()
  if(!is.null(threads)){
    threads <- suppressWarnings(as.numeric(threads))
    if(length(threads) != 1 || is.na(threads)){
      stop("threads must be one number")
    }
    if(threads < 1){
      stop("threads must be at least 1")
    }
    if(floor(threads) > cores){
      warning("threads=", threads, " is more than the ", cores, " available cores: using ", cores)
    }
    return(as.integer(min(floor(threads), cores)))
  }
  if(!parallel){
    return(1L)
  }
  cores <- min(cores, physical_cores(), na.rm=TRUE)

  return(as.integer(cores))
}

#' Count the physical cores
#'
#' `physical_cores` returns the number of physical cores of the machine, or NA when it cannot be read
#'
#' This function asks parallelly, which counts physical cores on Windows and macOS but gives the logical count on Linux; there it
#' counts the distinct (package, core) pairs of the CPUs in /sys/devices/system/cpu instead.
#' @returns An integer or NA
#' @noRd
physical_cores <- function(){
  if(Sys.info()[["sysname"]] != "Linux"){
    return(suppressWarnings(as.integer(parallelly::availableCores(logical=FALSE))))
  }
  cpus <- list.files("/sys/devices/system/cpu", pattern="^cpu[0-9]+$", full.names=TRUE)
  ids <- vapply(cpus, function(d){
    tryCatch(paste(readLines(file.path(d, "topology", "physical_package_id"), warn=FALSE)[1],
                   readLines(file.path(d, "topology", "core_id"), warn=FALSE)[1]),
             error=function(e) NA_character_, warning=function(w) NA_character_)
  }, character(1))
  ids <- ids[!is.na(ids)]
  if(length(ids) == 0){
    return(NA_integer_)
  }

  return(length(unique(ids)))
}

#' Check the output directory argument
#'
#' `check_outdir` stops unless outdir is one string, or NULL where it means writing no files
#'
#' This function takes the outdir argument of prioritize_lr or HyperCom and stops unless it is one string (or NULL when null.ok). A
#' cluster or a number there usually comes from a call by position written for the cl argument, which 0.2.0 removed.
#' @param outdir The outdir argument
#' @param null.ok A boolean of whether NULL (write no files) is allowed
#' @returns NULL, invisibly
#' @noRd
check_outdir <- function(outdir, null.ok=FALSE){
  if(is.null(outdir) && null.ok){
    return(invisible(NULL))
  }
  if(!(is.character(outdir) && length(outdir) == 1 && !is.na(outdir))){
    stop("outdir must be one string naming the output directory", if(null.ok) " (or NULL to write no files)",
         "; the cl argument before it was removed in 0.2.0, so drop the value given for cl from calls by position")
  }

  return(invisible(NULL))
}

#' Restore the state of the random number generator
#'
#' `restore_seed` puts back the random number state saved before a set.seed
#'
#' This function takes the value of .Random.seed saved before set.seed (NULL if there was none) and restores it, so that a seed
#' argument does not change the random numbers of the calling session.
#' @param old.seed An integer vector saved from .Random.seed, or NULL
#' @returns NULL, invisibly
#' @noRd
restore_seed <- function(old.seed){
  if(is.null(old.seed)){
    if(exists(".Random.seed", envir=globalenv(), inherits=FALSE)){
      rm(".Random.seed", envir=globalenv())
    }
  } else {
    assign(".Random.seed", old.seed, envir=globalenv())
  }

  return(invisible(NULL))
}
