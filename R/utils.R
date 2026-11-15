# Helpers shared by several steps: thread counts and node indexing

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
  cores <- min(cores, parallelly::availableCores(logical=FALSE), na.rm=TRUE)

  return(as.integer(cores))
}
