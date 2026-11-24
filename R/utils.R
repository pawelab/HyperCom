# Helpers shared by several steps: indexing of the node x node objects, thread counts, the output directory and the random number state

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

#' Read the dimensions of a node x node object
#'
#' `node_dim` returns the number of nodes of a hypercom_adjacency or hypercom_influence twice
#' @param x A hypercom_adjacency or a hypercom_influence
#' @returns An integer vector of two
#' @noRd
node_dim <- function(x){
  return(rep(nrow(x$H), 2L))
}

#' Read the node names of a node x node object
#'
#' `node_dimnames` returns the node names of a hypercom_adjacency or hypercom_influence twice, named node
#' @param x A hypercom_adjacency or a hypercom_influence
#' @returns A list of two vectors of strings
#' @noRd
node_dimnames <- function(x){
  return(list(node=rownames(x$H), node=rownames(x$H)))
}

#' Read the index of a node x node object
#'
#' `block_index` returns the rows and columns that x[i, j] selects
#'
#' This function takes the node names, the indices (missing when the call leaves them out, selecting every node), the number of
#' arguments of the call (without drop) and how the object is indexed, for the message when the call is not x[i, j].
#' @param nodes A vector of strings of the node names
#' @param i,j Indices of rows and columns
#' @param narg An integer number of arguments of the call
#' @param usage A string naming the object and how to index it
#' @returns A list with nodes, rows and cols
#' @noRd
block_index <- function(nodes, i, j, narg, usage){
  if(narg < 3){
    stop("index ", usage)
  }
  rows <- if(missing(i)) seq_along(nodes) else node_index(i, nodes)
  cols <- if(missing(j)) seq_along(nodes) else node_index(j, nodes)

  return(list(nodes=nodes, rows=rows, cols=cols))
}

#' Name a block of a node x node object
#'
#' `name_block` returns a block with its node names, dropping dimensions of extent one if asked
#' @param A A numeric matrix of the rows and columns of idx
#' @param idx A list from block_index
#' @param drop A boolean of whether to drop dimensions of extent one
#' @returns A matrix, or a vector when a dimension is dropped
#' @noRd
name_block <- function(A, idx, drop){
  dimnames(A) <- list(node=idx$nodes[idx$rows], node=idx$nodes[idx$cols])
  if(drop && any(dim(A) == 1)){
    A <- A[, , drop=TRUE]
  }

  return(A)
}

#' Choose the number of threads
#'
#' `get_threads` returns the number of threads used to score cell grids
#'
#' This function takes an optional thread count and returns it, or, when it is NULL, the available physical cores: hyperthreads and
#' big grids compete for the same caches. Never more than the available (logical) cores: a larger thread count warns. A thread count
#' below 1 stops.
#' @param threads An integer number of threads or NULL
#' @returns An integer number of threads
#' @noRd
get_threads <- function(threads=NULL){
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
#' `check_outdir` stops unless outdir is one string naming a directory, or NULL for writing no files
#'
#' This function takes the outdir argument of prioritize_lr or HyperCom and stops unless it is NULL or one string that does not name
#' an existing file.
#' @param outdir The outdir argument
#' @returns NULL, invisibly
#' @noRd
check_outdir <- function(outdir){
  if(is.null(outdir)){
    return(invisible(NULL))
  }
  if(!(is.character(outdir) && length(outdir) == 1 && !is.na(outdir))){
    stop("outdir must be one string naming the output directory, or NULL to write no files")
  }
  if(file.exists(outdir) && !dir.exists(outdir)){
    stop("outdir names a file, not a directory: ", outdir)
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
