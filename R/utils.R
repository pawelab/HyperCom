# Helpers shared by several steps: node indexing, the cluster, the output directory and the random number state

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

#' Check the output directory argument
#'
#' `check_outdir` stops unless outdir is one string, or NULL where it means writing no files
#'
#' This function takes the outdir argument of prioritize_lr or HyperCom and stops unless it is one string (or NULL when null.ok).
#' @param outdir The outdir argument
#' @param null.ok A boolean of whether NULL (write no files) is allowed
#' @returns NULL, invisibly
#' @noRd
check_outdir <- function(outdir, null.ok=FALSE){
  if(is.null(outdir) && null.ok){
    return(invisible(NULL))
  }
  if(!(is.character(outdir) && length(outdir) == 1 && !is.na(outdir))){
    stop("outdir must be one string naming the output directory", if(null.ok) " (or NULL to write no files)")
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

#' Choose the arguments of the default cluster
#'
#' `cluster_args` returns the arguments of parallelly::makeClusterPSOCK for the cluster prioritize_lr makes when none is given
#'
#' This function returns, inside an SGE parallel job (PE_HOSTFILE set), the job's workers started with qrsh, which can span several
#' nodes; elsewhere, one local worker per available core.
#' @returns A list of arguments of parallelly::makeClusterPSOCK
#' @noRd
cluster_args <- function(){
  if(nzchar(Sys.getenv("PE_HOSTFILE"))){
    return(list(parallelly::availableWorkers(), rshcmd="qrsh", rshopts=c("-inherit", "-nostdin", "-V"), autoStop=TRUE, outfile=""))
  }

  return(list(parallelly::availableCores(), autoStop=TRUE, outfile=""))
}
