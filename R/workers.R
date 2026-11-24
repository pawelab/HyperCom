# Worker processes of prioritize_lr: R processes (PSOCK, from parallelly) that score batches of cell grids, each holding its own copy
# of the within-sample blocks of the influence matrix. The main process makes each batch's inputs and the priority terms.

#' State of the session or of a worker process
#'
#' `session` holds what is read once per session (the physical memory) and, in a worker process, its blocks and the sample of each
#' cell
#' @noRd
session <- new.env(parent=emptyenv())

#' Start the scoring of a run
#'
#' `start_pool` returns how the grids of a run are scored: in the main process (one process) or by worker processes
#'
#' This function takes the number of processes and returns an environment for pool_grids. Nothing is made until the first grids are
#' scored, so a run that takes every pair from its cache starts no worker.
#' @param workers An integer number of processes, 1 for the main process alone
#' @returns An environment with workers, cl (the cluster, NULL until started) and blocks (the main process's blocks, NULL until made)
#' @noRd
start_pool <- function(workers=1){
  pool <- new.env(parent=emptyenv())
  pool$workers <- as.integer(workers)
  pool$cl <- NULL
  pool$blocks <- NULL

  return(pool)
}

#' Stop the worker processes of a run
#'
#' `stop_pool` stops the cluster of a pool, if started, and drops its blocks
#' @param pool An environment from start_pool
#' @returns NULL, invisibly
#' @noRd
stop_pool <- function(pool){
  if(!is.null(pool$cl)){
    try(parallel::stopCluster(pool$cl), silent=TRUE)
    pool$cl <- NULL
  }
  pool$blocks <- NULL

  return(invisible(NULL))
}

#' Make the blocks of a run
#'
#' `pool_ready` makes the within-sample blocks of the influence matrix once per run: kept in the pool for one process, else handed to
#' each worker process, which then keeps them for the run
#'
#' This function starts the worker processes (parallelly::makeClusterPSOCK, with this session's library paths so that they load this
#' HyperCom), gives each a copy of the blocks and the sample of each cell (worker_setup) and drops the main process's copy. The
#' cluster is stopped by stop_pool, not by parallelly's autoStop: that finalizer closes the cluster's connections by number when the
#' stopped cluster is collected, which can close the connections of a later cluster that reuses those numbers.
#' @param pool An environment from start_pool
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A data frame from generate_metadata
#' @returns The pool, invisibly
#' @noRd
pool_ready <- function(pool, infMat, metadata){
  if(!is.null(pool$blocks) || !is.null(pool$cl)){
    return(invisible(pool))
  }
  g <- grid_cells(rownames(infMat), metadata)
  blocks <- grid_blocks(infMat, g)
  if(pool$workers == 1){
    pool$blocks <- blocks
    return(invisible(pool))
  }
  pool$cl <- parallelly::makeClusterPSOCK(pool$workers, rscript_libs=.libPaths(), autoStop=FALSE)
  parallel::clusterCall(pool$cl, worker_setup, blocks, g$sample)
  rm(blocks)

  return(invisible(pool))
}

#' Set up a worker process
#'
#' `worker_setup` keeps the blocks and the cells' samples in a worker process and lowers its priority
#'
#' This function runs in each worker process: it keeps the blocks for worker_grids and lowers the process's priority (below normal on
#' Windows, nice 5 elsewhere), so the computer stays responsive while every core scores grids.
#' @param blocks A list from grid_blocks
#' @param sample An integer vector of the sample of each cell, in the order of grid_cells
#' @returns NULL, invisibly
#' @noRd
worker_setup <- function(blocks, sample){
  session$blocks <- blocks
  session$sample <- sample
  try(tools::psnice(value=if(.Platform$OS.type == "windows") 15L else 5L), silent=TRUE)

  return(invisible(NULL))
}

#' Score a batch of grids in a worker process
#'
#' `worker_grids` returns the cell scores of a batch of grids from the worker's blocks
#' @param gi A list with A, B, lc and rc from grid_inputs
#' @returns A list from score_blocks, cells in the order of grid_cells
#' @noRd
worker_grids <- function(gi){
  return(score_blocks(session$blocks, gi$A, gi$B, session$sample, gi$lc, gi$rc))
}

#' Score batches of grids
#'
#' `pool_grids` returns the cell scores of each batch of grids
#'
#' This function takes a pool, the influence matrix, the metadata and batches of grids, at most one per worker, and returns the cell
#' scores of each batch (score_grids): in the main process from the pool's blocks, or one batch per worker process, the main process
#' making each batch's inputs (grid_inputs) and putting the cells back in the order of the metadata.
#' @param pool An environment from start_pool
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A data frame from generate_metadata
#' @param batches A list of data frames from pair_grids, with columns ligand and receptor
#' @returns A list of lists of two matrices, sender and receiver, one per batch
#' @noRd
pool_grids <- function(pool, infMat, metadata, batches){
  pool_ready(pool, infMat, metadata)
  if(is.null(pool$cl)){
    return(lapply(batches, function(b) score_grids(infMat, metadata, b$ligand, b$receptor, pool$blocks)))
  }
  inputs <- lapply(batches, function(b) grid_inputs(infMat, metadata, b$ligand, b$receptor))
  scores <- parallel::clusterApply(pool$cl, lapply(inputs, function(gi) gi[c("A", "B", "lc", "rc")]), worker_grids)

  return(Map(function(sr, gi) metadata_order(sr, gi$g), scores, inputs))
}
