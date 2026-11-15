real.batch_pairs <- HyperCom:::batch_pairs
real.score_grids <- HyperCom:::score_grids

#' Stop a run after some batches of grids
#'
#' `stop_after` evaluates an expression in which grid scoring stops with an error after a number of batches
#'
#' This function takes a number of batches and an expression, scores every ligand in its own batch so that the fixture has many
#' batches, and lets score_grids stop with "run stopped" once the given number of batches is scored, as a time limit would.
#' @param batches An integer number of batches scored before the stop
#' @param expr An expression running prioritize_lr or HyperCom
#' @returns The value of the expression (never reached once the stop comes)
#' @noRd
stop_after <- function(batches, expr){
  calls <- 0

  return(testthat::with_mocked_bindings(expr, batch_pairs=function(ligands, receptors, size) real.batch_pairs(ligands, receptors, 1),
                                        score_grids=function(...){
                                          calls <<- calls + 1
                                          if(calls > batches){
                                            stop("run stopped")
                                          }
                                          return(real.score_grids(...))
                                        }, .package="HyperCom"))
}

#' Count the grids a run scores
#'
#' `count_grids` returns the value of an expression and the number of grids scored while it ran
#'
#' This function takes an expression running prioritize_lr and counts the grids passed to score_grids.
#' @param expr An expression running prioritize_lr or HyperCom
#' @returns A list with value and grids
#' @noRd
count_grids <- function(expr){
  grids <- 0
  value <- testthat::with_mocked_bindings(expr, score_grids=function(infMat, metadata, ligands, receptors, ...){
    grids <<- grids + length(ligands)
    return(real.score_grids(infMat, metadata, ligands, receptors, ...))
  }, .package="HyperCom")

  return(list(value=value, grids=grids))
}

test_that("a stopped run continues from the cache, scores only the rest and gives master's priorities and p-values", {
  skip_on_cran()
  outdir <- file.path(withr::local_tempdir(), "6.priority")
  cache.dir <- file.path(outdir, "cache")

  expect_error(stop_after(3, run_priority(significance=TRUE, seed=1, outdir=outdir)), "run stopped")
  expect_setequal(list.files(cache.dir), c("fingerprint.rds", "background.rds", "terms.rds"))
  cached <- readRDS(file.path(cache.dir, "terms.rds"))$terms
  expect_gt(nrow(cached), 0)
  expect_lt(nrow(cached), nrow(fx$output$lrs) + 1000)

  full <- count_grids(run_priority(significance=TRUE, seed=1))
  msgs <- character(0)
  resumed <- count_grids(withCallingHandlers(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata,
                                                           significance=TRUE, seed=1, threads=test_threads(2), outdir=outdir),
                                             message=function(m){
                                               msgs <<- c(msgs, conditionMessage(m))
                                               invokeRestart("muffleMessage")
                                             }))

  expect_master(resumed$value, fx$output$priority.sig)
  expect_lt(resumed$grids, full$grids)
  expect_equal(sum(grepl(paste0("Continuing from the cache in .*: [1-9][0-9]* of ", nrow(fx$output$lrs) + 1000, " pairs already scored"),
                         msgs)), 1)
  expect_false(dir.exists(cache.dir))
  expect_identical(list.files(outdir), "priority.1.csv")
})

test_that("without seed, a stopped run continues with the background pairs it drew", {
  skip_on_cran()
  outdir <- withr::local_tempdir()

  expect_error(withr::with_seed(7, stop_after(3, run_priority(significance=TRUE, outdir=outdir))), "run stopped")
  resumed <- withr::with_seed(99, run_priority(significance=TRUE, outdir=outdir))

  expect_identical(resumed, withr::with_seed(7, run_priority(significance=TRUE)))
})

test_that("weights and threads may change when a run continues, also on a factored influence matrix", {
  outdir <- withr::local_tempdir()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  expect_master(run_priority(weights=c(1, 0.5, 0), threads=test_threads(2), outdir=outdir), fx$output$priority)

  expect_error(stop_after(2, run_priority(adj=adj, infMat=infMat, threads=test_threads(2), outdir=outdir)), "run stopped")
  expect_identical(run_priority(adj=adj, infMat=infMat, weights=c(1, 0.5), outdir=outdir),
                   run_priority(adj=adj, infMat=infMat, weights=c(1, 0.5)))
})

test_that("a run stopped after scoring takes every pair from the cache and scores no grid", {
  outdir <- withr::local_tempdir()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))

  first <- testthat::with_mocked_bindings(run_priority(adj=adj, infMat=infMat, outdir=outdir), close_cache=function(dir) invisible(NULL),
                                          .package="HyperCom")
  expect_true(file.exists(file.path(outdir, "cache", "terms.rds")))
  again <- count_grids(run_priority(adj=adj, infMat=infMat, outdir=outdir))

  expect_identical(again$value, first)
  expect_equal(again$grids, 0)
  expect_false(dir.exists(file.path(outdir, "cache")))
})

test_that("a cache made from other inputs, or damaged, is moved to cache.old with a warning and can be moved back", {
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "cache")
  old.dir <- file.path(outdir, "cache.old")

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  expect_warning(priority <- run_priority(metadata=fx$output$metadata.none, outdir=outdir),
                 "made from other inputs, or is damaged, and is not used; it is moved to .*cache.old")
  expect_master(priority, fx$output$priority.none)
  expect_false(dir.exists(cache.dir))
  expect_true(file.rename(old.dir, cache.dir))
  resumed <- count_grids(expect_no_warning(run_priority(outdir=outdir)))
  expect_master(resumed$value, fx$output$priority["1"])
  expect_lt(resumed$grids, count_grids(run_priority())$grids)

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  writeLines("not a cache", file.path(cache.dir, "fingerprint.rds"))
  expect_warning(priority <- run_priority(outdir=outdir), "made from other inputs, or is damaged")
  expect_master(priority, fx$output$priority["1"])
  expect_true(file.exists(file.path(old.dir, "terms.rds")))
})

test_that("a half-written file is dropped silently, and metadata of other column types keep the cache", {
  outdir <- withr::local_tempdir()
  metadata <- fx$output$metadata
  metadata$sample <- factor(metadata$sample)
  metadata$group <- factor(metadata$group)

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  writeLines("half written", file.path(outdir, "cache", "terms.rds.tmp"))
  resumed <- count_grids(expect_no_warning(run_priority(metadata=metadata, outdir=outdir)))

  expect_master(resumed$value, fx$output$priority["1"])
  expect_lt(resumed$grids, count_grids(run_priority())$grids)
  expect_false(dir.exists(file.path(outdir, "cache")))
})

test_that("terms written by a run on other inputs into the same cache folder are not used", {
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "cache")

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  late <- readRDS(file.path(cache.dir, "terms.rds"))
  expect_warning(expect_error(stop_after(2, run_priority(metadata=fx$output$metadata.none, outdir=outdir)), "run stopped"),
                 "made from other inputs")
  saveRDS(late, file.path(cache.dir, "terms.rds"))
  resumed <- count_grids(expect_no_warning(run_priority(metadata=fx$output$metadata.none, outdir=outdir)))

  expect_master(resumed$value, fx$output$priority.none)
  expect_equal(resumed$grids, count_grids(run_priority(metadata=fx$output$metadata.none))$grids)
})

test_that("a cache that cannot be written warns and the run goes on; a call stopped before scoring leaves no cache", {
  outdir <- withr::local_tempdir()
  writeLines("a file where the cache folder would go", file.path(outdir, "cache"))

  expect_warning(priority <- run_priority(outdir=outdir), "could not write the cache file .*scoring goes on without the cache")
  expect_master(priority, fx$output$priority["1"])

  lrs <- fx$output$lrs
  lrs$ligand[1] <- "NotAGene"
  outdir <- withr::local_tempdir()
  expect_error(run_priority(lrs=lrs, significance=TRUE, seed=1, outdir=outdir), "not nodes of infMat: NotAGene")
  expect_length(list.files(outdir, recursive=TRUE, include.dirs=TRUE), 0)
})

test_that("cache=FALSE or outdir=NULL keep no cache, and cache must be TRUE or FALSE", {
  outdir <- withr::local_tempdir()

  expect_error(stop_after(2, run_priority(outdir=outdir, cache=FALSE)), "run stopped")
  expect_false(dir.exists(file.path(outdir, "cache")))
  expect_error(stop_after(2, run_priority()), "run stopped")
  expect_length(list.files(outdir, recursive=TRUE, include.dirs=TRUE), 0)

  expect_error(run_priority(cache=NA), "cache must be TRUE or FALSE")
  expect_error(run_priority(cache=c(TRUE, TRUE)), "cache must be TRUE or FALSE")
  expect_error(run_hypercom(file.path(outdir, "run"), cache="yes"), "cache must be TRUE or FALSE")
  expect_false(dir.exists(file.path(outdir, "run")))
})

test_that("HyperCom run again after a stop continues from the cache with the same results", {
  skip_on_cran()
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "6.priority", "cache")

  expect_error(stop_after(3, run_hypercom(outdir, significance=TRUE, seed=1)), "run stopped")
  expect_true(file.exists(file.path(cache.dir, "terms.rds")))
  res <- count_grids(run_hypercom(outdir, significance=TRUE, seed=1, threads=test_threads(2)))
  fresh <- count_grids(run_hypercom(withr::local_tempdir(), significance=TRUE, seed=1))

  expect_identical(res$value$priority, fresh$value$priority)
  expect_lt(res$grids, fresh$grids)
  expect_false(dir.exists(cache.dir))
})
