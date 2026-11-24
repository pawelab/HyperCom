# The cache of prioritize_lr (resume=TRUE): a stopped run continues where it stopped with the results of a run that never stopped, and
# a cache of other inputs or damaged is never used.

real.batch_pairs <- HyperCom:::batch_pairs
real.pool_grids <- HyperCom:::pool_grids

#' Stop a run after some rounds of batches of grids
#'
#' `stop_after` evaluates an expression in which grid scoring stops with an error after a number of rounds
#'
#' This function takes a number of rounds and an expression and lets pool_grids stop with "run stopped" once the given number of
#' rounds (one batch per process each) is done, as a time limit would. Every ligand is its own batch, so the fixture has many rounds.
#' @param rounds An integer number of rounds done before the stop
#' @param expr An expression running prioritize_lr or HyperCom
#' @returns The value of the expression (never reached once the stop comes)
#' @noRd
stop_after <- function(rounds, expr){
  calls <- 0
  stopper <- function(...){
    calls <<- calls + 1
    if(calls > rounds){
      stop("run stopped")
    }
    return(real.pool_grids(...))
  }
  testthat::local_mocked_bindings(batch_pairs=function(ligands, receptors, size) real.batch_pairs(ligands, receptors, 1),
                                  pool_grids=stopper, .package="HyperCom")

  return(expr)
}

#' Count the grids of a run that never stopped
#'
#' `full_grids` returns the number of grids prioritize_lr scores on the fixture without a cache, counted once
#' @noRd
full_grids <- local({
  grids <- NULL
  function(){
    if(is.null(grids)){
      grids <<- count_grids(run_priority())$grids
    }
    return(grids)
  }
})

test_that("a stopped run continues from the cache, scores only the rest and gives handwritten's priorities and p-values", {
  skip_on_cran()
  outdir <- file.path(withr::local_tempdir(), "6.priority")
  cache.dir <- file.path(outdir, "cache")

  expect_error(stop_after(3, run_priority(significance=TRUE, seed=1, outdir=outdir)), "run stopped")
  expect_setequal(list.files(cache.dir), c("fingerprint.rds", "background.rds", "terms.rds"))
  cached <- readRDS(file.path(cache.dir, "terms.rds"))$terms
  expect_gt(nrow(cached), 0)
  expect_lt(nrow(cached), nrow(fx$output$lrs))

  full <- count_grids(run_priority(significance=TRUE, seed=1))
  resumed <- count_grids(collect(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, significance=TRUE,
                                               seed=1, threads=test_threads(2), outdir=outdir)))

  expect_handwritten(resumed$value$value, fx$output$priority.sig)
  expect_lt(resumed$grids, full$grids)
  expect_equal(sum(grepl(paste0("Continuing from the cache in .*: [1-9][0-9]* of ", nrow(fx$output$lrs), " pairs already scored"),
                         resumed$value$messages)), 1)
  expect_false(dir.exists(cache.dir))
  expect_identical(list.files(outdir), "priority.1.csv")
})

test_that("a run stopped after its last batch takes every pair from the cache and scores no grid", {
  skip_on_cran()
  factored <- fx_factored()
  runs <- list(significance=function(outdir) run_priority(significance=TRUE, seed=1, outdir=outdir),
               factored=function(outdir) run_priority(adj=factored$adj, infMat=factored$infMat, outdir=outdir))
  for(run in runs){
    outdir <- withr::local_tempdir()
    first <- testthat::with_mocked_bindings(run(outdir), close_cache=function(dir) invisible(NULL), .package="HyperCom")
    again <- count_grids(run(outdir))

    expect_identical(again$value, first)
    expect_equal(again$grids, 0)
    expect_false(dir.exists(file.path(outdir, "cache")))
  }
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
  factored <- fx_factored()

  expect_error(stop_after(2, run_priority(threads=1, outdir=outdir)), "run stopped")
  expect_handwritten(run_priority(weights=c(1, 0.5, 0), threads=test_threads(2), outdir=outdir), fx$output$priority)

  expect_error(stop_after(2, run_priority(adj=factored$adj, infMat=factored$infMat, threads=test_threads(2), outdir=outdir)), "run stopped")
  expect_identical(run_priority(adj=factored$adj, infMat=factored$infMat, weights=c(1, 0.5), threads=1, outdir=outdir),
                   run_priority(adj=factored$adj, infMat=factored$infMat, weights=c(1, 0.5)))
})

test_that("a cache made from other inputs, or damaged, is moved to cache.old with a warning and can be moved back", {
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "cache")
  old.dir <- file.path(outdir, "cache.old")

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  expect_warning(priority <- run_priority(metadata=fx$output$metadata.none, outdir=outdir),
                 "made from other inputs, or is damaged, and is not used; it is moved to .*cache.old")
  expect_handwritten(priority, fx$output$priority.none)
  expect_false(dir.exists(cache.dir))
  expect_true(file.rename(old.dir, cache.dir))
  resumed <- count_grids(expect_no_warning(run_priority(outdir=outdir)))
  expect_handwritten(resumed$value, fx$output$priority["1"])
  expect_lt(resumed$grids, full_grids())

  expect_error(stop_after(2, run_priority(outdir=outdir)), "run stopped")
  writeLines("not a cache", file.path(cache.dir, "fingerprint.rds"))
  expect_warning(priority <- run_priority(outdir=outdir), "made from other inputs, or is damaged")
  expect_handwritten(priority, fx$output$priority["1"])
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

  expect_handwritten(resumed$value, fx$output$priority["1"])
  expect_lt(resumed$grids, full_grids())
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

  expect_handwritten(resumed$value, fx$output$priority.none)
  expect_equal(resumed$grids, full_grids())
})

test_that("a cache that cannot be written warns and the run goes on; a call stopped before scoring leaves no cache", {
  outdir <- withr::local_tempdir()
  writeLines("a file where the cache folder would go", file.path(outdir, "cache"))

  expect_warning(priority <- run_priority(outdir=outdir), "could not write the cache file .*scoring goes on without the cache")
  expect_handwritten(priority, fx$output$priority["1"])

  lrs <- fx$output$lrs
  lrs$ligand <- "NotAGene"
  outdir <- withr::local_tempdir()
  expect_error(run_priority(lrs=lrs, significance=TRUE, seed=1, outdir=outdir), "no interaction of lrs has all its genes in infMat")
  expect_length(list.files(outdir, recursive=TRUE, include.dirs=TRUE), 0)
})

test_that("resume=FALSE or outdir=NULL keep no cache", {
  outdir <- withr::local_tempdir()

  expect_error(stop_after(2, run_priority(outdir=outdir, resume=FALSE)), "run stopped")
  expect_false(dir.exists(file.path(outdir, "cache")))
  expect_error(stop_after(2, run_priority()), "run stopped")
  expect_length(list.files(outdir, recursive=TRUE, include.dirs=TRUE), 0)
})

test_that("HyperCom run again after a stop continues from the cache with the same results", {
  skip_on_cran()
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "6.priority", "cache")

  expect_error(stop_after(3, run_hypercom(outdir, significance=TRUE, seed=1)), "run stopped")
  expect_true(file.exists(file.path(cache.dir, "terms.rds")))
  res <- count_grids(run_hypercom(outdir, significance=TRUE, seed=1, threads=test_threads(2)))
  fresh <- count_grids(run_hypercom(NULL, significance=TRUE, seed=1))

  expect_identical(res$value$priority, fresh$value$priority)
  expect_lt(res$grids, fresh$grids)
  expect_false(dir.exists(cache.dir))
})
