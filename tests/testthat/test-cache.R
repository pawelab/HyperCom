# The cache of prioritize_lr (resume=TRUE): a stopped run continues where it stopped with the results of a run that never stopped, and
# a cache of other inputs, of the other kind of scoring or damaged is never used.

real.batch_pairs <- HyperCom:::batch_pairs
real.score_grids <- HyperCom:::score_grids
real.screen_grids <- HyperCom:::screen_grids

#' Stop a run after some batches of grids
#'
#' `stop_after` evaluates an expression in which grid scoring (or screening) stops with an error after a number of batches
#'
#' This function takes a number of batches and an expression and lets score_grids (or screen_grids) stop with "run stopped" once the
#' given number of batches is done, as a time limit would. Every ligand is its own batch, and screening batches hold 50 grids, so the
#' fixture has many batches.
#' @param batches An integer number of batches done before the stop
#' @param expr An expression running prioritize_lr or HyperCom
#' @param fun A string, "score_grids" or "screen_grids"
#' @returns The value of the expression (never reached once the stop comes)
#' @noRd
stop_after <- function(batches, expr, fun="score_grids"){
  calls <- 0
  real <- if(fun == "score_grids") real.score_grids else real.screen_grids
  stopper <- function(...){
    calls <<- calls + 1
    if(calls > batches){
      stop("run stopped")
    }
    return(real(...))
  }
  testthat::local_mocked_bindings(batch_pairs=function(ligands, receptors, size) real.batch_pairs(ligands, receptors, 1),
                                  .package="HyperCom")
  if(fun == "score_grids"){
    testthat::local_mocked_bindings(score_grids=stopper, .package="HyperCom")
  } else {
    testthat::local_mocked_bindings(batch_size=function(...) 50L, screen_grids=stopper, .package="HyperCom")
  }

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

  expect_error(stop_after(3, run_priority(significance=TRUE, seed=1, exact.priority=TRUE, outdir=outdir)), "run stopped")
  expect_setequal(list.files(cache.dir), c("fingerprint.rds", "background.rds", "terms.rds"))
  cached <- readRDS(file.path(cache.dir, "terms.rds"))$terms
  expect_gt(nrow(cached), 0)
  expect_lt(nrow(cached), nrow(fx$output$lrs))

  full <- count_grids(run_priority(significance=TRUE, seed=1, exact.priority=TRUE))
  resumed <- count_grids(collect(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, significance=TRUE,
                                               seed=1, exact.priority=TRUE, threads=test_threads(2), outdir=outdir)))

  expect_handwritten(resumed$value$value, fx$output$priority.sig)
  expect_lt(resumed$grids, full$grids)
  expect_equal(sum(grepl(paste0("Continuing from the cache in .*: [1-9][0-9]* of ", nrow(fx$output$lrs), " pairs already scored"),
                         resumed$value$messages)), 1)
  expect_false(dir.exists(cache.dir))
  expect_identical(list.files(outdir), "priority.1.csv")
})

test_that("a run stopped while screening continues from the cached bounds with the same results", {
  skip_on_cran()
  outdir <- withr::local_tempdir()
  cache.dir <- file.path(outdir, "cache")

  expect_error(stop_after(5, run_priority(weights=c(1, 0.5), significance=TRUE, seed=1, outdir=outdir), "screen_grids"), "run stopped")
  expect_setequal(list.files(cache.dir), c("fingerprint.rds", "background.rds", "bounds.rds"))
  expect_gt(nrow(readRDS(file.path(cache.dir, "bounds.rds"))$lo), 0)

  full <- count_grids(run_priority(weights=c(1, 0.5), significance=TRUE, seed=1), "screen_grids")
  resumed <- count_grids(collect(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5),
                                               significance=TRUE, seed=1, threads=test_threads(2), outdir=outdir)), "screen_grids")

  expect_identical(resumed$value$value, full$value)
  expect_lt(resumed$grids, full$grids)
  expect_equal(sum(grepl("Continuing from the cache in .*: [1-9][0-9]* of [0-9]+ (real|background) pairs already screened",
                         resumed$value$messages)), 1)
  expect_false(dir.exists(cache.dir))
})

test_that("a run stopped after its last batch takes every pair from the cache and neither scores nor screens a grid", {
  skip_on_cran()
  factored <- fx_factored()
  runs <- list(significance=function(outdir) run_priority(significance=TRUE, seed=1, outdir=outdir),
               factored=function(outdir) run_priority(adj=factored$adj, infMat=factored$infMat, outdir=outdir))
  for(run in runs){
    outdir <- withr::local_tempdir()
    first <- testthat::with_mocked_bindings(run(outdir), close_cache=function(dir) invisible(NULL), .package="HyperCom")
    again <- count_grids(count_grids(run(outdir)), "screen_grids")

    expect_identical(again$value$value, first)
    expect_equal(again$value$grids + again$grids, 0)
    expect_false(dir.exists(file.path(outdir, "cache")))
  }
})

test_that("without seed, a stopped run continues with the background pairs it drew", {
  skip_on_cran()
  outdir <- withr::local_tempdir()

  expect_error(withr::with_seed(7, stop_after(3, run_priority(significance=TRUE, exact.priority=TRUE, outdir=outdir))), "run stopped")
  resumed <- withr::with_seed(99, run_priority(significance=TRUE, exact.priority=TRUE, outdir=outdir))

  expect_identical(resumed, withr::with_seed(7, run_priority(significance=TRUE, exact.priority=TRUE)))
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

test_that("a memory mode run continues from its cache, and a cache of the other kind of scoring is moved aside, saying so", {
  outdir <- withr::local_tempdir()
  factored <- fx_factored()
  memory.gb <- memory_for(fx$output$metadata)
  run <- function(...) run_priority(adj=factored$adj, infMat=factored$infMat, outdir=outdir, ...)

  expect_error(stop_after(2, run(memory.gb=memory.gb)), "run stopped")
  resumed <- count_grids(run(memory.gb=memory.gb))
  expect_identical(resumed$value, run_priority(adj=factored$adj, infMat=factored$infMat, memory.gb=memory.gb))
  expect_lt(resumed$grids, full_grids())

  expect_error(stop_after(2, run(memory.gb=memory.gb)), "run stopped")
  expect_warning(run(memory.gb=Inf), "made with memory mode but this run uses the exact kernels")
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

  expect_error(stop_after(3, run_hypercom(outdir, significance=TRUE, seed=1, exact.priority=TRUE)), "run stopped")
  expect_true(file.exists(file.path(cache.dir, "terms.rds")))
  res <- count_grids(run_hypercom(outdir, significance=TRUE, seed=1, exact.priority=TRUE, threads=test_threads(2)))
  fresh <- count_grids(run_hypercom(NULL, significance=TRUE, seed=1, exact.priority=TRUE))

  expect_identical(res$value$priority, fresh$value$priority)
  expect_lt(res$grids, fresh$grids)
  expect_false(dir.exists(cache.dir))
})
