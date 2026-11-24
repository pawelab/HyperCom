test_that("prioritize_lr matches handwritten with a single group", {
  expect_handwritten(run_priority(metadata=fx$output$metadata.none), fx$output$priority.none)
})

test_that("priority tables are sorted, bounded and keep the LR columns", {
  priority <- run_priority(weights=c(1, 0))

  for(df in priority){
    expect_identical(colnames(df), c(colnames(fx$output$lrs), "priority"))
    expect_setequal(df$interaction, fx$output$lrs$interaction)
    expect_false(is.unsorted(rev(df$priority), na.rm=TRUE))
    expect_true(all(df$priority >= 0 & df$priority <= 1, na.rm=TRUE))
  }
})

test_that("prioritize_lr writes one priority file per weight, into a nested directory it creates (#9)", {
  outdir <- file.path(withr::local_tempdir(), "run", "6.priority")
  priority <- run_priority(weights=c(1, 0.5), outdir=outdir)

  expect_setequal(list.files(outdir), c("priority.1.csv", "priority.0.5.csv"))
  for(w in names(priority)){
    written <- utils::read.csv(file.path(outdir, paste0("priority.", w, ".csv")), check.names=FALSE)
    expect_equal(plain_df(written), plain_df(priority[[w]]), tolerance=tol)
  }
})

test_that("background pairs without a priority are left out of the p-values, with a warning", {
  expect_identical(HyperCom:::permutation_test(0.5, c(0.2, NaN, 0.7, NA, 0.5)), 2 / 3)
  expect_identical(HyperCom:::permutation_test(NaN, c(0.2, 0.7)), NA_real_)
  expect_identical(HyperCom:::permutation_test(0.5, c(NaN, NA)), NA_real_)

  keys <- HyperCom:::lrs_keys
  background <- withr::with_seed(1, HyperCom:::generate_background_lrs(fx$output$lrs, fx$output$metadata))
  broken <- keys(background)[999:1000]
  calculate <- HyperCom:::calculate_priority
  testthat::local_mocked_bindings(calculate_priority=function(lrs, ...){
    res <- calculate(lrs, ...)
    res[keys(lrs) %in% broken, "max"] <- NaN
    return(res)
  }, .package="HyperCom")
  expect_warning(priority <- run_priority(weights=c(1, 0.5), significance=TRUE, seed=1),
                 "2 of 1000 background pairs have no priority .* use the other 998")
  for(w in names(priority)){
    p <- priority[[w]]$p
    expect_false(anyNA(p))
    expect_equal(p * 998, round(p * 998), tolerance=1e-12)
  }
})

test_that("threads=NULL uses one process per physical core up to 4, never more than the available cores", {
  threads <- HyperCom:::get_threads()

  expect_gte(threads, 1)
  expect_lte(threads, 4)
  expect_equal(threads, min(4, parallelly::availableCores(), HyperCom:::physical_cores(), na.rm=TRUE))
  expect_lte(threads, parallelly::availableCores())
  expect_equal(HyperCom:::get_threads(1), 1)
  expect_lte(HyperCom:::physical_cores(), parallel::detectCores())
})

test_that("prioritize_lr reports progress as messages, not printed lines, and verbose=FALSE is silent", {
  dir <- withr::local_tempdir()
  run <- function(verbose){
    return(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5), threads=test_threads(2),
                         verbose=verbose, outdir=dir))
  }

  printed <- utils::capture.output(res <- collect(run(TRUE)))

  expect_length(printed, 0)
  expect_gte(sum(grepl("pairs scored", res$messages)), 1)
  expect_equal(sum(grepl("Weight 0.5: top pair", res$messages)), 1)
  expect_silent(run(FALSE))
})

test_that("format_duration gives seconds, minutes or hours", {
  expect_identical(HyperCom:::format_duration(30), "30 s")
  expect_identical(HyperCom:::format_duration(300), "5 min")
  expect_identical(HyperCom:::format_duration(9000), "2.5 h")
})

test_that("batches hold 256 grids and 64 MB of scores at most, a batch per process, and about 20 rounds or more when grids are big", {
  batch_size <- HyperCom:::batch_size
  big <- c(10000, 6667, 3333)

  expect_equal(batch_size(1, 1000, c(500, 500), 4090), 256)
  expect_equal(batch_size(16, 1000, c(500, 500), 4090), 256)
  expect_equal(batch_size(16, 1000, c(500, 500), 3), 1)
  expect_equal(batch_size(2, 1000, c(500, 500), 3), 2)
  expect_equal(batch_size(1, 20000, rep(10, 2000), 4090), floor(2^26 / (16 * 20000)))
  expect_equal(batch_size(1, 20000, big, 4090), ceiling(4090 / 20))
  expect_equal(batch_size(16, 20000, big, 4090), ceiling(4090 / (20 * 16)))
})

test_that("the processes that fit in memory are counted from the estimate, and one that does not fit stops", {
  skip_on_cran()
  need <- HyperCom:::run_memory(1e8)

  expect_equal(HyperCom:::fit_workers(1e8, 16, Inf), 16)
  expect_equal(HyperCom:::fit_workers(1e8, 16, (need[["handout"]] + 3.5 * need[["worker"]]) / 1e9), 3)
  expect_equal(HyperCom:::fit_workers(1e8, 16, (need[["handout"]] + 1.5 * need[["worker"]]) / 1e9), 1)
  expect_equal(HyperCom:::fit_workers(1e8, 1, Inf), 1)
  expect_error(HyperCom:::fit_workers(1e8, 16, need[["main"]] / 1e9 * 0.99), "scoring needs about 8.99 GB of memory but")
  res <- collect(prioritize_lr(fx$output$lrs[1:3, ], fx$output$adj, fx$output$infMat, fx$output$metadata, threads=test_threads(2),
                               memory.gb=1.1, outdir=NULL))
  expect_true(any(grepl("1 of the 2 processes asked for fit in 1.1 GB", res$messages)) || test_threads(2) == 1)
  expect_true(any(grepl("pairs in 1 R process", res$messages, fixed=TRUE)))
})

test_that("group purity is at least 0, so perfectly balanced groups give priority 0, not NaN", {
  balanced <- factor(paste0("t", rep(1:5, 4)))

  expect_lt(1 - 1 / 5 - HyperCom:::gini_impurity(balanced), 0)
  expect_identical(HyperCom:::group_purity(balanced, 1 - 1 / 5), 0)
  expect_identical(HyperCom:::group_purity(balanced[1:4], 0), 1)
  expect_identical(HyperCom:::group_purity(factor(c("a", "a", "a", "b")), 0.5), 0.5 - HyperCom:::gini_impurity(factor(c("a", "a", "a", "b"))))
})
