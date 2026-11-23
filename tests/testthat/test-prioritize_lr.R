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
  bound <- HyperCom:::bound_terms
  testthat::local_mocked_bindings(calculate_priority=function(lrs, ...){
    res <- calculate(lrs, ...)
    res[keys(lrs) %in% broken, "max"] <- NaN
    return(res)
  }, bound_terms=function(lrs, ...){
    res <- bound(lrs, ...)
    res$lo[keys(lrs) %in% broken, "max"] <- res$hi[keys(lrs) %in% broken, "max"] <- NA
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

test_that("threads=NULL uses one thread per physical core, never more than the available cores", {
  threads <- HyperCom:::get_threads()

  expect_gte(threads, 1)
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

test_that("batches hold a grid per thread at least and 64 MB of scores at most, and about 20 or more when grids are big", {
  batch_size <- HyperCom:::batch_size
  big <- c(10000, 6667, 3333)

  for(threads in c(1, 16)){
    for(bytes in c(16, 32)){
      small <- batch_size(threads, 1000, c(500, 500), 4090, bytes)
      expect_gte(small, threads)
      expect_lte(small, max(threads, 2^26 / (bytes * 1000)))
      expect_gte(small, min(max(256, 32 * threads), floor(2^26 / (bytes * 1000))))
      expect_lte(batch_size(threads, 20000, big, 4090, bytes), max(threads, ceiling(4090 / 20)))
      expect_equal(batch_size(threads, 20000, big, 40, bytes), max(threads, ceiling(40 / 20)))
    }
  }
})

test_that("group purity is at least 0, so perfectly balanced groups give priority 0, not NaN", {
  balanced <- factor(paste0("t", rep(1:5, 4)))

  expect_lt(1 - 1 / 5 - HyperCom:::gini_impurity(balanced), 0)
  expect_identical(HyperCom:::group_purity(balanced, 1 - 1 / 5), 0)
  expect_identical(HyperCom:::group_purity(balanced[1:4], 0), 1)
  expect_identical(HyperCom:::group_purity(factor(c("a", "a", "a", "b")), 0.5), 0.5 - HyperCom:::gini_impurity(factor(c("a", "a", "a", "b"))))
})
