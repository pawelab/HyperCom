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

test_that("prioritize_lr in parallel gives handwritten's priorities", {
  priority <- run_priority(weights=c(1, 0.5), parallel=TRUE)

  expect_handwritten(priority, fx$output$priority[c("1", "0.5")])
})

test_that("background pairs without a priority are left out of the p-values, with a warning", {
  expect_identical(HyperCom:::permutation_test(0.5, c(0.2, NaN, 0.7, NA, 0.5)), 2 / 3)
  expect_identical(HyperCom:::permutation_test(NaN, c(0.2, 0.7)), NA_real_)
  expect_identical(HyperCom:::permutation_test(0.5, c(NaN, NA)), NA_real_)

  keys <- function(lrs) HyperCom:::pair_keys(as.character(lrs$ligand), as.matrix(lrs[, -(1:2), drop=FALSE]))
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

test_that("seed does not change the random numbers of the session, and outdir must be one string or NULL", {
  withr::local_seed(42)
  before <- get(".Random.seed", envir=globalenv())
  priority <- run_priority(significance=TRUE, seed=1)

  expect_identical(get(".Random.seed", envir=globalenv()), before)
  expect_identical(run_priority(significance=TRUE, seed=1), priority)
  expect_error(run_priority(outdir=4), "outdir must be one string")
})

test_that("parallel uses one thread per physical core, and one thread without parallel", {
  physical <- min(parallelly::availableCores(), HyperCom:::physical_cores(), na.rm=TRUE)

  expect_equal(HyperCom:::get_threads(parallel=TRUE), physical)
  expect_equal(HyperCom:::get_threads(), 1)
  expect_lte(HyperCom:::physical_cores(), parallel::detectCores())
  if(Sys.info()[["sysname"]] != "Linux"){
    expect_equal(HyperCom:::physical_cores(), as.integer(parallelly::availableCores(logical=FALSE)))
  }
})

test_that("prioritize_lr reports progress as messages, not printed lines, and verbose=FALSE is silent", {
  dir <- withr::local_tempdir()
  run <- function(verbose){
    return(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5), threads=test_threads(2),
                         outdir=dir, verbose=verbose))
  }

  msgs <- character(0)
  printed <- utils::capture.output(x <- withCallingHandlers(run(TRUE), message=function(m){
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }))

  expect_length(printed, 0)
  expect_equal(sum(grepl("pairs scored", msgs)), 1)
  expect_equal(sum(grepl("Weight 0.5: top pair", msgs)), 1)
  expect_silent(run(FALSE))
})

test_that("format_duration gives seconds, minutes or hours", {
  expect_identical(HyperCom:::format_duration(30), "30 s")
  expect_identical(HyperCom:::format_duration(300), "5 min")
  expect_identical(HyperCom:::format_duration(9000), "2.5 h")
})

test_that("grids bigger than the kernel's cache budget are scored in about 20 batches or more, smaller ones in big batches", {
  batch_size <- HyperCom:::batch_size

  expect_equal(batch_size(16, 1000, c(500, 500), 4090), 512)
  expect_equal(batch_size(1, 1000, c(500, 500), 4090), 256)
  expect_equal(batch_size(16, 20000, c(10000, 6667, 3333), 400), 20)
  expect_equal(batch_size(16, 20000, c(10000, 6667, 3333), 40), 16)
  expect_equal(batch_size(1, 20000, c(10000, 6667, 3333), 4090), 205)
})

test_that("group purity is at least 0, so perfectly balanced groups give priority 0, not NaN", {
  balanced <- factor(paste0("t", rep(1:5, 4)))

  expect_lt(1 - 1 / 5 - HyperCom:::gini_impurity(balanced), 0)
  expect_identical(HyperCom:::group_purity(balanced, 1 - 1 / 5), 0)
  expect_identical(HyperCom:::group_purity(balanced[1:4], 0), 1)
  expect_identical(HyperCom:::group_purity(factor(c("a", "a", "a", "b")), 0.5), 0.5 - HyperCom:::gini_impurity(factor(c("a", "a", "a", "b"))))
})
