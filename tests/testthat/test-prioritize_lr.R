run_priority <- function(metadata=fx$output$metadata, ...){
  return(quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, metadata, ...)))
}

test_that("prioritize_lr matches master for several weights", {
  priority <- run_priority(weights=c(1, 0.5, 0), outdir=NULL)

  expect_identical(names(priority), names(fx$output$priority))
  for(w in names(priority)){
    expect_same_priority(priority[[w]], fx$output$priority[[w]])
  }
})

test_that("prioritize_lr matches master with a single group", {
  priority <- run_priority(metadata=fx$output$metadata.none, outdir=NULL)

  expect_same_priority(priority[["1"]], fx$output$priority.none[["1"]])
})

test_that("prioritize_lr p-values match master for the same seed and are BH adjusted", {
  skip_on_cran()
  set.seed(1)
  priority <- run_priority(significance=TRUE, outdir=NULL)
  df <- priority[["1"]]

  expect_same_priority(df, fx$output$priority.sig[["1"]])
  expect_true(all(df$p >= 0 & df$p <= 1, na.rm=TRUE))
  expect_equal(df$p.adj, stats::p.adjust(df$p, method="BH"))
})

test_that("priority tables are sorted, bounded and keep the LR columns", {
  priority <- run_priority(weights=c(1, 0), outdir=NULL)

  for(df in priority){
    expect_identical(colnames(df), c(colnames(fx$output$lrs), "priority"))
    expect_setequal(df$interaction, fx$output$lrs$interaction)
    expect_false(is.unsorted(rev(df$priority), na.rm=TRUE))
    expect_true(all(df$priority >= 0 & df$priority <= 1, na.rm=TRUE))
  }
})

test_that("prioritize_lr writes one priority file per weight", {
  outdir <- file.path(withr::local_tempdir(), "6.priority")
  priority <- run_priority(weights=c(1, 0.5), outdir=outdir)

  expect_setequal(list.files(outdir), c("priority.1.csv", "priority.0.5.csv"))
  for(w in names(priority)){
    written <- utils::read.csv(file.path(outdir, paste0("priority.", w, ".csv")), check.names=FALSE)
    expect_equal(plain_df(written), plain_df(priority[[w]]), tolerance=tol)
  }
})

test_that("prioritize_lr gives the same priorities in parallel", {
  skip_on_cran()
  cl <- parallel::makeCluster(2)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
  priority <- run_priority(weights=c(1, 0.5), parallel=TRUE, cl=cl, outdir=NULL)

  for(w in names(priority)){
    expect_same_priority(priority[[w]], fx$output$priority[[w]])
  }
})
