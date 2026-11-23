# Screening (R/screen.R, src/screen_grids.cpp): bounds that hold the exact cell scores and priority terms, the same at any number of
# threads and in both kernel modes, and p-values identical to scoring every background pair

#' Expect bounds to hold values
#'
#' `expect_bounds` checks that lower and upper bounds hold values, and are missing exactly where the values are
#' @param lo A numeric vector or matrix of lower bounds
#' @param hi A numeric vector or matrix of upper bounds
#' @param x A numeric vector or matrix of values
#' @param tolerance A double, the rounding allowed (bounds and values are divided differently)
#' @returns x, invisibly
#' @noRd
expect_bounds <- function(lo, hi, x, tolerance=1e-14){
  expect_identical(is.na(lo), is.na(x))
  expect_identical(is.na(hi), is.na(x))
  ok <- !is.na(x)
  expect_true(all(lo[ok] <= x[ok] + tolerance))
  expect_true(all(hi[ok] >= x[ok] - tolerance))

  return(invisible(x))
}

#' Draw the fixture's background pairs
#'
#' `fx_background` returns the background pairs prioritize_lr draws with seed 1 on the fixture
#' @returns A data frame of background pairs
#' @noRd
fx_background <- function(){
  return(withr::with_seed(1, HyperCom:::generate_background_lrs(fx$output$lrs, fx$output$metadata)))
}

test_that("screened bounds hold every cell's scores, the same at any number of threads and with ids kept between calls", {
  skip_on_cran()
  bg <- fx_background()
  grids <- HyperCom:::pair_grids(as.character(bg$ligand), as.matrix(bg[, -(1:2)]))
  md <- fx$output$metadata

  for(infMat in list(fx$output$infMat, fx_factored()$infMat)){
    blocks <- HyperCom:::grid_blocks(infMat, HyperCom:::grid_cells(rownames(infMat), md), 1)
    screen <- function(threads, keep.ids=FALSE) HyperCom:::screen_grids(infMat, md, grids$ligand, grids$receptor, threads, blocks, keep.ids)
    exact <- HyperCom:::score_grids(infMat, md, grids$ligand, grids$receptor, test_threads(2))
    b <- screen(1)
    expect_identical(screen(test_threads(3)), b)
    for(threads in c(test_threads(3), 1)){
      expect_identical(screen(threads, keep.ids=TRUE), b)
    }
    HyperCom:::screen_release_cpp()
    expect_identical(screen(test_threads(3)), b)
    expect_bounds(b$sender.lo, b$sender.hi, exact$sender)
    expect_bounds(b$receiver.lo, b$receiver.hi, exact$receiver)
    expect_lt(stats::median(b$sender.hi - b$sender.lo), 1e-3)
  }
})

test_that("screened bounds hold the scores of grids with ties, zeros, negative, infinite and missing values in interleaved samples", {
  set.seed(5)
  n <- 61
  sample <- c(rep(1:3, length.out=n - 1), 4L)
  sizes <- tabulate(sample)
  blocks <- signif(stats::runif(sum(sizes^2), 0.1, 10), 1)
  blocks[sample.int(length(blocks), 200)] <- 0
  blocks[sample.int(length(blocks), 30)] <- NaN
  blocks[sample.int(length(blocks), 30)] <- Inf
  blocks[sample.int(length(blocks), 100)] <- -blocks[sample.int(length(blocks), 100)]
  A <- matrix(signif(10^stats::runif(n * 4, -300, 300), 2), n, 4)
  A[1:3, 1] <- 0
  A[4, 2] <- NaN
  B <- matrix(signif(10^stats::runif(4 * n, -300, 300), 2), 4, n)
  B[1, 5] <- -1
  B[2, ] <- 1
  lig <- c(0L, 1L, 2L, 3L, 0L)
  rec <- c(0L, 1L, 2L, 3L, 1L)

  exact <- HyperCom:::score_blocks_cpp(blocks, A, B, sample, lig, rec, 1)
  for(threads in c(1, test_threads(3))){
    for(budget in c(2^28, 0)){
      b <- HyperCom:::screen_blocks_cpp(blocks, A, B, sample, lig, rec, threads, budget)
      expect_bounds(b$sender_lo, b$sender_hi, exact$sender)
      expect_bounds(b$receiver_lo, b$receiver_hi, exact$receiver)
    }
  }
  expect_true(all(is.nan(exact$sender[n, ])))
  equal <- HyperCom:::screen_blocks_cpp(rep(1, sum(sizes^2)), matrix(1, n, 1), matrix(1, 1, n), sample, 0L, 0L, 1, 2^28)
  expect_equal(equal$sender_lo[1:3, 1], rep(0, 3))
  expect_error(HyperCom:::screen_blocks_cpp(blocks[-1], A, B, sample, lig, rec, 1, 2^28), "blocks do not match the samples")
  expect_error(HyperCom:::screen_grids(fx$output$infMat, fx$output$metadata, NA, NA), "network score")
})

test_that("balance bounds are the smallest and largest balance over the possible numbers of receivers", {
  balance <- function(n, t) HyperCom:::gini_impurity(c(rep(FALSE, n - t), rep(TRUE, t))) / 0.5
  for(n in 0:9){
    for(t0 in 0:n){
      for(t1 in t0:n){
        all <- vapply(t0:t1, function(t) balance(n, t), numeric(1))
        expect_identical(HyperCom:::bound_balance(n, t0, t1), c(min(all), max(all)))
      }
    }
  }
})

test_that("a background pair needs its exact priority when its interval (L, U] holds a real priority or has no bound", {
  x <- c(0.2, 0.5, NA)
  L <- c(0.1, 0.1, 0.2, 0.3, NA, 0.6, 0.6)
  U <- c(0.15, 0.2, 0.3, 0.6, 0.3, Inf, 0.7)

  expect_identical(HyperCom:::ambiguous(L, U, x), c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that("term bounds hold the exact priority terms of background pairs", {
  skip_on_cran()
  bg <- fx_background()
  for(metadata in list(fx$output$metadata, fx$output$metadata.none)){
    exact <- HyperCom:::calculate_priority(bg, fx$output$adj, fx$output$infMat, metadata, test_threads(2))
    b <- HyperCom:::bound_terms(bg, fx$output$adj, fx$output$infMat, metadata, test_threads(2))
    ok <- !is.na(exact)
    expect_true(all(b$lo[ok] <= exact[ok] & exact[ok] <= b$hi[ok]))
    expect_true(all(is.na(b$lo[!ok]) & is.na(b$hi[!ok])))
    expect_identical(b$lo[, c("frac", "pct.min", "purity")], exact[, c("frac", "pct.min", "purity")])
  }
})

test_that("screened p-values are those of scoring every background pair, which few are, with screened or exact real pairs", {
  skip_on_cran()
  factored <- fx_factored()$infMat

  for(exact.priority in c(FALSE, TRUE)){
    screened <- count_grids(run_priority(weights=c(1, 0.5, 0), significance=TRUE, seed=1, exact.priority=exact.priority,
                                         threads=test_threads(2)))
    all <- count_grids(every_background(run_priority(weights=c(1, 0.5, 0), significance=TRUE, seed=1, exact.priority=exact.priority,
                                                     threads=test_threads(2))))
    expect_identical(screened$value, all$value)
    expect_lt(screened$grids, all$grids / 2)
  }
  expect_identical(run_priority(infMat=factored, metadata=fx$output$metadata.none, significance=TRUE, seed=3),
                   every_background(run_priority(infMat=factored, metadata=fx$output$metadata.none, significance=TRUE, seed=3)))
})

test_that("background pairs without a priority give the p-values and warning of scoring every background pair", {
  skip_on_cran()
  metadata <- fx$output$metadata
  metadata$sample <- as.character(metadata$sample)
  metadata$sample[match(metadata$cell[metadata$sample != "gene"][1], metadata$cell)] <- "alone"

  screened <- collect(run_priority(metadata=metadata, weights=c(1, 0.5), significance=TRUE, seed=1))
  all <- collect(every_background(run_priority(metadata=metadata, weights=c(1, 0.5), significance=TRUE, seed=1)))

  expect_identical(screened, all)
  expect_match(screened$warnings, "^[1-9][0-9]* of 1000 background pairs have no priority \\(NaN\\)", all=FALSE)
})
