# Background screening (R/screen.R; screen_grids and bound_ranks in R/grids.R): bounds that hold the exact ranks, cell scores and
# priority terms, the same from blocks made once and per batch, and p-values identical to scoring every background pair, on one process
# and on a cluster

real.score_grids <- HyperCom:::score_grids

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
#' @returns A dataframe of background pairs
#' @noRd
fx_background <- function(){
  return(withr::with_seed(1, HyperCom:::generate_background_lrs(fx$output$lrs, fx$output$metadata)))
}

#' Count the grids a run scores exactly
#'
#' `scored_grids` returns the value of an expression and the number of grids passed to score_grids while it ran
#' @param expr An expression running prioritize_lr
#' @returns A list with value and grids
#' @noRd
scored_grids <- function(expr){
  grids <- 0
  value <- testthat::with_mocked_bindings(expr, score_grids=function(infMat, metadata, ligands, receptors, ...){
    grids <<- grids + length(ligands)
    return(real.score_grids(infMat, metadata, ligands, receptors, ...))
  }, .package="HyperCom")

  return(list(value=value, grids=grids))
}

test_that("screened bounds hold every cell's scores and are the same from blocks made once, for a dense and a factored influence matrix", {
  skip_on_cran()
  bg <- fx_background()
  grids <- HyperCom:::pair_grids(as.character(bg$ligand), as.matrix(bg[, -(1:2)]))
  factored <- generate_infMat(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL), outfile=NULL)
  md <- fx$output$metadata

  for(infMat in list(fx$output$infMat, factored)){
    exact <- HyperCom:::score_grids(infMat, md, grids$ligand, grids$receptor)
    b <- HyperCom:::screen_grids(infMat, md, grids$ligand, grids$receptor)
    blocks <- HyperCom:::grid_blocks(infMat, HyperCom:::grid_cells(rownames(infMat), md))
    expect_identical(HyperCom:::screen_grids(infMat, md, grids$ligand, grids$receptor, blocks), b)
    expect_bounds(b$sender.lo, b$sender.hi, exact$sender)
    expect_bounds(b$receiver.lo, b$receiver.hi, exact$receiver)
    expect_lt(stats::median(b$sender.hi - b$sender.lo), 1e-3)
  }
})

test_that("screened bounds hold the scores of grids with ties, zeros, infinite and missing values in interleaved samples", {
  withr::local_seed(5)
  n <- 61
  cells <- paste0("c", seq_len(n))
  genes <- paste0("g", 1:4)
  nodes <- sample(c(cells, genes))
  F <- matrix(signif(stats::runif(length(nodes)^2, 0.1, 10), 1), length(nodes), length(nodes), dimnames=list(nodes, nodes))
  F[sample.int(length(F), 200)] <- 0
  F[sample.int(length(F), 30)] <- NaN
  F[sample.int(length(F), 30)] <- Inf
  F[cells, genes] <- signif(10^stats::runif(n * 4, -300, 300), 2)
  F[genes, cells] <- signif(10^stats::runif(4 * n, -300, 300), 2)
  F[cells[1:3], "g1"] <- 0
  F[cells[4], "g2"] <- NaN
  F["g2", cells] <- 1
  samples <- c(rep(1:3, length.out=n - 1), 4L)
  md <- data.frame(cell=c(cells, genes), sample=c(paste0("s", samples), rep("gene", 4)), group=c(rep(c("a", "b"), length.out=n), rep("gene", 4)))
  md <- md[rev(seq_len(nrow(md))), ]
  ligands <- c("g1", "g2", "g3", "g4", "g1")
  receptors <- c("g1", "g2", "g3", "g4", "g2")

  exact <- HyperCom:::score_grids(F, md, ligands, receptors)
  b <- HyperCom:::screen_grids(F, md, ligands, receptors)
  expect_bounds(b$sender.lo, b$sender.hi, exact$sender)
  expect_bounds(b$receiver.lo, b$receiver.hi, exact$receiver)
  expect_true(all(is.nan(exact$sender[md$cell[md$sample != "gene"] == cells[n], ])))
  ones <- F
  ones[] <- 1
  equal <- HyperCom:::screen_grids(ones, md, "g1", "g1")
  expect_true(all(equal$sender.lo == 0 & equal$sender.hi == 1, na.rm=TRUE))
  expect_equal(sum(!is.na(equal$sender.lo)), n - 1)
})

test_that("rank bounds hold the number of smaller entries of every entry, with ties, infinite and missing values", {
  withr::local_seed(7)
  vs <- list(signif(stats::runif(500, 0, 10), 1), c(10^stats::runif(300, -300, 300), Inf, Inf, NA, NaN, 0), rep(2.5, 50),
             c(1, 1 + 2^-52, 1 + 2^-51, NA), c(Inf, NA, Inf), c(NA_real_, NaN))

  for(v in vs){
    b <- HyperCom:::bound_ranks(v)
    ok <- !is.na(v)
    exact <- rank(v[ok], ties.method="min") - 1
    expect_identical(is.na(b$lo), !ok)
    expect_identical(is.na(b$hi), !ok)
    expect_true(all(b$lo[ok] <= exact & exact <= b$hi[ok]))
    expect_equal(b$k, sum(ok))
  }
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
    exact <- HyperCom:::calculate_priority(bg, fx$output$adj, fx$output$infMat, metadata)
    b <- HyperCom:::bound_terms(bg, fx$output$adj, fx$output$infMat, metadata)
    ok <- !is.na(exact)
    expect_true(all(b$lo[ok] <= exact[ok] & exact[ok] <= b$hi[ok]))
    expect_true(all(is.na(b$lo[!ok]) & is.na(b$hi[!ok])))
    expect_identical(b$lo[, c("frac", "pct.min", "purity")], exact[, c("frac", "pct.min", "purity")])
  }
})

test_that("screened p-values are those of scoring every background pair, which few are", {
  skip_on_cran()
  factored <- generate_infMat(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL), outfile=NULL)
  every <- function(expr) testthat::with_mocked_bindings(expr, ambiguous=function(L, U, x) rep(TRUE, length(L)), .package="HyperCom")

  for(seed in 1:2){
    screened <- scored_grids(run_priority(weights=c(1, 0.5, 0), significance=TRUE, seed=seed))
    all <- scored_grids(every(run_priority(weights=c(1, 0.5, 0), significance=TRUE, seed=seed)))
    expect_identical(screened$value, all$value)
    expect_lt(screened$grids, all$grids / 2)
  }
  expect_identical(run_priority(infMat=factored, metadata=fx$output$metadata.none, significance=TRUE, seed=3),
                   every(run_priority(infMat=factored, metadata=fx$output$metadata.none, significance=TRUE, seed=3)))
})

test_that("background pairs without a priority give the p-values and warning of scoring every background pair", {
  skip_on_cran()
  metadata <- fx$output$metadata
  metadata$sample <- as.character(metadata$sample)
  metadata$sample[match(metadata$cell[metadata$sample != "gene"][1], metadata$cell)] <- "alone"
  every <- function(expr) testthat::with_mocked_bindings(expr, ambiguous=function(L, U, x) rep(TRUE, length(L)), .package="HyperCom")
  run <- function(expr){
    warnings <- character(0)
    value <- withCallingHandlers(expr, warning=function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    return(list(value=value, warnings=warnings))
  }

  screened <- run(run_priority(metadata=metadata, weights=c(1, 0.5), significance=TRUE, seed=1))
  all <- run(every(run_priority(metadata=metadata, weights=c(1, 0.5), significance=TRUE, seed=1)))

  expect_identical(screened, all)
  expect_match(screened$warnings, "^[1-9][0-9]* of 1000 background pairs have no priority \\(NaN\\)", all=FALSE)
})

test_that("seeded p-values on a cluster are those of one process", {
  cl <- test_cluster()

  expect_identical(run_priority(weights=c(1, 0.5), significance=TRUE, seed=1, parallel=TRUE, cl=cl),
                   run_priority(weights=c(1, 0.5), significance=TRUE, seed=1))
})
