#' Read the handwritten fixture
#'
#' `read_fixture` returns the small input set and handwritten's output for every step
#'
#' This function reads fixtures/handwritten.rds, made by fixtures/make_fixtures.R with the handwritten branch installed.
#' @returns A list with about, input (counts, metadata, lrdb) and output (one element per pipeline step)
#' @noRd
read_fixture <- function(){
  return(readRDS(testthat::test_path("fixtures", "handwritten.rds")))
}

# testthat 3e also runs each test under C collation, but the fixture's node order depends on it (finding #15), so it is pinned
# explicitly for the whole run.
withr::local_collate("C", .local_envir=testthat::teardown_env())
fx <- read_fixture()
tol <- 1e-10
fx_pairs <- list(network=list(NULL, NULL), single=list("Thbs1", "Sdc4"), complex=list("Tgfb1", c("Tgfbr1", "Tgfbr2")),
                 triple=list("Tgfb1", c("Acvr1", "Tgfbr1", "Tgfbr2")), same=as.list(fx$output$same.pair))

#' Expect a score or priority result to equal handwritten's in the fixture
#'
#' `expect_handwritten` checks a score or priority result against handwritten's output in the fixture
#'
#' This function compares bit for bit: the base-R grids rank and average cell scores as handwritten did (dplyr::percent_rank, rowMeans
#' and colMeans), so from the same dense influence matrix the scores, priorities and seeded p-values are identical to handwritten's.
#' Results from a "woodbury" or "chol" influence matrix differ by rounding and are compared with a tolerance instead.
#' @param object A result from the build under test
#' @param expected Handwritten's result from the fixture
#' @noRd
expect_handwritten <- function(object, expected){
  testthat::expect_identical(object, expected)
}

test.env <- new.env()

#' Get the test cluster
#'
#' `test_cluster` returns a PSOCK cluster of two workers that load the HyperCom under test
#'
#' This function skips the test on CRAN, makes the cluster on first use (its workers get R_LIBS set to this session's library paths)
#' and returns the same cluster to every later test; it is stopped when the test run ends.
#' @returns A cluster object
#' @noRd
test_cluster <- function(){
  testthat::skip_on_cran()
  if(is.null(test.env$cl)){
    test.env$cl <- parallelly::makeClusterPSOCK(2, autoStop=TRUE, rscript_envs=c(R_LIBS=paste(.libPaths(), collapse=.Platform$path.sep)))
    withr::defer(parallel::stopCluster(test.env$cl), envir=testthat::teardown_env())
  }

  return(test.env$cl)
}

#' Run an expression without its printed output
#'
#' `quietly` returns the value of an expression and drops what it prints
#'
#' This function captures printed lines and messages so per-pair progress does not flood the test log. Warnings pass through.
#' @param expr An expression
#' @returns The value of the expression
#' @noRd
quietly <- function(expr){
  utils::capture.output(x <- suppressMessages(expr))
  return(x)
}

#' Run prioritize_lr on the fixture
#'
#' `run_priority` returns the priorities of prioritize_lr, by default on handwritten's inputs in the fixture, without writing them
#'
#' This function takes any input that replaces the fixture's and further arguments of prioritize_lr, and runs it quietly.
#' @param lrs A dataframe of ligand-receptor interactions
#' @param adj An adjacency
#' @param infMat An influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ... Further arguments of prioritize_lr
#' @param outdir A string naming the output directory, NULL to write nothing
#' @returns A list of dataframes named by weight
#' @noRd
run_priority <- function(lrs=fx$output$lrs, adj=fx$output$adj, infMat=fx$output$infMat, metadata=fx$output$metadata, ..., outdir=NULL){
  return(quietly(prioritize_lr(lrs, adj, infMat, metadata, ..., outdir=outdir)))
}

#' Run HyperCom on the fixture
#'
#' `run_hypercom` returns the results of HyperCom on the fixture's input, by sample and condition
#'
#' This function takes the output directory, any input that replaces the fixture's and further arguments of HyperCom, and runs it
#' quietly.
#' @param outdir A string naming the output directory
#' @param ... Further arguments of HyperCom
#' @param counts A matrix with row names being genes and column names being cells
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @returns The list HyperCom returns
#' @noRd
run_hypercom <- function(outdir, ..., counts=fx$input$counts, metadata=fx$input$metadata){
  return(quietly(HyperCom(fx$input$lrdb, counts, metadata, sample="sample", group="condition", outdir=outdir, ...)))
}

#' Expect a matrix to hold the same values as handwritten's
#'
#' `expect_same_matrix` checks dimnames exactly and values exactly or within a tolerance
#'
#' This function compares dense or sparse matrices by value, so storage class and explicit zeros do not matter.
#' @param object A matrix from the build under test
#' @param expected A matrix from handwritten
#' @param tolerance A numeric relative tolerance, 0 for exact
#' @noRd
expect_same_matrix <- function(object, expected, tolerance=0){
  testthat::expect_identical(dimnames(object), dimnames(expected))
  if(tolerance == 0){
    testthat::expect_identical(unname(as.matrix(object)), unname(as.matrix(expected)))
  } else {
    testthat::expect_equal(unname(as.matrix(object)), unname(as.matrix(expected)), tolerance=tolerance)
  }
}

#' Calculate handwritten's influence matrix on this machine
#'
#' `handwritten_infmat` returns the influence matrix as handwritten's generate_infMat makes it
#'
#' This function is handwritten b999a01's generate_infMat without the file write. The last bits of its result depend on the LAPACK and BLAS
#' R uses (the fixture was made with OpenBLAS on Windows; R's reference BLAS differs by ~3e-15), so a build is compared bit for bit
#' with this function run on the same machine.
#' @param adj An adjacency matrix with node names
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A dense influence matrix with node names
#' @noRd
handwritten_infmat <- function(adj, r=0.5){
  len <- dim(adj)[1]
  D <- as.matrix(Matrix::Diagonal(len))
  W <- adj/Matrix::rowSums(adj)
  W <- (1 - r) * W
  W <- as.matrix(W)
  M <- D - W
  infMat <- r * Matrix::solve(M)
  rownames(infMat) <- rownames(adj)
  colnames(infMat) <- rownames(adj)

  return(infMat)
}

#' Expect an influence matrix to be handwritten's
#'
#' `expect_handwritten_infmat` checks an influence matrix against handwritten's formula on this machine and against the fixture
#'
#' This function compares bit for bit with handwritten_infmat of the fixture's adjacency, and within 1e-12 with handwritten's matrix in the
#' fixture, whose last bits depend on the BLAS it was made with.
#' @param object An influence matrix from the build under test
#' @param r A double, 0.5 or 0.3 (the restart probabilities in the fixture)
#' @noRd
expect_handwritten_infmat <- function(object, r=0.5){
  expected <- switch(as.character(r), "0.5"=fx$output$infMat, "0.3"=fx$output$infMat.r03,
                     stop("the fixture has influence matrices for r = 0.5 and 0.3 only, not ", r))
  testthat::expect_identical(object, handwritten_infmat(fx$output$adj, r))
  expect_same_matrix(object, expected, tolerance=1e-12)
}

#' Make a dataframe comparable
#'
#' `plain_df` returns a plain dataframe without row names
#'
#' This function drops tibble classes and row names, which handwritten does not write to its output files.
#' @param df A dataframe or tibble
#' @returns A dataframe
#' @noRd
plain_df <- function(df){
  df <- as.data.frame(df, stringsAsFactors=FALSE)
  rownames(df) <- NULL
  return(df)
}

#' Expect HyperCom scores to match handwritten's
#'
#' `expect_same_scores` matches cells and compares every score column within a tolerance
#'
#' This function compares by cell, so the order of exactly tied scores does not matter, and checks the descending order separately.
#' @param object A dataframe from score_hypercom of the build under test
#' @param expected A dataframe from score_hypercom of handwritten
#' @param tolerance A numeric relative tolerance
#' @noRd
expect_same_scores <- function(object, expected, tolerance=tol){
  cols <- c("cell", "sender", "receiver", "dif", "avg", "max", "hypercom")
  testthat::expect_setequal(colnames(object), cols)
  testthat::expect_setequal(object$cell, expected$cell)
  a <- plain_df(object)[match(expected$cell, object$cell), cols]
  b <- plain_df(expected)[, cols]
  rownames(a) <- NULL
  testthat::expect_equal(a, b, tolerance=tolerance)
  testthat::expect_false(is.unsorted(rev(object$hypercom), na.rm=TRUE))
}

#' Expect LR priorities to match handwritten's
#'
#' `expect_same_priority` compares the ordered LR table, priorities and p-values
#'
#' This function checks the pair order exactly and the priority and p-value columns within a tolerance.
#' @param object A dataframe from prioritize_lr of the build under test (one weight)
#' @param expected A dataframe from prioritize_lr of handwritten (one weight)
#' @param tolerance A numeric relative tolerance
#' @noRd
expect_same_priority <- function(object, expected, tolerance=tol){
  testthat::expect_identical(colnames(object), colnames(expected))
  testthat::expect_identical(object$interaction, expected$interaction)
  testthat::expect_equal(plain_df(object), plain_df(expected), tolerance=tolerance)
}

#' Expect the result of a HyperCom run to be what it wrote and handwritten's
#'
#' `expect_run_result` checks the list HyperCom returns against the files of the run and against handwritten's output in the fixture
#'
#' This function checks the names of the result, that the adjacency and influence matrix returned are the ones written, and that the
#' priorities of weight 1, the metadata's cells and the interactions are handwritten's.
#' @param res The list HyperCom returned
#' @param rundir A string naming the output directory of the run
#' @noRd
expect_run_result <- function(res, rundir){
  testthat::expect_named(res, c("priority", "metadata", "infMat", "adj", "lrs"))
  testthat::expect_identical(res$adj, readRDS(file.path(rundir, "3.adj.rds")))
  testthat::expect_identical(res$infMat, readRDS(file.path(rundir, "4.infMat.rds")))
  expect_same_priority(res$priority[["1"]], fx$output$priority[["1"]])
  testthat::expect_identical(res$metadata$cell, fx$output$metadata$cell)
  testthat::expect_equal(plain_df(res$lrs), plain_df(fx$output$lrs))
}
