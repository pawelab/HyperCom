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

#' Cap a thread count at the available cores
#'
#' `test_threads` returns a thread count no larger than the available cores
#'
#' This function takes the number of threads a test asks for and returns it, or the available cores if fewer (2 under R CMD check),
#' so tests of several threads do not warn on small machines.
#' @param n An integer number of threads
#' @returns An integer number of threads
#' @noRd
test_threads <- function(n){
  return(min(n, parallelly::availableCores()))
}

#' Expect a score or priority result to equal handwritten's in the fixture
#'
#' `expect_handwritten` checks a result of the grid kernel against handwritten's output in the fixture
#'
#' This function compares within 1e-12. Since 0.2.1 the kernel divides each score's exact integer rank sum once (approved by J. Moy,
#' 2026-11-20), so cell scores lie within one ulp of handwritten's long double means on every computer; priorities, p-values and their
#' order were unchanged on the bundled and xenium test data.
#' @param object A result from the build under test
#' @param expected Handwritten's result from the fixture
#' @noRd
expect_handwritten <- function(object, expected){
  testthat::expect_equal(object, expected, tolerance=1e-12)
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
#' This function takes any input that replaces the fixture's and further arguments of prioritize_lr, and runs it quietly. Every real
#' pair is scored exactly unless exact.priority is FALSE (prioritize_lr's default screens them), so priorities compare with handwritten's.
#' @param lrs A dataframe of ligand-receptor interactions
#' @param adj An adjacency
#' @param infMat An influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ... Further arguments of prioritize_lr
#' @param outdir A string naming the output directory, NULL to write nothing
#' @param exact.priority A boolean of whether to score every real pair exactly
#' @returns A list of dataframes named by weight
#' @noRd
run_priority <- function(lrs=fx$output$lrs, adj=fx$output$adj, infMat=fx$output$infMat, metadata=fx$output$metadata, ..., outdir=NULL,
                         exact.priority=TRUE){
  return(quietly(prioritize_lr(lrs, adj, infMat, metadata, ..., outdir=outdir, exact.priority=exact.priority)))
}

#' Run HyperCom on the fixture
#'
#' `run_hypercom` returns the results of HyperCom on the fixture's input, by sample and condition
#'
#' This function takes the output directory, any input that replaces the fixture's and further arguments of HyperCom, and runs it
#' quietly, every real pair scored exactly unless exact.priority is FALSE.
#' @param outdir A string naming the output directory
#' @param ... Further arguments of HyperCom
#' @param counts A matrix with row names being genes and column names being cells
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param exact.priority A boolean of whether to score every real pair exactly
#' @returns The list HyperCom returns
#' @noRd
run_hypercom <- function(outdir, ..., counts=fx$input$counts, metadata=fx$input$metadata, exact.priority=TRUE){
  return(quietly(HyperCom(fx$input$lrdb, counts, metadata, sample="sample", group="condition", outdir=outdir, exact.priority=exact.priority,
                          ...)))
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
