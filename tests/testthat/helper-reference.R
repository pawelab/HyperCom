#' Read the master fixture
#'
#' `read_fixture` returns the small input set and master's output for every step
#'
#' This function reads fixtures/master.rds, made by fixtures/make_fixtures.R with the master branch installed.
#' @returns A list with about, input (counts, metadata, lrdb) and output (one element per pipeline step)
#' @noRd
read_fixture <- function(){
  return(readRDS(testthat::test_path("fixtures", "master.rds")))
}

withr::local_collate("C", .local_envir=testthat::teardown_env())
fx <- read_fixture()
tol <- 1e-10
x87 <- identical(.Machine$longdouble.digits, 64L)

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

#' Expect a score or priority result to equal master's in the fixture
#'
#' `expect_master` checks a result of the grid kernel against master's output in the fixture
#'
#' This function compares bit for bit where long double is the 80-bit x87 type of the x86-64 machine the fixture was made on, and
#' within 1e-12 elsewhere: R's means and the kernel's sums use the platform's long double (double on arm64 macOS and Windows, 128 bits
#' on arm64 Linux, absent on builds without long double), so their last bits differ from the fixture's there.
#' @param object A result from the build under test
#' @param expected Master's result from the fixture
#' @noRd
expect_master <- function(object, expected){
  if(x87){
    testthat::expect_identical(object, expected)
  } else {
    testthat::expect_equal(object, expected, tolerance=1e-12)
  }
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

#' Get master's adjacency object from any build's adjacency
#'
#' `as_master_adjacency` returns the adjacency as master's dsCMatrix
#'
#' This function takes an adjacency from generate_adjacency_matrix: master's dsCMatrix is returned as is, the incidence kept by
#' builds without node x node matrices is expanded to it.
#' @param adj An adjacency from generate_adjacency_matrix
#' @returns A dsCMatrix
#' @noRd
as_master_adjacency <- function(adj){
  if(inherits(adj, "hypercom_adjacency")){
    return(HyperCom:::adjacency_sparse(adj))
  }

  return(adj)
}

#' Expect a matrix to hold the same values as master's
#'
#' `expect_same_matrix` checks dimnames exactly and values exactly or within a tolerance
#'
#' This function compares dense or sparse matrices by value, so storage class and explicit zeros do not matter.
#' @param object A matrix from the build under test
#' @param expected A matrix from master
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

#' Calculate master's influence matrix on this machine
#'
#' `master_infmat` returns the influence matrix as master's generate_infMat makes it
#'
#' This function is master b999a01's generate_infMat without the file write. The last bits of its result depend on the LAPACK and BLAS
#' R uses (the fixture was made with OpenBLAS on Windows; R's reference BLAS differs by ~3e-15), so a build is compared bit for bit
#' with this function run on the same machine.
#' @param adj An adjacency matrix with node names
#' @param r A double between 0 and 1 representing the restart probability
#' @returns A dense influence matrix with node names
#' @noRd
master_infmat <- function(adj, r=0.5){
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

#' Expect an influence matrix to be master's
#'
#' `expect_master_infmat` checks an influence matrix against master's formula on this machine and against the fixture
#'
#' This function compares bit for bit with master_infmat of the fixture's adjacency, and within 1e-12 with master's matrix in the
#' fixture, whose last bits depend on the BLAS it was made with.
#' @param object An influence matrix from the build under test
#' @param r A double, 0.5 or 0.3 (the restart probabilities in the fixture)
#' @noRd
expect_master_infmat <- function(object, r=0.5){
  expected <- if(r == 0.5) fx$output$infMat else fx$output$infMat.r03
  testthat::expect_identical(object, master_infmat(fx$output$adj, r))
  expect_same_matrix(object, expected, tolerance=1e-12)
}

#' Make a dataframe comparable
#'
#' `plain_df` returns a plain dataframe without row names
#'
#' This function drops tibble classes and row names, which master does not write to its output files.
#' @param df A dataframe or tibble
#' @returns A dataframe
#' @noRd
plain_df <- function(df){
  df <- as.data.frame(df, stringsAsFactors=FALSE)
  rownames(df) <- NULL
  return(df)
}

#' Expect HyperCom scores to match master's
#'
#' `expect_same_scores` matches cells and compares every score column within a tolerance
#'
#' This function compares by cell, so the order of exactly tied scores does not matter, and checks the descending order separately.
#' @param object A dataframe from score_hypercom of the build under test
#' @param expected A dataframe from score_hypercom of master
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

#' Expect LR priorities to match master's
#'
#' `expect_same_priority` compares the ordered LR table, priorities and p-values
#'
#' This function checks the pair order exactly and the priority and p-value columns within a tolerance.
#' @param object A dataframe from prioritize_lr of the build under test (one weight)
#' @param expected A dataframe from prioritize_lr of master (one weight)
#' @param tolerance A numeric relative tolerance
#' @noRd
expect_same_priority <- function(object, expected, tolerance=tol){
  testthat::expect_identical(colnames(object), colnames(expected))
  testthat::expect_identical(object$interaction, expected$interaction)
  testthat::expect_equal(plain_df(object), plain_df(expected), tolerance=tolerance)
}

#' Skip a test of a known master bug unless asked
#'
#' `skip_known_bug` skips tests that state the intended behaviour where master is known to differ
#'
#' This function runs the test only when the environment variable HYPERCOM_KNOWN_BUGS is true, so the suite passes on master.
#' @param id A string naming the finding (see the project findings list)
#' @noRd
skip_known_bug <- function(id){
  if(!identical(tolower(Sys.getenv("HYPERCOM_KNOWN_BUGS")), "true")){
    testthat::skip(paste0("known master bug ", id, " (set HYPERCOM_KNOWN_BUGS=true to run)"))
  }
}
