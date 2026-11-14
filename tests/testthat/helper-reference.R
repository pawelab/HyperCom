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
