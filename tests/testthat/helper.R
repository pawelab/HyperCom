# The fixture (fixtures/master.rds, made by fixtures/make_fixtures.R with 0.0.0.9000 installed) holds a 160-cell input set and the
# output of every step of 0.0.0.9000 on it. Its node order is that of the C collation, so the collation is pinned for the whole run.
withr::local_collate("C", .local_envir=testthat::teardown_env())
fx <- readRDS(testthat::test_path("fixtures", "master.rds"))

#' Run an expression without its messages
#'
#' `quietly` returns the value of an expression and drops its messages
#'
#' This function suppresses progress messages so they do not flood the test log. Warnings pass through.
#' @param expr An expression
#' @returns The value of the expression
#' @noRd
quietly <- function(expr){
  return(suppressMessages(expr))
}

#' Drop the row names and tibble class of a dataframe
#'
#' `plain_df` returns a dataframe without row names or tibble class
#'
#' This function takes a dataframe or tibble and returns a plain dataframe with default row names, so tables made in different ways
#' can be compared by value.
#' @param df A dataframe or tibble
#' @returns A dataframe
#' @noRd
plain_df <- function(df){
  df <- as.data.frame(df)
  rownames(df) <- NULL

  return(df)
}

#' Expect cell scores to match by cell
#'
#' `expect_same_scores` checks two score_hypercom results cell by cell
#'
#' This function sorts both results by cell before comparing, as exact ties in the HyperCom score may order cells differently.
#' @param object A dataframe from score_hypercom
#' @param expected A dataframe from score_hypercom
#' @noRd
expect_same_scores <- function(object, expected){
  testthat::expect_identical(plain_df(object[order(object$cell), ]), plain_df(expected[order(expected$cell), ]))
}

#' Run prioritize_lr on the fixture
#'
#' `run_priority` returns prioritize_lr's result on the fixture's outputs, without files or messages
#'
#' This function runs prioritize_lr on the fixture's LR table, adjacency, influence matrix and metadata; any of them and any other
#' argument can be replaced.
#' @param lrs,adj,infMat,metadata Inputs of prioritize_lr
#' @param ... Other arguments of prioritize_lr
#' @returns A list of dataframes named by weight
#' @noRd
run_priority <- function(lrs=fx$output$lrs, adj=fx$output$adj, infMat=fx$output$infMat, metadata=fx$output$metadata, ...){
  return(quietly(prioritize_lr(lrs, adj, infMat, metadata, outdir=NULL, ...)))
}

#' Run HyperCom on the fixture
#'
#' `run_hypercom` returns HyperCom's result on the fixture's inputs, without messages
#'
#' This function runs HyperCom on the fixture's lrdb, counts and metadata with sample and group set; any of them and any other
#' argument can be replaced.
#' @param outdir A string naming the output directory
#' @param lrdb,counts,metadata Inputs of HyperCom
#' @param ... Other arguments of HyperCom
#' @returns The list HyperCom returns
#' @noRd
run_hypercom <- function(outdir, lrdb=fx$input$lrdb, counts=fx$input$counts, metadata=fx$input$metadata, ...){
  return(quietly(HyperCom(lrdb, counts, metadata, sample="sample", group="condition", outdir=outdir, ...)))
}
