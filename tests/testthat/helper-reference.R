# Shared test helpers: the handwritten fixture (b999a01's inputs and outputs), runs of prioritize_lr and HyperCom on it, counters of
# the grids a run scores, and expectations against handwritten's outputs.

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

# The fixture's node order is that of the C collation, so it is pinned for the whole run.
withr::local_collate("C", .local_envir=testthat::teardown_env())
fx <- read_fixture()
tol <- 1e-10
fx_pairs <- list(network=list(NULL, NULL), single=list("Thbs1", "Sdc4"), complex=list("Tgfb1", c("Tgfbr1", "Tgfbr2")),
                 triple=list("Tgfb1", c("Acvr1", "Tgfbr1", "Tgfbr2")), same=as.list(fx$output$same.pair))

#' Cap a process count at the available cores
#'
#' `test_threads` returns a process count (the threads argument) no larger than the available cores
#'
#' This function takes the number of processes a test asks for and returns it, or the available cores if fewer (2 under R CMD check),
#' so tests of several processes do not warn on small machines.
#' @param n An integer number of processes
#' @returns An integer number of processes
#' @noRd
test_threads <- function(n){
  return(min(n, parallelly::availableCores()))
}

#' Make the factored influence matrix of the fixture once
#'
#' `fx_factored` returns the fixture's adjacency and its influence matrix by method "woodbury", made on the first call and kept
#' @returns A list with adj (a hypercom_adjacency) and infMat (a hypercom_influence)
#' @noRd
fx_factored <- local({
  made <- NULL
  function(){
    if(is.null(made)){
      adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
      made <<- list(adj=adj, infMat=generate_infMat(adj, outfile=NULL))
    }
    return(made)
  }
})

#' Expect a score or priority result to equal handwritten's in the fixture
#'
#' `expect_handwritten` checks a score or priority result against handwritten's output in the fixture
#'
#' This function compares bit for bit: the R grids rank and average cell scores as handwritten did (dplyr::percent_rank, rowMeans and
#' colMeans), so from the same dense influence matrix the scores, priorities and seeded p-values are identical to handwritten's.
#' Results from a "woodbury" or "chol" influence matrix differ by rounding and are compared with a tolerance instead.
#' @param object A result from the build under test
#' @param expected Handwritten's result from the fixture
#' @noRd
expect_handwritten <- function(object, expected){
  testthat::expect_identical(object, expected)
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

#' Run an expression and collect its messages and warnings
#'
#' `collect` returns the value of an expression with the messages and warnings it gave, which are not shown
#' @param expr An expression
#' @returns A list with value, messages and warnings
#' @noRd
collect <- function(expr){
  messages <- character(0)
  warnings <- character(0)
  value <- withCallingHandlers(expr, message=function(m){
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  }, warning=function(w){
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  return(list(value=value, messages=messages, warnings=warnings))
}

#' Count the grids a run scores
#'
#' `count_grids` returns the value of an expression and the number of grids passed to pool_grids while it ran, in one process or in
#' worker processes
#' @param expr An expression running prioritize_lr or HyperCom
#' @returns A list with value and grids
#' @noRd
count_grids <- function(expr){
  real <- HyperCom:::pool_grids
  grids <- 0
  counter <- function(pool, infMat, metadata, batches){
    grids <<- grids + sum(vapply(batches, nrow, integer(1)))
    return(real(pool, infMat, metadata, batches))
  }
  testthat::local_mocked_bindings(pool_grids=counter, .package="HyperCom")
  value <- expr

  return(list(value=value, grids=grids))
}

#' Run prioritize_lr on the fixture
#'
#' `run_priority` returns the priorities of prioritize_lr, by default on handwritten's inputs in the fixture, without writing them
#'
#' This function takes any input that replaces the fixture's and further arguments of prioritize_lr, and runs it quietly with the
#' package defaults, except one process unless threads is given (worker processes take seconds to start).
#' @param lrs A data frame of ligand-receptor interactions
#' @param adj An adjacency
#' @param infMat An influence matrix
#' @param metadata A data frame containing columns named cell, sample and group
#' @param ... Further arguments of prioritize_lr
#' @param threads The threads argument of prioritize_lr, 1 by default
#' @param outdir A string naming the output directory, NULL to write nothing
#' @returns A list of data frames named by weight
#' @noRd
run_priority <- function(lrs=fx$output$lrs, adj=fx$output$adj, infMat=fx$output$infMat, metadata=fx$output$metadata, ..., threads=1,
                         outdir=NULL){
  return(quietly(prioritize_lr(lrs, adj, infMat, metadata, ..., threads=threads, outdir=outdir)))
}

#' Run HyperCom on the fixture
#'
#' `run_hypercom` returns the results of HyperCom on the fixture's input, by sample and condition
#'
#' This function takes the output directory, any input that replaces the fixture's and further arguments of HyperCom, and runs it
#' quietly with the package defaults, except one process unless threads is given.
#' @param outdir A string naming the output directory, or NULL
#' @param ... Further arguments of HyperCom
#' @param threads The threads argument of HyperCom, 1 by default
#' @param counts A matrix with row names being genes and column names being cells
#' @param metadata A data frame containing cell metadata with row names being cells
#' @returns The list HyperCom returns
#' @noRd
run_hypercom <- function(outdir, ..., threads=1, counts=fx$input$counts, metadata=fx$input$metadata){
  return(quietly(HyperCom(fx$input$lrdb, counts, metadata, sample="sample", group="condition", outdir=outdir, threads=threads, ...)))
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
  testthat::expect_identical(unname(dimnames(object)), unname(dimnames(expected)))
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
#' This function is handwritten b999a01's generate_infMat without the file write. The last bits of its result depend on the LAPACK and
#' BLAS R uses, so a build is compared bit for bit with this function run on the same machine.
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
#' This function compares bit for bit with handwritten_infmat of the fixture's adjacency, and within 1e-12 with handwritten's matrix in
#' the fixture, whose last bits depend on the BLAS it was made with.
#' @param object An influence matrix from the build under test
#' @param r A double, 0.5 or 0.3 (the restart probabilities in the fixture)
#' @noRd
expect_handwritten_infmat <- function(object, r=0.5){
  expected <- switch(as.character(r), "0.5"=fx$output$infMat, "0.3"=fx$output$infMat.r03,
                     stop("the fixture has influence matrices for r = 0.5 and 0.3 only, not ", r))
  testthat::expect_identical(object, handwritten_infmat(fx$output$adj, r))
  expect_same_matrix(object, expected, tolerance=1e-12)
}

#' Make a data frame comparable
#'
#' `plain_df` returns a plain data frame without row names
#'
#' This function drops tibble classes and row names, which handwritten does not write to its output files.
#' @param df A data frame or tibble
#' @returns A data frame
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
#' @param object A data frame from score_hypercom of the build under test
#' @param expected A data frame from score_hypercom of handwritten
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
#' @param object A data frame from prioritize_lr of the build under test (one weight)
#' @param expected A data frame from prioritize_lr of handwritten (one weight)
#' @param tolerance A numeric relative tolerance
#' @noRd
expect_same_priority <- function(object, expected, tolerance=tol){
  testthat::expect_identical(colnames(object), colnames(expected))
  testthat::expect_identical(object$interaction, expected$interaction)
  testthat::expect_equal(plain_df(object), plain_df(expected), tolerance=tolerance)
}

#' Expect the result of a HyperCom run to be what it wrote and handwritten's
#'
#' `expect_run_result` checks the list HyperCom returns, and the files of the run, against handwritten's output in the fixture
#'
#' This function checks the names of the result, the priorities of weight 1, the metadata's cells and the interactions against
#' handwritten's; with a run directory also that every step's file exists, that the adjacency and influence matrix returned are the
#' ones written, and that the counts, adjacency, influence matrix and priority files are handwritten's.
#' @param res The list HyperCom returned
#' @param rundir A string naming the output directory of the run, or NULL for a run without files
#' @noRd
expect_run_result <- function(res, rundir=NULL){
  testthat::expect_named(res, c("priority", "metadata", "infMat", "adj", "lrs"))
  expect_same_priority(res$priority[["1"]], fx$output$priority[["1"]])
  testthat::expect_identical(res$metadata$cell, fx$output$metadata$cell)
  testthat::expect_equal(plain_df(res$lrs), plain_df(fx$output$lrs))
  if(is.null(rundir)){
    return(invisible(res))
  }

  files <- c("0.counts.rds", "1.lrs.csv", "2.hyperedges.csv", "3.adj.rds", "4.infMat.rds", "5.metadata.csv", "6.priority/priority.1.csv")
  testthat::expect_true(all(file.exists(file.path(rundir, files))))
  testthat::expect_identical(res$adj, readRDS(file.path(rundir, "3.adj.rds")))
  testthat::expect_identical(res$infMat, readRDS(file.path(rundir, "4.infMat.rds")))
  expect_same_matrix(readRDS(file.path(rundir, "0.counts.rds")), fx$output$counts)
  expect_same_matrix(res$adj, fx$output$adj)
  expect_same_matrix(res$infMat, fx$output$infMat, tolerance=tol)
  expect_same_priority(utils::read.csv(file.path(rundir, "6.priority", "priority.1.csv"), check.names=FALSE), fx$output$priority[["1"]])

  return(invisible(res))
}
