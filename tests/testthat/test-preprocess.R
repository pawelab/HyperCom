test_that("filter_counts matches master", {
  counts <- filter_counts(fx$input$counts, fx$input$lrdb, outfile=NULL)

  expect_same_matrix(counts, fx$output$counts)
})

test_that("filter_counts returns cells x expressed ligand and receptor genes", {
  counts <- filter_counts(fx$input$counts, fx$input$lrdb, outfile=NULL)
  lr.genes <- unique(stats::na.omit(unlist(fx$input$lrdb[, -1])))
  expressed <- rownames(fx$input$counts)[Matrix::rowSums(fx$input$counts) > 0]

  expect_identical(rownames(counts), colnames(fx$input$counts))
  expect_setequal(colnames(counts), intersect(expressed, lr.genes))
  expect_true(all(Matrix::colSums(counts) > 0))
})

test_that("filter_counts gives the same values for dense and sparse input", {
  counts <- filter_counts(as.matrix(fx$input$counts), fx$input$lrdb, outfile=NULL)

  expect_same_matrix(counts, fx$output$counts)
})

test_that("generate_lrs_table matches master", {
  lrs <- generate_lrs_table(fx$output$counts, fx$input$lrdb, outfile=NULL)

  expect_s3_class(lrs, "data.frame")
  expect_equal(plain_df(lrs), plain_df(fx$output$lrs))
})

test_that("generate_lrs_table keeps interactions whose genes are all expressed", {
  lrs <- generate_lrs_table(fx$output$counts, fx$input$lrdb, outfile=NULL)
  genes <- colnames(fx$output$counts)
  all_present <- function(df){
    apply(df[, -1, drop=FALSE], 1, function(x) all(stats::na.omit(x) %in% genes))
  }
  dropped <- fx$input$lrdb[!fx$input$lrdb$interaction %in% lrs$interaction, ]

  expect_true(all(all_present(lrs)))
  expect_false(any(all_present(dropped)))
  expect_false(any(vapply(lrs, function(x) all(is.na(x)), logical(1))))
  expect_identical(colnames(lrs)[1:2], c("interaction", "ligand"))
})
