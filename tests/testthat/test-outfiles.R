test_that("each step writes what it returns to outfile", {
  dir <- withr::local_tempdir()
  out <- function(name){
    file.path(dir, name)
  }
  counts <- filter_counts(fx$input$counts, fx$input$lrdb, outfile=out("0.counts.rds"))
  lrs <- generate_lrs_table(counts, fx$input$lrdb, outfile=out("1.lrs.csv"))
  hyperedges <- generate_hyperedges(counts, lrs, outfile=out("2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=out("3.adj.rds"))
  infMat <- generate_infMat(adj, outfile=out("4.infMat.rds"))
  md <- generate_metadata(infMat, fx$input$metadata, sample="sample", group="condition", outfile=out("5.metadata.csv"))
  read_csv <- function(name){
    plain_df(utils::read.csv(out(name), check.names=FALSE))
  }

  expect_identical(readRDS(out("0.counts.rds")), counts)
  expect_equal(read_csv("1.lrs.csv"), plain_df(lrs))
  expect_equal(read_csv("2.hyperedges.csv"), plain_df(hyperedges))
  expect_identical(readRDS(out("3.adj.rds")), adj)
  expect_identical(readRDS(out("4.infMat.rds")), infMat)
  expect_equal(read_csv("5.metadata.csv"), plain_df(data.frame(lapply(md, as.character), check.names=FALSE)))
})

test_that("outfile=NULL writes nothing", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  counts <- filter_counts(fx$input$counts, fx$input$lrdb, outfile=NULL)
  lrs <- generate_lrs_table(counts, fx$input$lrdb, outfile=NULL)
  hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
  adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL)
  md <- generate_metadata(infMat, fx$input$metadata, sample="sample", group="condition", outfile=NULL)
  priority <- quietly(prioritize_lr(lrs, adj, infMat, md, outdir=NULL))

  expect_identical(list.files(dir, recursive=TRUE, all.files=TRUE), character(0))
})
