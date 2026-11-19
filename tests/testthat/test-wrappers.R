expect_run_matches_handwritten <- function(rundir){
  files <- c("0.counts.rds", "1.lrs.csv", "2.hyperedges.csv", "3.adj.rds", "4.infMat.rds", "5.metadata.csv", "6.priority/priority.1.csv")
  testthat::expect_true(all(file.exists(file.path(rundir, files))))
  expect_same_matrix(readRDS(file.path(rundir, "0.counts.rds")), fx$output$counts)
  expect_same_matrix(readRDS(file.path(rundir, "3.adj.rds")), fx$output$adj)
  expect_same_matrix(readRDS(file.path(rundir, "4.infMat.rds")), fx$output$infMat, tolerance=tol)
  priority <- utils::read.csv(file.path(rundir, "6.priority", "priority.1.csv"), check.names=FALSE)
  expect_same_priority(priority, fx$output$priority[["1"]])
}

test_that("HyperCom runs the whole pipeline with defaults and reproduces handwritten", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  res <- run_hypercom("run")

  expect_run_result(res, "run")
  expect_run_matches_handwritten(file.path(dir, "run"))
  expect_invisible(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition", outdir="run2",
                            verbose=FALSE))
})

test_that("HyperCom_seurat reproduces handwritten from a Seurat object", {
  skip_if_not_installed("SeuratObject")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  seurat <- suppressWarnings(SeuratObject::CreateSeuratObject(counts=fx$input$counts, meta.data=fx$input$metadata))
  res <- quietly(HyperCom_seurat(seurat, fx$input$lrdb, sample="sample", group="condition", outdir="run"))

  expect_run_result(res, "run")
  expect_run_matches_handwritten(file.path(dir, "run"))
})

test_that("get_seurat_layer returns the counts of v3 and v5 assays", {
  skip_if_not_installed("SeuratObject")
  for(version in c("v3", "v5")){
    withr::local_options(Seurat.object.assay.version=version)
    seurat <- suppressWarnings(SeuratObject::CreateSeuratObject(counts=fx$input$counts))
    counts <- HyperCom:::get_seurat_layer(seurat, assay="RNA", layer="counts")

    expect_same_matrix(counts, fx$input$counts)
  }
})
