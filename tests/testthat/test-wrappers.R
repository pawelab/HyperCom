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

test_that("HyperCom_sce reproduces handwritten from a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=fx$input$counts), colData=fx$input$metadata)
  res <- quietly(HyperCom_sce(sce, fx$input$lrdb, sample="sample", group="condition", outdir="run"))

  expect_run_result(res, "run")
  expect_run_matches_handwritten(file.path(dir, "run"))
})

test_that("get_sce_assay reads delayed assays for the lrdb genes and stops on bad input", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("DelayedArray")
  counts <- fx$input$counts
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=DelayedArray::DelayedArray(as.matrix(counts))))
  delayed <- HyperCom:::get_sce_assay(sce, fx$input$lrdb, assay="counts")
  genes <- unique(unlist(lapply(fx$input$lrdb[, -1, drop=FALSE], as.character), use.names=FALSE))

  expect_s4_class(delayed, "dgCMatrix")
  expect_equal(colnames(delayed), colnames(counts))
  expect_true(all(rownames(delayed) %in% genes))
  expect_same_matrix(filter_counts(delayed, fx$input$lrdb, outfile=NULL), filter_counts(counts, fx$input$lrdb, outfile=NULL))
  expect_error(HyperCom:::get_sce_assay(sce, fx$input$lrdb, assay="logcounts"), "no assay 'logcounts'")
  expect_error(HyperCom:::get_sce_assay(unname(sce), fx$input$lrdb), "no cell names")
  expect_error(HyperCom:::get_sce_assay(counts, fx$input$lrdb), "must be a SingleCellExperiment")
  ensembl <- sce
  rownames(ensembl) <- paste0("ENSMUSG", seq_len(nrow(ensembl)))
  expect_error(HyperCom:::get_sce_assay(ensembl, fx$input$lrdb), "must be gene symbols")
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
