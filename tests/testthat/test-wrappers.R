# HyperCom and its wrappers for Seurat and SingleCellExperiment objects run the whole pipeline and reproduce handwritten.

test_that("HyperCom runs the whole pipeline with defaults, returns its results invisibly and keeps the working directory (#6)", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)

  expect_invisible(res <- quietly(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition",
                                           outdir="run")))
  expect_run_result(res, "run")
  expect_identical(normalizePath(getwd()), normalizePath(dir))
})

test_that("HyperCom_seurat and HyperCom_sce reproduce handwritten from a Seurat object and a SingleCellExperiment", {
  wrappers <- list(SeuratObject=function(){
    seurat <- suppressWarnings(SeuratObject::CreateSeuratObject(counts=fx$input$counts, meta.data=fx$input$metadata))
    return(HyperCom_seurat(seurat, fx$input$lrdb, sample="sample", group="condition", outdir=rundir))
  }, SingleCellExperiment=function(){
    sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=fx$input$counts), colData=fx$input$metadata)
    return(HyperCom_sce(sce, fx$input$lrdb, sample="sample", group="condition", outdir=rundir))
  })
  for(pkg in names(wrappers)){
    if(!requireNamespace(pkg, quietly=TRUE)){
      next
    }
    rundir <- file.path(withr::local_tempdir(), "run")

    expect_run_result(quietly(wrappers[[pkg]]()), rundir)
  }
})

test_that("the wrappers pass their other arguments on to HyperCom, and unknown ones stop", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=fx$input$counts), colData=fx$input$metadata)
  res <- quietly(HyperCom_sce(sce, fx$input$lrdb, sample="sample", group="condition", r=0.3, infmat.method="lu", outdir=NULL))

  expect_handwritten_infmat(res$infMat, r=0.3)
  expect_error(HyperCom_sce(sce, fx$input$lrdb, sample="sample", parallel=TRUE, outdir=NULL), "unused argument")
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
  expect_error(HyperCom:::get_sce_assay(sce, fx$input$lrdb, assay="logcounts"), "no assay \"logcounts\"")
  expect_error(HyperCom:::get_sce_assay(unname(sce), fx$input$lrdb), "no cell names")
  expect_error(HyperCom:::get_sce_assay(counts, fx$input$lrdb), "must be a SingleCellExperiment")
  ensembl <- sce
  rownames(ensembl) <- paste0("ENSMUSG", seq_len(nrow(ensembl)))
  expect_error(HyperCom:::get_sce_assay(ensembl, fx$input$lrdb), "must be gene symbols")
})

test_that("get_seurat_layer returns the counts of v3 and v5 assays and names the layers when one is missing", {
  skip_if_not_installed("SeuratObject")
  missing <- c(v3="is not a layer of assay \"RNA\"", v5="no layer matching \"spliced\" in assay \"RNA\"")
  for(version in c("v3", "v5")){
    withr::local_options(Seurat.object.assay.version=version)
    seurat <- suppressWarnings(SeuratObject::CreateSeuratObject(counts=fx$input$counts))

    expect_same_matrix(HyperCom:::get_seurat_layer(seurat, assay="RNA", layer="counts"), fx$input$counts)
    expect_error(HyperCom:::get_seurat_layer(seurat, assay="RNA", layer="spliced"), missing[[version]], fixed=TRUE)
  }
})
