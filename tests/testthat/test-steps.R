# Every step on the fixture's input gives the output of 0.0.0.9000. Each step is fed 0.0.0.9000's own input from the fixture, so a
# difference points at the step that made it.

test_that("filter_counts and generate_lrs_table match 0.0.0.9000", {
  counts <- filter_counts(fx$input$counts, fx$input$lrdb, outfile=NULL)
  lrs <- generate_lrs_table(fx$output$counts, fx$input$lrdb, outfile=NULL)

  expect_identical(counts, fx$output$counts)
  expect_identical(lrs, fx$output$lrs)
})

test_that("generate_hyperedges and generate_adjacency_matrix match 0.0.0.9000", {
  hyperedges <- generate_hyperedges(fx$output$counts, fx$output$lrs, outfile=NULL)
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)

  expect_identical(hyperedges, fx$output$hyperedges)
  expect_identical(adj, fx$output$adj)
})

test_that("generate_infMat matches 0.0.0.9000", {
  expect_equal(generate_infMat(fx$output$adj, outfile=NULL), fx$output$infMat, tolerance=1e-12)
  expect_equal(generate_infMat(fx$output$adj, r=0.3, outfile=NULL), fx$output$infMat.r03, tolerance=1e-12)
})

test_that("generate_metadata matches 0.0.0.9000 with a group, without and with kept columns", {
  md <- fx$input$metadata
  infMat <- fx$output$infMat

  expect_identical(generate_metadata(infMat, md, sample="sample", group="condition", outfile=NULL), fx$output$metadata)
  expect_identical(generate_metadata(infMat, md, sample="sample", outfile=NULL), fx$output$metadata.none)
  expect_identical(generate_metadata(infMat, md, sample="sample", group="condition", keep=c("celltype", "UMAP_1"), outfile=NULL),
                   fx$output$metadata.keep)
})

test_that("score_hypercom matches 0.0.0.9000 for the network, a pair, complexes and a same-gene pair", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata

  expect_identical(score_hypercom(infMat, md), fx$output$scores$network)
  expect_identical(score_hypercom(infMat, md, ligand="Thbs1", receptors="Sdc4"), fx$output$scores$single)
  expect_identical(score_hypercom(infMat, md, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2")), fx$output$scores$complex)
  expect_identical(score_hypercom(infMat, md, ligand="Tgfb1", receptors=c("Acvr1", "Tgfbr1", "Tgfbr2")), fx$output$scores$triple)
  expect_identical(score_hypercom(infMat, md, ligand=fx$output$same.pair[1], receptors=fx$output$same.pair[2]), fx$output$scores$same)
})

test_that("prioritize_lr matches 0.0.0.9000 for three weights and without groups", {
  expect_identical(run_priority(weights=c(1, 0.5, 0)), fx$output$priority)
  expect_identical(run_priority(metadata=fx$output$metadata.none), fx$output$priority.none)
})

test_that("#10 seeded p-values match 0.0.0.9000 with set.seed and with seed, and seed restores the random state", {
  skip_on_cran()
  set.seed(1)
  by.set.seed <- run_priority(significance=TRUE)
  set.seed(42)
  before <- .Random.seed
  by.seed <- run_priority(significance=TRUE, seed=1)

  expect_identical(by.set.seed, fx$output$priority.sig)
  expect_identical(by.seed, fx$output$priority.sig)
  expect_identical(.Random.seed, before)
})

test_that("HyperCom matches the steps and writes what it returns", {
  dir <- withr::local_tempdir()
  res <- run_hypercom(file.path(dir, "run"), weights=c(1, 0.5, 0))

  expect_identical(res$adj, fx$output$adj)
  expect_identical(res$metadata, fx$output$metadata)
  expect_identical(res$priority, fx$output$priority)
  expect_identical(readRDS(file.path(dir, "run", "3.adj.rds")), res$adj)
  expect_identical(readRDS(file.path(dir, "run", "4.infMat.rds")), res$infMat)
  expect_equal(utils::read.csv(file.path(dir, "run", "6.priority", "priority.0.5.csv")), plain_df(res$priority[["0.5"]]))
})

test_that("HyperCom_sce matches HyperCom on a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  dir <- withr::local_tempdir()
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=fx$input$counts), colData=fx$input$metadata)
  res <- quietly(HyperCom_sce(sce, fx$input$lrdb, sample="sample", group="condition", outdir=file.path(dir, "run"),
                              weights=c(1, 0.5, 0)))

  expect_identical(res$adj, fx$output$adj)
  expect_identical(res$metadata, fx$output$metadata)
  expect_identical(res$priority, fx$output$priority)
})

test_that("get_sce_assay reads delayed assays for the lrdb genes and stops on bad input", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("DelayedArray")
  counts <- fx$input$counts
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=DelayedArray::DelayedArray(as.matrix(counts))))
  delayed <- get_sce_assay(sce, fx$input$lrdb, assay="counts")
  genes <- unique(unlist(lapply(fx$input$lrdb[, -1, drop=FALSE], as.character), use.names=FALSE))

  expect_s4_class(delayed, "dgCMatrix")
  expect_identical(colnames(delayed), colnames(counts))
  expect_true(all(rownames(delayed) %in% genes))
  expect_equal(as.matrix(filter_counts(delayed, fx$input$lrdb, outfile=NULL)),
               as.matrix(filter_counts(counts, fx$input$lrdb, outfile=NULL)))
  expect_error(get_sce_assay(sce, fx$input$lrdb, assay="logcounts"), "no assay 'logcounts'")
  expect_error(get_sce_assay(unname(sce), fx$input$lrdb), "no cell names")
  expect_error(get_sce_assay(counts, fx$input$lrdb), "must be a SingleCellExperiment")
  ensembl <- sce
  rownames(ensembl) <- paste0("ENSMUSG", seq_len(nrow(ensembl)))
  expect_error(get_sce_assay(ensembl, fx$input$lrdb), "must be gene symbols")
})
