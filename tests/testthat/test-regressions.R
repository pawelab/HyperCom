# One test per bug of 0.0.0.9000 (b999a01) fixed in 0.0.1. "#n" are the findings of the review of b999a01, "r8 #n" the later review's
# findings that b999a01 shares. #10 (seeded p-values) is in test-steps.R; #13 has no behaviour to assert.

test_that("#1 one kept column keeps its name", {
  md <- generate_metadata(fx$output$infMat, fx$input$metadata, sample="sample", group="condition", keep="celltype", outfile=NULL)

  expect_identical(colnames(md), c("cell", "sample", "group", "celltype"))
})

test_that("#2 a cell with a missing sample stops instead of becoming a gene", {
  md <- fx$input$metadata
  md$sample[1] <- NA

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", outfile=NULL), "1 cells miss a sample")
})

test_that("#3 a group named gene stops with a clear message", {
  md <- fx$input$metadata
  md$condition[md$condition == "TAC"] <- "gene"

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", outfile=NULL), "cannot name a sample or group")
})

test_that("#4 the documented lrdb format with an edges column works", {
  lrdb <- fx$input$lrdb
  colnames(lrdb)[1] <- "edges"
  lrs <- generate_lrs_table(fx$output$counts, lrdb, outfile=NULL)

  expect_equal(plain_df(generate_hyperedges(fx$output$counts, lrs, outfile=NULL)), plain_df(fx$output$hyperedges))
})

test_that("#5 scoring a gene missing from the graph fails with a clear message", {
  expect_error(score_hypercom(fx$output$infMat, fx$output$metadata, ligand="NotAGene", receptors="Sdc4"), "NotAGene")
})

test_that("#6 HyperCom leaves the working directory alone and #9 makes a nested outdir", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  run_hypercom(file.path("a", "b", "run"))

  expect_identical(normalizePath(getwd()), normalizePath(dir))
  expect_true(file.exists(file.path(dir, "a", "b", "run", "6.priority", "priority.1.csv")))
})

test_that("#7 a cluster passed as cl is left running", {
  skip_on_cran()
  cl <- parallelly::makeClusterPSOCK(2, autoStop=TRUE, rscript_envs=c(R_LIBS=paste(.libPaths(), collapse=.Platform$path.sep)))
  on.exit(parallel::stopCluster(cl))
  priority <- run_priority(parallel=TRUE, cl=cl)

  expect_identical(priority, fx$output$priority["1"])
  expect_identical(unlist(parallel::clusterEvalQ(cl, 1)), c(1, 1))
})

test_that("#8 the qrsh cluster is made only inside an SGE parallel job", {
  withr::local_envvar(PE_HOSTFILE="")
  local.args <- cluster_args()
  withr::local_envvar(PE_HOSTFILE="pe_hostfile")
  sge.args <- suppressWarnings(cluster_args())

  expect_null(local.args$rshcmd)
  expect_identical(local.args[[1]], parallelly::availableCores())
  expect_identical(sge.args$rshcmd, "qrsh")
})

test_that("#11 the bundled counts can be subset right after library(HyperCom)", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, c("-e", shQuote("library(HyperCom); cat(ncol(alexanian2021.counts[, 1:10]))")),
                                  stdout=TRUE, stderr=TRUE))

  expect_identical(utils::tail(out, 1), "10")
})

test_that("#12 prioritize_lr prints nothing per pair and is silent with verbose=FALSE", {
  expect_silent(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=NULL, verbose=FALSE))
  msgs <- testthat::capture_messages(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=NULL))

  expect_length(msgs, 2)
  expect_match(msgs[2], "top pair")
})

test_that("#14 an lrdb with no interaction present gives an empty table", {
  lrs <- generate_lrs_table(fx$output$counts, fx$input$lrdb[fx$input$lrdb$ligand == "NotAGene", ], outfile=NULL)

  expect_identical(nrow(lrs), 0L)
})

test_that("#15 node order does not depend on the collation locale", {
  english <- c("en_US.UTF-8", "English_United States.utf8", "English_United States.1252", "en_US")
  usable <- vapply(english, function(loc) suppressWarnings(withr::with_collate(loc, Sys.getlocale("LC_COLLATE"))) != "C", logical(1))
  skip_if_not(any(usable), "no English collation locale on this system")
  adj.en <- withr::with_collate(english[usable][1], generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))

  expect_identical(adj.en, fx$output$adj)
})

test_that("#16 a factor kept column labels gene rows as gene", {
  md <- fx$input$metadata
  md$celltype <- factor(md$celltype)
  clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep="celltype", outfile=NULL)
  is.gene <- clean.df$sample == "gene"

  expect_true(all(clean.df$celltype[is.gene] == "gene"))
  expect_false(anyNA(clean.df$celltype))
})

test_that("#17 ligand-receptor columns stored as factors give the same priorities as characters", {
  lrs <- fx$output$lrs
  lrs[] <- lapply(lrs, factor)
  priority <- run_priority(lrs=lrs)

  expect_identical(priority[["1"]]$priority, fx$output$priority[["1"]]$priority)
  expect_identical(as.character(priority[["1"]]$interaction), fx$output$priority[["1"]]$interaction)
})

test_that("#18 an lrdb read back with stringsAsFactors keeps its interactions", {
  f <- withr::local_tempfile(fileext=".csv")
  utils::write.csv(fx$input$lrdb, f, row.names=FALSE)
  lrdb <- utils::read.csv(f, stringsAsFactors=TRUE)
  counts <- filter_counts(fx$input$counts, lrdb, outfile=NULL)
  lrs <- generate_lrs_table(counts, lrdb, outfile=NULL)

  expect_identical(counts, fx$output$counts)
  expect_equal(nrow(lrs), nrow(fx$output$lrs))
  expect_identical(generate_adjacency_matrix(generate_hyperedges(counts, lrs, outfile=NULL), outfile=NULL), fx$output$adj)
})

test_that("#19 two cells with one name stop instead of merging into one node", {
  counts <- fx$output$counts
  rownames(counts)[2] <- rownames(counts)[1]

  expect_error(generate_hyperedges(counts, fx$output$lrs, outfile=NULL), "cell names")
})

test_that("#20 kept columns named like the output's columns are renamed", {
  md <- fx$input$metadata
  md$cell <- "x"
  clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep=c("cell", "sample"), outfile=NULL)

  expect_identical(colnames(clean.df), c("cell", "sample", "group", "cell.1", "sample.1"))
})

test_that("#21 metadata whose row names are not the cell names stops", {
  md <- fx$input$metadata
  rownames(md) <- NULL

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", outfile=NULL), "row names")
})

test_that("#22 factor hyperedge columns give the node order of strings", {
  he <- as.data.frame(fx$output$hyperedges)
  he$node <- factor(he$node, levels=rev(unique(he$node)))
  he$edge <- factor(he$edge)

  expect_identical(generate_adjacency_matrix(he, outfile=NULL), fx$output$adj)
})

test_that("#23 metadata rows in another order give the same scores and priorities", {
  md <- fx$output$metadata[rev(seq_len(nrow(fx$output$metadata))), ]

  expect_same_scores(score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors="Sdc4"), fx$output$scores$single)
  expect_same_scores(score_hypercom(fx$output$infMat, md, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2")), fx$output$scores$complex)
  expect_identical(run_priority(metadata=md)[["1"]]$priority, fx$output$priority[["1"]]$priority)
})

test_that("#24 cells of the counts missing from the metadata are left out, with a warning, and never drawn as background genes", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  missing <- rownames(fx$input$metadata)[1:5]
  expect_warning(res <- run_hypercom(file.path(dir, "run"), metadata=fx$input$metadata[-(1:5), ], significance=TRUE, seed=1),
                 "5 cells of counts are not in metadata and are left out")

  expect_false(any(missing %in% rownames(res$adj)))
  expect_false(any(missing %in% res$metadata$cell))
  expect_false(anyNA(res$priority[["1"]]$p))
})

test_that("#25 a missing count does not stop filter_counts", {
  counts <- fx$input$counts
  counts[which(rownames(counts) == "Thbs1"), 1] <- NA

  expect_true("Thbs1" %in% colnames(filter_counts(counts, fx$input$lrdb, outfile=NULL)))
})

test_that("#26 a date column can be kept, as strings", {
  md <- fx$input$metadata
  md$day <- as.Date("2026-01-01") + seq_len(nrow(md))
  clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep=c("day", "celltype"), outfile=NULL)
  is.gene <- !clean.df$cell %in% rownames(md)

  expect_true(all(clean.df$day[is.gene] == "gene"))
  expect_identical(clean.df$day[!is.gene], as.character(md[clean.df$cell[!is.gene], "day"]))
})

test_that("#27 r8 #17 weights that print alike, are negative or are not numbers stop", {
  expect_error(run_priority(weights=c(0.3, 0.1 + 0.2)), "weights must differ")
  expect_error(run_priority(weights=-4), "at least 0")
  expect_error(run_priority(weights="1"), "at least 0")
})

test_that("#28 gene names used by several rows of the counts give a warning", {
  counts <- fx$input$counts
  counts <- rbind(counts, counts[rownames(counts) == "Tgfbr1", , drop=FALSE] * 2)

  expect_warning(filter_counts(counts, fx$input$lrdb, outfile=NULL), "Tgfbr1")
})

test_that("r8 #1 an lrdb with no interaction expressed stops HyperCom", {
  dir <- withr::local_tempdir()
  lrdb <- fx$input$lrdb
  lrdb[, -1] <- lapply(lrdb[, -1], toupper)

  expect_error(run_hypercom(file.path(dir, "run"), lrdb=lrdb), "species")
})

test_that("r8 #2 a sample with one cell stops, and background pairs without a priority are left out of the p-values", {
  md <- fx$input$metadata
  md$sample[1] <- "lonely"

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", outfile=NULL), "at least 2 cells")
  expect_identical(permutation_test(0.5, c(0.1, NaN, 0.7, NA)), 0.5)
  expect_identical(permutation_test(0.5, c(NaN, NA)), NA_real_)
})

test_that("r8 #3 cells balanced over the groups give purity 0, not a NaN priority", {
  md <- fx$output$metadata
  cells <- md$sample != "gene"
  md$group <- factor(as.character(md$group), levels=c("A", "B", "gene"))
  md$group[cells] <- rep(c("A", "B"), length.out=sum(cells))
  priority <- run_priority(metadata=md)

  expect_false(anyNA(priority[["1"]]$priority))
})

test_that("r8 #4 score_hypercom with only a ligand or only receptors stops, and NA receptors are left out", {
  expect_identical(score_hypercom(fx$output$infMat, fx$output$metadata, ligand="Tgfb1", receptors=c("Tgfbr1", NA, "Tgfbr2")),
                   fx$output$scores$complex)
  expect_error(score_hypercom(fx$output$infMat, fx$output$metadata, ligand="Thbs1", receptors=NA), "at least one gene")
  expect_error(score_hypercom(fx$output$infMat, fx$output$metadata, ligand="Thbs1"), "both ligand and receptors")
  expect_error(score_hypercom(fx$output$infMat, fx$output$metadata, receptors="Sdc4"), "both ligand and receptors")
})

test_that("r8 #5 r8 #10 r8 #11 scaled, transposed or swapped inputs stop", {
  dir <- withr::local_tempdir()
  counts <- fx$input$counts
  counts[which(rownames(counts) == "Thbs1"), 1] <- -1

  expect_error(filter_counts(counts, fx$input$lrdb, outfile=NULL), "negative")
  expect_no_warning(filter_counts(fx$input$counts * NA, fx$input$lrdb, outfile=NULL))
  expect_error(filter_counts(Matrix::t(fx$input$counts), fx$input$lrdb, outfile=NULL), "genes as column names")
  expect_error(filter_counts(fx$input$lrdb, fx$input$counts, outfile=NULL), "lrdb must be")
  expect_error(run_hypercom(file.path(dir, "run"), counts=Matrix::t(fx$input$counts)), "cells as row names")
  expect_error(generate_metadata(fx$output$infMat, fx$input$metadata, sample="Sample", outfile=NULL), "not columns of metadata: Sample")
})

test_that("r8 #9 a SeuratObject v3 assay is read from the layer asked for", {
  skip_if_not_installed("SeuratObject")
  assay <- SeuratObject::CreateAssayObject(counts=fx$input$counts)
  seurat <- SeuratObject::CreateSeuratObject(assay)

  expect_identical(get_seurat_layer(seurat, layer="counts"), methods::slot(assay, "counts"))
})

test_that("r8 #13 r8 #14 HyperCom returns its results and warns about priority files of other weights", {
  dir <- withr::local_tempdir()
  res <- run_hypercom(file.path(dir, "run"), weights=c(1, 0.5))

  expect_named(res, c("priority", "metadata", "infMat", "adj", "lrs"))
  expect_warning(run_hypercom(file.path(dir, "run")), "priority.0.5.csv")
})
