# Intended behaviour where master b999a01 is known to differ. Skipped unless HYPERCOM_KNOWN_BUGS=true; a failure here means the bug
# is still present. Numbers follow the project findings list.

test_that("#1 generate_metadata keeps the name of a single keep column", {
  skip_known_bug("#1")
  md <- generate_metadata(fx$output$infMat, fx$input$metadata, sample="sample", group="condition", keep="celltype", outfile=NULL)

  expect_true("celltype" %in% colnames(md))
})

test_that("#2 a cell with a missing sample is not labelled as a gene", {
  skip_known_bug("#2")
  md.in <- fx$input$metadata
  md.in$sample[1] <- NA
  res <- tryCatch(generate_metadata(fx$output$infMat, md.in, sample="sample", group="condition", outfile=NULL), error=function(e) e)

  if(inherits(res, "error")){
    expect_match(conditionMessage(res), "sample")
  } else {
    expect_false(identical(as.character(res$sample[res$cell == rownames(md.in)[1]]), "gene"))
  }
})

test_that("#3 a group named gene is scored as cells or refused clearly", {
  skip_known_bug("#3")
  md.in <- fx$input$metadata
  md.in$condition[md.in$condition == "TAC"] <- "gene"
  res <- tryCatch(generate_metadata(fx$output$infMat, md.in, sample="sample", group="condition", outfile=NULL), error=function(e) e)

  if(inherits(res, "error")){
    expect_match(conditionMessage(res), "gene")
  } else {
    expect_equal(nrow(score_hypercom(fx$output$infMat, res, ligand="Thbs1", receptors="Sdc4")), nrow(md.in))
  }
})

test_that("#4 the documented lrdb format with an edges column works", {
  skip_known_bug("#4")
  lrdb <- fx$input$lrdb
  colnames(lrdb)[1] <- "edges"
  lrs <- generate_lrs_table(fx$output$counts, lrdb, outfile=NULL)
  res <- tryCatch(generate_hyperedges(fx$output$counts, lrs, outfile=NULL), error=function(e) e)

  expect_false(inherits(res, "error"))
  if(!inherits(res, "error")){
    expect_equal(nrow(res), nrow(fx$output$hyperedges))
  }
})

test_that("#5 scoring a gene missing from the graph fails with a clear message", {
  skip_known_bug("#5")

  expect_error(score_hypercom(fx$output$infMat, fx$output$metadata, ligand="NotAGene", receptors="Sdc4"), "NotAGene")
})

test_that("#6 HyperCom restores the working directory", {
  skip_known_bug("#6")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  status <- quietly(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition", outdir="run"))

  expect_identical(normalizePath(getwd()), normalizePath(dir))
})

test_that("#7 prioritize_lr leaves a cluster it was given running", {
  skip_known_bug("#7")
  skip_on_cran()
  cl <- parallel::makeCluster(2)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
  priority <- quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, parallel=TRUE, cl=cl, outdir=NULL))
  alive <- tryCatch(unlist(parallel::clusterEvalQ(cl, 1)), error=function(e) NULL)

  expect_identical(alive, c(1, 1))
})

test_that("#9 prioritize_lr creates a nested output directory", {
  skip_known_bug("#9")
  outdir <- file.path(withr::local_tempdir(), "a", "b")
  res <- tryCatch(quietly(suppressWarnings(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=outdir))),
                  error=function(e) e)

  expect_false(inherits(res, "error"))
  expect_true(file.exists(file.path(outdir, "priority.1.csv")))
})

test_that("#15 node order does not depend on the collation locale", {
  skip_known_bug("#15")
  english <- c("en_US.UTF-8", "English_United States.utf8", "English_United States.1252", "en_US")
  usable <- vapply(english, function(loc) suppressWarnings(withr::with_collate(loc, Sys.getlocale("LC_COLLATE"))) != "C", logical(1))
  skip_if_not(any(usable), "no English collation locale on this system")
  adj.en <- withr::with_collate(english[usable][1], generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))
  adj.c <- withr::with_collate("C", generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))

  expect_identical(rownames(adj.en), rownames(adj.c))
})

test_that("#11 the bundled counts can be subset right after library(HyperCom)", {
  skip_known_bug("#11")
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, c("-e", shQuote("library(HyperCom); cat(ncol(alexanian2021.counts[, 1:10]))")),
                                  stdout=TRUE, stderr=TRUE))

  expect_identical(utils::tail(out, 1), "10")
})

test_that("#8 prioritize_lr in parallel without a cluster runs on this machine", {
  skip_known_bug("#8")
  skip_on_cran()
  res <- tryCatch(quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, parallel=TRUE, outdir=NULL)),
                  error=function(e) e)

  expect_false(inherits(res, "error"))
})

test_that("#10 prioritize_lr takes a seed for reproducible p-values", {
  skip_known_bug("#10")
  skip_on_cran()
  run <- function(){
    return(quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, significance=TRUE, seed=1,
                                 outdir=NULL)))
  }

  expect_identical(run(), run())
})

test_that("#12 prioritize_lr can run without printing", {
  skip_known_bug("#12")
  out <- utils::capture.output(priority <- prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata,
                                                          verbose=FALSE, outdir=NULL))

  expect_length(out, 0)
})

test_that("#14 generate_lrs_table returns an empty table when no interaction is present", {
  skip_known_bug("#14")
  lrs <- generate_lrs_table(fx$output$counts, fx$input$lrdb[0, ], outfile=NULL)

  expect_equal(nrow(lrs), 0)
  expect_identical(colnames(lrs), colnames(fx$input$lrdb))
})

test_that("#16 a factor column kept in the metadata labels genes as gene", {
  skip_known_bug("#16")
  md <- fx$input$metadata
  md$celltype <- factor(md$celltype)
  clean.df <- suppressWarnings(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep="celltype", outfile=NULL))
  is.gene <- !clean.df$cell %in% rownames(md)

  expect_true(all(as.character(clean.df[is.gene, ncol(clean.df)]) == "gene"))
})

test_that("#17 ligand-receptor columns stored as factors give the same priorities as characters", {
  skip_known_bug("#17")
  lrs <- fx$output$lrs
  lrs[] <- lapply(lrs, factor)
  priority <- quietly(prioritize_lr(lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=NULL))

  expect_equal(priority[["1"]]$priority, fx$output$priority[["1"]]$priority)
  expect_identical(as.character(priority[["1"]]$interaction), fx$output$priority[["1"]]$interaction)
})

test_that("#18 an lrdb read back with stringsAsFactors keeps its interactions", {
  skip_known_bug("#18")
  f <- withr::local_tempfile(fileext=".csv")
  utils::write.csv(fx$input$lrdb, f, row.names=FALSE)
  lrdb <- utils::read.csv(f, stringsAsFactors=TRUE)
  counts <- filter_counts(fx$input$counts, lrdb, outfile=NULL)
  lrs <- generate_lrs_table(counts, lrdb, outfile=NULL)

  expect_identical(counts, fx$output$counts)
  expect_equal(nrow(lrs), nrow(fx$output$lrs))
  expect_identical(rownames(generate_adjacency_matrix(generate_hyperedges(counts, lrs, outfile=NULL), outfile=NULL)), rownames(fx$output$adj))
})

test_that("#19 two cells with one name stop instead of merging into one node", {
  skip_known_bug("#19")
  counts <- fx$output$counts
  rownames(counts)[2] <- rownames(counts)[1]

  expect_error(generate_hyperedges(counts, fx$output$lrs, outfile=NULL), "cell names")
})

test_that("#20 a kept column named sample is kept under a new name", {
  skip_known_bug("#20")
  md <- fx$input$metadata
  md$orig <- md$sample
  clean.df <- generate_metadata(fx$output$infMat, md, sample="orig", group="condition", keep=c("sample", "celltype"), outfile=NULL)

  expect_true(all(c("sample.1", "celltype") %in% colnames(clean.df)))
  expect_identical(clean.df$sample, fx$output$metadata$sample)
})

test_that("#21 metadata whose row names are not the cell names stops", {
  skip_known_bug("#21")
  md <- fx$input$metadata
  rownames(md) <- NULL

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", outfile=NULL), "row names")
})

test_that("#22 factor hyperedge columns give the node order of strings", {
  skip_known_bug("#22")
  he <- as.data.frame(fx$output$hyperedges)
  he$node <- factor(he$node, levels=rev(unique(he$node)))
  he$edge <- factor(he$edge)
  adj <- generate_adjacency_matrix(he, outfile=NULL)

  expect_identical(rownames(adj), rownames(fx$output$adj))
})

test_that("#23 metadata rows in another order give the same scores", {
  skip_known_bug("#23")
  md <- fx$output$metadata[rev(seq_len(nrow(fx$output$metadata))), ]
  score <- score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors="Sdc4")

  expect_same_scores(score, fx$output$scores$single)
})

test_that("#24 HyperCom warns about cells of the counts missing from the metadata", {
  skip_known_bug("#24")
  dir <- withr::local_tempdir()
  md <- fx$input$metadata[-(1:5), ]

  expect_warning(quietly(HyperCom(fx$input$lrdb, fx$input$counts, md, sample="sample", group="condition", outdir=file.path(dir, "run"))),
                 "not in metadata")
})

test_that("#25 a missing count does not stop filter_counts", {
  skip_known_bug("#25")
  counts <- fx$input$counts
  counts[which(rownames(counts) == "Thbs1"), 1] <- NA
  res <- tryCatch(filter_counts(counts, fx$input$lrdb, outfile=NULL), error=function(e) e)

  expect_false(inherits(res, "error"))
  if(!inherits(res, "error")){
    expect_true("Thbs1" %in% colnames(res))
  }
})

test_that("#26 a date column can be kept", {
  skip_known_bug("#26")
  md <- fx$input$metadata
  md$day <- as.Date("2026-01-01") + seq_len(nrow(md))
  clean.df <- tryCatch(generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep=c("day", "celltype"), outfile=NULL),
                       error=function(e) e)

  expect_false(inherits(clean.df, "error"))
  if(!inherits(clean.df, "error")){
    expect_true(all(clean.df$day[!clean.df$cell %in% rownames(md)] == "gene"))
  }
})

test_that("#27 weights that would share a name stop", {
  skip_known_bug("#27")
  dir <- withr::local_tempdir()

  expect_error(quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(0.3, 0.1 + 0.2),
                                     outdir=file.path(dir, "p"))), "weights")
})

test_that("#28 gene names used by several rows of the counts give a warning", {
  skip_known_bug("#28")
  counts <- fx$input$counts
  counts <- rbind(counts, counts[rownames(counts) == "Tgfbr1", , drop=FALSE] * 2)

  expect_warning(filter_counts(counts, fx$input$lrdb, outfile=NULL), "Tgfbr1")
})
