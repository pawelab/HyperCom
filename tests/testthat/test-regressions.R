# Each test is a bug of handwritten b999a01, numbered as in the project's findings list (wiki/findings.md), fixed on this branch. The fixes
# of #1, #14, #16, #20 and #27 are asserted in test-inputs.R, of #9 and #12 in test-prioritize_lr.R, of #10 (seeded p-values) in
# test-exact.R; #13 has no behaviour to assert; #7 and #8 (the cluster) went with the cl argument.

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

test_that("#6 HyperCom restores the working directory", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  run_hypercom("run")

  expect_identical(normalizePath(getwd()), normalizePath(dir))
})

test_that("#11 the bundled counts can be subset right after library(HyperCom)", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, c("-e", shQuote("library(HyperCom); cat(ncol(alexanian2021.counts[, 1:10]))")),
                                  stdout=TRUE, stderr=TRUE))

  expect_identical(utils::tail(out, 1), "10")
})

test_that("#15 node order does not depend on the collation locale", {
  english <- c("en_US.UTF-8", "English_United States.utf8", "English_United States.1252", "en_US")
  usable <- vapply(english, function(loc) suppressWarnings(withr::with_collate(loc, Sys.getlocale("LC_COLLATE"))) != "C", logical(1))
  skip_if_not(any(usable), "no English collation locale on this system")
  adj.en <- withr::with_collate(english[usable][1], generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))
  adj.c <- withr::with_collate("C", generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))

  expect_identical(rownames(adj.en), rownames(adj.c))
})

test_that("#17 ligand-receptor columns stored as factors give the same priorities as characters", {
  lrs <- fx$output$lrs
  lrs[] <- lapply(lrs, factor)
  priority <- run_priority(lrs=lrs)

  expect_equal(priority[["1"]]$priority, fx$output$priority[["1"]]$priority)
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
  expect_identical(rownames(generate_adjacency_matrix(generate_hyperedges(counts, lrs, outfile=NULL), outfile=NULL)), rownames(fx$output$adj))
})

test_that("#19 two cells with one name stop instead of merging into one node", {
  counts <- fx$output$counts
  rownames(counts)[2] <- rownames(counts)[1]

  expect_error(generate_hyperedges(counts, fx$output$lrs, outfile=NULL), "cell names")
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
  adj <- generate_adjacency_matrix(he, outfile=NULL)

  expect_identical(rownames(adj), rownames(fx$output$adj))
})

test_that("#23 metadata rows in another order give the same scores", {
  md <- fx$output$metadata[rev(seq_len(nrow(fx$output$metadata))), ]
  score <- score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors="Sdc4")

  expect_same_scores(score, fx$output$scores$single)
})

test_that("#24 cells of the counts missing from the metadata are left out, with a warning, and never drawn as background genes", {
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

test_that("#28 gene names used by several rows of the counts give a warning", {
  counts <- fx$input$counts
  counts <- rbind(counts, counts[rownames(counts) == "Tgfbr1", , drop=FALSE] * 2)

  expect_warning(filter_counts(counts, fx$input$lrdb, outfile=NULL), "Tgfbr1")
})
