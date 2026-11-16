# Argument and input validation: what the steps and HyperCom accept, what they stop on (and with which message) and what they warn
# about.

test_that("threads must be one number of at least 1 and is capped at the available cores with a warning", {
  expect_error(run_priority(threads=NA), "threads")
  expect_error(run_priority(threads=c(2, 4)), "threads")
  expect_error(run_priority(threads="two"), "threads")
  expect_error(run_priority(threads=0), "at least 1")
  expect_error(run_priority(threads=-3), "at least 1")
  expect_warning(priority <- run_priority(threads=1e10), paste0("available cores: using ", parallelly::availableCores()))
  expect_master(priority, fx$output$priority["1"])
})

test_that("a kept factor column, also listed twice, labels genes as gene without a warning (#16)", {
  md <- fx$input$metadata
  md$celltype <- factor(md$celltype)

  expect_no_warning(clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep=c("celltype", "celltype"),
                                                  outfile=NULL))
  is.gene <- !clean.df$cell %in% rownames(md)
  expect_true(all(as.character(clean.df$celltype[is.gene]) == "gene"))
  expect_true(all(as.character(clean.df$celltype.1[is.gene]) == "gene"))
})

test_that("pairs need a receptor and metadata needs cells", {
  md <- fx$output$metadata

  expect_error(score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors=character(0)), "receptors")
  expect_error(score_hypercom(fx$output$infMat, md, ligand="Tgfb1", receptors=c("Tgfbr1", NA)), "ligand and a receptor")
  md$sample <- factor(rep("gene", nrow(md)))
  expect_error(score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors="Sdc4"), "no cells")
})

test_that("metadata cells must be distinct nodes of the influence matrix", {
  md <- fx$output$metadata
  cell <- which(md$sample != "gene")[1]

  expect_error(score_hypercom(fx$output$infMat, rbind(md, md[cell, ]), ligand="Thbs1", receptors="Sdc4"), "distinct nodes")
  md$cell[cell] <- "NotACell"
  expect_error(score_hypercom(fx$output$infMat, md, ligand="Thbs1", receptors="Sdc4"), "distinct nodes")
})

test_that("an empty LR table keeps the lrdb columns, and counts without an expressed ligand or receptor gene stop clearly (#14)", {
  counts <- fx$input$counts
  counts[rownames(counts) %in% unlist(fx$input$lrdb[, -1]), ] <- 0
  cc <- filter_counts(counts, fx$input$lrdb, outfile=NULL)
  lrs <- generate_lrs_table(cc, fx$input$lrdb, outfile=NULL)
  empty <- generate_lrs_table(fx$output$counts, fx$input$lrdb[0, ], outfile=NULL)
  he <- generate_hyperedges(cc, lrs, outfile=NULL)
  adj <- generate_adjacency_matrix(he, outfile=NULL)

  for(df in list(lrs, empty)){
    expect_equal(nrow(df), 0)
    expect_identical(colnames(df), colnames(fx$input$lrdb))
  }
  expect_equal(nrow(adj), 0)
  expect_error(generate_infMat(adj, outfile=NULL), "no nodes")
})

test_that("kept columns keep their names, and kept columns named cell or sample are renamed and keep their values (#1, #20)", {
  md <- fx$input$metadata
  md$cell <- paste0("x", seq_len(nrow(md)))
  md$orig <- md$sample
  one <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep="celltype", outfile=NULL)
  clean.df <- generate_metadata(fx$output$infMat, md, sample="orig", group="condition", keep=c("cell", "sample", "celltype"), outfile=NULL)
  is.cell <- clean.df$cell %in% rownames(md)

  expect_identical(colnames(one), c("cell", "sample", "group", "celltype"))
  expect_identical(colnames(clean.df), c("cell", "sample", "group", "cell.1", "sample.1", "celltype"))
  expect_identical(clean.df$cell.1[is.cell], md$cell[match(clean.df$cell[is.cell], rownames(md))])
  expect_identical(clean.df$cell, fx$output$metadata$cell)
  expect_identical(clean.df$sample, fx$output$metadata$sample)
})

test_that("numeric node names are ordered as numbers and hyperedges need node and edge columns", {
  he <- data.frame(edge=c("a", "a", "b", "b", "c", "c"), node=c(10, 2, 2, 100, 100, 10))
  adj <- generate_adjacency_matrix(he, outfile=NULL)

  expect_identical(rownames(adj), c("2", "10", "100"))
  expect_error(generate_adjacency_matrix(data.frame(Edge=he$edge, Node=he$node), outfile=NULL), "node and edge")
})

test_that("generate_infMat takes a dataframe and keeps an unnamed adjacency unnamed", {
  dense <- as.matrix(fx$output$adj)

  expect_master_infmat(generate_infMat(as.data.frame(dense), outfile=NULL))
  expect_null(dimnames(generate_infMat(unname(dense), outfile=NULL)))
})

test_that("an adjacency with one-sided dimnames gives the same priorities", {
  adj <- fx$output$adj
  adj@Dimnames <- list(NULL, rownames(adj))

  expect_master(run_priority(adj=adj), fx$output$priority["1"])
})

test_that("HyperCom checks its arguments before any step runs, and weights that print alike stop (#27)", {
  dir <- withr::local_tempdir()
  args <- list(lrdb=fx$input$lrdb, counts=fx$input$counts, metadata=fx$input$metadata, sample="sample", group="condition",
               outdir=file.path(dir, "run"), threads=test_threads(2), verbose=FALSE)
  run <- function(...){
    return(do.call(HyperCom, utils::modifyList(args, list(...))))
  }

  expect_error(run(lrdb=fx$input$counts, counts=fx$input$lrdb), "lrdb must be")
  expect_error(run(sample="smaple"), "not columns of metadata: smaple")
  expect_error(run(group="conditon", keep="celltyp"), "conditon, celltyp")
  expect_error(run(sample=fx$input$metadata$sample), "one column name")
  expect_error(run(metadata=as.matrix(fx$input$metadata)), "dataframe")
  expect_error(run(weights=-1), "at least 0")
  expect_error(run(outdir=4), "cl argument before it was removed")
  expect_error(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", outdir=NULL), "outdir must be one string")
  expect_error(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata[0, ], sample="sample", outdir=file.path(dir, "run")),
               "no cell of counts is a row name of metadata")
  expect_error(run(weights="1"), "at least 0")
  expect_error(run(weights=c(0.3, 0.1 + 0.2)), "weights must differ")
  expect_error(run(infmat.method="qr"), "should be one of")
  expect_false(dir.exists(file.path(dir, "run")))
  expect_error(run_priority(weights=c(0.3, 0.1 + 0.2)), "weights must differ")
})

test_that("counts of another species, transposed counts and scaled data stop clearly", {
  dir <- withr::local_tempdir()
  counts <- fx$input$counts
  upper <- counts
  rownames(upper) <- toupper(rownames(upper))

  expect_error(filter_counts(Matrix::t(counts), fx$input$lrdb, outfile=NULL), "genes as column names")
  expect_error(filter_counts(counts - 1, fx$input$lrdb, outfile=NULL), "negative values")
  expect_error(run_hypercom(dir, counts=upper), "species")
})

test_that("samples need 2 cells in the graph, and metadata may be a tibble with row names", {
  md <- fx$input$metadata
  md$sample <- as.character(md$sample)
  md[intersect(rownames(md), rownames(fx$output$infMat))[1], "sample"] <- "lonely"

  expect_error(generate_metadata(fx$output$infMat, md, sample="sample", outfile=NULL), "at least 2 cells.*lonely")
  tb <- suppressWarnings({x <- dplyr::as_tibble(fx$input$metadata); rownames(x) <- rownames(fx$input$metadata); x})
  expect_identical(generate_metadata(fx$output$infMat, tb, sample="sample", group="condition", outfile=NULL),
                   generate_metadata(fx$output$infMat, fx$input$metadata, sample="sample", group="condition", outfile=NULL))
})

test_that("score_hypercom needs both or neither of ligand and receptors, and one ligand", {
  md <- fx$output$metadata

  expect_error(score_hypercom(fx$output$infMat, md, receptors="Sdc4"), "both ligand and receptors")
  expect_error(score_hypercom(fx$output$infMat, md, ligand="Thbs1"), "both ligand and receptors")
  expect_error(score_hypercom(fx$output$infMat, md, ligand=c("Thbs1", "Tgfb1"), receptors="Sdc4"), "one gene")
})

test_that("HyperCom warns about priority files of other weights already in outdir", {
  dir <- withr::local_tempdir()

  expect_no_warning(run_hypercom(dir, weights=c(1, 0.5), threads=test_threads(2)))
  expect_no_warning(run_hypercom(dir, weights=c(1, 0.5), threads=test_threads(2)))
  expect_warning(run_hypercom(dir, weights=1, threads=test_threads(2)), "priority.0.5.csv")
})

test_that("within-sample cell pairs may pass 2^32 and stop at 2^52, before any block is made", {
  skip_on_cran()
  past32 <- tryCatch(HyperCom:::influence_blocks(NULL, list(sample=rep(1L, 65537))), error=conditionMessage)
  big <- list(sample=rep(1L, 2^26 + 1))

  expect_false(grepl("sample sizes", past32))
  expect_error(HyperCom:::influence_blocks(NULL, big), "below 2^52", fixed=TRUE)
})
