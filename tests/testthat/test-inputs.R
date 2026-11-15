# Input checks of the arguments and paths master does not have (threads, name mapping, batches) and of inputs master rejects too,
# with clearer messages.

test_that("threads must be one number of at least 1 and is capped at the available cores with a warning", {
  skip_if_not("threads" %in% names(formals(prioritize_lr)), "build without a threads argument")
  args <- list(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=NULL, verbose=FALSE)

  expect_error(do.call(prioritize_lr, c(args, threads=NA)), "threads")
  expect_error(do.call(prioritize_lr, c(args, list(threads=c(2, 4)))), "threads")
  expect_error(do.call(prioritize_lr, c(args, threads="two")), "threads")
  expect_error(do.call(prioritize_lr, c(args, threads=0)), "at least 1")
  expect_error(do.call(prioritize_lr, c(args, threads=-3)), "at least 1")
  expect_warning(priority <- do.call(prioritize_lr, c(args, threads=1e10)), "available cores")
  expect_master(priority, fx$output$priority["1"])
})

test_that("prioritize_lr reports progress as messages, not printed lines, and verbose=FALSE is silent", {
  dir <- withr::local_tempdir()
  run <- function(verbose){
    return(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5), threads=2,
                         outdir=dir, verbose=verbose))
  }

  expect_message(run(TRUE), "pairs scored")
  expect_message(run(TRUE), "Weight 0.5: top pair")
  expect_length(utils::capture.output(x <- suppressMessages(run(TRUE))), 0)
  expect_silent(run(FALSE))
  msgs <- character(0)
  withCallingHandlers(run(TRUE), message=function(m){
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_equal(sum(grepl("pairs scored", msgs)), 1)
  expect_identical(HyperCom:::format_duration(c(30)), "30 s")
  expect_identical(HyperCom:::format_duration(300), "5 min")
  expect_identical(HyperCom:::format_duration(9000), "2.5 h")
})

test_that("grids bigger than the kernel's cache budget are scored in about 20 batches or more, smaller ones in big batches", {
  batch_size <- HyperCom:::batch_size

  expect_equal(batch_size(16, 1000, c(500, 500), 4090), 512)
  expect_equal(batch_size(1, 1000, c(500, 500), 4090), 256)
  expect_equal(batch_size(16, 20000, c(10000, 6667, 3333), 400), 20)
  expect_equal(batch_size(16, 20000, c(10000, 6667, 3333), 40), 16)
  expect_equal(batch_size(1, 20000, c(10000, 6667, 3333), 4090), 205)
})

test_that("score_hypercom gives the same scores on several threads, from a dense or factored influence matrix", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=2)
  md <- fx$output$metadata

  for(p in list(list(NULL, NULL), list("Thbs1", "Sdc4"), list("Tgfb1", c("Tgfbr1", "Tgfbr2")))){
    expect_identical(score_hypercom(fx$output$infMat, md, p[[1]], p[[2]], threads=test_threads(3)),
                     score_hypercom(fx$output$infMat, md, p[[1]], p[[2]]))
    expect_identical(score_hypercom(infMat, md, p[[1]], p[[2]], parallel=TRUE), score_hypercom(infMat, md, p[[1]], p[[2]]))
  }
  expect_error(score_hypercom(infMat, md, "Thbs1", "Sdc4", threads=0), "at least 1")
})

test_that("group purity is at least 0, so perfectly balanced groups give priority 0, not NaN", {
  balanced <- factor(paste0("t", rep(1:5, 4)))

  expect_lt(1 - 1 / 5 - HyperCom:::gini_impurity(balanced), 0)
  expect_identical(HyperCom:::group_purity(balanced, 1 - 1 / 5), 0)
  expect_identical(HyperCom:::group_purity(balanced[1:4], 0), 1)
  expect_identical(HyperCom:::group_purity(factor(c("a", "a", "a", "b")), 0.5), 0.5 - HyperCom:::gini_impurity(factor(c("a", "a", "a", "b"))))
})

test_that("a kept factor column listed twice labels genes in both copies", {
  md <- fx$input$metadata
  md$celltype <- factor(md$celltype)
  clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep=c("celltype", "celltype"), outfile=NULL)
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

test_that("counts without an expressed ligand or receptor gene stop clearly", {
  counts <- fx$input$counts
  counts[rownames(counts) %in% unlist(fx$input$lrdb[, -1]), ] <- 0
  cc <- filter_counts(counts, fx$input$lrdb, outfile=NULL)
  lrs <- generate_lrs_table(cc, fx$input$lrdb, outfile=NULL)
  he <- generate_hyperedges(cc, lrs, outfile=NULL)
  adj <- generate_adjacency_matrix(he, outfile=NULL)

  expect_equal(nrow(adj), 0)
  expect_error(generate_infMat(adj, outfile=NULL), "no nodes")
})

test_that("a kept column named cell keeps the metadata's own values", {
  md <- fx$input$metadata
  md$cell <- paste0("x", seq_len(nrow(md)))
  clean.df <- generate_metadata(fx$output$infMat, md, sample="sample", group="condition", keep="cell", outfile=NULL)
  is.cell <- clean.df$cell %in% rownames(md)

  expect_identical(clean.df$cell.1[is.cell], md$cell[match(clean.df$cell[is.cell], rownames(md))])
  expect_identical(clean.df$cell, fx$output$metadata$cell)
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
  priority <- quietly(prioritize_lr(fx$output$lrs, adj, fx$output$infMat, fx$output$metadata, outdir=NULL))

  expect_master(priority, fx$output$priority["1"])
})

test_that("adjacency and influence matrix without node x node matrices index like the matrices they stand for", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  skip_if_not(inherits(adj, "hypercom_adjacency"), "build with node x node matrices")
  infMat <- generate_infMat(adj, outfile=NULL, threads=2)
  A <- as.matrix(adj)
  F <- as.matrix(infMat)
  nodes <- rownames(A)
  index <- list(list(5, nodes[c(3, 9)]), list(c(9, 1, 9), -1), list(nodes[1:4] %in% nodes[c(2, 4)], 7), list(3, 3),
                list(factor(nodes[c(10, 2)]), c(2, 10, 2)), list(nodes[200], nodes))

  expect_s3_class(infMat, "hypercom_influence")
  expect_identical(dim(adj), dim(A))
  expect_identical(dimnames(infMat), dimnames(fx$output$infMat))
  for(k in seq_along(index)){
    i <- index[[k]][[1]]
    j <- index[[k]][[2]]
    if(is.factor(i)){
      i.base <- as.character(i)
    } else {
      i.base <- i
    }
    expect_identical(adj[i, j], A[i.base, j], info=k)
    expect_identical(infMat[i, j], F[i.base, j], info=k)
    expect_identical(infMat[i, j, drop=FALSE], F[i.base, j, drop=FALSE], info=k)
  }
  expect_identical(adj[, 2], A[, 2])
  expect_identical(infMat[2, ], F[2, ])
  expect_error(adj[5], "adj[i, j]", fixed=TRUE)
  expect_error(infMat[5], "infMat[i, j]", fixed=TRUE)
  expect_error(infMat["NotANode", 1], "subscript out of bounds")
  expect_error(adj[1, 1000], "subscript out of bounds")
  expect_error(adj[c(1, NA), 1], "NA")
  expect_error(infMat[c(TRUE, logical(nrow(A))), 1], "too long")
  expect_identical(as.matrix(infMat, threads=test_threads(3)), F)
  expect_output(print(adj), "227 nodes")
  expect_output(print(infMat), "227 nodes")
})

test_that("parallel uses one thread per physical core and an explicit count up to the logical cores", {
  skip_if_not("threads" %in% names(formals(prioritize_lr)), "build without a threads argument")
  physical <- min(parallelly::availableCores(), parallelly::availableCores(logical=FALSE), na.rm=TRUE)

  expect_equal(HyperCom:::get_threads(parallel=TRUE), physical)
  expect_equal(HyperCom:::get_threads(parallel=TRUE, cl=2), min(2, physical))
  expect_warning(n <- HyperCom:::get_threads(threads=1e10), "available cores")
  expect_equal(n, as.integer(parallelly::availableCores()))
  expect_equal(HyperCom:::get_threads(), 1)
})

test_that("HyperCom checks its arguments before any step runs", {
  dir <- withr::local_tempdir()
  args <- list(lrdb=fx$input$lrdb, counts=fx$input$counts, metadata=fx$input$metadata, sample="sample", group="condition",
               outdir=file.path(dir, "run"), threads=2, verbose=FALSE)
  run <- function(...){
    return(do.call(HyperCom, utils::modifyList(args, list(...))))
  }

  expect_error(run(lrdb=fx$input$counts, counts=fx$input$lrdb), "lrdb must be")
  expect_error(run(sample="smaple"), "not columns of metadata: smaple")
  expect_error(run(group="conditon", keep="celltyp"), "conditon, celltyp")
  expect_error(run(sample=fx$input$metadata$sample), "one column name")
  expect_error(run(metadata=as.matrix(fx$input$metadata)), "dataframe")
  expect_error(run(weights=-1), "at least 0")
  expect_error(run(weights="1"), "at least 0")
  expect_error(run(infmat.method="hypergraph"), "now called \"woodbury\"")
  expect_false(dir.exists(file.path(dir, "run")))
})

test_that("counts of another species, transposed counts and scaled data stop clearly", {
  dir <- withr::local_tempdir()
  counts <- fx$input$counts
  upper <- counts
  rownames(upper) <- toupper(rownames(upper))

  expect_error(filter_counts(Matrix::t(counts), fx$input$lrdb, outfile=NULL), "genes as column names")
  expect_error(filter_counts(counts - 1, fx$input$lrdb, outfile=NULL), "negative values")
  expect_error(quietly(HyperCom(fx$input$lrdb, upper, fx$input$metadata, sample="sample", outdir=dir, verbose=FALSE)), "species")
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
  run <- function(weights){
    return(quietly(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition", weights=weights,
                            outdir=dir, threads=2, verbose=FALSE)))
  }

  expect_no_warning(run(c(1, 0.5)))
  expect_no_warning(run(c(1, 0.5)))
  expect_warning(run(1), "priority.0.5.csv")
})
