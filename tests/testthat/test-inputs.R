# Input checks of the arguments and paths master does not have (threads, name mapping, batches) and of inputs master rejects too,
# with clearer messages.

test_that("threads must be one number and is capped at the available cores", {
  skip_if_not("threads" %in% names(formals(prioritize_lr)), "build without a threads argument")
  args <- list(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, outdir=NULL, verbose=FALSE)

  expect_error(do.call(prioritize_lr, c(args, threads=NA)), "threads")
  expect_error(do.call(prioritize_lr, c(args, list(threads=c(2, 4)))), "threads")
  expect_error(do.call(prioritize_lr, c(args, threads="two")), "threads")
  expect_identical(do.call(prioritize_lr, c(args, threads=1e10)), fx$output$priority["1"])
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

  expect_identical(generate_infMat(as.data.frame(dense), outfile=NULL), fx$output$infMat)
  expect_null(dimnames(generate_infMat(unname(dense), outfile=NULL)))
})

test_that("an adjacency with one-sided dimnames gives the same priorities", {
  adj <- fx$output$adj
  adj@Dimnames <- list(NULL, rownames(adj))
  priority <- quietly(prioritize_lr(fx$output$lrs, adj, fx$output$infMat, fx$output$metadata, outdir=NULL))

  expect_identical(priority, fx$output$priority["1"])
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
  expect_identical(as.matrix(infMat, threads=3), F)
  expect_output(print(adj), "227 nodes")
  expect_output(print(infMat), "227 nodes")
})

test_that("parallel uses one thread per physical core and an explicit count up to the logical cores", {
  skip_if_not("threads" %in% names(formals(prioritize_lr)), "build without a threads argument")
  physical <- min(parallelly::availableCores(), parallelly::availableCores(logical=FALSE), na.rm=TRUE)

  expect_equal(HyperCom:::get_threads(parallel=TRUE), physical)
  expect_equal(HyperCom:::get_threads(parallel=TRUE, cl=2), min(2, physical))
  expect_equal(HyperCom:::get_threads(threads=1e10), as.integer(parallelly::availableCores()))
  expect_equal(HyperCom:::get_threads(), 1)
})
