test_that("generate_hyperedges matches handwritten", {
  hyperedges <- generate_hyperedges(fx$output$counts, fx$output$lrs, outfile=NULL)

  expect_s3_class(hyperedges, "data.frame")
  expect_equal(plain_df(hyperedges), plain_df(fx$output$hyperedges))
})

test_that("hyperedges hold interaction, gene and expression edges", {
  hyperedges <- generate_hyperedges(fx$output$counts, fx$output$lrs, outfile=NULL)
  counts <- fx$output$counts
  cells <- rownames(counts)
  genes <- colnames(counts)

  expect_identical(colnames(hyperedges), c("edge", "node", "weight"))
  expect_true(all(hyperedges$weight == 1))
  expect_equal(sum(hyperedges$node %in% cells), sum(counts != 0))
  expect_true(all(hyperedges$edge[hyperedges$node %in% cells] %in% genes))
  expect_setequal(hyperedges$node[hyperedges$edge == hyperedges$node], genes)
  expect_true(all(fx$output$lrs$interaction %in% hyperedges$edge))
})

test_that("generate_hyperedges keeps the largest connected part of the hypergraph, with a warning, and later steps its interactions", {
  lrdb <- data.frame(interaction=c("IA", "IB"), ligand=c("A1", "B1"), r1=c("A2", "B2"))
  counts <- Matrix::sparseMatrix(i=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4), j=c(1:4, 3:6, 7, 8, 8, 9), x=1, dims=c(4, 9),
                                 dimnames=list(c("A1", "A2", "B1", "B2"), paste0("c", 1:9)))
  md <- data.frame(sample=rep(c("s1", "s2", "s3"), each=3), row.names=colnames(counts))
  cc <- filter_counts(counts, lrdb, outfile=NULL)
  lrs <- generate_lrs_table(cc, lrdb, outfile=NULL)

  expect_warning(he <- generate_hyperedges(cc, lrs, outfile=NULL), "3 cells and 2 genes are not connected to the rest of the hypergraph")
  expect_setequal(rownames(generate_adjacency_matrix(he, outfile=NULL)), c("A1", "A2", paste0("c", 1:6)))
  expect_warning(res <- quietly(HyperCom(lrdb, counts, md, sample="sample", outdir=NULL)), "not connected")
  expect_identical(res$lrs$interaction, "IA")
  expect_identical(res$priority[["1"]]$interaction, "IA")

  inf <- generate_infMat(generate_adjacency_matrix(he, outfile=NULL), outfile=NULL, threads=1)
  meta <- generate_metadata(inf, md, sample="sample", outfile=NULL)
  expect_warning(pr <- prioritize_lr(lrs, generate_adjacency_matrix(he, outfile=NULL), inf, meta, outdir=NULL, verbose=FALSE, threads=1),
                 "1 interactions of lrs have genes that are not nodes of infMat")
  expect_identical(pr[["1"]]$interaction, "IA")
  expect_error(score_hypercom(inf, meta, ligand="B1", receptors="B2", threads=1), "not connected to the rest of the hypergraph")
})

test_that("generate_adjacency_matrix matches handwritten, by value and as handwritten's dsCMatrix", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)

  expect_same_matrix(adj, fx$output$adj)
  expect_identical(HyperCom:::adjacency_sparse(adj), fx$output$adj)
})

test_that("the adjacency matrix is symmetric with an empty diagonal and no isolated nodes", {
  A <- as.matrix(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))

  expect_identical(rownames(A), colnames(A))
  expect_true(isSymmetric(unname(A)))
  expect_true(all(diag(A) == 0))
  expect_true(all(rowSums(A) > 0))
  expect_true(all(A >= 0 & A == round(A)))
})

test_that("generate_infMat by lu (the default without hyperedges) is handwritten's bit for bit from every kind of adjacency", {
  adj <- fx$output$adj
  dense <- as.matrix(adj)
  dense.copy <- dense + 0
  inputs <- list(adj, methods::as(adj, "generalMatrix"), dense, as.data.frame(dense), fx_factored()$adj)

  for(r in c(0.5, 0.3)){
    for(x in inputs){
      expect_handwritten_infmat(generate_infMat(x, r=r, method=if(inherits(x, "hypercom_adjacency")) "lu", outfile=NULL), r=r)
    }
  }
  expect_identical(dense, dense.copy)
  expect_null(dimnames(generate_infMat(unname(dense), outfile=NULL)))
})

test_that("influence matrix rows are random walk with restart distributions", {
  for(r in c(0.5, 0.3)){
    infMat <- as.matrix(generate_infMat(fx$output$adj, r=r, outfile=NULL))

    expect_identical(rownames(infMat), rownames(fx$output$adj))
    expect_equal(unname(rowSums(infMat)), rep(1, nrow(infMat)), tolerance=1e-8)
    expect_true(all(infMat > -1e-12))
    expect_true(all(diag(infMat) >= r - 1e-12))
  }
})
