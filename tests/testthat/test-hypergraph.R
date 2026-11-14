test_that("generate_hyperedges matches master", {
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

test_that("generate_adjacency_matrix matches master", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)

  expect_same_matrix(adj, fx$output$adj)
})

test_that("the adjacency matrix is symmetric with an empty diagonal and no isolated nodes", {
  A <- as.matrix(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL))

  expect_identical(rownames(A), colnames(A))
  expect_true(isSymmetric(unname(A)))
  expect_true(all(diag(A) == 0))
  expect_true(all(rowSums(A) > 0))
  expect_true(all(A >= 0 & A == round(A)))
})

test_that("generate_infMat matches master", {
  infMat <- generate_infMat(fx$output$adj, outfile=NULL)
  infMat.r03 <- generate_infMat(fx$output$adj, r=0.3, outfile=NULL)

  expect_same_matrix(infMat, fx$output$infMat, tolerance=tol)
  expect_same_matrix(infMat.r03, fx$output$infMat.r03, tolerance=tol)
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
