# generate_infMat methods "chol" (v1-fast's Cholesky inverse) and "woodbury" (Woodbury through the incidence, kept factored): the
# same matrix as "lu" (handwritten, bit for bit) within ~1e-15 relative; scoring the factored matrix equals scoring its dense matrix bit for
# bit. An adjacency and an influence matrix without node x node matrices index like the matrices they stand for.

max_rel <- function(a, b){
  return(max(abs(as.matrix(a) - as.matrix(b)) / pmax(abs(as.matrix(b)), .Machine$double.xmin)))
}

test_that("methods chol and woodbury give handwritten's influence matrix within 1e-13", {
  adj <- fx$output$adj
  for(r in c(0.5, 0.3)){
    expected <- if(r == 0.5) fx$output$infMat else fx$output$infMat.r03
    chol <- generate_infMat(adj, r=r, outfile=NULL, method="chol")
    wood <- generate_infMat(adj, r=r, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges, threads=1)

    expect_identical(dimnames(chol), dimnames(expected))
    expect_identical(dimnames(wood), dimnames(expected))
    expect_lt(max_rel(chol, expected), 1e-13)
    expect_lt(max_rel(wood, expected), 1e-13)
  }
})

test_that("method chol reads sparse, general, lower and dense adjacency alike and leaves a dense input unchanged", {
  adj <- fx$output$adj
  dense <- as.matrix(adj)
  dense.copy <- dense + 0
  chol <- generate_infMat(adj, outfile=NULL, method="chol")

  expect_identical(generate_infMat(methods::as(adj, "generalMatrix"), outfile=NULL, method="chol"), chol)
  expect_identical(generate_infMat(Matrix::t(adj), outfile=NULL, method="chol"), chol)
  expect_identical(generate_infMat(dense, outfile=NULL, method="chol"), chol)
  expect_identical(dense, dense.copy)
})

test_that("method woodbury is the same at any thread count, from factor hyperedges and from a dense adjacency", {
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  wood <- generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=1)
  factors <- as.data.frame(hyperedges)
  factors$edge <- factor(factors$edge)
  factors$node <- factor(factors$node)

  expect_identical(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=test_threads(3)), wood)
  expect_identical(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=factors, threads=test_threads(2)), wood)
  expect_identical(generate_infMat(as.matrix(adj), outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=1), wood)
})

test_that("method woodbury counts nodes listed twice and ignores hyperedges of one distinct node", {
  hyperedges <- data.frame(edge=c("e1", "e1", "e1", "e2", "e2", "e3", "e3", "e3", "e4", "e4", "e5", "e6"),
                           node=c("a", "b", "b", "b", "c", "c", "d", "e", "a", "a", "f", "d"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)

  expect_identical(rownames(adj), c("a", "b", "c", "d", "e"))
  for(r in c(0.5, 0.3, 1)){
    lu <- generate_infMat(adj, r=r, outfile=NULL, method="lu")

    expect_lt(max_rel(generate_infMat(adj, r=r, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=test_threads(2)), lu),
              1e-13)
    expect_lt(max_rel(generate_infMat(adj, r=r, outfile=NULL, method="chol"), lu), 1e-13)
  }
})

test_that("methods chol and woodbury stop on input they cannot invert", {
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  dense <- as.matrix(adj)
  changed <- dense
  changed[1, 2] <- changed[2, 1] <- changed[1, 2] + 1
  asymmetric <- dense
  asymmetric[1, 2] <- asymmetric[1, 2] + 1
  isolated <- rbind(cbind(dense, z=0), z=0)
  extra <- rbind(as.data.frame(hyperedges)[1, ], as.data.frame(hyperedges))

  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury"), "needs the hyperedges")
  expect_error(generate_infMat(adj, r=0, outfile=NULL, method="chol"), "r must be one number")
  expect_error(generate_infMat(adj, r=1.5, outfile=NULL, method="woodbury", hyperedges=hyperedges), "r must be one number")
  expect_error(generate_infMat(asymmetric, outfile=NULL, method="chol"), "symmetric")
  expect_error(generate_infMat(isolated, outfile=NULL, method="chol"), "positive row sums")
  expect_error(generate_infMat(isolated, outfile=NULL, method="woodbury", hyperedges=hyperedges), "positive row sums")
  expect_error(generate_infMat(changed, outfile=NULL, method="woodbury", hyperedges=hyperedges), "do not match adj")
  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=extra), "do not match adj")
  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges[-1, ]), "do not match adj")
  expect_error(generate_infMat(unname(dense), outfile=NULL, method="woodbury", hyperedges=hyperedges), "do not match adj")
  expect_error(generate_infMat(adj, outfile=NULL, method="qr"), "should be one of")
})

test_that("method woodbury takes the nodes of adj in any order and falls back to chol when hyperedges are at least as many as nodes", {
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  o <- rev(seq_len(nrow(adj)))
  dense <- as.matrix(adj)[o, o]
  few <- data.frame(edge=c("e1", "e1", "e2", "e2", "e3", "e3", "e4", "e4", "e5", "e5", "e5"),
                    node=c("a", "b", "a", "b", "b", "c", "a", "c", "a", "b", "c"))
  adj.few <- generate_adjacency_matrix(few, outfile=NULL)

  expect_lt(ncol(HyperCom:::incidence_matrix(hyperedges, adj)), nrow(adj))
  expect_lt(max_rel(generate_infMat(dense, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=test_threads(2)),
                    generate_infMat(dense, outfile=NULL)), 1e-13)
  expect_identical(generate_infMat(adj.few, outfile=NULL, method="woodbury", hyperedges=few),
                   generate_infMat(adj.few, outfile=NULL, method="chol"))
})

test_that("method woodbury rejects hyperedges with the same nodes and row sums but other neighbours", {
  given <- data.frame(edge=c("e1", "e1", "e2", "e2", "e3", "e3", "e4", "e4", "e5", "e5", "e5", "e5", "e5"),
                      node=c("a", "c", "b", "d", "a", "b", "c", "d", "a", "b", "c", "d", "x"))
  other <- given
  other$node[c(2, 4)] <- c("d", "c")
  adj <- HyperCom:::adjacency_sparse(generate_adjacency_matrix(other, outfile=NULL))
  adj.given <- HyperCom:::adjacency_sparse(generate_adjacency_matrix(given, outfile=NULL))

  expect_identical(rownames(adj), rownames(adj.given))
  expect_identical(Matrix::rowSums(adj), Matrix::rowSums(adj.given))
  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=given), "other entries")
})

test_that("methods chol and woodbury stop on a computationally singular system and on non-finite weights", {
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  nan <- as.matrix(adj)
  nan[1, 2] <- nan[2, 1] <- NaN

  expect_error(generate_infMat(adj, r=1e-300, outfile=NULL, method="chol"), "singular|not positive")
  expect_error(generate_infMat(adj, r=1e-300, outfile=NULL, method="woodbury", hyperedges=hyperedges), "singular|not positive")
  expect_error(generate_infMat(nan, outfile=NULL, method="chol"), "finite")
  expect_error(generate_infMat(nan, outfile=NULL, method="woodbury", hyperedges=hyperedges), "finite")
})

test_that("HyperCom with infmat.method lu writes handwritten's adjacency and influence matrix bit for bit", {
  rundir <- file.path(withr::local_tempdir(), "run")
  res <- run_hypercom(rundir, threads=test_threads(2), infmat.method="lu")

  expect_run_result(res, rundir)
  expect_identical(HyperCom:::adjacency_sparse(res$adj), fx$output$adj)
  expect_handwritten_infmat(res$infMat)
})

test_that("generate_infMat on an adjacency without node x node matrix: factored by default, lu as handwritten, chol as from the matrix", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))

  expect_s3_class(adj, "hypercom_adjacency")
  expect_s3_class(infMat, "hypercom_influence")
  expect_identical(infMat, generate_infMat(fx$output$adj, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges,
                                           threads=test_threads(3)))
  expect_handwritten_infmat(generate_infMat(adj, outfile=NULL, method="lu"))
  expect_handwritten_infmat(generate_infMat(adj, r=0.3, outfile=NULL, method="lu"), r=0.3)
  expect_identical(generate_infMat(adj, outfile=NULL, method="chol"), generate_infMat(fx$output$adj, outfile=NULL, method="chol"))
})

test_that("scoring a factored influence matrix equals scoring its dense matrix bit for bit", {
  skip_on_cran()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges, threads=test_threads(2))
  dense <- as.matrix(infMat)
  md <- fx$output$metadata
  md.rev <- md[rev(seq_len(nrow(md))), ]

  expect_s3_class(infMat, "hypercom_influence")
  for(threads in c(1, test_threads(3))){
    expect_identical(run_priority(adj=adj, infMat=infMat, weights=c(1, 0.5, 0), threads=threads),
                     run_priority(infMat=dense, weights=c(1, 0.5, 0), threads=threads))
  }
  expect_identical(run_priority(adj=adj, infMat=infMat, significance=TRUE, seed=1, threads=test_threads(2)),
                   run_priority(infMat=dense, significance=TRUE, seed=1, threads=test_threads(2)))
  for(p in fx_pairs){
    expect_identical(score_hypercom(infMat, md, p[[1]], p[[2]]), score_hypercom(dense, md, p[[1]], p[[2]]))
    expect_identical(score_hypercom(infMat, md.rev, p[[1]], p[[2]]), score_hypercom(dense, md.rev, p[[1]], p[[2]]))
  }
})

test_that("generate_infMat defaults to woodbury given the hyperedges, and an unknown method stops", {
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  infMat <- generate_infMat(adj, outfile=NULL, hyperedges=hyperedges, threads=test_threads(2))

  expect_s3_class(infMat, "hypercom_influence")
  expect_identical(infMat, generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=test_threads(2)))
  expect_error(generate_infMat(adj, outfile=NULL, method="qr", hyperedges=hyperedges), "should be one of")
})

test_that("adjacency and influence matrix without node x node matrices index like the matrices they stand for", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))
  A <- as.matrix(adj)
  F <- as.matrix(infMat)
  nodes <- rownames(A)
  index <- list(list(5, nodes[c(3, 9)]), list(c(9, 1, 9), -1), list(nodes[1:4] %in% nodes[c(2, 4)], 7), list(3, 3),
                list(factor(nodes[c(10, 2)]), c(2, 10, 2)), list(nodes[200], nodes))

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
