# generate_infMat methods "chol" (v1-fast's Cholesky inverse) and "woodbury" (Woodbury through the incidence, kept factored): the
# same matrix as "lu" (master, bit for bit) within ~1e-15 relative; scoring the factored matrix equals scoring its dense matrix bit for
# bit. Builds without the method argument (master) skip these tests.

skip_if_no_methods <- function(){
  testthat::skip_if_not("method" %in% names(formals(generate_infMat)), "build without generate_infMat methods")
}

max_rel <- function(a, b){
  return(max(abs(as.matrix(a) - as.matrix(b)) / pmax(abs(as.matrix(b)), .Machine$double.xmin)))
}

test_that("methods chol and woodbury give master's influence matrix within 1e-13", {
  skip_if_no_methods()
  adj <- fx$output$adj
  for(r in c(0.5, 0.3)){
    expected <- if(r == 0.5) fx$output$infMat else fx$output$infMat.r03
    chol <- generate_infMat(adj, r=r, outfile=NULL, method="chol")
    hyper <- generate_infMat(adj, r=r, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges, threads=1)

    expect_identical(dimnames(chol), dimnames(expected))
    expect_identical(dimnames(hyper), dimnames(expected))
    expect_lt(max_rel(chol, expected), 1e-13)
    expect_lt(max_rel(hyper, expected), 1e-13)
  }
})

test_that("method chol reads sparse, general, lower and dense adjacency alike and leaves a dense input unchanged", {
  skip_if_no_methods()
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
  skip_if_no_methods()
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  hyper <- generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=1)
  factors <- as.data.frame(hyperedges)
  factors$edge <- factor(factors$edge)
  factors$node <- factor(factors$node)

  expect_identical(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=test_threads(3)), hyper)
  expect_identical(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=factors, threads=2), hyper)
  expect_identical(generate_infMat(as.matrix(adj), outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=1), hyper)
})

test_that("method woodbury counts nodes listed twice and ignores hyperedges of one distinct node", {
  skip_if_no_methods()
  hyperedges <- data.frame(edge=c("e1", "e1", "e1", "e2", "e2", "e3", "e3", "e3", "e4", "e4", "e5", "e6"),
                           node=c("a", "b", "b", "b", "c", "c", "d", "e", "a", "a", "f", "d"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)

  expect_identical(rownames(adj), c("a", "b", "c", "d", "e"))
  for(r in c(0.5, 0.3, 1)){
    lu <- generate_infMat(adj, r=r, outfile=NULL, method="lu")

    expect_lt(max_rel(generate_infMat(adj, r=r, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=2), lu), 1e-13)
    expect_lt(max_rel(generate_infMat(adj, r=r, outfile=NULL, method="chol"), lu), 1e-13)
  }
})

test_that("methods chol and woodbury stop on input they cannot invert", {
  skip_if_no_methods()
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

test_that("method woodbury takes the nodes of adj in any order and falls back to chol without fewer hyperedges than nodes", {
  skip_if_no_methods()
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  o <- rev(seq_len(nrow(adj)))
  dense <- as.matrix(adj)[o, o]
  few <- data.frame(edge=c("e1", "e1", "e2", "e2", "e3", "e3", "e4", "e4", "e5", "e5", "e5"),
                    node=c("a", "b", "a", "b", "b", "c", "a", "c", "a", "b", "c"))
  adj.few <- generate_adjacency_matrix(few, outfile=NULL)

  expect_lt(ncol(HyperCom:::incidence_matrix(hyperedges, adj)), nrow(adj))
  expect_lt(max_rel(generate_infMat(dense, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=2),
                    generate_infMat(dense, outfile=NULL)), 1e-13)
  expect_identical(generate_infMat(adj.few, outfile=NULL, method="woodbury", hyperedges=few),
                   generate_infMat(adj.few, outfile=NULL, method="chol"))
})

test_that("method woodbury rejects hyperedges with the same nodes and row sums but other neighbours", {
  skip_if_no_methods()
  given <- data.frame(edge=c("e1", "e1", "e2", "e2", "e3", "e3", "e4", "e4", "e5", "e5", "e5", "e5", "e5"),
                      node=c("a", "c", "b", "d", "a", "b", "c", "d", "a", "b", "c", "d", "x"))
  other <- given
  other$node[c(2, 4)] <- c("d", "c")
  adj <- as_master_adjacency(generate_adjacency_matrix(other, outfile=NULL))
  adj.given <- as_master_adjacency(generate_adjacency_matrix(given, outfile=NULL))

  expect_identical(rownames(adj), rownames(adj.given))
  expect_identical(Matrix::rowSums(adj), Matrix::rowSums(adj.given))
  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=given), "other entries")
})

test_that("methods chol and woodbury stop on a computationally singular system and on non-finite weights", {
  skip_if_no_methods()
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  nan <- as.matrix(adj)
  nan[1, 2] <- nan[2, 1] <- NaN

  expect_error(generate_infMat(adj, r=1e-300, outfile=NULL, method="chol"), "singular|not positive")
  expect_error(generate_infMat(adj, r=1e-300, outfile=NULL, method="woodbury", hyperedges=hyperedges), "singular|not positive")
  expect_error(generate_infMat(nan, outfile=NULL, method="chol"), "finite")
  expect_error(generate_infMat(nan, outfile=NULL, method="woodbury", hyperedges=hyperedges), "finite")
})

test_that("HyperCom with infmat.method woodbury ranks the pairs as master", {
  skip_if_no_methods()
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  res <- quietly(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition", outdir="run",
                          threads=2, infmat.method="woodbury"))
  priority <- utils::read.csv(file.path("run", "6.priority", "priority.1.csv"), check.names=FALSE)

  expect_identical(res$infMat, readRDS(file.path("run", "4.infMat.rds")))
  expect_lt(max_rel(readRDS(file.path("run", "4.infMat.rds")), fx$output$infMat), 1e-13)
  expect_same_priority(priority, fx$output$priority[["1"]])
})

test_that("HyperCom with infmat.method lu writes master's adjacency and influence matrix bit for bit", {
  skip_if_no_methods()
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  res <- quietly(HyperCom(fx$input$lrdb, fx$input$counts, fx$input$metadata, sample="sample", group="condition", outdir="run",
                          threads=2, infmat.method="lu"))
  priority <- utils::read.csv(file.path("run", "6.priority", "priority.1.csv"), check.names=FALSE)

  expect_master_infmat(res$infMat)
  expect_identical(as_master_adjacency(readRDS(file.path("run", "3.adj.rds"))), fx$output$adj)
  expect_master_infmat(readRDS(file.path("run", "4.infMat.rds")))
  expect_same_priority(priority, fx$output$priority[["1"]])
})

test_that("generate_infMat on an adjacency without node x node matrix: factored by default, lu as master, chol as from the matrix", {
  skip_if_no_methods()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  skip_if_not(inherits(adj, "hypercom_adjacency"), "build with node x node matrices")
  infMat <- generate_infMat(adj, outfile=NULL, threads=2)

  expect_s3_class(infMat, "hypercom_influence")
  expect_identical(infMat, generate_infMat(fx$output$adj, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges, threads=test_threads(3)))
  expect_lt(max_rel(infMat, fx$output$infMat), 1e-13)
  expect_master_infmat(generate_infMat(adj, outfile=NULL, method="lu"))
  expect_master_infmat(generate_infMat(adj, r=0.3, outfile=NULL, method="lu"), r=0.3)
  expect_identical(generate_infMat(adj, outfile=NULL, method="chol"), generate_infMat(fx$output$adj, outfile=NULL, method="chol"))
})

test_that("scoring a factored influence matrix equals scoring its dense matrix bit for bit", {
  skip_if_no_methods()
  skip_on_cran()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=fx$output$hyperedges, threads=2)
  skip_if_not(inherits(infMat, "hypercom_influence"), "build without a factored influence matrix")
  dense <- as.matrix(infMat)
  md <- fx$output$metadata
  md.rev <- md[rev(seq_len(nrow(md))), ]
  lrs <- fx$output$lrs
  pairs <- list(list(NULL, NULL), list("Thbs1", "Sdc4"), list("Tgfb1", c("Tgfbr1", "Tgfbr2")))

  for(threads in c(1, test_threads(3))){
    expect_identical(quietly(prioritize_lr(lrs, adj, infMat, md, weights=c(1, 0.5, 0), threads=threads, outdir=NULL)),
                     quietly(prioritize_lr(lrs, fx$output$adj, dense, md, weights=c(1, 0.5, 0), threads=threads, outdir=NULL)))
  }
  expect_identical(quietly(prioritize_lr(lrs, adj, infMat, md, significance=TRUE, seed=1, threads=2, outdir=NULL)),
                   quietly(prioritize_lr(lrs, fx$output$adj, dense, md, significance=TRUE, seed=1, threads=2, outdir=NULL)))
  for(p in pairs){
    expect_identical(score_hypercom(infMat, md, p[[1]], p[[2]]), score_hypercom(dense, md, p[[1]], p[[2]]))
    expect_identical(score_hypercom(infMat, md.rev, p[[1]], p[[2]]), score_hypercom(dense, md.rev, p[[1]], p[[2]]))
  }
})

test_that("generate_infMat defaults to woodbury, and to lu for an adjacency matrix given without hyperedges", {
  skip_if_no_methods()
  skip_if_not(identical(eval(formals(generate_infMat)$method)[1], "woodbury"), "build without the woodbury default")
  adj <- fx$output$adj
  hyperedges <- fx$output$hyperedges
  infMat <- generate_infMat(adj, outfile=NULL, hyperedges=hyperedges, threads=2)

  expect_s3_class(infMat, "hypercom_influence")
  expect_identical(infMat, generate_infMat(adj, outfile=NULL, method="woodbury", hyperedges=hyperedges, threads=2))
  expect_master_infmat(generate_infMat(adj, outfile=NULL))
  expect_error(generate_infMat(adj, outfile=NULL, method="woodbury"), "needs the hyperedges")
  expect_error(generate_infMat(adj, outfile=NULL, method="hypergraph", hyperedges=hyperedges), "now called \"woodbury\"")
})
