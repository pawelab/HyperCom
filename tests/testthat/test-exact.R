# Scores and priorities equal master's bit for bit (identical(), not within a tolerance), with one thread or several, where long
# double is the x87 type the fixture was made with (x86-64; within 1e-12 elsewhere, see expect_master). Builds without the threads,
# seed or verbose arguments (master) skip those tests. The influence matrix is compared bit for bit with master's formula run on
# this machine (its last bits depend on the BLAS) and within 1e-12 with the fixture.

has_args <- function(...){
  return(all(c(...) %in% names(formals(prioritize_lr))))
}

test_that("score_hypercom equals master bit for bit", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata
  pairs <- list(network=list(NULL, NULL), single=list("Thbs1", "Sdc4"), complex=list("Tgfb1", c("Tgfbr1", "Tgfbr2")),
                triple=list("Tgfb1", c("Acvr1", "Tgfbr1", "Tgfbr2")), same=as.list(fx$output$same.pair))

  for(n in names(pairs)){
    expect_master(score_hypercom(infMat, md, ligand=pairs[[n]][[1]], receptors=pairs[[n]][[2]]), fx$output$scores[[n]])
  }
})

test_that("prioritize_lr equals master bit for bit with one thread and with several", {
  skip_if_not(has_args("threads"), "build without a threads argument")
  for(threads in c(1, test_threads(4))){
    priority <- quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5, 0),
                                      threads=threads, outdir=NULL))

    expect_master(priority, fx$output$priority)
  }
})

test_that("seeded p-values equal master's after set.seed", {
  skip_if_not(has_args("seed", "threads"), "build without seed and threads arguments")
  skip_on_cran()
  priority <- quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, significance=TRUE, seed=1,
                                    threads=2, outdir=NULL))

  expect_master(priority, fx$output$priority.sig)
})

test_that("prioritize_lr reads an adjacency without node x node matrix as master's matrix, bit for bit", {
  skip_if_not(has_args("seed", "threads"), "build without seed and threads arguments")
  skip_on_cran()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  priority <- quietly(prioritize_lr(fx$output$lrs, adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5, 0), threads=2,
                                    outdir=NULL))
  sig <- quietly(prioritize_lr(fx$output$lrs, adj, fx$output$infMat, fx$output$metadata, significance=TRUE, seed=1, threads=2,
                               outdir=NULL))

  expect_master(priority, fx$output$priority)
  expect_master(sig, fx$output$priority.sig)
})

test_that("a factor ligand is read by name, as in master", {
  score <- score_hypercom(fx$output$infMat, fx$output$metadata, ligand=factor("Thbs1"), receptors="Sdc4")

  expect_master(score, fx$output$scores$single)
})

test_that("the adjacency matrix is master's object, not only its values", {
  expect_identical(as_master_adjacency(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)), fx$output$adj)
})

test_that("the influence matrix equals master's from sparse, general and dense adjacency, and leaves a dense input unchanged", {
  adj <- fx$output$adj
  dense <- as.matrix(adj)
  dense.copy <- dense + 0

  expect_master_infmat(generate_infMat(adj, outfile=NULL))
  expect_master_infmat(generate_infMat(adj, r=0.3, outfile=NULL), r=0.3)
  expect_master_infmat(generate_infMat(methods::as(adj, "generalMatrix"), outfile=NULL))
  expect_master_infmat(generate_infMat(dense, outfile=NULL))
  expect_identical(dense, dense.copy)
})

test_that("priorities do not depend on how the pairs fall into batches", {
  skip_if_not(has_args("threads"), "build without a threads argument")
  lrs <- fx$output$lrs[rev(seq_len(nrow(fx$output$lrs))), ]
  priority <- quietly(prioritize_lr(lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, threads=1, outdir=NULL))

  a <- priority[["1"]][order(priority[["1"]]$interaction), ]
  b <- fx$output$priority[["1"]][order(fx$output$priority[["1"]]$interaction), ]
  rownames(a) <- NULL
  rownames(b) <- NULL

  expect_master(a, b)
})

test_that("priorities and seeded p-values are the same when every ligand is its own batch", {
  skip_if_not(has_args("seed", "threads"), "build without seed and threads arguments")
  skip_if_not(exists("batch_pairs", envir=asNamespace("HyperCom"), inherits=FALSE), "build without batches")
  skip_on_cran()
  lrs <- fx$output$lrs
  md <- fx$output$metadata
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=2)
  factored <- quietly(prioritize_lr(lrs, adj, infMat, md, weights=c(1, 0.5, 0), threads=2, outdir=NULL))
  batch <- HyperCom:::batch_pairs
  testthat::local_mocked_bindings(batch_pairs=function(ligands, receptors, size) batch(ligands, receptors, 1), .package="HyperCom")

  expect_gt(length(HyperCom:::batch_pairs(lrs$ligand, as.matrix(lrs[, 3:ncol(lrs)]), 256)), 1)
  expect_master(quietly(prioritize_lr(lrs, fx$output$adj, fx$output$infMat, md, weights=c(1, 0.5, 0), threads=2, outdir=NULL)),
                fx$output$priority)
  expect_master(quietly(prioritize_lr(lrs, fx$output$adj, fx$output$infMat, md, significance=TRUE, seed=1, threads=2, outdir=NULL)),
                fx$output$priority.sig)
  expect_identical(quietly(prioritize_lr(lrs, adj, infMat, md, weights=c(1, 0.5, 0), threads=2, outdir=NULL)), factored)
})
