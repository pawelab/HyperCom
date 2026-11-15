# Scores and priorities equal master's bit for bit (identical(), not within a tolerance), with one thread or several, where long
# double is the x87 type the fixture was made with (x86-64; within 1e-12 elsewhere, see expect_master). The influence matrix is
# compared bit for bit with master's formula run on this machine (its last bits depend on the BLAS) and within 1e-12 with the fixture.

test_that("score_hypercom equals master bit for bit", {
  for(n in names(fx_pairs)){
    expect_master(score_hypercom(fx$output$infMat, fx$output$metadata, ligand=fx_pairs[[n]][[1]], receptors=fx_pairs[[n]][[2]]),
                  fx$output$scores[[n]])
  }
})

test_that("prioritize_lr equals master bit for bit with one thread and with several", {
  for(threads in c(1, test_threads(4))){
    expect_master(run_priority(weights=c(1, 0.5, 0), threads=threads), fx$output$priority)
  }
})

test_that("seeded p-values equal master's with seed=1 and after set.seed(1), lie in [0, 1] and are BH adjusted", {
  skip_on_cran()
  priority <- run_priority(significance=TRUE, seed=1, threads=test_threads(2))
  df <- priority[["1"]]

  expect_master(priority, fx$output$priority.sig)
  expect_identical(withr::with_seed(1, run_priority(significance=TRUE)), priority)
  expect_true(all(df$p >= 0 & df$p <= 1, na.rm=TRUE))
  expect_identical(df$p.adj, stats::p.adjust(df$p, method="BH"))
})

test_that("prioritize_lr reads an adjacency without node x node matrix as master's matrix, bit for bit", {
  skip_on_cran()
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)

  expect_master(run_priority(adj=adj, weights=c(1, 0.5, 0), threads=test_threads(2)), fx$output$priority)
  expect_master(run_priority(adj=adj, significance=TRUE, seed=1, threads=test_threads(2)), fx$output$priority.sig)
})

test_that("a factor ligand is read by name, as in master", {
  score <- score_hypercom(fx$output$infMat, fx$output$metadata, ligand=factor("Thbs1"), receptors="Sdc4")

  expect_master(score, fx$output$scores$single)
})

test_that("the influence matrix equals master's from a general and a dense adjacency, and leaves a dense input unchanged", {
  adj <- fx$output$adj
  dense <- as.matrix(adj)
  dense.copy <- dense + 0

  expect_master_infmat(generate_infMat(methods::as(adj, "generalMatrix"), outfile=NULL))
  expect_master_infmat(generate_infMat(dense, outfile=NULL))
  expect_identical(dense, dense.copy)
})

test_that("priorities do not depend on how the pairs fall into batches", {
  lrs <- fx$output$lrs[rev(seq_len(nrow(fx$output$lrs))), ]
  priority <- run_priority(lrs=lrs, threads=1)

  a <- priority[["1"]][order(priority[["1"]]$interaction), ]
  b <- fx$output$priority[["1"]][order(fx$output$priority[["1"]]$interaction), ]
  rownames(a) <- NULL
  rownames(b) <- NULL

  expect_master(a, b)
})

test_that("priorities and seeded p-values are the same when every ligand is its own batch", {
  skip_on_cran()
  lrs <- fx$output$lrs
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))
  factored <- run_priority(lrs, adj, infMat, weights=c(1, 0.5, 0), threads=test_threads(2))
  batch <- HyperCom:::batch_pairs
  testthat::local_mocked_bindings(batch_pairs=function(ligands, receptors, size) batch(ligands, receptors, 1), .package="HyperCom")

  expect_gt(length(HyperCom:::batch_pairs(lrs$ligand, as.matrix(lrs[, 3:ncol(lrs)]), 256)), 1)
  expect_master(run_priority(weights=c(1, 0.5, 0), threads=test_threads(2)), fx$output$priority)
  expect_master(run_priority(significance=TRUE, seed=1, threads=test_threads(2)), fx$output$priority.sig)
  expect_identical(run_priority(lrs, adj, infMat, weights=c(1, 0.5, 0), threads=test_threads(2)), factored)
})
