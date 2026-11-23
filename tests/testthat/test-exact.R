# Scores, priorities and seeded p-values equal handwritten's within 1e-12, at any thread count and however the pairs fall into batches.

test_that("score_hypercom equals handwritten within one ulp, also for a factor ligand", {
  for(n in names(fx_pairs)){
    expect_handwritten(score_hypercom(fx$output$infMat, fx$output$metadata, ligand=fx_pairs[[n]][[1]], receptors=fx_pairs[[n]][[2]]),
                       fx$output$scores[[n]])
  }
  expect_handwritten(score_hypercom(fx$output$infMat, fx$output$metadata, ligand=factor("Thbs1"), receptors="Sdc4"), fx$output$scores$single)
})

test_that("prioritize_lr equals handwritten with one thread and with several", {
  for(threads in c(1, test_threads(4))){
    expect_handwritten(run_priority(weights=c(1, 0.5, 0), threads=threads), fx$output$priority)
  }
})

test_that("seeded p-values equal handwritten's, repeat with set.seed, lie in [0, 1], are BH adjusted and leave the session's seed", {
  skip_on_cran()
  withr::local_seed(42)
  before <- get(".Random.seed", envir=globalenv())
  priority <- run_priority(significance=TRUE, seed=1, exact.priority=TRUE, threads=test_threads(2))
  df <- priority[["1"]]

  expect_identical(get(".Random.seed", envir=globalenv()), before)
  expect_handwritten(priority, fx$output$priority.sig)
  expect_identical(withr::with_seed(1, run_priority(significance=TRUE, exact.priority=TRUE)), priority)
  expect_true(all(df$p >= 0 & df$p <= 1, na.rm=TRUE))
  expect_identical(df$p.adj, stats::p.adjust(df$p, method="BH"))
})

test_that("prioritize_lr reads an adjacency without node x node matrix as handwritten's matrix", {
  skip_on_cran()
  adj <- fx_factored()$adj

  expect_handwritten(run_priority(adj=adj, weights=c(1, 0.5, 0), threads=test_threads(2)), fx$output$priority)
  expect_handwritten(run_priority(adj=adj, significance=TRUE, seed=1, exact.priority=TRUE, threads=test_threads(2)), fx$output$priority.sig)
})

test_that("priorities and seeded p-values do not depend on how the pairs fall into batches", {
  skip_on_cran()
  lrs <- fx$output$lrs[rev(seq_len(nrow(fx$output$lrs))), ]
  reversed <- run_priority(lrs=lrs, threads=1)[["1"]]
  batch <- HyperCom:::batch_pairs
  testthat::local_mocked_bindings(batch_pairs=function(ligands, receptors, size) batch(ligands, receptors, 1), .package="HyperCom")

  expect_gt(length(HyperCom:::batch_pairs(lrs$ligand, as.matrix(lrs[, 3:ncol(lrs)]), 256)), 1)
  expect_handwritten(plain_df(reversed[order(reversed$interaction), ]),
                     plain_df(fx$output$priority[["1"]][order(fx$output$priority[["1"]]$interaction), ]))
  expect_handwritten(run_priority(weights=c(1, 0.5, 0), threads=test_threads(2)), fx$output$priority)
  expect_handwritten(run_priority(significance=TRUE, seed=1, exact.priority=TRUE, threads=test_threads(2)), fx$output$priority.sig)
})
