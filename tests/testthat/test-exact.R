# Scores and priorities equal master's bit for bit (identical(), not within a tolerance), with one thread or several.
# Builds without the threads, seed or verbose arguments (master) skip those tests.

has_args <- function(...){
  return(all(c(...) %in% names(formals(prioritize_lr))))
}

test_that("score_hypercom equals master bit for bit", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata
  pairs <- list(network=list(NULL, NULL), single=list("Thbs1", "Sdc4"), complex=list("Tgfb1", c("Tgfbr1", "Tgfbr2")),
                triple=list("Tgfb1", c("Acvr1", "Tgfbr1", "Tgfbr2")), same=as.list(fx$output$same.pair))

  for(n in names(pairs)){
    expect_identical(score_hypercom(infMat, md, ligand=pairs[[n]][[1]], receptors=pairs[[n]][[2]]), fx$output$scores[[n]])
  }
})

test_that("prioritize_lr equals master bit for bit with one thread and with several", {
  skip_if_not(has_args("threads"), "build without a threads argument")
  for(threads in c(1, 4)){
    priority <- quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5, 0),
                                      threads=threads, outdir=NULL))

    expect_identical(priority, fx$output$priority)
  }
})

test_that("seeded p-values equal master's after set.seed", {
  skip_if_not(has_args("seed", "threads"), "build without seed and threads arguments")
  skip_on_cran()
  priority <- quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, significance=TRUE, seed=1,
                                    threads=2, outdir=NULL))

  expect_identical(priority, fx$output$priority.sig)
})

test_that("a factor ligand is read by name, as in master", {
  score <- score_hypercom(fx$output$infMat, fx$output$metadata, ligand=factor("Thbs1"), receptors="Sdc4")

  expect_identical(score, fx$output$scores$single)
})

test_that("a network score asked with several receptors averages copies of the network score, as in master", {
  net <- score_hypercom(fx$output$infMat, fx$output$metadata)
  three <- score_hypercom(fx$output$infMat, fx$output$metadata, receptors=c("Sdc4", "Tgfbr1", "Tgfbr2"))
  three <- three[match(net$cell, three$cell), ]

  expect_identical(three$sender, (net$sender + net$sender + net$sender) / 3)
  expect_identical(three$receiver, (net$receiver + net$receiver + net$receiver) / 3)
})

test_that("the adjacency matrix is master's object, not only its values", {
  expect_identical(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL), fx$output$adj)
})

test_that("the influence matrix equals master's from sparse, general and dense adjacency, and leaves a dense input unchanged", {
  adj <- fx$output$adj
  dense <- as.matrix(adj)
  dense.copy <- dense + 0

  expect_identical(generate_infMat(adj, outfile=NULL), fx$output$infMat)
  expect_identical(generate_infMat(adj, r=0.3, outfile=NULL), fx$output$infMat.r03)
  expect_identical(generate_infMat(methods::as(adj, "generalMatrix"), outfile=NULL), fx$output$infMat)
  expect_identical(generate_infMat(dense, outfile=NULL), fx$output$infMat)
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

  expect_identical(a, b)
})
