# Memory mode (log-domain tiles, passes) and screened real pairs: both depart from exact scoring by approved, measured amounts
# (decisions 2026-11-20 and 2026-11-22), so they are checked against the exact kernels within those, and against themselves exactly.

factored_infmat <- function(){
  return(generate_infMat(generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL), outfile=NULL, threads=test_threads(2)))
}

# memory (GB) that leaves the exact kernels out and fits memory mode
memory_for <- function(metadata){
  cells.idx <- metadata$sample != "gene"
  need <- HyperCom:::kernel_memory(sum(as.numeric(table(metadata$sample[cells.idx]))^2))
  return((need[["tiles"]] + 1) / 1e9)
}

test_that("the kernel is exact when it fits, memory mode in 3 passes otherwise, and a stop with the estimate", {
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  expect_identical(HyperCom:::choose_kernel(factored, metadata, Inf)$type, "exact")
  kernel <- HyperCom:::choose_kernel(factored, metadata, memory_for(metadata))
  expect_identical(c(kernel$type, kernel$passes), c("tiles", "3"))
  expect_error(HyperCom:::choose_kernel(factored, metadata, 1e-6), "scoring needs about .* GB of memory")
  expect_error(HyperCom:::choose_kernel(fx$output$infMat, metadata, memory_for(metadata)), "memory mode needs an influence matrix")
})

test_that("memory mode scores within the log-domain keys' rounding of exact, the same for any number of passes", {
  skip_on_cran()
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  g <- HyperCom:::grid_cells(rownames(factored), metadata)
  grids <- HyperCom:::pair_grids(as.character(fx$output$lrs$ligand), as.matrix(fx$output$lrs[, -(1:2)]))
  exact <- HyperCom:::score_grids(factored, metadata, grids$ligand, grids$receptor, test_threads(2))
  tiles <- lapply(1:4, function(passes){
    lb <- HyperCom:::lambda_blocks(factored, g, test_threads(2), passes)
    HyperCom:::score_grids(factored, metadata, grids$ligand, grids$receptor, test_threads(2), lb)
  })

  for(passes in 2:4){
    expect_identical(tiles[[passes]], tiles[[1]])
  }
  # a same-gene grid's entries (i, j) and (j, i) are equal: the keys keep them tied, the exact kernel orders them by rounding, one
  # rank of a block apart (here 1 / 80^2)
  same <- grids$ligand == grids$receptor
  n <- min(table(metadata$sample[metadata$sample != "gene"]))
  for(k in c("sender", "receiver")){
    expect_equal(tiles[[1]][[k]][, !same], exact[[k]][, !same], tolerance=1e-9)
    expect_equal(tiles[[1]][[k]][, same], exact[[k]][, same], tolerance=2 / n^2)
  }
})

test_that("memory mode screening bounds memory mode scores", {
  skip_on_cran()
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  g <- HyperCom:::grid_cells(rownames(factored), metadata)
  lb <- HyperCom:::lambda_blocks(factored, g, test_threads(2), 3L)
  grids <- HyperCom:::pair_grids(as.character(fx$output$lrs$ligand), as.matrix(fx$output$lrs[, -(1:2)]))
  sr <- HyperCom:::score_grids(factored, metadata, grids$ligand, grids$receptor, test_threads(2), lb)
  for(threads in c(1, 3)){
    sb <- HyperCom:::screen_grids(factored, metadata, grids$ligand, grids$receptor, test_threads(threads), lb)
    ok <- !is.na(sr$sender)
    expect_true(all(sb$sender.lo[ok] <= sr$sender[ok] & sr$sender[ok] <= sb$sender.hi[ok]))
    expect_true(all(sb$receiver.lo[ok] <= sr$receiver[ok] & sr$receiver[ok] <= sb$receiver.hi[ok]))
  }
})

test_that("prioritize_lr in memory mode gives the p-values and order of exact scoring here, priorities within rounding", {
  skip_on_cran()
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  exact <- run_priority(infMat=factored, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), memory=Inf)
  tiles <- run_priority(infMat=factored, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2),
                        memory=memory_for(metadata))

  for(w in names(exact)){
    expect_identical(tiles[[w]]$interaction, exact[[w]]$interaction)
    expect_equal(tiles[[w]]$priority, exact[[w]]$priority, tolerance=1e-9)
    expect_identical(tiles[[w]][c("p", "p.adj")], exact[[w]][c("p", "p.adj")])
  }
})

test_that("screened real pairs give the p-values and order of exact scoring, priorities near exact", {
  skip_on_cran()
  factored <- factored_infmat()
  for(infMat in list(fx$output$infMat, factored)){
    exact <- run_priority(infMat=infMat, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), exact.priority=TRUE)
    screened <- run_priority(infMat=infMat, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), exact.priority=FALSE)
    for(w in names(exact)){
      expect_identical(screened[[w]]$interaction, exact[[w]]$interaction)
      expect_identical(screened[[w]][c("p", "p.adj")], exact[[w]][c("p", "p.adj")])
      expect_equal(screened[[w]]$priority, exact[[w]]$priority, tolerance=1e-6)
    }
  }
  expect_identical(quietly(prioritize_lr(fx$output$lrs, fx$output$adj, fx$output$infMat, fx$output$metadata, weights=c(1, 0.5),
                                         significance=TRUE, seed=1, outdir=NULL)),
                   run_priority(weights=c(1, 0.5), significance=TRUE, seed=1, exact.priority=FALSE))
  expect_identical(run_priority(weights=c(1, 0.5), exact.priority=FALSE), run_priority(weights=c(1, 0.5)))
})

test_that("intervals that meet are found", {
  expect_identical(HyperCom:::meets(c(0, 2, 5), c(1, 3, 6), c(1, 4), c(1.5, 4.5)), c(TRUE, FALSE, FALSE))
  expect_identical(HyperCom:::meets(c(0, 2), c(1, 3), numeric(0), numeric(0)), c(FALSE, FALSE))
  expect_identical(HyperCom:::meets_other(c(0, 0.5, 2, 3), c(1, 0.7, 2.5, 3)), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(HyperCom:::meets_other(c(0, 1), c(1, 2)), c(TRUE, TRUE))
})

test_that("a cache from the other grid kernels is not used", {
  factored <- factored_infmat()
  exact <- HyperCom:::cache_fingerprint(fx$output$adj, factored, fx$output$metadata, "exact")
  tiles <- HyperCom:::cache_fingerprint(fx$output$adj, factored, fx$output$metadata, "tiles")
  expect_true(HyperCom:::same_fingerprint(c(exact, token="x"), exact))
  expect_false(HyperCom:::same_fingerprint(c(exact, token="x"), tiles))
})

test_that("memory must be NULL or one positive number and exact.priority TRUE or FALSE", {
  for(memory in list(-1, 0, NA_real_, c(1, 2), "8")){
    expect_error(run_priority(memory=memory), "memory must be NULL or one positive number of gigabytes")
  }
  for(exact.priority in list(NULL, NA, c(TRUE, FALSE), "yes")){
    expect_error(run_priority(exact.priority=exact.priority), "exact.priority must be TRUE or FALSE")
  }
  expect_error(run_hypercom(tempfile(), memory=-1), "memory must be NULL")
  expect_error(run_hypercom(tempfile(), exact.priority="no"), "exact.priority must be TRUE or FALSE")
})
