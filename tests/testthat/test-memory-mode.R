# Memory mode (log-domain tiles, passes) and screened real pairs: both depart from exact scoring by known amounts (memory mode ~1e-11,
# screened real pairs ~1e-9 typical, ~1e-7 at most on the bundled data), so they are checked against the exact kernels within those, and against themselves
# exactly. Also the memory budget and the checks memory mode makes before it starts.

factored_infmat <- function(){
  return(fx_factored()$infMat)
}

test_that("the memory budget is the smaller of the computer's memory and the control group's limit", {
  root <- withr::local_tempdir()
  proc <- file.path(root, "cgroup")
  put <- function(path, value){
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    writeLines(value, path)
  }
  writeLines("0::/job/step", proc)
  put(file.path(root, "job", "memory.max"), "8000000000")
  put(file.path(root, "job", "step", "memory.max"), "max")
  expect_identical(HyperCom:::cgroup_memory(proc, root), 8e9)
  writeLines(c("4:memory:/slurm/job", "3:cpu,cpuacct:/x"), proc)
  put(file.path(root, "memory", "slurm", "job", "memory.limit_in_bytes"), "9223372036854771712")
  put(file.path(root, "memory", "slurm", "memory.limit_in_bytes"), "4000000000")
  expect_identical(HyperCom:::cgroup_memory(proc, root), 4e9)
  expect_identical(HyperCom:::cgroup_memory(file.path(root, "none"), file.path(root, "none")), Inf)

  testthat::local_mocked_bindings(physical_memory_cpp=function() 16e9, cgroup_memory=function(...) 8e9, .package="HyperCom")
  expect_identical(HyperCom:::memory_budget(), 0.9 * 8e9)
  expect_identical(HyperCom:::memory_budget(2), 2e9)
})

test_that("memory mode stops at the start on a hypergraph that is not connected", {
  he <- data.frame(edge=c("IA", "IA", "A1", "A1", "A1", "A2", "A2", "A2", "IB", "IB", "B1", "B1", "B2", "B2"),
                   node=c("A1", "A2", "A1", "c1", "c2", "A2", "c2", "c3", "B1", "B2", "B1", "c4", "B2", "c5"))
  infMat <- generate_infMat(generate_adjacency_matrix(he, outfile=NULL), outfile=NULL, threads=1)
  md <- generate_metadata(infMat, data.frame(sample=c(1, 1, 1, 2, 2), row.names=paste0("c", 1:5)), sample="sample", outfile=NULL)

  expect_s3_class(infMat, "hypercom_influence")
  expect_error(score_hypercom(infMat, md, "A1", "A2", memory.gb=memory_for(md)), "memory mode needs a connected hypergraph")
  expect_no_error(score_hypercom(infMat, md, "A1", "A2", memory.gb=Inf))
})

test_that("score_hypercom scores in memory mode when the exact kernels do not fit", {
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  for(p in fx_pairs[c("single", "complex")]){
    exact <- score_hypercom(factored, metadata, p[[1]], p[[2]], memory.gb=Inf)
    tiles <- score_hypercom(factored, metadata, p[[1]], p[[2]], memory.gb=memory_for(metadata))
    tiles <- tiles[match(exact$cell, tiles$cell), ]

    expect_equal(tiles$sender, exact$sender, tolerance=1e-9)
    expect_equal(tiles$receiver, exact$receiver, tolerance=1e-9)
  }
})

test_that("the kernel is exact when it fits (also in the computer's memory), memory mode otherwise, else a stop with the estimate", {
  factored <- factored_infmat()
  metadata <- fx$output$metadata
  expect_identical(HyperCom:::choose_kernel(factored, metadata, Inf)$type, "exact")
  expect_identical(HyperCom:::choose_kernel(factored, metadata, NULL)$type, "exact")
  kernel <- HyperCom:::choose_kernel(factored, metadata, memory_for(metadata))
  expect_identical(kernel$type, "memory")
  expect_error(HyperCom:::choose_kernel(factored, metadata, 0.5), "scoring needs about 1.07 GB of memory")
  expect_error(HyperCom:::choose_kernel(fx$output$infMat, metadata, memory_for(metadata)), "memory mode needs an influence matrix")
})

test_that("memory mode scores samples of any size, also not multiples of 8 and of 2 cells, within rounding of exact", {
  factored <- factored_infmat()
  md <- fx$output$metadata
  md$sample <- as.character(md$sample)
  cells <- which(md$sample != "gene")
  md$sample[cells] <- c(rep("two", 2), rep("odd", 13), rep("rest", length(cells) - 15))
  g <- HyperCom:::grid_cells(rownames(factored), md)
  grids <- HyperCom:::pair_grids(as.character(fx$output$lrs$ligand), as.matrix(fx$output$lrs[, -(1:2)]))
  grids <- grids[grids$ligand != grids$receptor, ]
  lb <- HyperCom:::lambda_blocks(factored, g, 1)
  exact <- HyperCom:::score_grids(factored, md, grids$ligand, grids$receptor, 1)
  tiles <- HyperCom:::score_grids(factored, md, grids$ligand, grids$receptor, test_threads(2), lb)
  sb <- HyperCom:::screen_grids(factored, md, grids$ligand, grids$receptor, test_threads(2), lb)
  ok <- !is.na(tiles$sender)

  expect_equal(tiles, exact, tolerance=1e-9)
  expect_true(all(sb$sender.lo[ok] <= tiles$sender[ok] & tiles$sender[ok] <= sb$sender.hi[ok]))
  expect_true(all(sb$receiver.lo[ok] <= tiles$receiver[ok] & tiles$receiver[ok] <= sb$receiver.hi[ok]))
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

test_that("HyperCom in memory mode, real pairs screened, gives the p-values and order of exact scoring", {
  skip_on_cran()
  exact <- run_hypercom(NULL, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), memory.gb=Inf,
                        exact.priority=TRUE)$priority
  tiles <- run_hypercom(NULL, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2),
                        memory.gb=memory_for(fx$output$metadata))$priority

  for(w in names(exact)){
    expect_identical(tiles[[w]]$interaction, exact[[w]]$interaction)
    expect_equal(tiles[[w]]$priority, exact[[w]]$priority, tolerance=1e-6)
    expect_identical(tiles[[w]][c("p", "p.adj")], exact[[w]][c("p", "p.adj")])
  }
})

test_that("screened real pairs give the p-values and order of exact scoring, priorities within 1e-6 relative of exact", {
  skip_on_cran()
  factored <- factored_infmat()
  for(infMat in list(fx$output$infMat, factored)){
    exact <- run_priority(infMat=infMat, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), exact.priority=TRUE)
    screened <- run_priority(infMat=infMat, weights=c(1, 0.5), significance=TRUE, seed=1, threads=test_threads(2), exact.priority=FALSE)
    for(w in names(exact)){
      expect_identical(screened[[w]]$interaction, exact[[w]]$interaction)
      expect_identical(screened[[w]][c("p", "p.adj")], exact[[w]][c("p", "p.adj")])
      expect_lte(max(abs(screened[[w]]$priority - exact[[w]]$priority) / abs(exact[[w]]$priority), na.rm=TRUE), 1e-6)
    }
  }
  expect_identical(run_priority(weights=c(1, 0.5), exact.priority=TRUE), run_priority(weights=c(1, 0.5)))
})

test_that("intervals that meet are found", {
  expect_identical(HyperCom:::meets(c(0, 2, 5), c(1, 3, 6), c(1, 4), c(1.5, 4.5)), c(TRUE, FALSE, FALSE))
  expect_identical(HyperCom:::meets(c(0, 2), c(1, 3), numeric(0), numeric(0)), c(FALSE, FALSE))
  expect_identical(HyperCom:::meets_other(c(0, 0.5, 2, 3), c(1, 0.7, 2.5, 3)), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(HyperCom:::meets_other(c(0, 1), c(1, 2)), c(TRUE, TRUE))
})

test_that("memory.gb must be NULL or one positive number, and the switches TRUE or FALSE, in prioritize_lr and HyperCom", {
  for(memory.gb in list(-1, 0, NA_real_, c(1, 2), "8")){
    expect_error(run_priority(memory.gb=memory.gb), "memory.gb must be NULL or one positive number of gigabytes")
  }
  for(flag in c("significance", "exact.priority", "resume", "verbose")){
    for(value in list(NULL, NA, c(TRUE, FALSE), "yes")){
      expect_error(do.call(run_priority, stats::setNames(list(value), flag)), paste(flag, "must be TRUE or FALSE"))
    }
  }
  dir <- file.path(withr::local_tempdir(), "run")
  expect_error(run_hypercom(dir, memory.gb=-1), "memory.gb must be NULL")
  expect_error(run_hypercom(dir, resume="no"), "resume must be TRUE or FALSE")
  expect_false(dir.exists(dir))
})
