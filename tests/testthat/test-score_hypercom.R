test_that("a receptor complex scores the mean of its subunits", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata
  a <- score_hypercom(infMat, md, ligand="Tgfb1", receptors="Tgfbr1")
  b <- score_hypercom(infMat, md, ligand="Tgfb1", receptors="Tgfbr2")
  ab <- score_hypercom(infMat, md, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2"))
  b <- b[match(a$cell, b$cell), ]
  ab <- ab[match(a$cell, ab$cell), ]

  expect_equal(ab$sender, (a$sender + b$sender) / 2, tolerance=tol)
  expect_equal(ab$receiver, (a$receiver + b$receiver) / 2, tolerance=tol)
})

test_that("scores are bounded and consistent", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata
  cells <- md$cell[md$sample != "gene"]
  sr <- score_hypercom(infMat, md, ligand="Thbs1", receptors="Sdc4")
  same <- score_hypercom(infMat, md, ligand=fx$output$same.pair[1], receptors=fx$output$same.pair[2])

  expect_setequal(sr$cell, cells)
  expect_true(all(sr$sender >= 0 & sr$sender <= 1 & sr$receiver >= 0 & sr$receiver <= 1, na.rm=TRUE))
  expect_equal(sr$dif, sr$receiver - sr$sender, tolerance=tol)
  expect_equal(sr$avg, (sr$receiver + sr$sender) / 2, tolerance=tol)
  expect_equal(sr$max, pmax(sr$sender, sr$receiver), tolerance=tol)
  expect_true(all(abs(sr$hypercom) <= sr$avg + 1e-12, na.rm=TRUE))
  expect_true(all(same$dif == 0 & same$hypercom == 0))
})

test_that("one sample is scored as handwritten scored its whole grid, also the network grid", {
  infMat <- fx$output$infMat
  md <- fx$output$metadata
  md$sample <- ifelse(md$sample == "gene", "gene", "one")
  cells <- as.character(md$cell[md$sample != "gene"])
  F <- infMat[cells, cells]
  whole <- function(trans){
    diag(trans) <- 0
    trans[trans == 0] <- NA
    R <- matrix(dplyr::percent_rank(as.vector(trans)), nrow(trans))
    return(list(sender=rowMeans(R, na.rm=TRUE), receiver=colMeans(R, na.rm=TRUE)))
  }
  grid <- whole(outer(infMat[cells, "Thbs1"], infMat["Sdc4", cells]) / F)
  network <- whole(F)

  sr <- HyperCom:::score_grids(infMat, md, c("Thbs1", NA), c("Sdc4", NA))
  expect_identical(sr$sender[, 1], unname(grid$sender))
  expect_identical(sr$receiver[, 1], unname(grid$receiver))
  expect_identical(sr$sender[, 2], unname(network$sender))
  expect_identical(sr$receiver[, 2], unname(network$receiver))
})

test_that("score_hypercom scores in this R session: no threads argument, and memory.gb stops a run that cannot fit", {
  md <- fx$output$metadata

  expect_error(score_hypercom(fx$output$infMat, md, "Thbs1", "Sdc4", threads=2), "unused argument")
  expect_error(score_hypercom(fx$output$infMat, md, "Thbs1", "Sdc4", memory.gb=0.5), "scoring needs about 1.07 GB of memory")
  expect_message(score_hypercom(fx$output$infMat, md, "Thbs1", "Sdc4", verbose=TRUE), "Scoring needs at most about 1.07 GB")
})
