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

test_that("score_hypercom gives the same scores on several threads, for a dense and for a factored influence matrix", {
  adj <- generate_adjacency_matrix(fx$output$hyperedges, outfile=NULL)
  infMat <- generate_infMat(adj, outfile=NULL, threads=test_threads(2))
  md <- fx$output$metadata

  for(p in fx_pairs){
    expect_identical(score_hypercom(fx$output$infMat, md, p[[1]], p[[2]], threads=test_threads(3)),
                     score_hypercom(fx$output$infMat, md, p[[1]], p[[2]]))
    expect_identical(score_hypercom(infMat, md, p[[1]], p[[2]], parallel=TRUE), score_hypercom(infMat, md, p[[1]], p[[2]]))
  }
  expect_error(score_hypercom(infMat, md, "Thbs1", "Sdc4", threads=0), "at least 1")
})
