test_that("generate_metadata matches handwritten, except that kept numeric and logical columns keep their type (NA for genes)", {
  infMat <- fx$output$infMat
  md <- fx$input$metadata
  md$flag <- md$UMAP_2 > 0
  keep <- fx$output$metadata.keep
  is.gene <- keep$sample == "gene"
  keep$UMAP_1 <- ifelse(is.gene, NA, suppressWarnings(as.numeric(keep$UMAP_1)))
  keep$flag <- md$flag[match(keep$cell, rownames(md))]
  kept <- generate_metadata(infMat, md, sample="sample", group="condition", keep=c("celltype", "UMAP_1", "flag"), outfile=NULL)

  expect_equal(plain_df(generate_metadata(infMat, md, sample="sample", group="condition", outfile=NULL)),
               plain_df(fx$output$metadata))
  expect_equal(plain_df(generate_metadata(infMat, md, sample="sample", outfile=NULL)),
               plain_df(fx$output$metadata.none))
  expect_equal(plain_df(kept), plain_df(keep))
  expect_type(kept$flag, "logical")
})

test_that("metadata rows follow the influence matrix and label gene nodes", {
  infMat <- fx$output$infMat
  md <- fx$input$metadata
  clean.df <- generate_metadata(infMat, md, sample="sample", group="condition", outfile=NULL)
  none.df <- generate_metadata(infMat, md, sample="sample", outfile=NULL)
  is.gene <- !clean.df$cell %in% rownames(md)
  cells <- clean.df$cell[!is.gene]

  expect_identical(colnames(clean.df), c("cell", "sample", "group"))
  expect_identical(clean.df$cell, rownames(infMat))
  expect_true(all(clean.df$sample[is.gene] == "gene" & clean.df$group[is.gene] == "gene"))
  expect_identical(as.character(clean.df$sample[!is.gene]), as.character(md[cells, "sample"]))
  expect_identical(as.character(clean.df$group[!is.gene]), md[cells, "condition"])
  expect_true(all(none.df$group[!is.gene] == "none"))
})

test_that("metadata rows outside the graph may have missing values", {
  md <- fx$input$metadata
  extra <- md[1, ]
  extra$sample <- NA
  extra$condition <- NA
  rownames(extra) <- "NOT_IN_GRAPH"
  clean.df <- generate_metadata(fx$output$infMat, rbind(md, extra), sample="sample", group="condition", outfile=NULL)

  expect_equal(plain_df(clean.df), plain_df(fx$output$metadata))
})
