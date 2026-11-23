# Input checks of 0.3.0: the lrdb and lrs tables, metadata made by generate_metadata, counts held another way than a matrix, and the
# stops HyperCom makes before it writes a file.

test_that("lrdb blanks, spaces and factors are cleaned, giving the same interactions", {
  lrdb <- fx$input$lrdb
  blank <- lrdb
  blank[is.na(blank)] <- ""
  padded <- lrdb
  padded$ligand <- paste0(" ", padded$ligand, " ")
  factors <- lrdb
  factors[] <- lapply(factors, factor)
  expected <- generate_lrs_table(fx$output$counts, lrdb, outfile=NULL)

  for(x in list(blank, padded, factors)){
    expect_identical(generate_lrs_table(fx$output$counts, x, outfile=NULL), expected)
  }
  for(db in list(cellchat.mm.lrdb, cellchat.hs.lrdb)){
    expect_identical(expect_no_warning(HyperCom:::clean_lrdb(db)), db)
  }
})

test_that("lrdb with other columns, repeated or gene-like interaction names, or no interactions stops naming the problem", {
  lrdb <- fx$input$lrdb
  dup <- lrdb
  dup$interaction[1:2] <- "TGFb"
  gene <- lrdb
  gene$interaction[1] <- lrdb$r1[2]
  bad <- list("columns not expected: pathway"=cbind(lrdb, pathway="x"),
              "columns not expected: interaction"=lrdb[, c(2, 1, 3:6)],
              "CellChat's interaction table"=data.frame(interaction_name="TGFB1_TGFBR1", ligand="Tgfb1", receptor="TGFbR1_R2"),
              "used more than once: TGFb"=dup,
              "named like a gene"=gene,
              "lrdb has no interactions"=lrdb[0, ],
              "must be a data frame"=as.matrix(lrdb))

  for(msg in names(bad)){
    expect_error(filter_counts(fx$input$counts, bad[[msg]], outfile=NULL), msg, fixed=TRUE)
  }
  expect_error(run_priority(lrs=dup), "interaction names of lrs must be distinct")
})

test_that("interactions without a ligand or receptor are left out with a warning, and NA gene names never match", {
  lrdb <- fx$input$lrdb
  lrdb$ligand[1] <- NA
  lrdb[2, paste0("r", 1:4)] <- NA
  counts <- fx$input$counts
  rownames(counts)[1] <- NA

  expect_warning(lrs <- generate_lrs_table(fx$output$counts, lrdb, outfile=NULL), "2 interactions of lrdb have no ligand or no receptor")
  expect_false(any(fx$input$lrdb$interaction[1:2] %in% lrs$interaction))
  expect_false(anyNA(colnames(suppressWarnings(filter_counts(counts, lrdb, outfile=NULL)))))
})

test_that("metadata must come from generate_metadata, and needs its gene rows only for p-values", {
  cells <- fx$output$metadata[fx$output$metadata$sample != "gene", ]

  expect_handwritten(run_priority(metadata=cells), fx$output$priority["1"])
  expect_error(run_priority(metadata=cells, significance=TRUE, seed=1), "with its gene rows")
  expect_error(run_priority(metadata=fx$input$metadata), "metadata must come from generate_metadata")
  expect_error(score_hypercom(fx$output$infMat, fx$input$metadata, "Thbs1", "Sdc4"), "metadata must come from generate_metadata")
})

test_that("HyperCom stops before writing a file on bad labels of expressed cells and on memory that fits no kernel", {
  dir <- file.path(withr::local_tempdir(), "run")
  cell <- rownames(fx$output$counts)[Matrix::rowSums(fx$output$counts) > 0][1]
  na <- fx$input$metadata
  na[cell, "sample"] <- NA
  gene <- fx$input$metadata
  gene[cell, "condition"] <- "gene"

  expect_error(run_hypercom(dir, metadata=na), "1 cells miss a sample")
  expect_error(run_hypercom(dir, metadata=gene), "cannot name a sample or group")
  expect_error(run_hypercom(dir, memory.gb=0.5), "scoring needs about 1.07 GB of memory")
  expect_false(dir.exists(dir))
})

test_that("HyperCom with outdir=NULL writes no file and returns the results", {
  withr::local_dir(withr::local_tempdir())
  res <- run_hypercom(NULL)

  expect_length(list.files(".", recursive=TRUE, all.files=TRUE, include.dirs=TRUE), 0)
  expect_handwritten(res$priority["1"], fx$output$priority["1"])
})

test_that("HyperCom reads counts held as a DelayedArray", {
  skip_if_not_installed("DelayedArray")
  delayed <- DelayedArray::DelayedArray(as.matrix(fx$input$counts))

  expect_handwritten(run_hypercom(NULL, counts=delayed)$priority["1"], fx$output$priority["1"])
})
