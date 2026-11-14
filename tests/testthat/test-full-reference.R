# Every step on the full bundled data against master's golden outputs (fixtures/make_reference.R). Skipped unless HYPERCOM_REFERENCE
# names the golden folder. Each step gets master's own input, so a difference points at the step that made it. Several minutes.

ref_path <- function(name){
  return(file.path(Sys.getenv("HYPERCOM_REFERENCE"), name))
}

skip_without_reference <- function(){
  testthat::skip_if(Sys.getenv("HYPERCOM_REFERENCE") == "", "set HYPERCOM_REFERENCE to a golden output folder")
  testthat::skip_if_not(file.exists(ref_path("6.priority.rds")), "golden output folder is incomplete")
}

test_that("full data: preprocessing and hypergraph match master", {
  skip_without_reference()
  counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
  lrs <- generate_lrs_table(readRDS(ref_path("0.counts.rds")), cellchat.mm.lrdb, outfile=NULL)
  ref.lrs <- utils::read.csv(ref_path("1.lrs.csv"), check.names=FALSE)
  hyperedges <- generate_hyperedges(readRDS(ref_path("0.counts.rds")), ref.lrs, outfile=NULL)
  ref.he <- utils::read.csv(ref_path("2.hyperedges.csv"), check.names=FALSE)
  adj <- generate_adjacency_matrix(ref.he, outfile=NULL)

  expect_same_matrix(counts, readRDS(ref_path("0.counts.rds")))
  expect_equal(plain_df(lrs), plain_df(ref.lrs))
  expect_equal(plain_df(hyperedges), plain_df(ref.he))
  expect_same_matrix(adj, readRDS(ref_path("3.adj.rds")))
})

test_that("full data: influence matrix and metadata match master", {
  skip_without_reference()
  infMat <- generate_infMat(readRDS(ref_path("3.adj.rds")), outfile=NULL)
  ref.inf <- readRDS(ref_path("4.infMat.rds"))
  md <- generate_metadata(ref.inf, alexanian2021.metadata, sample="sample", group="condition", outfile=NULL)
  ref.md <- utils::read.csv(ref_path("5.metadata.csv"), check.names=FALSE)

  expect_same_matrix(infMat, ref.inf, tolerance=tol)
  expect_equal(plain_df(data.frame(lapply(md, as.character))), plain_df(data.frame(lapply(ref.md, as.character))))
})

test_that("full data: scores and priorities match master", {
  skip_without_reference()
  skip_on_cran()
  ref.inf <- readRDS(ref_path("4.infMat.rds"))
  adj <- readRDS(ref_path("3.adj.rds"))
  lrs <- generate_lrs_table(readRDS(ref_path("0.counts.rds")), cellchat.mm.lrdb, outfile=NULL)
  md <- generate_metadata(ref.inf, alexanian2021.metadata, sample="sample", group="condition", outfile=NULL)
  ref.scores <- readRDS(ref_path("scores.rds"))
  scores <- list(network=score_hypercom(ref.inf, md),
                 Thbs1_Sdc4=score_hypercom(ref.inf, md, ligand="Thbs1", receptors="Sdc4"),
                 Tgfb1_Tgfbr1_Tgfbr2=score_hypercom(ref.inf, md, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2")))
  for(n in names(ref.scores)){
    expect_same_scores(scores[[n]], ref.scores[[n]])
  }

  workers <- as.integer(Sys.getenv("HYPERCOM_TEST_WORKERS", "4"))
  set.seed(1)
  if("threads" %in% names(formals(prioritize_lr))){
    priority <- quietly(prioritize_lr(lrs, adj, ref.inf, md, weights=c(1, 0.5), significance=TRUE, threads=workers, outdir=NULL))
  } else if("cores" %in% names(formals(prioritize_lr))){
    priority <- quietly(prioritize_lr(lrs, adj, ref.inf, md, weights=c(1, 0.5), significance=TRUE, parallel=TRUE, cores=workers,
                                      outdir=NULL))
  } else {
    cl <- parallel::makeCluster(workers)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
    priority <- quietly(prioritize_lr(lrs, adj, ref.inf, md, weights=c(1, 0.5), significance=TRUE, parallel=TRUE, cl=cl, outdir=NULL))
  }
  ref.priority <- readRDS(ref_path("6.priority.rds"))
  for(w in names(ref.priority)){
    expect_same_priority(priority[[w]], ref.priority[[w]])
  }
})
