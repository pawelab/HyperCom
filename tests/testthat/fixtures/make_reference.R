# Golden outputs of HyperCom master (b999a01) on the bundled data, for test-full-reference.R (about 23 MB, not in git)
# Run with R_LIBS=<library holding master> so PSOCK workers load master too (about 8 min on 12 workers):
#   Rscript tests/testthat/fixtures/make_reference.R <outdir> <workers>
library(HyperCom)

args <- commandArgs(trailingOnly=TRUE)
outdir <- if(length(args) > 0) args[1] else "reference/master_b999a01"
workers <- if(length(args) > 1) as.integer(args[2]) else 20

stopifnot(packageVersion("HyperCom") == "0.0.0.9000")
Sys.setlocale("LC_COLLATE", "C")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
timings <- data.frame(step=character(), seconds=numeric())
time_step <- function(step, expr){
  t0 <- Sys.time()
  x <- expr
  timings[nrow(timings) + 1, ] <<- list(step, as.numeric(difftime(Sys.time(), t0, units="secs")))
  print(paste0(step, ",", round(timings$seconds[nrow(timings)], 2)))
  return(x)
}

counts <- time_step("filter_counts", filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=file.path(outdir, "0.counts.rds")))
lrs <- time_step("generate_lrs_table", generate_lrs_table(counts, cellchat.mm.lrdb, outfile=file.path(outdir, "1.lrs.csv")))
hyperedges <- time_step("generate_hyperedges", generate_hyperedges(counts, lrs, outfile=file.path(outdir, "2.hyperedges.csv")))
adj <- time_step("generate_adjacency_matrix", generate_adjacency_matrix(hyperedges, outfile=file.path(outdir, "3.adj.rds")))
infMat <- time_step("generate_infMat", generate_infMat(adj, outfile=file.path(outdir, "4.infMat.rds")))
metadata <- time_step("generate_metadata", generate_metadata(infMat, alexanian2021.metadata, sample="sample", group="condition",
                                                             outfile=file.path(outdir, "5.metadata.csv")))

scores <- list(network=score_hypercom(infMat, metadata),
               Thbs1_Sdc4=score_hypercom(infMat, metadata, ligand="Thbs1", receptors="Sdc4"),
               Tgfb1_Tgfbr1_Tgfbr2=score_hypercom(infMat, metadata, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2")))
saveRDS(scores, file.path(outdir, "scores.rds"))

cl <- parallelly::makeClusterPSOCK(workers, autoStop=TRUE)
worker.versions <- unique(unlist(parallel::clusterEvalQ(cl, as.character(utils::packageVersion("HyperCom")))))
stopifnot(identical(worker.versions, "0.0.0.9000"))
set.seed(1)
priority <- time_step("prioritize_lr", prioritize_lr(lrs, adj, infMat, metadata, weights=c(1, 0.5), significance=TRUE, parallel=TRUE,
                                                     cl=cl, outdir=file.path(outdir, "6.priority")))
saveRDS(priority, file.path(outdir, "6.priority.rds"))

utils::write.csv(timings, file.path(outdir, "timings.csv"), row.names=FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
print(head(priority[["1"]]))
