# Build tests/testthat/fixtures/handwritten.rds: a small input set and the handwritten branch's output for every step
# Handwritten (b999a01) is the truth for outputs. Run from the package root with branch handwritten installed:
#   R_LIBS=<library holding handwritten> Rscript tests/testthat/fixtures/make_fixtures.R
library(Matrix)
library(HyperCom)

stopifnot(packageVersion("HyperCom") == "0.0.0.9000")
Sys.setlocale("LC_COLLATE", "C")

set.seed(1)
md <- alexanian2021.metadata
cells <- c(sample(which(md$sample == 2), 80), sample(which(md$sample == 3), 80))
counts.all <- alexanian2021.counts[, cells]
lrdb.all <- cellchat.mm.lrdb

lrs.all <- generate_lrs_table(filter_counts(counts.all, lrdb.all, outfile=NULL), lrdb.all, outfile=NULL)
num.receptors <- rowSums(!is.na(lrs.all[, 3:ncol(lrs.all)]))
same.gene <- lrs.all$ligand == lrs.all$r1 & num.receptors == 1
picked <- c("THBS1_SDC4", "TGFB1_TGFBR1_TGFBR2", "TGFB1_ACVR1_TGFBR1", "SEMA3C_NRP1_NRP2", "MMP2_INTEGRIN_AVB3_COMPLEX")
single <- setdiff(lrs.all$interaction[num.receptors == 1 & !same.gene], picked)
double <- setdiff(lrs.all$interaction[num.receptors == 2], picked)
picked <- c(picked, sample(single, 15), sample(double, 10), sample(lrs.all$interaction[same.gene], 2))
absent <- sample(setdiff(lrdb.all$interaction, lrs.all$interaction), 4)

lrdb <- lrdb.all[lrdb.all$interaction %in% c(picked, absent), ]
stopifnot(!anyDuplicated(lrdb$interaction))
lr.genes <- intersect(unique(stats::na.omit(unlist(lrdb[, -1]))), rownames(counts.all))
other.genes <- sample(setdiff(rownames(counts.all), unique(stats::na.omit(unlist(lrdb.all[, -1])))), 40)
counts <- counts.all[sort(c(lr.genes, other.genes)), ]
metadata <- md[cells, ]

input <- list(counts=counts, metadata=metadata, lrdb=lrdb)
output <- list()
output$counts <- filter_counts(counts, lrdb, outfile=NULL)
output$lrs <- generate_lrs_table(output$counts, lrdb, outfile=NULL)
output$hyperedges <- generate_hyperedges(output$counts, output$lrs, outfile=NULL)
output$adj <- generate_adjacency_matrix(output$hyperedges, outfile=NULL)
output$infMat <- generate_infMat(output$adj, outfile=NULL)
output$infMat.r03 <- generate_infMat(output$adj, r=0.3, outfile=NULL)
output$metadata <- generate_metadata(output$infMat, metadata, sample="sample", group="condition", outfile=NULL)
output$metadata.none <- generate_metadata(output$infMat, metadata, sample="sample", outfile=NULL)
output$metadata.keep <- generate_metadata(output$infMat, metadata, sample="sample", group="condition", keep=c("celltype", "UMAP_1"),
                                          outfile=NULL)

same.pair <- output$lrs[output$lrs$ligand == output$lrs$r1 & is.na(output$lrs$r2), ][1, ]
output$scores <- list(network=score_hypercom(output$infMat, output$metadata),
                      single=score_hypercom(output$infMat, output$metadata, ligand="Thbs1", receptors="Sdc4"),
                      complex=score_hypercom(output$infMat, output$metadata, ligand="Tgfb1", receptors=c("Tgfbr1", "Tgfbr2")),
                      triple=score_hypercom(output$infMat, output$metadata, ligand="Tgfb1", receptors=c("Acvr1", "Tgfbr1", "Tgfbr2")),
                      same=score_hypercom(output$infMat, output$metadata, ligand=same.pair$ligand, receptors=same.pair$r1))
output$same.pair <- c(same.pair$ligand, same.pair$r1)

output$priority <- prioritize_lr(output$lrs, output$adj, output$infMat, output$metadata, weights=c(1, 0.5, 0), outdir=NULL)
output$priority.none <- prioritize_lr(output$lrs, output$adj, output$infMat, output$metadata.none, outdir=NULL)
set.seed(1)
output$priority.sig <- prioritize_lr(output$lrs, output$adj, output$infMat, output$metadata, significance=TRUE, outdir=NULL)

about <- list(commit="b999a01", version=as.character(packageVersion("HyperCom")), made=format(Sys.time()), R=R.version.string,
              platform=R.version$platform, collate=Sys.getlocale("LC_COLLATE"))
saveRDS(list(about=about, input=input, output=output), "tests/testthat/fixtures/handwritten.rds", compress="xz")
print(paste0("pairs ", nrow(output$lrs), ", nodes ", nrow(output$infMat), ", cells ", ncol(counts),
             ", NaN priorities ", sum(is.na(output$priority[["1"]]$priority))))
