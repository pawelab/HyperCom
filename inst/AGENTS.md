# AGENTS.md: HyperCom

This file is for AI agents that help a user run HyperCom. The installed copy is at `system.file("AGENTS.md", package="HyperCom")`.

## First: check the version

Run `packageVersion("HyperCom")`. Several versions share the package name but have different arguments:

- **0.3.1.1** (branch `base-R`): this file. HyperCom 0.3.1 in R only (no compiled code): R processes, resume after a stop, no memory
  mode, no `exact.priority`, no `threads` in `generate_infMat()` or `score_hypercom()`. For up to about 10,000 cells per sample.
- **0.3.0, 0.3.1** (branch `master`): compiled code, threads, memory mode; up to about 100,000 cells per sample. Same results.
- **0.2.x:** has `parallel=` and `cache=`, and uses one thread by default; no `memory.gb`, `exact.priority` or `r` in `HyperCom()`.
- **0.1.x** (branch `spatial`, formerly `v1-fast`): has `cores=`, `infMat.cell=`, `cell.info=`, `method="auto"` and spatial scores.
  0.3.x has none of these.
- **0.0.1** (branch `base-R` before 0.3.1.1): R code only, with the same bug fixes. Has `parallel=`, `cl=`, `cache=`, `seed=` and
  `verbose=`; no `threads`, `resume` or `infmat.method`.
- **0.0.0.9000** (branch `handwritten`): the original code of the paper, with bugs fixed in 0.0.1. Needs mltools (archived on
  CRAN).

For this version, read `?HyperCom`, `?prioritize_lr`, `news(package="HyperCom")` and the vignette (the PDF linked from the README
of branch `base-R`, or `vignettes/articles/` in the source repository).

## Run

```r
library(HyperCom)
res <- HyperCom(cellchat.mm.lrdb, counts, metadata, sample="sample", group="condition",
                significance=TRUE, seed=1, outdir="HyperCom")
res$priority[["1"]]
s <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")
```

- **Seurat:** `HyperCom_seurat(seurat, lrdb, sample=, assay="RNA", layer="counts", ...)`. Needs SeuratObject.
- **SingleCellExperiment:** `HyperCom_sce(sce, lrdb, sample=, assay="counts", ...)`. Columns come from `colData`, cell names from
  `colnames(sce)`. Needs SummarizedExperiment (Bioconductor; installed with SingleCellExperiment).
- The wrappers pass every other argument on to `HyperCom()`. A misspelled or removed argument (`parallel=`, `cache=`, `cl=`,
  `exact.priority=`) stops with "unused argument".
- **Install check:** `example(HyperCom)` takes seconds.
- **Processes:** `threads` is the number of R processes that score. The default starts one per physical core, at most 4 (more gain
  little), fewer when they would not fit in `memory.gb` (with a message). `threads=1` scores in the R session alone. Results do not depend on it.
- **Output folder:** set `outdir` yourself; the default is `HyperCom/` in the working directory. `outdir=NULL` writes no files.

## Check the input before a run

Most failed runs come from the input.

- `counts` is genes x cells: gene symbols as row names, cell names as column names. Sparse matrices are fine.
- Raw or normalized counts are fine (HyperCom only uses whether a value is above 0). Scaled data stops the run (negative values).
- Cell names must be unique. When samples are merged, prefix each barcode with its sample.
- `metadata` has one row per cell, with the cell names as row names. Cells of counts that are not in metadata are left out, with a
  warning.
- `sample` is the library or replicate column, not the cell type. Cells are compared only within their sample, and each sample needs
  at least 2 cells that express a ligand or receptor gene. One sample? Add `metadata$sample <- "all"`; `sample` is required.
- `group` (optional) is the condition column. It has no NA in cells that express a ligand or receptor gene, and no value "gene". It
  feeds the purity term of the priority (pairs whose cells concentrate in some groups rank higher); there is no ranking per group.
- `lrdb` must match the species of the counts: `cellchat.mm.lrdb` for mouse (`Tgfb1`), `cellchat.hs.lrdb` for human (`TGFB1`). The
  wrong species stops the run.
- A custom `lrdb`: column 1 is the interaction name, column 2 is `ligand`, the rest are `r1`, `r2`, ... (receptor subunits, NA when
  a receptor has fewer).
- A warning "... cells and ... genes are not connected to the rest of the hypergraph and are left out" is normal when the numbers
  are small. These cells get no scores. Many such cells point to the wrong species or gene names.

## Estimate size before a run

- Memory of one R process: about 79 bytes x (sum over samples of the squared number of cells), plus about 1.1 GB. In R:
  `sum(as.numeric(table(metadata$sample))^2) * 79 / 1e9 + 1.1` GB (an upper bound). Each further process needs about as much
  again (0.4 GB instead of 1.1 GB on top).
- HyperCom starts only as many processes as fit in `memory.gb` (default: 90% of the computer's memory, or of a container's or
  cluster job's limit on Linux). If one process does not fit, the run stops at the start with the estimate. On a cluster, pass the
  job's request, e.g. `memory.gb=64`.
- More than about 10,000 cells in one sample: suggest the compiled version (branch `master`, same results, about 4 times less
  memory, memory mode up to about 100,000 cells in 64 GB, much faster).
- More samples of the same total size are cheap. Do not split a sample only to save memory without asking the user: it changes which
  cells are compared, and so the results.
- Measured, one sample (Intel Core i9-12900K), the bundled mouse database (about 1,000 pairs expressed), with
  p-values, 4 processes (the default): 1,000 cells 1.4 min and 1.8 GB; 3,000 cells 20 min and 4.0 GB; 10,000 cells 3.1 h and 32 GB.
  Time also grows with the number of pairs. Past about 4 processes runs gain little (memory-bound).
- With `significance=TRUE`, 1,000 background pairs are scored too. Every pair is scored exactly.

## Long runs and clusters

- If a run stops (time limit, crash), run the same call again. It continues from `<outdir>/6.priority/cache` and loses at most one
  batch of pairs.
- The cache is not reused, but moved to `cache.old` with a warning, when the new call differs in counts, lrdb, the cells' sample or
  group, `r`, `infmat.method` or the HyperCom version. Keep these the same, and keep `seed`. `weights` and `threads` may change.
- Use the same node type for the same last digits: another BLAS can change them by about 1e-13.
- Do not delete the cache during a run. HyperCom removes it at the end.

## Output

- `HyperCom()` returns its list invisibly, so assign it: `priority` (a list named by weight, e.g. `"1"`), `metadata`, `infMat`,
  `adj` and `lrs`.
- Files: `0.counts.rds`, `1.lrs.csv`, `2.hyperedges.csv`, `3.adj.rds`, `4.infMat.rds`, `5.metadata.csv`,
  `6.priority/priority.<weight>.csv`.
- Call `library(HyperCom)` before `readRDS` of `3.adj.rds` or `4.infMat.rds`.
- `adj` and `infMat` are stored compactly, not as matrices. `dim`, `dimnames` and `x[i, j]` work. `as.matrix()` builds the full
  nodes x nodes matrix, which is huge on large data: take slices instead.
- Priority table: the columns of `lrs` (those of `lrdb`, without receptor columns that are empty for every pair) plus `priority`,
  best first. With significance it also has `p` and `p.adj` (Benjamini-Hochberg).
- `p` moves in steps of 0.001, and `p = 0` means below 0.001: report it as `p < 0.001`, not 0. Its `p.adj` is then 0 too and is
  also only a bound.
- Score table: `cell`, `sender` and `receiver` (0 to 1), `dif` = receiver - sender, `avg`, `max`, and `hypercom` (negative for
  senders, positive for receivers).
- A NaN score means the cell has no score for that pair (rare). Cells that express no ligand or receptor gene are not in the graph
  and have no row.
- For a receptor complex, give all subunits: `receptors=c("Tgfbr1", "Tgfbr2")`. Their scores are averaged.
- `weights` is the weight of the cell-fraction term. Values below 1 favor rare pairs; a vector gives one table per weight.
- A few interactions of the bundled databases have the same genes as another under another name (e.g. ADM_CALCR_RAMP1 and
  ADM_CALCRL_RAMP1): they get the same priority.

## Reproduce

- `seed=` makes p-values reproducible. The session's random state is restored afterwards.
- `infmat.method="woodbury"` (the default) is fast. `"lu"` reproduces the numbers of the original release (0.0.0.9000, branch
  `handwritten`) bit for bit, scores and seeded p-values included, but holds the full matrix; the two differ by about 1e-15.

## Change the package (source repository)

- Outputs must match the original code (branch `handwritten`; fixture `tests/testthat/fixtures/handwritten.rds`) bit for bit from
  the same dense influence matrix. Never loosen a check to make a test pass. R code only: branch `master` has the compiled version;
  `R/grids.R` stands in for its `src/`, so the other files stay close to master's.
- Tests use testthat. `NOT_CRAN=true` runs all of them. `HYPERCOM_REFERENCE=<golden output folder>` adds the full-data test, which is
  skipped otherwise; `HYPERCOM_TEST_WORKERS=<n>` sets its processes (default 4).
- Documentation is roxygen in `R/`. Regenerate `man/` and `NAMESPACE`; never edit them by hand.
- `R CMD check --as-cran` must end with Status OK.
- Add every user-visible change to `NEWS.md`.
