# HyperCom 0.3.1.1

The pure R version now follows HyperCom 0.3.1 (branch `master`): the same functions, arguments, checks and fixes, with cell grids
scored in R instead of compiled code. It still installs from CRAN packages only.

* **Much faster than 0.0.1.** Cells are compared only within their sample (0.0.1 built the whole cell x cell grid and dropped the
  pairs across samples), each grid of a ligand and a receptor subunit is scored once for all the pairs and weights that use it, and
  the influence matrix is kept as its factors (`infmat.method="woodbury"`, the new default). The bundled data with p-values take
  about a minute on 8 cores instead of about 8 minutes on 12.
* **Same results.** With `infmat.method="lu"`, scores, priorities and seeded p-values are identical to 0.0.1's (and to those of the
  paper's code, branch `handwritten`). The default `"woodbury"` differs from them by rounding (about 1e-15), which can reorder pairs
  with exactly tied priorities.
* **P-values are exact.** Every real and background pair is scored exactly. The main version screens pairs to save time; in R,
  screening a grid takes longer than scoring it, so this version does not, and has no `exact.priority`.
* **Memory is checked before the run.** Each R process that scores grids needs about 79 bytes per pair of cells in the same sample
  (one sample of 10,000 cells: 7.9 GB). A run that cannot fit stops at the start, giving the estimate.

## Breaking changes

* **`parallel` and `cl` are gone; `threads` gives the number of R processes.** The default (`threads=NULL`) starts one process per
  physical core, at a lower priority than other programs, and fewer when they would not fit in memory (with a message). `threads=1`
  runs in the R session alone. Results do not depend on it. Cluster jobs and containers use only their allocated cores;
  `options(mc.cores=)` caps the default. Clusters across several machines (`cl` on SGE) are no longer supported.
* **`cache` is now `resume`.** A run stopped under 0.0.1 is not resumed (the version changed).
* **R 4.2 or later** and Matrix 1.5 or later are needed.
* **New argument orders**, as in the main version. The leading arguments are unchanged, so calls that give the first arguments by
  position still work; give the others by name:
  * `HyperCom(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=1, significance=FALSE, seed=NULL, r=0.5,
    infmat.method="woodbury", threads=NULL, memory.gb=NULL, resume=TRUE, verbose=TRUE, outdir="HyperCom")`
  * `prioritize_lr(lrs, adj, infMat, metadata, weights=1, significance=FALSE, seed=NULL, threads=NULL, memory.gb=NULL, resume=TRUE,
    verbose=TRUE, outdir="6.priority")`
  * `score_hypercom(infMat, metadata, ligand=NULL, receptors=NULL, memory.gb=NULL, verbose=FALSE)`
  * `generate_infMat(adj, r=0.5, method=NULL, hyperedges=NULL, outfile="4.infMat.rds")`
* **Wrappers pass arguments on.** `HyperCom_seurat(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, ...)`
  and `HyperCom_sce(sce, lrdb, sample, assay="counts", group=NULL, keep=NULL, ...)` pass every other argument on to `HyperCom()`; a
  misspelled or removed argument (`parallel=`, `cl=`, `cache=`, `exact.priority=`) stops with "unused argument".
* **The adjacency and the influence matrix are stored compactly** (`generate_adjacency_matrix()` keeps which cells and genes share a
  hyperedge, `"woodbury"` keeps the factors of the influence matrix). They work like matrices: `dim`, `dimnames`, `x[i, j]` and
  `as.matrix()`; call `library(HyperCom)` before `readRDS()` of `3.adj.rds` or `4.infMat.rds`.
* The bundled databases: `cellchat.mm.lrdb` drops 3 interactions whose receptor has no mouse gene (never used);
  `cellchat.hs.lrdb` fixes CXCL8_CXCR1 (it scored CXCL1) and CD99_PILRB (it had mouse symbols), so their results change.

## New arguments

* `memory.gb` in `HyperCom()`, the wrappers, `prioritize_lr()` and `score_hypercom()`: the gigabytes (10^9 bytes) a run may use. The
  default is 90% of the computer's memory, or on Linux of a container's or cluster job's memory limit when that is smaller. On a
  cluster, give the job's memory request.
* `r` in `HyperCom()` and the wrappers: the restart probability of `generate_infMat()` (default 0.5, as before).
* `infmat.method` in `HyperCom()` and the wrappers: `"woodbury"` (default), `"chol"` or `"lu"`.
* `verbose` in `score_hypercom()`: reports the memory the scoring needs.
* `HyperCom(outdir=NULL)` runs without writing files (and keeps no cache to resume from).

## Fixes and checks

Every fix and check of the main version up to 0.3.1 (see its [NEWS](https://github.com/pawelab/HyperCom/blob/master/NEWS.md)),
among them:

* Cells and genes that no path of the hypergraph connects to the rest are left out by `generate_hyperedges()`, with a warning; the
  interactions of their genes are left out too.
* Without gene rows in `metadata`, priorities counted one group too few.
* `HyperCom()` reads counts held on disk or delayed (DelayedArray, HDF5) for the genes of the database only.
* The interaction table is checked (column names, interaction names, empty entries) and cleaned; interactions without a ligand or a
  receptor are left out with a warning.
* `HyperCom()` stops before it writes a file on cells with a missing sample or group (or one named "gene"), on a run that cannot fit
  in memory, and on `outdir` naming a file.
* `significance`, `resume` and `verbose` must each be `TRUE` or `FALSE`; `threads` one number of at least 1.

# HyperCom 0.0.1

Bug fixes for the R-only version (base-R). On valid input every output is unchanged from 0.0.0.9000 (run under `LC_COLLATE=C`),
except where 0.0.0.9000 was wrong, as listed below.

## New function

* `HyperCom_sce()` runs the pipeline on a SingleCellExperiment, as `HyperCom_seurat()` does on a Seurat object: one assay
  (`assay="counts"` by default) and the cell metadata (`colData`). Assays held on disk (DelayedArray, HDF5Array) are read for the
  genes of `lrdb` only. Needs SingleCellExperiment (Bioconductor).

## New arguments and return values

* `prioritize_lr()`, `HyperCom()` and `HyperCom_seurat()` gain `seed` (reproducible p-values; the session's random state is restored
  afterwards) and `verbose`. They come after `outdir`, so calls by position keep working.
* `HyperCom()` and `HyperCom_seurat()` return their results invisibly (`priority`, `metadata`, `infMat`, `adj`, `lrs`) instead of 0.
* Progress is reported with a few `message()` lines instead of printing every pair's score.

## Parallel runs

* A cluster passed as `cl` is no longer stopped by `prioritize_lr()`.
* With `parallel=TRUE` and no `cl`, the qrsh cluster over `availableWorkers()` is made only inside an SGE parallel job (`PE_HOSTFILE`
  set); elsewhere a local cluster of `availableCores()` workers. A cluster made by `prioritize_lr()` is stopped also after an error.

## Results that were wrong

* Ligand-receptor columns stored as factors indexed the adjacency by factor codes (wrong priorities); an lrdb read with
  `stringsAsFactors=TRUE` lost its interactions. Genes are now always read by name.
* Metadata rows in another order than the influence matrix gave wrong scores: cells are now matched by name.
* Node order followed the session's collation (and factor codes for factor hyperedge columns), so the same data gave different
  node orders, and seeded p-values, on different machines. Nodes are now ordered by name in the C locale.
* A pair whose cells are perfectly balanced over the groups got priority NaN (purity rounded below 0); the purity is now at least 0.
* Background pairs without a priority (NaN) made every p-value NA; they are now left out, with a warning.
* `keep` with one column named it `additional`; a factor `keep` column got NA on gene rows; date columns stopped. Kept columns
  keep their names (renamed `sample.1` etc. when they collide with `cell`, `sample` or `group`), factors get a `"gene"` level, dates
  become strings.
* The lrdb's first column no longer has to be named `interaction` (the documented `edges` format errored).

## Input that now stops or warns instead of giving nonsense

* Stops: an NA sample or group of a cell in the graph (it became a gene); a sample or group named `"gene"`; a sample with fewer than
  2 cells in the graph; two cells with one name; metadata whose row names are not the cell names; counts with negative values
  (scaled data); transposed counts; `lrdb` and `counts` swapped; misspelt metadata columns; no interaction of lrdb expressed (e.g.
  lrdb of another species); weights that are negative, not finite or print alike; `score_hypercom()` with only a ligand or only
  receptors, or with a gene that is not in the graph.
* Warns: cells of counts missing from the metadata (left out; they used to become background genes); gene names used by several
  rows of counts (the first is kept); priority files of other weights already in `outdir`.
* Missing (NA) counts are ignored instead of stopping `filter_counts()`.

## Other

* `HyperCom()` no longer changes the working directory; `outdir` may be a nested path.
* `library(HyperCom)` loads Matrix, so the bundled counts can be subset at once.
* SeuratObject 4 assays are read from the layer asked for (`GetAssayData(layer=)` ignored it there).
* mltools (archived on CRAN) is no longer needed: the Gini impurity is computed by HyperCom itself, with the same numbers.
* A test suite (testthat) compares every step with the outputs of 0.0.0.9000 and has one test per fixed bug.
