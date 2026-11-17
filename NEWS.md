# HyperCom 0.0.1

Bug fixes for the R-only version (base-R). On valid input every output is unchanged from 0.0.0.9000 (run under `LC_COLLATE=C`),
except where 0.0.0.9000 was wrong, as listed below.

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
* A test suite (testthat) compares every step with the outputs of 0.0.0.9000 and has one test per fixed bug.
