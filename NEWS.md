# HyperCom 0.3.2

* **`sample` is optional**, like `group`: without it all cells are one sample (`"all"`), in `HyperCom()`, `HyperCom_seurat()`,
  `HyperCom_sce()` and `generate_metadata()`. `HyperCom()` also needs no `metadata` when there is neither sample, group nor
  `keep`: `HyperCom(cellchat.hs.lrdb, counts)` runs on every cell of counts. Results equal those of `metadata$sample <- "all"`.
* **`prioritize=FALSE`** in `HyperCom()` and the wrappers skips ranking the pairs (`prioritize_lr`), which takes most of the run
  time: the run stops after the influence matrix and metadata, for `score_hypercom()` on chosen pairs. `res$priority` is then
  `NULL` and `6.priority` is not written.
* **Faster exact scoring of very large samples.** The biggest value buckets (one sample of more than ~80,000 cells) are now cut
  into in-cache pieces before sorting instead of taking a slower general sort: grid scoring takes ~18% less time at the bucket size
  of 100,000 cells in one sample and ~40% less at 200,000-400,000 cells (measured on 50,000 cells with fewer buckets). Results are
  identical to the last bit. Memory mode and screening are unchanged.

# HyperCom 0.3.1

* **`cellchat.hs.lrdb`: two interactions had the wrong genes**, copied from CellChat's gene symbol columns: CXCL8_CXCR1 scored
  CXCL1 (the same genes as CXCL1_CXCR1), and CD99_PILRB had the mouse symbols Cd99 / Pilrb, so it never matched human data. Both now
  have CellChat's genes (CXCL8; CD99 / PILRB). Results for these two interactions change; on human data where they are expressed,
  the background pairs (drawn from the interactions) and p.adj of other pairs can change too. `cellchat.mm.lrdb` matches CellChat on
  every row.
* The help of both tables now says why the 12 `*_CALCRL_RAMP1-3` interactions have the same genes as the `*_CALCR_RAMP1-3`
  ones: CellChat's complex table gives CALCRL_RAMP1-3 the subunit CALCR, not CALCRL. The tables keep CellChat's genes.

# HyperCom 0.3.0

Bigger samples in less memory, faster p-values and simpler arguments.

* **Memory mode.** A sample whose exact scoring does not fit in memory is scored in a mode that needs about a third of the memory,
  with priorities within ~1e-11 of exact scoring. About 100,000 cells in one sample fit in 64 GB (8.2 h on 16 threads, with
  p-values).
* **Faster p-values.** With `significance=TRUE`, the real pairs are screened like the background pairs, and only those whose rank
  or p-value is in doubt are scored fully. The other priorities are estimates: typically within ~1e-9 of exact (0.2.2's), up to
  ~1e-7 on the bundled data (1.5e-8 measured on the full bundled run). P-values, p.adj and order are exact. `exact.priority=TRUE`
  gives exact values.
* Runs without p-values that fit in memory give 0.2.2's results.

## Breaking changes

* **All cores by default.** `parallel` is gone from every function, and every function now uses one thread per physical core by
  default (`threads=NULL`); before, the default was one thread. Give `threads=1` to run on one thread, or another number. Results
  do not depend on the number of threads. The threads run at a lower priority than other programs, so the computer stays
  responsive; cluster jobs and containers use only their allocated cores, and `options(mc.cores=)` caps the default.
* **`cache` is now `resume`**, in `HyperCom()`, `prioritize_lr()` and the wrappers.
* **R 4.2 or later** (64-bit) and Matrix 1.5 or later are needed.
* **Esc / Ctrl-C now stops a running step** within about a second; a run stopped while pairs are scored continues from its cache
  when run again.
* **New argument orders.** The leading arguments are unchanged, so calls that give the first arguments by position still work; give
  the others by name:
  * `HyperCom(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=1, significance=FALSE, seed=NULL, r=0.5,
    infmat.method="woodbury", exact.priority=FALSE, threads=NULL, memory.gb=NULL, resume=TRUE, verbose=TRUE, outdir="HyperCom")`
  * `prioritize_lr(lrs, adj, infMat, metadata, weights=1, significance=FALSE, seed=NULL, exact.priority=FALSE, threads=NULL,
    memory.gb=NULL, resume=TRUE, verbose=TRUE, outdir="6.priority")`
  * `score_hypercom(infMat, metadata, ligand=NULL, receptors=NULL, threads=NULL, memory.gb=NULL, verbose=FALSE)`
  * `generate_infMat(adj, r=0.5, method=NULL, hyperedges=NULL, threads=NULL, outfile="4.infMat.rds")`: `outfile` moved last;
    `method=NULL` picks `"woodbury"` when the hyperedges are known and `"lu"` otherwise, as leaving out `method` did before.
* **Wrappers pass arguments on.** `HyperCom_seurat(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, ...)`
  and `HyperCom_sce(sce, lrdb, sample, assay="counts", group=NULL, keep=NULL, ...)` pass every other argument on to `HyperCom()`;
  a misspelled argument stops with "unused argument".
* The bundled mouse database `cellchat.mm.lrdb` has 2,374 interactions: the 3 whose receptor has no mouse gene (RARRES2_GPR1,
  H2-D_KLRA, H2-K1_KLRA) are left out. They were never used.
* Kept numeric and logical metadata columns (`keep`) stay numeric or logical: gene rows get `NA` in them instead of the text "gene"
  (text and factor columns still get "gene").
* The rows and columns of an influence matrix made by `"woodbury"` are named `node`, like those of the adjacency
  (`dimnames(infMat)`, `as.matrix(infMat)`, `infMat[i, j]`).
* The cache of a run stopped under 0.2.x is not reused (the version changed).

## New arguments

* `memory.gb` in `HyperCom()`, the wrappers, `prioritize_lr()` and `score_hypercom()`: the gigabytes (10^9 bytes) a run may use.
  The default is 90% of the computer's memory, or on Linux of a container's or cluster job's memory limit when that is smaller.
  HyperCom estimates the memory of scoring before it starts (about 17.6 bytes per pair of cells in the same sample, plus about
  1.1 GB) and switches to memory mode (about 5.5 bytes per pair) when that does not fit; memory mode needs the default
  `infmat.method="woodbury"`. A run that cannot fit even so stops at once, giving the estimate. On a cluster, give the job's memory
  request. `score_hypercom()` switches the same way, so it works after a `HyperCom()` run that needed memory mode.
* `exact.priority` in `HyperCom()`, the wrappers and `prioritize_lr()` (default `FALSE`): `TRUE` scores every real pair exactly
  when `significance=TRUE`.
* `r` in `HyperCom()` and the wrappers: the restart probability of `generate_infMat()` (default 0.5, as before).
* `verbose` in `score_hypercom()`: reports the memory the scoring needs.
* `HyperCom(outdir=NULL)` runs without writing files (and keeps no cache to resume from).

## Fixes

* Cells and genes that no path of the hypergraph connects to the rest are left out by `generate_hyperedges()`, with a warning that
  counts them (they could neither send nor receive anything); `HyperCom()` and `prioritize_lr()` then leave out the interactions of
  their genes, `prioritize_lr()` with a warning.
  Connected data are unchanged.
* `print()` and `as.matrix()` of the adjacency and of the influence matrix work without `library(HyperCom)` (e.g.
  `HyperCom::HyperCom()`), and `as.matrix()` gives a matrix when the Matrix package is attached.
* Without gene rows in `metadata`, priorities counted one group too few; groups are now counted over the cells.
* `HyperCom()` on counts held on disk or delayed (e.g. a DelayedArray or HDF5 matrix) reads the genes of the database only, as
  `HyperCom_sce()` did, instead of failing at the hypergraph step.
* `generate_infMat(method="lu")` checks `r` like the other methods.

## New checks

* The interaction table (`lrdb`, and `lrs` in `generate_hyperedges()` and `prioritize_lr()`) is checked: the second column must be
  `ligand` and the others `r1`, `r2`, ... (another column stops, naming it, with a hint for CellChat's own interaction table);
  interaction names must be present, distinct and differ from gene names; factors and spaces are cleaned and empty names (`""`) read
  as missing. Interactions without a ligand or a receptor are left out with a warning; an empty table stops.
* Missing gene names (`NA` row names of counts) are never matched.
* `prioritize_lr()` and `score_hypercom()` stop at once on metadata that does not come from `generate_metadata()` (columns cell,
  sample and group), and `prioritize_lr(significance=TRUE)` when it has no gene rows.
* `HyperCom()` stops before it writes a file when a cell that expresses a ligand or receptor gene has a missing sample or group, or
  one named "gene", and when the scoring cannot fit in memory; `outdir` naming an existing file stops too.
* Memory mode stops at the start, giving the memory of exact scoring, when the hypergraph is not connected (possible only with
  hyperedges not made by `generate_hyperedges()`).
* A cache that is not reused because the run switched between exact scoring and memory mode (which follows `memory.gb` and the
  computer's memory) now says so; caches made from other inputs, by another version of HyperCom, or damaged are moved aside as before.
* `significance`, `exact.priority`, `resume` and `verbose` must each be `TRUE` or `FALSE`.
* A missing `sample` argument stops with how to give one (one sample: `metadata$sample <- "all"`).

# HyperCom 0.2.2

Faster background screening on big samples (13% less time per screening batch on one sample of 30,000 cells), same results.

* Screening keeps each cell's two rank sums and skipped count side by side instead of in three arrays, so each entry touches one cache
  line. Priorities, p-values and order are unchanged. A cache made by an earlier version is not reused (the version changed).

# HyperCom 0.2.1

Faster scoring and screening (whole `prioritize_lr` with p-values 1.11-1.29x faster on 3,000-20,000 cells, 1.5x on one sample of
30,000 cells), same priorities, p-values and order on the bundled data.

* Cell scores are each sample's exact integer rank sum divided once, correctly rounded, instead of a long double running mean. They
  are within one ulp (~1e-16) of 0.2.0's and the same on every computer. A cache made by an earlier version is not reused.
* The grid kernel samples big grids more sparsely, writes keys to their buckets with non-temporal stores (x86), ranks with one wide
  radix digit and, for big buckets, sorts and ranks each cache-sized piece at once; background screening keeps its bucket ids from
  one batch to the next and sums with fewer operations. None of these changes a result.

# HyperCom 0.2.0

The same results as 0.0.0.9000, much faster and with far less memory. On the bundled data every step, the cell scores, the
priorities and the seeded p-values are identical to 0.0.0.9000's. Two exceptions: the default influence matrix method `"woodbury"`
differs from 0.0.0.9000's by about 1e-15 (`infmat.method="lu"` gives its numbers exactly; on the bundled data the priorities and
p-values are the same either way), and pairs whose cells are perfectly balanced over groups get priority 0 where 0.0.0.9000 gave NaN.

## What changes for you

* `HyperCom()` and `HyperCom_seurat()` return, invisibly, a list of their results instead of 0: `priority` (one dataframe per weight),
  `metadata`, `infMat`, `adj` and `lrs`, the objects also written to `outdir` (the files are as before). Results need not be read back:
  `res <- HyperCom(...)`, then `res$priority[["1"]]` or `score_hypercom(res$infMat, res$metadata, "Thbs1", "Sdc4")`.
* Progress is reported with `message()` instead of `print()`, so `suppressMessages()` or `verbose=FALSE` silences it: what each step
  found, the pairs scored about every 5 s with the time left, and the top pair of each weight. The per-pair score lines are gone; the
  scores are in the returned tables and the CSV files.
* The adjacency and influence matrices are no longer node x node matrices by default (see below). They work like matrices for `dim`,
  `dimnames`, `rownames`, `x[i, j]` and `as.matrix()`; code that used them as base matrices (`rowSums(adj)`, heatmaps) needs
  `as.matrix()` first. `.rds` files of them need `library(HyperCom)` to be read. Plain matrices, e.g. files of earlier versions, work
  everywhere as before.
* Parallel runs use threads instead of a cluster of R workers: `parallel=TRUE` uses one thread per physical core and `threads` sets
  the number. The `cl` argument of `prioritize_lr`, `HyperCom` and `HyperCom_seurat` is gone: remove `cl=` from calls (it now stops
  with "unused argument"), and in calls by position drop the value given for `cl`, which would otherwise be read as `outdir`.
* The package has compiled code: installing from source needs build tools (see README).
* Coming from branch spatial (was v1-fast, 0.1.x): 0.2.0 has none of its additions (spatial scores and adjacency, `generate_cell_info`,
  `auto_grid_mem_gb`, `available_memory_gb`, `prioritize_lr(cores=, infMat.cell=, cell.info=)`, `method="auto"`). Its version is
  higher, so `update.packages()` moves spatial (v1-fast) installs to 0.2.0.

## New function

* `HyperCom_sce()` runs the pipeline on a SingleCellExperiment, as `HyperCom_seurat()` does on a Seurat object: one assay
  (`assay="counts"` by default) and the cell metadata (`colData`); `sample`, `group` and `keep` name columns of `colData`. Assays held
  on disk (DelayedArray, HDF5Array) are read for the genes of `lrdb` only. Needs SummarizedExperiment (Bioconductor).

## New arguments

* `prioritize_lr`, `HyperCom`, `HyperCom_seurat`: `threads` (number of threads; below 1 stops, more than the available cores warns and
  uses them all), `seed` (set before drawing the background pairs, so p-values are reproducible; ignored without `significance`; the session's
  random numbers are restored afterwards) and `verbose`, after `outdir`.
* `score_hypercom`: `parallel` and `threads` as `prioritize_lr` (one thread by default). On 20,000 cells one pair took 4.6 s on 1 thread
  and 0.9 s on 16. Scores are the same at any number of threads.
* `generate_infMat`: `method` (`"woodbury"`, the default, `"lu"` or `"chol"`), `hyperedges` and `threads`; `HyperCom` and
  `HyperCom_seurat`: `infmat.method`. On a plain adjacency matrix given without its `hyperedges` the default is `"lu"`.
* `prioritize_lr`, `HyperCom`, `HyperCom_seurat`: `cache` (default `TRUE`). With an output directory, the pairs scored so far are
  kept in `<outdir>/cache` (`<outdir>/6.priority/cache` for `HyperCom`) while pairs are scored, and the folder is removed once the
  priority files are written. A run that stops (time limit of a cluster job, crash) continues where it stopped when the same call
  is run again, losing at most the last batch of pairs, with the same results as a run that never stopped on the same type of
  computer (another CPU can change the last digits, about 1e-13); without `seed` it keeps
  the background pairs drawn before the stop. A cache made from other inputs is not used and is moved to `cache.old`, with a
  warning; a cache that cannot be written is skipped, with a warning. The cache holds 5
  numbers per scored pair and 10 per screened background pair (about 200 KB for the bundled data with significance) and the names
  of the cells (0.8 MB at 20,000 cells); on 20,000
  cells its writes took 0.1 s of a 77 s scoring step, with the same peak memory.
* `prioritize_lr` documents its p-values: the fraction of 1000 background pairs scoring at least as high, in steps of 0.001; 0 means
  none did.

## Speed and memory

On 20,000 synthetic cells in 3 samples (16 threads, a tenth of the CellChat pairs, no p-values) the whole pipeline takes about 1.3 min
and peaks at 3.5 GB; 0.0.0.9000 ran out of memory at this size. 50,000 cells take about 10 min and 17 GB. The bundled data take
seconds.

* Cells are scored by compiled code on threads that share one influence matrix, only within their sample (every cross-sample entry
  was NA), once per ligand and receptor subunit for all weights and pairs. Weights and significance no longer multiply the run time.
* With `significance`, the 1000 background pairs are screened instead of scored: one pass over a background pair's cell grids bounds
  its priority without ranking the cells, and only the few whose bounds hold a real pair's priority (about 1%) are scored in full.
  The p-values are the same as with every background pair scored. On 16 threads a run with p-values took 1.7 s instead of 2.9 s on
  the bundled data (12 s instead of 20 s on 1 thread), 31 s instead of 90 s on 3,000 cells, 85 s instead of 240 s on 5,000 cells
  (one sample) and 8.4 min instead of 23 min on 20,000 cells in 3 samples, with no more memory. A stopped run keeps the screened
  bounds in its cache too.
* The adjacency and influence matrices are kept as their factors through the hyperedges (genes and ligand-receptor interactions),
  which grow with the number of cells rather than its square; only the within-sample cell blocks the scores need are made, while
  pairs are scored. Small data with at least as many hyperedges as nodes (with the full CellChat database, roughly under 1000 cells)
  get the dense influence matrix as before, which is then no bigger than the factors.
* Memory grows with the sum of squared sample sizes (8 bytes per within-sample cell pair for the blocks, 8 to 10 more while a grid is
  ranked) plus 8 bytes per hyperedge and node, so more samples of the same size cost little.
* No limit on cells per sample but memory: earlier builds of this version stopped when the sum of squared sample sizes reached 2^32
  (one sample of 65,535 cells). One sample of 65,535 cells needs about 70 GB, 100,000 cells about 160 GB. The run time grows with the
  within-sample cell pairs times the grids: one sample of 65,540 cells took about 45 s per grid on 16 threads of a desktop.
* Influence matrix methods: `"woodbury"` 0.4 s at 20,000 cells, `"chol"` 46 s with one dense node x node matrix, `"lu"` 105-121 s with
  two (0.0.0.9000 held five). `"chol"` and `"woodbury"` differ from `"lu"` by about 1e-15 relative, which can reorder exactly tied
  grid entries: on 20,000 cells the pairs keep their order and priorities move by at most 1.3e-13 (`"woodbury"`) and 1.9e-14
  (`"chol"`).
* `generate_infMat` saves the influence matrix uncompressed (gzip saved 7% of the file and took 50 times longer than writing it).

## Input checks

Mistakes that ran to meaningless results, or failed late with unclear messages, now stop at once with a message that says what to
fix. Results of valid inputs are unchanged.

* `HyperCom` checks `lrdb`, the `sample` / `group` / `keep` columns, `weights` and `infmat.method` before any step runs (a misspelt
  column used to fail after the influence matrix; `lrdb` and `counts` swapped took 27 s and 6 GB), and stops when no interaction of
  `lrdb` has all its genes expressed, e.g. a human `lrdb` on mouse counts (it used to finish with an empty priority table).
* `filter_counts` stops on counts with genes as column names (transposed) and on negative values (scaled data, which ran as if every
  value were expression).
* `generate_metadata` stops when a cell of the graph has no sample or group (they were turned into genes), when a sample or group is
  named `"gene"`, when a sample has fewer than 2 cells in the graph (its cells get no scores, and with `significance` every p-value
  became NA), and when no node is a row name of the metadata; `sample` and `group` must be one column name each.
* `generate_hyperedges` stops on cells sharing a barcode (they were mangled with `make.names` and merged or lost).
* `score_hypercom` stops when given only one of `ligand` and `receptors` (a ligand alone gave the network score, receptors alone
  all-zero scores), when `ligand` names more than one gene and when `receptors` is empty; `score_hypercom` and `prioritize_lr` name
  genes missing from the graph instead of failing with a subscript error.
* `prioritize_lr` stops on weights that are not numbers of at least 0 (negative weights gave infinite priorities) or that would share
  a name and an output file, on a `threads` that is not one number, and on an `outdir` that is not one string (or NULL).
* `HyperCom` leaves out cells of the counts that are not in the metadata, with a warning (they used to become gene nodes and be
  drawn as background genes, which skewed the p-values), stops when no cell is, and warns when `outdir` already holds priority
  files of other weights, which it keeps. `filter_counts` warns about gene names used by several rows (the first row is kept).

## Bug fixes

* Cells perfectly balanced over 5, 10 or 11 groups no longer give NaN priorities: the group purity term rounded to -1.1e-16, whose log
  is NaN, and every p-value then became NA. Purity below 0 now counts as 0, so such pairs get priority 0.
* A background pair without a priority (NaN) no longer makes every p-value NA: it is left out of the p-values, with a warning.
* Metadata rows are matched to the influence matrix by name, so metadata in another order gives the same scores (they were matched by
  position, giving wrong scores).
* Ligand and receptor columns stored as factors (an lrdb read with `stringsAsFactors=TRUE`) are read by name; the adjacency matrix was
  indexed by factor codes and gave wrong priorities without a warning, and an all-NA receptor column dropped every interaction.
* Nodes are ordered in the C locale and factor `node` / `edge` columns by name, so results no longer depend on the session's
  collation.
* `generate_metadata`: a single `keep` column keeps its name; kept columns named `cell`, `sample` or `group` are renamed `cell.1`, ...
  instead of breaking the join; only gene nodes are labelled `"gene"` (cells keep their missing `keep` values); factor `keep` columns
  get a `"gene"` level instead of NA, also when kept twice; dates and other classed columns are kept as strings.
* `generate_hyperedges` takes the interaction names from the first column, whatever its name.
* `filter_counts` ignores missing counts in the gene means (it stopped with "subscript out of bounds") and works with a single gene;
  `generate_lrs_table` keeps its columns when no interaction is present.
* `HyperCom` no longer changes the working directory; it writes into `outdir`, created recursively. `prioritize_lr` creates nested
  output directories. With the cluster gone, it no longer stops a cluster it was given or defaults to an SGE `qrsh` one.
* `HyperCom_seurat` reads the layer asked for with SeuratObject 4 as well (its `GetAssayData` ignored `layer` and returned `data`).
* `library(HyperCom)` loads Matrix, so the bundled sparse counts can be subset straight away.
* When R runs out of memory inside the compiled code (very large data), the inputs of that step are released again; they used to
  stay in memory until R was restarted.

## Package

* Dropped dependencies: collapse, doParallel, foreach, mltools, Rfast; SeuratObject moved to Suggests (only `HyperCom_seurat` needs
  it). Minimum versions: dplyr 1.1.0, parallelly 1.23.0 (older versions failed with cryptic errors).
* README: prerequisites and installation per operating system, the common installation errors, memory and run time, and how
  to install the R-only version (0.0.0.9000) from the `base-R` branch. DESCRIPTION describes the method and has URL and BugReports; `citation("HyperCom")`
  gives the bioRxiv preprint. The bundled datasets are documented.
* Every help page has examples on the bundled data (each runs in under a second). The overview vignette is part of the package
  (`vignettes/HyperCom_Overview.Rmd`, R Markdown instead of Quarto), so R CMD build renders it and R CMD check runs it; it needs
  knitr, rmarkdown, ggplot2 and ggridges (Suggests). The rendered HTML and the outdated PDF are no longer in the repository.
* Docs: `filter_counts` returns a cell x gene matrix; `generate_metadata` returns a dataframe; the lrdb's first column may have any
  name; vignette `HyperCom_seurat` example fixed.

## Developer notes

* Grid kernel (`src/score_grids.cpp`): an exact integer sort instead of `percent_rank` on n x n values; means summed in long double in
  0.0.0.9000's order, so scores are bit for bit; about 10 bytes per within-sample entry for the grids in flight. Grids too big for the
  cache are scored one at a time by all threads, with 8 bytes per entry (6-byte keys) and value buckets of about equal counts (bounds
  from a histogram of a sample of the grid).
* `generate_infMat`: `"lu"` has LAPACK's `dgesv` (the routine `solve()` calls) factor the system in place; `"chol"` is HyperCom
  0.1.0's Cholesky inverse (`dpotrf` + `dpotri`), bit for bit; `"woodbury"` uses the Woodbury identity on adj = H H' (H the node x
  hyperedge incidence) and checks exactly that given `hyperedges` make the adjacency. The package links R's LAPACK (on Windows through
  an import library for `Rlapack.dll`) and is compiled as C++17 on every R version.
* `generate_adjacency_matrix` keeps the incidence H; `as.matrix(adj)` accumulates the upper triangle in C++ (0.0.0.9000's matrix
  exactly, 4x faster than `xtabs` + `tcrossprod`). A `hypercom_influence` keeps `H`, the row sums, `delta` and the hyperedge x node
  matrix `Z`; its entries equal the dense `"woodbury"` matrix's bit for bit.
* Tests compare with fixtures made by 0.0.0.9000: the LU influence matrix bit for bit with 0.0.0.9000's formula run on the same
  machine and within 1e-12 with the fixture (whose last bits depend on the BLAS it was made with); scores and priorities bit for bit
  where long double is x87 (x86-64), within 1e-12 elsewhere (arm64).
* The R code is split by pipeline stage (`R/preprocess.R`, `hypergraph.R`, `influence.R`, `metadata.R`, `score.R`, `prioritize.R`,
  `pipeline.R`, `utils.R`) instead of one `R/functions.R`; the code itself is unchanged.
* No more `:::` to the package's own functions; `weighted.mean` is called as `stats::weighted.mean`.
* Large R allocations in the compiled code go through `Rcpp::unwindProtect` (`src/r_alloc.h`), so an allocation error unwinds the
  C++ stack (destructors run, inputs are unprotected) before R's error is raised. The grid kernel packs its 6-byte keys byte by
  byte, so big grids rank the same on big-endian machines (a 6-byte `memcpy` of a `uint64_t` kept the wrong bytes there).
