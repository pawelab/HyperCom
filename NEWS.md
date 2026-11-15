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
* Coming from branch v1-fast (0.1.x): 0.2.0 has none of its additions (spatial scores and adjacency, `generate_cell_info`,
  `auto_grid_mem_gb`, `available_memory_gb`, `prioritize_lr(cores=, infMat.cell=, cell.info=)`, `method="auto"`). Its version is
  higher, so `update.packages()` moves v1-fast installs to 0.2.0.

## New arguments

* `prioritize_lr`, `HyperCom`, `HyperCom_seurat`: `threads` (number of threads; below 1 stops, more than the available cores warns and
  uses them all), `seed` (set before drawing the background pairs, so p-values are reproducible; ignored without `significance`; the session's
  random numbers are restored afterwards) and `verbose`, after `outdir`.
* `score_hypercom`: `parallel` and `threads` as `prioritize_lr` (one thread by default). On 20,000 cells one pair took 4.6 s on 1 thread
  and 0.9 s on 16. Scores are the same at any number of threads.
* `generate_infMat`: `method` (`"woodbury"`, the default, `"lu"` or `"chol"`), `hyperedges` and `threads`; `HyperCom` and
  `HyperCom_seurat`: `infmat.method`. On a plain adjacency matrix given without its `hyperedges` the default is `"lu"`.
* `prioritize_lr` documents its p-values: the fraction of 1000 background pairs scoring at least as high, in steps of 0.001; 0 means
  none did.

## Speed and memory

On 20,000 synthetic cells in 3 samples (16 threads) the whole pipeline takes about 1.5 min and peaks at 3.5 GB; 0.0.0.9000 ran out of
memory at this size. 50,000 cells take about 28 min and 19 GB. The bundled data take seconds.

* Cells are scored by compiled code on threads that share one influence matrix, only within their sample (every cross-sample entry
  was NA), once per ligand and receptor subunit for all weights and pairs. Weights and significance no longer multiply the run time.
* The adjacency and influence matrices are kept as their factors through the hyperedges (genes and ligand-receptor interactions),
  which grow with the number of cells rather than its square; only the within-sample cell blocks the scores need are made, while
  pairs are scored. Small data with at least as many hyperedges as nodes (with the full CellChat database, roughly under 1000 cells)
  get the dense influence matrix as before, which is then no bigger than the factors.
* Memory grows with the sum of squared sample sizes (8 bytes per within-sample cell pair for the blocks, about 10 more while a grid is
  ranked) plus 8 bytes per hyperedge and node, so more samples of the same size cost little.
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
* Every help page has examples on the bundled data (each runs in under a second). The overview vignette's HTML is
  rendered from the current source; its outdated PDF is gone.
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
