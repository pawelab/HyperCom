# HyperCom 0.2.0

Same results as 0.0.0.9000 (master b999a01) bit for bit on the bundled data: counts, LR table, hyperedges, adjacency and influence
matrices, metadata, cell scores, priorities and seeded p-values are `identical()` to master's under the C collation. One exception
by default: the influence matrix is now made with method `"woodbury"`, about 1e-15 from master's (on the bundled data `HyperCom()`
still writes master's priorities and seeded p-values); method `"lu"` gives master's numbers bit for bit. Where master gave NaN
priorities for cells perfectly balanced over groups, priorities are now 0 (see Bug fixes). The adjacency and influence matrices are
no longer stored as node x node matrices by default (see below).

## Speed

* Cell grids are scored in C++ (`src/score_grids.cpp`): only within-sample blocks (every cross-sample entry is NA), an exact integer
  sort instead of `percent_rank` on n x n values, means summed in master's order, no n x n R matrices.
* `prioritize_lr` scores each ligand / receptor-subunit grid once, for all weights and all pairs sharing it, and the background pairs
  together with the real ones. Weights no longer multiply the run time.
* Parallel runs use C++ threads sharing one influence matrix instead of PSOCK workers that each copy it.
* `generate_hyperedges` reads the sparse counts directly instead of a dense data frame of every cell x gene.
* `generate_infMat` saves the influence matrix uncompressed: gzip saved 7% of the file and took 50 times longer than writing it.

## Memory

On 3000 synthetic cells (3 samples, 3648 nodes) the whole pipeline with significance peaks at 0.68 GB instead of 1.64 GB (0.42 GB
above its inputs instead of 1.37 GB), and every step is as fast or faster; outputs are unchanged bit for bit.

* The grid kernel keeps about 10 bytes per within-sample entry for the grids in flight: grids too big for the cache are scored one at
  a time by all threads (the row ranges of a big sample summed as a pipeline, in master's order) instead of one per thread, each value
  bucket is one exact run of memory instead of linked slabs, and the influence matrix is copied only for samples whose cells are not
  consecutive, while those copies fit the cache (else read in place). Kernel memory on 3000 cells at 16 threads: 46 MB instead of
  673 MB; two samples of 5000 cells: 0.49 GB instead of 2.35 GB. Faster at 8 or more threads.
* Grids too big for the cache, on several threads, keep 8 bytes per entry instead of 10 (each value bucket keeps at most 48 key
  bits, in 6 bytes) and their value buckets hold about equal numbers of entries (bounds from a histogram of a sample of the grid)
  instead of equal value ranges, so every bucket is sorted within the core's cache. On 20,000 synthetic cells (16 threads) the whole
  pipeline takes 4-7% less time and peaks 9% lower (3.29 instead of 3.63 GB); the gain grows with the grids. Ranks, and so every
  score, are the same bit for bit; smaller grids and one thread are scored as before.
* `generate_infMat` writes the system matrix straight from the sparse adjacency and has LAPACK's `dgesv` (the routine `solve()`
  calls) factor it in place into the identity that becomes the result: at most 2 dense node x node matrices instead of 5, and faster.
  The package now links R's LAPACK (on Windows through an import library for `Rlapack.dll`).
* When the adjacency matrix itself is made (`as.matrix(adj)`, methods `"lu"` and `"chol"`), its upper triangle is accumulated directly
  in C++ instead of `xtabs` + `tcrossprod` (whose full product, diagonal reset and node removal each made a copy): allocated once,
  master's object exactly, 4x faster.
* `prioritize_lr` scores pairs in batches of whole ligands and keeps only one batch of cell scores (at most about 64 MB), instead of
  the scores of every grid at once (16 bytes x cells x grids); no grid is scored twice. Which cells express each gene is kept as a list
  of cells, not a cells x genes table. The gene x cell counts are read from the hyperedge incidence (or
  the stored triangle of an adjacency matrix) instead of a full symmetric expansion.
* `filter_counts` keeps the ligand and receptor genes before transposing, instead of transposing the whole counts matrix.

## No node x node matrices

The adjacency and influence matrices have as many rows and columns as nodes (cells + genes), so they grew with the square of the
number of cells. They are now kept as their factors through the hyperedges (genes and ligand-receptor interactions), which grow with
the number of cells; only the within-sample cell blocks the scores need are made, while pairs are scored. Small data with at least as
many hyperedges as nodes (with the full CellChat database, roughly under 1000 cells) get the dense influence matrix as before, which is
then no bigger than the factors. On 20,000 synthetic cells
in 3 samples (20,653 nodes, 1641 hyperedges, 16 threads) the whole pipeline peaks at 3.65 GB instead of 9.1 GB, the adjacency takes
0.06 s instead of 5.4 s and the influence matrix 0.4 s instead of 105 s; scores are as before.

* `generate_adjacency_matrix` returns a `hypercom_adjacency`: the node x hyperedge incidence `H` (how often each node is in each
  hyperedge); the adjacency is H H' with a zero diagonal. It works like a matrix for `dim`, `dimnames`, `rownames`, `adj[i, j]` (those
  entries) and `as.matrix` (master's matrix exactly). `3.adj.rds` holds it.
* `generate_infMat` on such an adjacency returns a `hypercom_influence` by default (method `"woodbury"`): the factors of the
  influence matrix (`H`, the row sums, `delta` and the hyperedge x node matrix `Z`), with the same methods; its entries equal the dense
  `"woodbury"` matrix's bit for bit. `4.infMat.rds` holds it. Methods `"lu"` and `"chol"` return the dense matrix as before; on a plain
  adjacency matrix given without its `hyperedges` the default is still `"lu"`.
* `prioritize_lr` and `score_hypercom` take either kind; with a factored influence matrix they make the within-sample cell blocks (the
  sum of squared sample sizes of doubles) and the genes' rows and columns, and give the same scores as from its dense matrix, bit for
  bit. `generate_metadata` reads only the node names.
* Plain matrices, e.g. files saved by earlier versions, work everywhere as before. Code that used the adjacency or influence matrix
  as a base matrix (`rowSums(adj)`, `Matrix` functions, heatmaps) needs `as.matrix()` first, which makes the node x node matrix
  (`as.matrix(infMat, threads=)` makes it on several threads; indexing uses one).
* Memory now grows with the sum of squared sample sizes (8 bytes per within-sample cell pair for the blocks, about 10 for the pair
  being ranked) plus 8 bytes per hyperedge x node for `Z`, not with the number of nodes squared: more samples of the same size cost
  little. On 50,000 synthetic cells in 3 samples (16 threads) the pipeline peaks at 19.0 GB; before, it needed about 46 GB.

## Influence matrix methods

`generate_infMat(method=)` chooses how the inverse is taken: `"woodbury"` (the default), `"lu"` (master's, bit for bit; still what
an adjacency matrix gets by default when it comes without its `hyperedges`) or `"chol"`. On 20,000 synthetic
cells the dense influence matrix took 121 s with `"lu"` and 46 s with `"chol"`; the factors take 0.4 s (`"woodbury"`).

* `"chol"`: v1-fast's Cholesky inverse of the symmetric form of the system, in place: one dense node x node matrix instead of two,
  about 3/8 of the arithmetic. The same numbers as v1-fast's `method="chol"` bit for bit.
* `"woodbury"`, the hypergraph method: the adjacency is the hyperedge incidence times its transpose, so the inverse is a diagonal plus a term of rank at
  most the number of hyperedges (Woodbury identity), kept as its factors. Work grows with nodes x incidences instead of nodes cubed,
  and it is the same at any `threads`. From an adjacency matrix it needs `hyperedges`, the table the matrix was made from, checks
  exactly that they give the adjacency (in any node order) and stops otherwise. With at least as many hyperedges as nodes it gains
  nothing and returns `"chol"`'s dense matrix.
* Both stop like `"lu"` when the system is computationally singular (`r` near 0) and on non-finite weights.
* Both differ from `"lu"` by about 1e-15 relative (6e-15 at 20,000 cells), which can break exact ties when grids are ranked; on 20,000
  cells the pairs keep master's order and priorities move by at most 1.3e-13 (`"woodbury"`) and 1.9e-14 (`"chol"`).
* `HyperCom` and `HyperCom_seurat` take `infmat.method` (default `"woodbury"`; `"lu"` for master's numbers bit for bit) and pass
  their thread count (`parallel`, `cl`, `threads`: one thread without `parallel`, as for scoring).

## New arguments and changed behaviour

* `score_hypercom` takes `parallel` and `threads` like `prioritize_lr` (one thread by default); scores are the same at any number of
  threads. On 20,000 cells one pair took 4.6 s on 1 thread and 0.9 s on 16.

* `prioritize_lr`, `HyperCom`, `HyperCom_seurat`: `threads` (number of threads: below 1 stops, more than the available cores warns and uses them all), `seed` (set before drawing the background
  interactions, so p-values are reproducible; ignored without `significance`), `verbose` (`FALSE` reports nothing). They come after
  `outdir`, so positional calls written for master still work. `parallel=TRUE` uses one thread per physical core (`parallelly::availableCores(logical=FALSE)`); `cl` is
  accepted but only caps the thread count and is never stopped.
* `HyperCom` and `HyperCom_seurat` return, invisibly, a list of their results instead of 0: `priority` (one dataframe per weight),
  `metadata`, `infMat`, `adj` and `lrs`, the objects also written to `outdir`. Results need not be read back from the files:
  `res <- HyperCom(...)`, then `res$priority[["1"]]` or `score_hypercom(res$infMat, res$metadata, "Thbs1", "Sdc4")`. Code that compared
  the value with 0 needs changing; the files are as before.
* Progress is reported with `message()` instead of `print()`, so `suppressMessages()` silences it: what each `HyperCom` step
  found (genes, cells, interactions, nodes), the pairs scored about every 5 s with an estimate of the time left, and the top pair
  and output file of each weight. On data whose grids exceed the kernel's 24 MB cache budget (about 1600 cells in one sample), pairs
  are scored in about 20 batches or more, so that progress comes often; smaller data keep the larger batches. Each batch sets up the
  kernel's buffers again, so a run of few pairs on big data takes a little longer (20,000 cells, 104 pairs: 90 instead of 86 s); with
  the full CellChat database the number of batches is about the same as before. The
  per-pair score lines (one per pair and background pair and weight, 800-3600 lines) are gone; the scores are in the returned tables
  and the CSV files.
* `prioritize_lr` documents its p-values: the fraction of 1000 background pairs scoring at least as high, in steps of 0.001; 0 means
  none did.
* `generate_metadata` stops when a cell of the graph has no sample or group (master turned such cells into genes); metadata rows
  that are not in the graph may still have missing values.
* The package has compiled code: installing from source needs a C++ compiler (Rtools on Windows).

## Input checks

Mistakes that ran to meaningless results, or failed late with unclear messages, now stop at once with a message that says what to
fix. Results of valid inputs are unchanged.

* `HyperCom` checks `lrdb`, the `sample` / `group` / `keep` columns, `weights` and `infmat.method` before any step runs (a misspelt
  column used to fail after the influence matrix with "undefined columns selected"; `lrdb` and `counts` swapped took 27 s and 6 GB).
* `HyperCom` stops when no interaction of `lrdb` has all its genes expressed, e.g. a human `lrdb` on mouse counts; it used to finish
  with an empty priority table.
* `filter_counts` stops on counts with genes as column names (transposed) and on negative values (scaled data, which ran as if every
  value were expression).
* `generate_metadata` stops when a sample has fewer than 2 cells in the graph (its cells get no scores, and with `significance` every
  p-value became NA); `sample` and `group` must be one column name each; tibbles are read by column name.
* `score_hypercom` stops when given only one of `ligand` and `receptors` (a ligand alone gave the network score, receptors alone
  all-zero scores) and when `ligand` names more than one gene.
* `prioritize_lr` stops on weights that are not numbers of at least 0 (negative weights gave infinite priorities).
* The old method name `"hypergraph"` stops with its new name, `"woodbury"`.
* `HyperCom` warns when `outdir` already holds priority files of other weights, which it keeps.
* `HyperCom_seurat` reads v3 assays from the layer asked for with SeuratObject 4 as well (SeuratObject 4's `GetAssayData` ignored
  `layer` and returned the `data` slot).

## Bug fixes

* Cells perfectly balanced over 5, 10 or 11 groups no longer give NaN priorities: the group purity term rounded to -1.1e-16, whose
  log is NaN, and every p-value then became NA (also in master). Purity below 0 now counts as 0, so such pairs get priority 0. Results
  change only where master gave NaN.

* Cells sharing a barcode stop with a message in `generate_hyperedges` (master mangled them with `make.names` and merged or lost
  cells); metadata whose row names match no node stops in `generate_metadata`.
* Metadata rows are matched to the influence matrix by name, so metadata in another order gives the same scores (master matched by
  position and gave wrong scores).
* Factor `node` / `edge` columns of the hyperedges are ordered by name, as strings are (master sorted factor codes, so the node
  order and seeded p-values followed the collation the levels were made in); hyperedges without `node` and `edge` columns stop.
* `generate_metadata`: kept columns named `cell`, `sample` or `group` are renamed `cell.1`, ... instead of breaking the join (and
  a kept `cell` column holds the metadata's own values, not the row names); a factor column kept twice gets the `"gene"` level in both copies; dates and other classed columns are kept as strings.
* `filter_counts` reads an lrdb of factors by name and leaves missing counts out of the gene means (master stopped with "subscript out of bounds") and warns about gene
  names used by several rows (the first row is kept, as in master).
* `HyperCom` warns when cells of the counts are not in the metadata (they count as genes, as in master).
* `prioritize_lr` stops on weights that would share a name and an output file; `threads` must be one number; a number given as `cl`
  caps the threads; an lrdb read with `stringsAsFactors=TRUE` with an all-NA receptor column keeps its interactions (master kept none).
* Clear errors for a receptor complex with no gene, a pair with only one of ligand and receptor, metadata without cells, and counts
  where no ligand or receptor gene is expressed.

* `generate_metadata`: a single `keep` column keeps its name; a missing sample or group, or one named `"gene"`, stops with a clear
  message instead of turning cells into genes or failing on duplicated factor levels; only gene nodes are labelled `"gene"` (cells
  keep their missing `keep` values); factor `keep` columns get a `"gene"` level instead of NA.
* `generate_hyperedges` takes the interaction names from the first column, whatever its name.
* `generate_adjacency_matrix` orders nodes in the C locale, so results no longer depend on the session's collation.
* `score_hypercom` and `prioritize_lr` name genes missing from the graph instead of failing with a subscript error.
* `HyperCom` no longer changes the working directory; it writes into `outdir` (created recursively).
* `prioritize_lr` creates nested output directories, does not stop a cluster it was given, and no longer defaults to an SGE `qrsh`
  cluster.
* `generate_lrs_table` keeps its columns when no interaction is present (it returned a table with no columns); `filter_counts`
  works with a single gene.
* Ligand and receptor columns stored as factors are read by name; master indexed the adjacency matrix by factor codes and gave wrong
  priorities without a warning.
* `library(HyperCom)` loads Matrix, so the bundled sparse counts can be subset straight away.

## Package

* Dropped dependencies: collapse, doParallel, foreach, mltools, Rfast; SeuratObject moved to Suggests (only `HyperCom_seurat` needs it).
* Minimum versions: dplyr 1.1.0 (`join_by`), parallelly 1.23.0 (`availableCores(logical=)`); older versions failed with cryptic
  errors. C++ is compiled as C++17 on every R version.
* README: installation per operating system (build tools for the C++ code) with `pak`.
* Tests compare the LU influence matrix bit for bit with master's formula run on the same machine and within 1e-12 with the fixture,
  whose last bits depend on the BLAS it was made with (they failed on Linux and with R's reference BLAS). Scores and priorities are
  compared bit for bit with the fixture where long double is the x87 type it was made with (x86-64), within 1e-12 elsewhere (arm64).
* Version 0.2.0, above the 0.1.x of branch v1-fast, so installs of that branch update to it. DESCRIPTION has a description of the
  method, URL and BugReports; `citation("HyperCom")` gives the bioRxiv preprint.
* No more `:::` to the package's own functions; `weighted.mean` called as `stats::weighted.mean`.
* Docs: `filter_counts` returns a cell x gene matrix; `generate_metadata` returns a data frame; the lrdb's first column may have any
  name; vignette `HyperCom_seurat` example fixed.
