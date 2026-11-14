# HyperCom 0.0.0.9001 (development)

Same results as 0.0.0.9000 (master b999a01) bit for bit on the bundled data: counts, LR table, hyperedges, adjacency and influence
matrices, metadata, cell scores, priorities and seeded p-values are `identical()` to master's under the C collation.

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
* `generate_infMat` writes the system matrix straight from the sparse adjacency and has LAPACK's `dgesv` (the routine `solve()`
  calls) factor it in place into the identity that becomes the result: at most 2 dense node x node matrices instead of 5, and faster.
  The package now links R's LAPACK (on Windows through an import library for `Rlapack.dll`).
* `generate_adjacency_matrix` accumulates the upper triangle directly in C++ instead of `xtabs` + `tcrossprod` (whose full product,
  diagonal reset and node removal each made a copy): the result is allocated once, master's object exactly, 4x faster.
* `prioritize_lr` scores pairs in batches of whole ligands and keeps only one batch of cell scores, instead of the scores of every
  grid at once (16 bytes x cells x grids); no grid is scored twice. The gene x cell counts are read from the stored triangle of the
  adjacency instead of a full symmetric expansion.
* `filter_counts` keeps the ligand and receptor genes before transposing, instead of transposing the whole counts matrix.

## Influence matrix methods

`generate_infMat(method=)` chooses how the inverse is taken; the default `"lu"` is master's, bit for bit. On 20,000 synthetic cells
(20,653 nodes, 1641 hyperedges, 16 threads) the influence matrix takes 121 s with `"lu"`, 46 s with `"chol"` and 1.4 s with
`"hypergraph"`, and the whole pipeline peaks at 7.7 GB instead of 9.1 GB with either of the new methods.

* `"chol"`: v1-fast's Cholesky inverse of the symmetric form of the system, in place: one dense node x node matrix instead of two,
  about 3/8 of the arithmetic. The same numbers as v1-fast's `method="chol"` bit for bit.
* `"hypergraph"` (needs `hyperedges`, the table the adjacency was made from): the adjacency is the hyperedge incidence times its
  transpose, so the inverse is a diagonal plus a term of rank at most the number of hyperedges (Woodbury identity). Work grows with
  nodes x incidences instead of nodes cubed, the result is the only dense matrix, and it is the same at any `threads`. It checks
  exactly that the hyperedges give the adjacency (in any node order) and stops otherwise; with at least as many hyperedges as nodes it
  computes `"chol"`.
* Both stop like `"lu"` when the system is computationally singular (`r` near 0) and on non-finite weights.
* Both differ from `"lu"` by about 1e-15 relative (6e-15 at 20,000 cells), which can break exact ties when grids are ranked; on 20,000
  cells the pairs keep master's order and priorities move by at most 1.3e-13 (`"hypergraph"`) and 1.9e-14 (`"chol"`).
* `HyperCom` and `HyperCom_seurat` take `infmat.method` (default `"lu"`) and pass their hyperedges and `threads`.

## New arguments and changed behaviour

* `prioritize_lr`, `HyperCom`, `HyperCom_seurat`: `threads` (number of threads, capped at the available cores), `seed` (set before drawing the background
  interactions, so p-values are reproducible; ignored without `significance`), `verbose` (`FALSE` prints nothing). They come after
  `outdir`, so positional calls written for master still work. `parallel=TRUE` uses `parallelly::availableCores()` threads; `cl` is
  accepted but only caps the thread count and is never stopped.
* The per-pair score lines are printed once all pairs are scored, not as progress.
* `generate_metadata` stops when a cell of the graph has no sample or group (master turned such cells into genes); metadata rows
  that are not in the graph may still have missing values.
* The package has compiled code: installing from source needs a C++ compiler (Rtools on Windows).

## Bug fixes

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
* No more `:::` to the package's own functions; `weighted.mean` called as `stats::weighted.mean`.
* Docs: `filter_counts` returns a cell x gene matrix; `generate_metadata` returns a data frame; the lrdb's first column may have any
  name; vignette `HyperCom_seurat` example fixed.
