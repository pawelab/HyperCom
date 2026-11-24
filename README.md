# HyperCom (pure R version)

**See which single cells send and receive each ligand-receptor signal, with no clustering.**

HyperCom gives every cell a sender and receiver score for each ligand-receptor pair, then ranks the pairs that drive communication
in your data. With a condition column, pairs whose cells concentrate in some conditions rank higher. Start from a Seurat object, a
SingleCellExperiment or a counts matrix.

Paper: [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1)
(bioRxiv, 2026).

![HyperCom method overview](https://raw.githubusercontent.com/pawelab/HyperCom/master/images/MethodOverview.jpg)

## Which version?

This branch (`base-R`, HyperCom 0.3.1.1) is HyperCom 0.3.1 written in R only. It needs no build tools to install. The main version,
on [`master`](https://github.com/pawelab/HyperCom), is faster and handles bigger samples.

| | pure R (this branch) | main version (`master`) |
|---|---|---|
| Install | any computer with R 4.2 or later | ready-built on Windows and macOS with current R; build tools otherwise |
| Bundled data with p-values | about a minute (4 processes) | 2 s (8 threads) |
| Cells per sample | up to about 10,000 | up to about 100,000 (64 GB) |
| Memory per pair of cells in a sample | 79 bytes per R process | 17.6 bytes (5.5 in memory mode) |
| Results | the same, to about 1e-15 | |

Both are the package `HyperCom`: installing one replaces the other. Functions and arguments are the same, except that this version
has no `exact.priority` (every pair is scored exactly) and no `threads` in `generate_infMat()` and `score_hypercom()`.

## Install

```r
install.packages("pak")
pak::pak("pawelab/HyperCom@base-R")
```

Check that it works (takes a few seconds):

```r
library(HyperCom)
example(HyperCom)
```

## Run

**Check these four things first:**

1. **Database matches your species:** `cellchat.hs.lrdb` for human, `cellchat.mm.lrdb` for mouse. Genes are gene symbols (`TGFB1`), not
   Ensembl IDs.
2. **`sample`** names the metadata column that holds the sample (library or replicate). Cells are compared only within their sample.
   One sample? Add `metadata$sample <- "all"`. **`group`** (optional) names the condition column: pairs whose cells concentrate in
   some conditions rank higher. There is one ranking for the whole dataset, not one per condition.
3. **Counts are raw or normalized, not scaled.**
4. **Cell names are unique.** When you merge samples, add the sample name to each barcode.

Then pick the line that matches your data:

```r
library(HyperCom)

# Seurat object (needs install.packages("SeuratObject"))
res <- HyperCom_seurat(seurat, cellchat.hs.lrdb, sample="sample", group="condition")

# SingleCellExperiment (needs pak::pak("SingleCellExperiment"), which installs SummarizedExperiment too)
res <- HyperCom_sce(sce, cellchat.hs.lrdb, sample="sample", group="condition")

# Counts matrix (genes x cells) and a metadata table (one row per cell)
res <- HyperCom(cellchat.hs.lrdb, counts, metadata, sample="sample", group="condition")
```

Then look at the results:

```r
res$priority[["1"]]    # ligand-receptor pairs, ranked (one table per weight; the default weight is 1)
scores <- score_hypercom(res$infMat, res$metadata, ligand="TGFB1", receptors=c("TGFBR1", "TGFBR2"))    # score per cell
```

- **P-values:** add `significance=TRUE, seed=1`. `p = 0` means p < 0.001 (no background pair ranked as high): report it as
  p < 0.001. Its `p.adj` is then 0 too and is also only a bound.
- **Output files:** every result is also saved to `outdir` (default `HyperCom/` in the working directory; `outdir=NULL` writes none).
  Running again with the same `outdir` overwrites its files.
- **Cores:** HyperCom starts up to 4 R processes (one per physical core), at a lower priority than your other programs, so the
  computer stays responsive. More gain little in R. Each process holds its own copy of the data, so HyperCom starts fewer when
  they would not fit in memory, and says so. `threads=8` asks for 8, `threads=1` runs in your R session alone. On a cluster it uses only the job's cores; give the job's memory
  with `memory.gb=`.
- **Run stopped?** Run the same call again. It continues where it stopped.

## Run time and memory

Time and memory depend on the number of cells **per sample** (they grow with its square). Splitting the same cells into more
samples reduces both, but changes which cells are compared. One sample, with p-values, 4 R processes (the default):

| Cells in the sample | Time | Peak memory |
|---|---|---|
| 1,000 | 1.4 min | 1.8 GB |
| 3,000 | 20 min | 4.0 GB |
| 10,000 | 3.1 h | 32 GB |

- **Memory:** each R process needs about 79 bytes per pair of cells in the same sample (one sample of 10,000 cells: 10^8 pairs,
  7.9 GB), plus 0.4 GB; the run needs about 1.1 GB more. A run that cannot fit even in one process stops at the start, giving the
  estimate: use the main version.
- **Pairs:** time also grows with the number of ligand-receptor pairs expressed.

<sub>Measured with `HyperCom()` on one sample of cells resampled from the bundled mouse heart data (counts thinned to 80%), with the bundled mouse database (about 1,000 pairs expressed), `significance=TRUE`, one run each, on a desktop (Intel Core i9-12900K, 128 GB). Memory is the peak of all processes together. GB = 10<sup>9</sup> bytes.</sub>

## Learn more

- **Tutorial:** the <a href="https://raw.githubusercontent.com/pawelab/HyperCom/base-R/vignettes/articles/HyperCom_Overview.pdf"
  download>overview vignette</a> (source: [vignettes/articles/HyperCom_Overview.qmd](vignettes/articles/HyperCom_Overview.qmd))
  explains every step, the scores and how to plot them.
- **Help:** `?HyperCom`, `?HyperCom_seurat`, `?HyperCom_sce`, `?score_hypercom`, `?prioritize_lr` (how pairs are ranked).
- **What changed:** [NEWS.md](NEWS.md).
- **Questions or bugs:** [open an issue](https://github.com/pawelab/HyperCom/issues).
- **Cite:** `citation("HyperCom")`.
