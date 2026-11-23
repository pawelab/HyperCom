# HyperCom

[![r-universe](https://pawelab.r-universe.dev/badges/HyperCom)](https://pawelab.r-universe.dev/HyperCom)

**See which single cells send and receive each ligand-receptor signal, with no clustering.**

HyperCom gives every cell a sender and receiver score for each ligand-receptor pair, then ranks the pairs that drive communication
in your data. With a condition column, pairs whose cells concentrate in some conditions rank higher. Start from a Seurat object, a
SingleCellExperiment or a counts matrix.

Paper: [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1)
(bioRxiv, 2026).

![HyperCom method overview](https://raw.githubusercontent.com/pawelab/HyperCom/master/images/MethodOverview.jpg)

## Install

```r
install.packages("HyperCom", repos = c("https://pawelab.r-universe.dev", "https://cloud.r-project.org"))
```

Check that it works (takes a few seconds):

```r
library(HyperCom)
example(HyperCom)
```

| Your setup | What happens |
|---|---|
| Windows or macOS, current or previous R release (4.6 or 4.5 as of November 2026) | Installs a ready-built binary in seconds. Nothing else needed. |
| **Older R (4.2 to 4.4), or Linux** | **Builds from source: set up build tools first (below).** |
| R 4.1 or older | Not supported: update R. |

Not sure which R you have? Run `R.version.string`.

### Older R or Linux: set up build tools first

HyperCom compiles C++ when it installs from source. Set up build tools once, then run the same `install.packages()` line above.
On older R, some dependencies also build from source, so the first install can take 10 minutes or more. Installing straight from
GitHub (`pak::pak("pawelab/HyperCom")`) also builds from source.

| System | Build tools |
|---|---|
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version (Rtools 4.4 for R 4.4), then restart R |
| macOS | `xcode-select --install` in the Terminal, then GNU Fortran for your R version from [mac.r-project.org/tools](https://mac.r-project.org/tools/) |
| Linux | `sudo apt install r-base-dev` (Ubuntu, Debian) or `sudo dnf install R-devel` (Fedora, RHEL) |

**Install failed?**
- `'make' not found` (Windows): Rtools is missing or doesn't match your R version.
- `library 'quadmath' not found` or `gfortran` missing (macOS): install GNU Fortran (see the table).
- Not sure? Run `install.packages("pkgbuild")`, then `pkgbuild::check_build_tools(debug=TRUE)`.

<details>
<summary><b>Pure R version (no build tools)</b></summary>

Branch [`base-R`](https://github.com/pawelab/HyperCom/tree/base-R) (HyperCom 0.0.1) is the paper's original R code (branch
[`handwritten`](https://github.com/pawelab/HyperCom/tree/handwritten)) with the same bug fixes as this version. It installs from
CRAN packages only, so Windows and macOS need no build tools:

```r
pak::pak("pawelab/HyperCom@base-R")
```

- **Much slower:** about 8 minutes for the bundled data on 12 cores, instead of seconds.
- **Much more memory:** it runs out of memory at about 20,000 cells.
- **Same results** as this version, to about 1e-15 without p-values. With p-values, this version screens the real pairs:
  priorities typically within ~1e-9 of exact, up to ~1e-7 on the bundled data (1.5e-8 measured on the full bundled run); p-values
  and order exact; `exact.priority=TRUE` for exact values.
- **Uses `cl` instead of `threads`** for parallel runs. On machines with many cores, pass a small cluster, e.g.
  `cl=parallelly::makeClusterPSOCK(8)`.
- **It replaces this version:** both are the package `HyperCom`, so installing one removes the other.

</details>

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
- **Threads:** HyperCom uses one thread per physical core; `threads=4` uses 4.
- **Run stopped?** Run the same call again. It continues where it stopped.

## Run time and memory

Time and memory depend on the number of cells **per sample** (they grow with its square). Splitting the same cells into more
samples reduces both, but changes which cells are compared. One sample, with p-values:

| Cells in the sample | Time | Peak memory | Threads | HyperCom |
|---|---|---|---|---|
| 1,000 | 6 s | 0.5 GB | 8 | 0.2.x |
| 3,000 | 42 s | 0.8 GB | 8 | 0.2.x |
| 10,000 | 8 min | 2.6 GB | 8 | 0.2.x |
| 30,000 | 1.3 h | 16 GB | 8 | 0.2.x |
| 60,000 | 5.9 h | 61 GB | 8 | 0.2.x |
| 100,000 | 8.2 h | 52 GB | 16 | 0.3.0, memory mode |

- **Memory mode:** when exact scoring does not fit in memory, HyperCom 0.3.0 switches to a mode that needs about a third of the
  memory (priorities within ~1e-11 of exact scoring). About 100,000 cells in one sample fit in 64 GB. On a cluster, give the job's
  memory: `memory.gb=64`.
- **Pairs:** time also grows with the number of ligand-receptor pairs expressed (about 700 here).
- **P-values** took about 40% of the time in the 0.2.x runs. With p-values, 0.3.0 also screens the real pairs, in the same time or
  less: priorities typically within ~1e-9 of exact, up to ~1e-7 on the bundled data (1.5e-8 measured on the full bundled run);
  p-values and order exact; `exact.priority=TRUE` for exact values.

<sub>Measured with `HyperCom_seurat()` on one sample of a Xenium 5K human ovarian cancer slide, with the CellChat human database
(about 700 pairs expressed) and `significance=TRUE`. Up to 60,000 cells: HyperCom 0.2.x on a cluster node (Intel Xeon Gold 6426Y,
64 GB job), median of 5 runs (3 at 60,000 cells). 100,000 cells: HyperCom 0.3.0 on a desktop (Intel Core i9-12900K, 128 GB,
`memory.gb=64`), one run. GB = 10<sup>9</sup> bytes.</sub>

## Learn more

- **Tutorial:** the [overview vignette](https://pawelab.r-universe.dev/articles/HyperCom/HyperCom_Overview.html) (source:
  [vignettes/HyperCom_Overview.Rmd](vignettes/HyperCom_Overview.Rmd)) explains every step, the scores and how to plot them. After
  installing: `vignette("HyperCom_Overview", package="HyperCom")`.
- **Help:** `?HyperCom`, `?HyperCom_seurat`, `?HyperCom_sce`, `?score_hypercom`, `?prioritize_lr` (how pairs are ranked).
- **Questions or bugs:** [open an issue](https://github.com/pawelab/HyperCom/issues).
- **Cite:** `citation("HyperCom")`.
