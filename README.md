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

Branch [`base-R`](https://github.com/pawelab/HyperCom/tree/base-R) (HyperCom 0.3.1.1) is this version written in R only: the same
functions, arguments and checks. It installs from CRAN packages only, so it needs no build tools on any system:

```r
pak::pak("pawelab/HyperCom@base-R")
```

- **Slower:** about a minute for the bundled data with p-values, instead of seconds. One sample, with its default 4 processes:
  20 minutes and 4 GB for 3,000 cells, 3.1 hours and 32 GB for 10,000.
- **More memory:** about 79 bytes per pair of cells in the same sample for each of its R processes, instead of 17.6; practical up
  to about 10,000 cells per sample.
- **Same results**, to about 1e-15 (bit for bit with the paper's code under `infmat.method="lu"`). It scores every pair exactly,
  so it has no `exact.priority`; with p-values this version's screened priorities differ from it by typically ~1e-9, up to ~1e-7
  on the bundled data (p-values and order are the same).
- **`threads` counts R processes**, at most 4 by default; `generate_infMat()` and `score_hypercom()` take no `threads`.
- **It replaces this version:** both are the package `HyperCom`, so installing one removes the other.

</details>

## Run

**Check these four things first:**

1. **Database matches your species:** `cellchat.hs.lrdb` for human, `cellchat.mm.lrdb` for mouse. Genes are gene symbols (`TGFB1`), not
   Ensembl IDs.
2. **`sample`** (optional) names the metadata column that holds the sample (library or replicate). Cells are compared only within
   their sample; without `sample`, all cells are one sample. **`group`** (optional) names the condition column: pairs whose cells
   concentrate in some conditions rank higher. There is one ranking for the whole dataset, not one per condition.
3. **No negative expression values:** raw or normalized counts are fine; scaled data is not.
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

# One sample, no conditions: counts alone
res <- HyperCom(cellchat.hs.lrdb, counts)
```

Then look at the results:

```r
res$priority[["1"]]    # ligand-receptor pairs, ranked (one table per weight; the default weight is 1)
scores <- score_hypercom(res$infMat, res$metadata, ligand="TGFB1", receptors=c("TGFBR1", "TGFBR2"))    # score per cell
```

- **Cell scores only:** `prioritize=FALSE` skips ranking the pairs, which takes most of the run time; then `score_hypercom()`
  scores the pairs you choose.
- **P-values:** add `significance=TRUE, seed=1`. `p = 0` means p < 0.001 (no background pair ranked as high): report it as
  p < 0.001. Its `p.adj` is then 0 too and is also only a bound.
- **Output files:** every result is also saved to `outdir` (default `HyperCom/` in the working directory; `outdir=NULL` writes none).
  Running again with the same `outdir` overwrites its files.
- **Threads:** HyperCom uses one thread per physical core, at a lower priority than your other programs, so the computer stays
  responsive; `threads=4` uses 4, `threads=1` one. On a cluster it uses only the job's cores.
- **Run stopped?** Run the same call again. It continues where it stopped.

## Run time and memory

Time and memory depend on the number of cells **per sample** (they grow with its square). Splitting the same cells into more
samples reduces both, but changes which cells are compared. One sample, with p-values, 8 threads:

| Cells in the sample | Time | Memory needed |
|---|---|---|
| 1,000 | 4 s | 0.5 GB |
| 3,000 | 34 s | 0.7 GB |
| 10,000 | 6 min | 1.5 GB |
| 30,000 | 53 min | 6.4 GB |
| 60,000 | 3.6 h | 21 GB |
| 100,000 | 10.5 h | 53 GB |

- **More memory, faster:** HyperCom fits its scoring to the memory it has and is up to 1.4× faster when it has more: 10,000
  cells in 4.4 min with 2.7 GB, 30,000 in 38 min with 15 GB, 60,000 in 2.8 h with 57 GB. It uses 90% of the computer's memory,
  or on Linux of the cluster job's limit; on a cluster, give the job's memory: `memory.gb=64`.
- **Pairs:** time also grows with the number of ligand-receptor pairs expressed (about 700 here).
- **P-values** (1,000 background pairs) take about two thirds of the time. With p-values, HyperCom also screens the real pairs:
  priorities typically within ~1e-9 of exact, up to ~1e-7 on the bundled data (1.5e-8 measured on the full bundled run); p-values
  and order exact; `exact.priority=TRUE` for exact values.

<sub>Measured with `HyperCom_seurat()` on one sample of a Xenium 5K human ovarian cancer slide, with the CellChat human database
(about 700 pairs expressed) and `significance=TRUE`: HyperCom 0.3.1 on a cluster node (Intel Xeon Gold 6326, 8 cores), median of 5
runs (3 at 60,000 and 100,000 cells). "Memory needed": peak memory with `memory.gb` just above HyperCom's own estimate; "more
memory": with no `memory.gb` limit, in a 64 GiB job. GB = 10<sup>9</sup> bytes.</sub>

## Learn more

- **Tutorial:** the [overview vignette](https://pawelab.r-universe.dev/articles/HyperCom/HyperCom_Overview.html) (source:
  [vignettes/HyperCom_Overview.Rmd](vignettes/HyperCom_Overview.Rmd)) explains every step, the scores and how to plot them. After
  installing: `vignette("HyperCom_Overview", package="HyperCom")`.
- **Help:** `?HyperCom`, `?HyperCom_seurat`, `?HyperCom_sce`, `?score_hypercom`, `?prioritize_lr` (how pairs are ranked).
- **Questions or bugs:** [open an issue](https://github.com/pawelab/HyperCom/issues).
- **Cite:** `citation("HyperCom")`.
