# HyperCom

HyperCom is the R package implementation of the method described in [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1) on bioRxiv. It is a hypergraph-based method to understand cell-cell communication at a cluster-free single-cell resolution.

HyperCom provides the following hypergraph-based features described in the [overview vignette](vignettes/HyperCom_Overview.Rmd):

1. Per Cell Ligand-Receptor Pair Scoring

2. Per Dataset Ligand-Receptor Pair Prioritization

## Method Overview
![HyperCom Overview](images/MethodOverview.jpg)

## Installation

### Prerequisites

- **R 4.1 or newer.**
- **Build tools.** HyperCom contains C++ code that is compiled when it is installed from GitHub. Set the tools up once:
  - **Windows:** install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your version of R (Rtools 4.4 for R 4.4,
    Rtools 4.5 for R 4.5) with its default options, then restart R.
  - **macOS:** run `xcode-select --install` in the Terminal, then install GNU Fortran for your version of R from
    [mac.r-project.org/tools](https://mac.r-project.org/tools/) (HyperCom uses R's linear algebra library, which needs it).
  - **Linux:** install R's development package: `sudo apt install r-base-dev` (Ubuntu, Debian) or `sudo dnf install R-devel`
    (Fedora, RHEL).

  To check the tools from R: `install.packages("pkgbuild")`, then `pkgbuild::check_build_tools(debug=TRUE)`.
- **R packages** are installed automatically: dplyr (1.1.0 or newer), Matrix, parallelly (1.23.0 or newer), Rcpp and tidyr.
  `HyperCom_seurat()` also needs SeuratObject: `install.packages("SeuratObject")`.

### Install

This version lives on the `master` branch:

```r
install.packages("pak")
pak::pak("pawelab/HyperCom")
```

To check the installation, `library(HyperCom)` then `example(HyperCom)` runs a small analysis of the bundled data in a few seconds.

If the installation stops with **`'make' not found`** (Windows), Rtools is missing or does not match your version of R. If it
stops with **`library 'quadmath' not found`** or a missing **`gfortran`** (macOS), install GNU Fortran as above.

### Memory and run time

Cells are compared within their sample, so memory and time grow with the number of cells per sample squared: memory is about 16
bytes per pair of cells in the same sample, and the same cells split into more samples need less of both. Time also grows with the
number of ligand-receptor pairs. Measured with `HyperCom_seurat()` on one sample of a Xenium 5K human ovarian cancer slide, the
CellChat human database (about 700 ligand-receptor pairs expressed) and p-values (`significance=TRUE`), `threads=8`, on 8 cores of
an Intel Xeon Gold 6426Y node with a 64 GiB job (median of 5 runs; 3 at 60,000 cells):

| Cells in the sample | Time | Peak memory |
|---|---|---|
| 1,000 | 6 s | 0.5 GiB |
| 3,000 | 42 s | 0.8 GiB |
| 10,000 | 8 min | 2.4 GiB |
| 30,000 | 1.3 h | 15 GiB |
| 60,000 | 5.9 h | 57 GiB |

60,000 cells is about the largest single sample that fits in 64 GiB. The p-values (1000 background pairs) take about 40% of the
time; without them a run takes about 60% as long. A run that stops while pairs are scored (e.g. the time limit of a cluster job)
continues where it stopped when the same call is run again (see the vignette).

### Pure R version

The R-only implementation is kept on the [`base-R` branch](https://github.com/pawelab/HyperCom/tree/base-R) as HyperCom 0.0.1: the
original code of the paper (branch [`handwritten`](https://github.com/pawelab/HyperCom/tree/handwritten)) with the same bug fixes as
this version (bad input stops with a message, factor inputs and reordered metadata give correct results, `seed` for reproducible
p-values, `verbose`, results returned by `HyperCom()`; its own NEWS lists them). It no longer needs mltools, which was archived on
CRAN, so it installs from CRAN packages only (collapse, doParallel, foreach, Rfast and SeuratObject besides the ones above) and
needs no build tools on Windows and macOS. It is much slower and needs far more memory: the bundled data takes about 8 minutes on
12 cores instead of seconds, and since it holds a dense matrix over all cells and genes of all samples (one copy per worker) it runs out
of memory around 20,000 cells. It has no `threads` or `cache`. Its results agree with this version to about 1e-15 (exactly with
`HyperCom(..., infmat.method="lu")`), and on valid input they equal the original code's except where that was wrong.

```r
pak::pak("pawelab/HyperCom@base-R")
```

Both versions are the package `HyperCom`, so installing one replaces the other. The R-only version keeps the `cl` argument (a PSOCK
cluster, left running when given) instead of `threads`. With `parallel=TRUE` and no `cl` it starts one R process per core (inside an
SGE parallel job, the job's workers through qrsh): on machines with many cores, give it a smaller cluster, e.g.
`cl=parallelly::makeClusterPSOCK(8)`.

## Usage
Counts are a gene x cell matrix (raw or normalised, not scaled) and the metadata has one row per cell, named by its barcode:
```r
library(HyperCom)
res <- HyperCom(cellchat.mm.lrdb, alexanian2021.counts, alexanian2021.metadata, sample="sample", group="condition", parallel=TRUE)
head(res$priority[["1"]])
scores <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")
```
Every step is also written to `outdir` (default `HyperCom/`). See `?HyperCom` and the vignette; cite with `citation("HyperCom")`.
