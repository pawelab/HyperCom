# HyperCom

HyperCom is the R package implementation of the method described in [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1) on bioRxiv. It is a hypergraph-based method to understand cell-cell communication at a cluster-free single-cell resolution.

HyperCom provides the following hypergraph-based features described in the [overview vignette](vignettes/articles/HyperCom_Overview.qmd):

1. Per Cell Ligand-Receptor Pair Scoring

2. Per Dataset Ligand-Receptor Pair Prioritization

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

This version lives on the `speed-fixes` branch:

```r
install.packages("pak")
pak::pak("pawelab/HyperCom@speed-fixes")
```

Without `@speed-fixes`, pak installs the default branch, which holds the R-only version described below. To check the
installation, `library(HyperCom)` then `example(HyperCom)` runs a small analysis of the bundled data in a few seconds.

If the installation stops with **`'make' not found`** (Windows), Rtools is missing or does not match your version of R. If it
stops with **`library 'quadmath' not found`** or a missing **`gfortran`** (macOS), install GNU Fortran as above.

### Memory and run time

Memory grows with the number of cells per sample squared, since cells are compared within their sample. With `parallel=TRUE`
(one thread per physical core), 1000 cells take seconds; 20,000 cells in 3 samples took about 1.5 minutes and 3.5 GB on 16
threads, and 50,000 cells in 3 samples about 28 minutes and 19 GB.

### Pure R version

The original, R-only implementation (HyperCom 0.0.0.9000) is kept on the [`base-R` branch](https://github.com/pawelab/HyperCom/tree/base-R).
It needs no build tools on Windows and macOS, but it is much slower, needs far more memory (it runs out of memory around 20,000
cells) and lacks the fixes listed in [NEWS](NEWS.md). Its results agree with this version to about 1e-15 (exactly with
`HyperCom(..., infmat.method="lu")`).

```r
pak::pak("pawelab/HyperCom@base-R")
```

Both versions are the package `HyperCom`, so installing one replaces the other. In the R-only version, `HyperCom()` changes the
working directory to `outdir`, and `parallel=TRUE` starts one R process per core: on machines with many cores, give it a smaller
cluster, e.g. `cl=parallelly::makeClusterPSOCK(8)`.

## Usage
Counts are a gene x cell matrix (raw or normalised, not scaled) and the metadata has one row per cell, named by its barcode:
```r
library(HyperCom)
res <- HyperCom(cellchat.mm.lrdb, alexanian2021.counts, alexanian2021.metadata, sample="sample", group="condition", parallel=TRUE)
head(res$priority[["1"]])
scores <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")
```
Every step is also written to `outdir` (default `HyperCom/`). See `?HyperCom` and the vignette; cite with `citation("HyperCom")`.

## Method Overview
![HyperCom Overview](images/MethodOverview.jpg)
