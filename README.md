# HyperCom

HyperCom is the R package implementation of the method described in [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1) on biorXiv. It is a hypergraph based method to understand cell-cell communication at a cluster free single-cell resolution.

HyperCom provides the following hypergraph based features described in the <a href="https://raw.githubusercontent.com/pawelab/HyperCom/master/vignettes/articles/HyperCom_Overview.pdf" download>overview vignette</a>:

1. Per Cell Ligand-Receptor Pair Scoring

2. Per Dataset Ligand-Receptor Pair Prioritization

## Installation
HyperCom needs R 4.1 or newer. It contains C++ code, so installing it from GitHub compiles that code, which needs build tools.
Set them up once for your system:

- **Windows:** install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your version of R (e.g. Rtools 4.4 for R 4.4).
- **macOS:** run `xcode-select --install` in the Terminal, then install GNU Fortran for your version of R from
  [mac.r-project.org/tools](https://mac.r-project.org/tools/).
- **Linux:** install R's development package (Ubuntu/Debian: `sudo apt install r-base-dev`).

Then, in R:
```r
install.packages("pak")
pak::pak("pawelab/HyperCom")
```

`HyperCom_seurat` also needs the SeuratObject package (`install.packages("SeuratObject")`).

## Method Overview
![HyperCom Overview](images/MethodOverview.jpg)
