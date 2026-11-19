# HyperCom

**See which single cells send and receive each ligand-receptor signal, with no clustering.**

HyperCom gives every cell a sender and receiver score for each ligand-receptor pair, then ranks the pairs that drive communication
in your data, across the whole dataset or by condition. Start from a Seurat object, a SingleCellExperiment or a counts matrix.

Paper: [HyperCom: a hypergraph-based method to infer cell-resolved cell-cell communication](https://www.biorxiv.org/content/10.64898/2026.09.30.755742v1)
(bioRxiv, 2026).

![HyperCom method overview](images/MethodOverview.jpg)

## Install

```r
install.packages("pak")
pak::pak("pawelab/HyperCom")
```

Check that it works (takes a few seconds):

```r
library(HyperCom)
example(HyperCom)
```

HyperCom compiles C++ when it installs, so it needs R 4.1 or newer and build tools. Set them up once:

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
- **Same results** to about 1e-15 (exactly the same with `HyperCom(..., infmat.method="lu")` in this version).
- **Uses `cl` instead of `threads`** for parallel runs. On machines with many cores, pass a small cluster, e.g.
  `cl=parallelly::makeClusterPSOCK(8)`.
- **It replaces this version:** both are the package `HyperCom`, so installing one removes the other.

</details>

## Run

Pick the line that matches your data:

```r
library(HyperCom)

# Seurat object (needs install.packages("SeuratObject"))
res <- HyperCom_seurat(seurat, cellchat.hs.lrdb, sample="sample", group="condition", parallel=TRUE)

# SingleCellExperiment (needs pak::pak("SingleCellExperiment"))
res <- HyperCom_sce(sce, cellchat.hs.lrdb, sample="sample", group="condition", parallel=TRUE)

# Counts matrix (genes x cells) and a metadata table (one row per cell)
res <- HyperCom(cellchat.hs.lrdb, counts, metadata, sample="sample", group="condition", parallel=TRUE)
```

Then look at the results:

```r
res$priority[["1"]]    # ligand-receptor pairs, ranked
scores <- score_hypercom(res$infMat, res$metadata, ligand="TGFB1", receptors=c("TGFBR1", "TGFBR2"))    # score per cell
```

**Check these four things before you run:**

1. **Database matches your species:** `cellchat.hs.lrdb` for human, `cellchat.mm.lrdb` for mouse. Genes are gene symbols (`TGFB1`), not
   Ensembl IDs.
2. **`sample`** names the metadata column that holds the sample (library or replicate). Cells are compared only within their sample.
   **`group`** (optional) names the condition column; use it to rank pairs by condition.
3. **Counts are raw or normalized, not scaled.**
4. **Cell names are unique.** When you merge samples, add the sample name to each barcode.

- **P-values:** add `significance=TRUE, seed=1`.
- **Output files:** every result is also saved to `outdir` (default `HyperCom/`).
- **Run stopped?** Run the same call again. It continues where it stopped.

## Run time and memory

Time and memory depend on the number of cells **per sample** (they grow with its square). Splitting the same cells into more
samples reduces both. One sample, with p-values, 8 cores:

| Cells in the sample | Time | Peak memory |
|---|---|---|
| 1,000 | 6 s | 0.5 GiB |
| 3,000 | 42 s | 0.8 GiB |
| 10,000 | 8 min | 2.4 GiB |
| 30,000 | 1.3 h | 15 GiB |
| 60,000 | 5.9 h | 57 GiB |

Time also grows with the number of ligand-receptor pairs expressed (about 700 here). P-values took about 40% of the time at
every size from 3,000 cells up, so runs without them take about 60% as long; with fewer pairs, p-values are a larger share.
60,000 cells is about the largest sample that fits in 64 GiB.

<sub>Measured with `HyperCom_seurat()` on one sample of a Xenium 5K human ovarian cancer slide: CellChat human database (about 700
pairs expressed), `significance=TRUE`, `threads=8` on an Intel Xeon Gold 6426Y node with a 64 GiB job; median of 5 runs (3 at
60,000 cells).</sub>

## Learn more

- **Tutorial:** the [overview vignette](vignettes/HyperCom_Overview.Rmd) explains every step, the scores and how to plot them.
- **Help:** `?HyperCom`, `?HyperCom_seurat`, `?HyperCom_sce`, `?score_hypercom`.
- **Questions or bugs:** [open an issue](https://github.com/pawelab/HyperCom/issues).
- **Cite:** `citation("HyperCom")`.
