# AGENTS.md: HyperCom

For AI agent. Use HyperCom right. Installed copy: `system.file("AGENTS.md", package="HyperCom")`.

## First: which version

`packageVersion("HyperCom")`. Same name, different API:

- 0.2.x: this file. Threads, cache, compiled code.
- 0.1.x (branch `spatial`, was `v1-fast`): has `cores=`, `infMat.cell=`, `cell.info=`, `method="auto"`, spatial. 0.2.x has none.
- 0.0.1 (branch `base-R`): R only, same bug fixes. Has `cl=`, `seed=`, `verbose=`; no `threads`, `cache`, `infmat.method`.
- 0.0.0.9000 (branch `handwritten`, b999a01): original code of paper. Bugs fixed in 0.0.1 still there. Needs mltools (archived on CRAN).

Read for this version: `?HyperCom`, `news(package="HyperCom")`, vignette `HyperCom_Overview` (in repo `vignettes/` if not built).

## Run

```r
library(HyperCom)
res <- HyperCom(cellchat.mm.lrdb, counts, metadata, sample="sample", group="condition",
                significance=TRUE, seed=1, parallel=TRUE, outdir="HyperCom")
res$priority[["1"]]
s <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4", parallel=TRUE)
```

- Seurat: `HyperCom_seurat(seurat, lrdb, sample=, assay="RNA", layer="counts", ...)`. Need SeuratObject.
- SingleCellExperiment: `HyperCom_sce(sce, lrdb, sample=, assay="counts", ...)`. Columns from `colData`. `colnames(sce)` = cell names. Need SummarizedExperiment (Bioconductor).
- Install check: `example(HyperCom)`. Seconds.
- Default 1 thread. Slow. Set `parallel=TRUE` (1 per physical core) or `threads=n`.
- Set `outdir` yourself. Default `HyperCom/` in working dir.
- No `cl=`. Gone. Stop with "unused argument".

## Input: check before run

Most failed run = bad input.

- `counts`: gene x cell. Genes = row names. Cells = column names. Sparse fine.
- Raw or normalised fine (only > 0 used). Never scaled: negative value stop run.
- Cell names unique. Merge samples: prefix barcode with sample.
- `metadata`: one row per cell. Row names = cell names of counts. Cell not in metadata: dropped, warn.
- `sample` = library or replicate, not cell type. Cell only compared inside own sample. Each sample need 2+ cells in graph.
- `group` = condition. Optional. No NA in graph cells. Value never `"gene"`.
- `lrdb` species = counts species. `cellchat.mm.lrdb` mouse (`Tgfb1`). `cellchat.hs.lrdb` human (`TGFB1`). Wrong species: stop.
- Own `lrdb`: column 1 = interaction name, 2 = ligand, rest = receptor subunits (NA pad).

## Size: estimate before run

- Memory about 16 bytes x sum over samples of (cells in sample)^2. More samples, same size: cheap.
- Split sample to save memory? Change which cells compared, so change results. Ask user first.
- Measured, one sample, 8 cores (Xeon Gold 6426Y), all CellChat human pairs (~700 expressed), with p-values: 1k 6 s 0.5 GiB. 10k 8 min
  2.4 GiB. 30k 1.3 h 15 GiB. 60k 5.9 h 57 GiB (about max one sample for 64 GiB). Time grows with LR pairs too.
- `significance=TRUE` add 1000 background pairs: about 40% of time. Without: about 60% as long.

## Long run, cluster

- Run stop (time limit, crash)? Run same call again. Resume from `<outdir>/6.priority/cache`. Lose one batch max.
- Keep `seed`, inputs, `infmat.method` same. Other inputs: old cache moved to `cache.old`, warn.
- Same node type for same last digits (about 1e-13 drift across CPU).
- Do not delete cache mid-run. Removed by itself at end.

## Output

- `HyperCom()` return list invisible. Assign it: `priority` (list, name = weight, e.g. `"1"`), `metadata`, `infMat`, `adj`, `lrs`.
- Files: `0.counts.rds`, `1.lrs.csv`, `2.hyperedges.csv`, `3.adj.rds`, `4.infMat.rds`, `5.metadata.csv`, `6.priority/priority.<weight>.csv`.
- `.rds` of `adj` / `infMat`: `library(HyperCom)` before `readRDS`.
- `adj`, `infMat` = factored objects, not matrix. `dim`, `x[i, j]`, `dimnames` work. `as.matrix()` build full nodes^2 matrix: huge. Slice instead.
- Priority table: `lrdb` columns + `priority`, best first. With significance: `p`, `p.adj` (BH).
- `p` step 0.001. `p = 0` mean below 0.001. Report `p < 0.001`, not 0.
- Score table: `cell`, `sender`, `receiver` (0 to 1), `dif` = receiver - sender, `avg`, `max`, `hypercom` (< 0 sender, > 0 receiver).
- NaN score = cell no path to pair. Normal. Cell with no LR gene: not in graph, no row.
- Receptor complex: give all subunits. `receptors=c("Tgfbr1", "Tgfbr2")`. Scores averaged.
- `weights` = exponent of cell fraction term. Below 1 favor rare pairs. Vector give one table each.

## Reproduce

- `seed=` make p-values same. Session RNG restored after.
- `infmat.method="woodbury"` default, fast. `"lu"` = exact numbers of the R-only version (branch `base-R`), full matrix memory. Differ about 1e-15.

## Change package (source repo)

- Outputs must match the original code (branch `handwritten`, b999a01; fixture `tests/testthat/fixtures/handwritten.rds`). Tests check bit for bit. Never loosen tolerance to pass.
- Tests: testthat. `NOT_CRAN=true` run all. `HYPERCOM_REFERENCE=<golden dir>` add full-data test, else skip.
- Docs: roxygen in `R/`. Regenerate `man/`, `NAMESPACE`. Never hand edit. C++ signature change: `Rcpp::compileAttributes()`.
- `R CMD check --as-cran`: Status OK.
- User-visible change: add to `NEWS.md`.
