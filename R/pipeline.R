# The whole pipeline: HyperCom() on a counts matrix and cell metadata, HyperCom_seurat() on a Seurat object, HyperCom_sce() on a
# SingleCellExperiment

#' Run the HyperCom pipeline
#'
#' `HyperCom` returns the ligand-receptor pairs ranked by priority and the objects to score cells with
#'
#' This function takes a gene x cell counts matrix, its cell metadata (optional) and a ligand-receptor database and runs every step:
#' `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`, `generate_metadata`
#' and `prioritize_lr`. The results are returned and written to outdir, so they need not be read back: `res <- HyperCom(...)`, then
#' `res$priority[["1"]]` or `score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")`. Ranking the pairs
#' (`prioritize_lr`) takes most of the time; `prioritize=FALSE` skips it, for scoring chosen pairs with `score_hypercom` only.
#'
#' The input is checked first. Before the first file is written, the run stops when a cell that expresses a ligand or receptor gene
#' has a missing sample or group (or one named "gene"), and when the scoring cannot fit in memory (see memory.gb); it also stops when
#' no interaction of lrdb has all its genes expressed (e.g. lrdb of another species). Cells of counts that are not row names of
#' metadata are left out, with a warning, and so are cells and genes not connected to the rest of the hypergraph (see
#' `generate_hyperedges`). Priority files of other weights already in outdir are kept, with a warning.
#'
#' Time and memory grow with the square of the number of cells per sample (see the Memory section of `prioritize_lr`). If a run
#' stops while pairs are scored, run the same call again: it repeats the earlier steps (seconds) and continues the scoring from the
#' cache in `<outdir>/6.priority/cache` (see the Stopped runs section of `prioritize_lr`).
#' @inheritParams filter_counts
#' @inheritParams generate_metadata
#' @inheritParams prioritize_lr
#' @inheritParams generate_infMat
#' @param metadata A data frame containing cell metadata with row names being cells, or NULL (the default) when there is neither a
#' sample, group nor keep column: then every cell of counts is used, in one sample and one group
#' @param infmat.method A string naming the method of `generate_infMat`: "woodbury" (the default: fastest, least memory, and the only
#' one memory mode can use), "chol" or "lu" (a dense matrix of all cells and genes squared). The methods agree to about 1e-15
#' @param resume A boolean of whether to keep the pairs scored so far in `<outdir>/6.priority/cache`, so that a stopped run continues
#' where it stopped when the same call is run again; TRUE by default, ignored when outdir is NULL
#' @param prioritize A boolean of whether to rank the ligand-receptor pairs with `prioritize_lr`, the step that takes most of the time;
#' TRUE by default. FALSE stops after `generate_metadata`, with everything `score_hypercom` needs, and ignores weights, significance,
#' seed, exact.priority and resume
#' @param verbose A boolean of whether to report progress as messages: what each step found, the memory needed, the pairs scored with
#' the time left, and the top pair of each weight
#' @param outdir A string naming the output directory of the run ("HyperCom" in the working directory by default), or NULL to write
#' no files (and keep no cache to resume from)
#' @returns Invisibly, a list:
#'
#' - priority: the ranked pairs from `prioritize_lr`, one data frame per weight, named by weight (e.g. `res$priority[["1"]]`); NULL
#'   with `prioritize=FALSE`;
#' - metadata: cell (and gene) metadata from `generate_metadata`;
#' - infMat: the influence matrix from `generate_infMat`;
#' - adj: the adjacency from `generate_adjacency_matrix`;
#' - lrs: the interactions whose genes are all in the hypergraph.
#'
#' With outdir, the same results are written there: 0.counts.rds (the counts of the ligand and receptor genes), 1.lrs.csv,
#' 2.hyperedges.csv, 3.adj.rds, 4.infMat.rds, 5.metadata.csv and 6.priority/priority.<weight>.csv (not with `prioritize=FALSE`).
#' Read the .rds files after
#' `library(HyperCom)`.
#' @examples
#' # a few ligands keep the example quick; use the whole database for real data. A smaller
#' # database makes a smaller hypergraph (984 cells instead of 1,000), so the scores differ too
#' lrdb <- cellchat.mm.lrdb[cellchat.mm.lrdb$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' res <- HyperCom(lrdb, alexanian2021.counts, alexanian2021.metadata, sample="sample",
#'                 group="condition", outdir=file.path(tempdir(), "HyperCom"))
#' head(res$priority[["1"]])
#' list.files(file.path(tempdir(), "HyperCom"), recursive=TRUE)
#' scores <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")
#' @export
HyperCom <- function(lrdb, counts, metadata=NULL, sample=NULL, group=NULL, keep=NULL, weights=1, significance=FALSE, seed=NULL,
                     r=0.5, infmat.method="woodbury", exact.priority=FALSE, threads=NULL, memory.gb=NULL, resume=TRUE,
                     prioritize=TRUE, verbose=TRUE, outdir="HyperCom"){
  lrdb <- clean_lrdb(lrdb)
  if(!is.null(metadata) || !is.null(c(sample, group, keep))){
    check_columns(metadata, sample, group, keep)
  }
  check_weights(weights)
  check_flag(significance, "significance")
  check_flag(exact.priority, "exact.priority")
  check_flag(resume, "resume")
  check_flag(prioritize, "prioritize")
  check_flag(verbose, "verbose")
  check_memory(memory.gb)
  check_outdir(outdir)
  check_restart(r)
  infmat.method <- match_method(infmat.method)
  threads <- get_threads(threads)
  counts <- read_counts(counts, lrdb_genes(lrdb))
  if(is.null(metadata)){
    if(is.null(colnames(counts)) || anyDuplicated(colnames(counts))){
      stop("counts needs distinct cell names as column names")
    }
    metadata <- data.frame(row.names=colnames(counts))
  }
  in.metadata <- colnames(counts) %in% rownames(metadata)
  if(!any(in.metadata)){
    if(any(rownames(counts) %in% rownames(metadata))){
      stop("counts has cells as row names: give a gene x cell matrix, e.g. t(counts)")
    }
    stop("no cell of counts is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  if(prioritize && !is.null(outdir)){
    stale <- setdiff(list.files(file.path(outdir, "6.priority"), pattern="^priority\\..*\\.csv$"), priority_filename(weights))
    if(length(stale) > 0){
      warning("outdir already holds priorities of other weights, which are kept: ", paste(stale, collapse=", "))
    }
  }
  if(!all(in.metadata)){
    warning(sum(!in.metadata), " cells of counts are not in metadata and are left out")
    counts <- counts[, in.metadata, drop=FALSE]
  }

  counts <- filter_counts(counts, lrdb, outfile=NULL)
  lrs <- generate_lrs_table(counts, lrdb, outfile=NULL)
  if(nrow(lrs) == 0){
    stop("no interaction of lrdb has all its genes expressed in counts (", ncol(counts), " of its genes are): check that lrdb is of the ",
         "species of counts (cellchat.mm.lrdb mouse, cellchat.hs.lrdb human) and that counts is gene x cell")
  }
  check_run_cells(counts, metadata, sample, group, memory.gb, infmat.method, interaction_edges(lrs))
  if(!is.null(outdir)){
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    saveRDS(counts, step_file(outdir, "0.counts.rds"))
  }
  if(verbose){
    message(ncol(counts), " ligand and receptor genes expressed; ", nrow(lrs), " interactions of lrdb have all their genes expressed")
  }
  hyperedges <- generate_hyperedges(counts, lrs, outfile=step_file(outdir, "2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=step_file(outdir, "3.adj.rds"))
  lrs <- lrs[lrs_in_graph(lrs, rownames(adj)), , drop=FALSE]
  if(!is.null(outdir)){
    utils::write.csv(lrs, step_file(outdir, "1.lrs.csv"), row.names=FALSE)
  }
  cells <- sum(rownames(adj) %in% rownames(counts))
  rm(counts, hyperedges)
  infMat <- generate_infMat(adj, r=r, method=infmat.method, threads=threads, outfile=step_file(outdir, "4.infMat.rds"))
  if(verbose){
    message("Hypergraph of ", cells, " cells and ", nrow(adj) - cells, " genes; influence matrix made (method \"", infmat.method, "\")")
  }
  metadata <- generate_metadata(infMat, metadata, sample, group, keep, outfile=step_file(outdir, "5.metadata.csv"))
  priority <- NULL
  if(prioritize){
    priority <- prioritize_lr(lrs, adj, infMat, metadata, weights=weights, significance=significance, seed=seed,
                              exact.priority=exact.priority, threads=threads, memory.gb=memory.gb, resume=resume, verbose=verbose,
                              outdir=step_file(outdir, "6.priority"))
  }
  if(verbose){
    message("Done", if(!prioritize) " (pairs not ranked: prioritize=FALSE)", if(!is.null(outdir)) paste0(": results in ", outdir))
  }

  return(invisible(list(priority=priority, metadata=metadata, infMat=infMat, adj=adj, lrs=lrs)))
}

#' Name the file of a pipeline step
#'
#' `step_file` returns the path of a step's output file in the output directory, or NULL when no files are written
#' @param outdir A string naming the output directory, or NULL
#' @param name A string naming the file
#' @returns A string, or NULL
#' @noRd
step_file <- function(outdir, name){
  if(is.null(outdir)){
    return(NULL)
  }

  return(file.path(outdir, name))
}

#' Find the interactions whose genes are all in the hypergraph
#'
#' `lrs_in_graph` returns which interactions have all their genes among the nodes of the hypergraph
#' @param lrs A data frame of ligand-receptor interactions
#' @param nodes A vector of strings of the node names
#' @returns A logical vector, one per interaction
#' @noRd
lrs_in_graph <- function(lrs, nodes){
  genes <- as.matrix(lrs[, -1, drop=FALSE])

  return(rowSums(!(genes %in% nodes) & !is.na(genes)) == 0)
}

#' Count the interactions that are hyperedges of the adjacency
#'
#' `interaction_edges` returns how many interactions have at least two distinct genes, the only ones that join nodes
#' @param lrs A data frame of ligand-receptor interactions in the format of lrdb
#' @returns An integer
#' @noRd
interaction_edges <- function(lrs){
  genes <- as.matrix(lrs[, -1, drop=FALSE])
  distinct <- apply(genes, 1, function(g) length(unique(g[!is.na(g)])))

  return(sum(distinct >= 2))
}

#' Check the cells of a run before it starts
#'
#' `check_run_cells` stops early on cell labels and memory that would stop the run only after its influence matrix
#'
#' This function takes the filtered counts, the metadata and its sample and group columns, the memory.gb argument, the influence
#' method and the number of interactions (method "woodbury" factors the matrix unless the data have no more cells than interactions), and stops when a cell that expresses a ligand or receptor gene has a missing sample or group or
#' one named "gene" (check_labels), or when the grid kernels of these cells cannot fit in memory (kernel_type).
#' @param counts A matrix with row names being cells and column names being genes, from filter_counts
#' @param metadata A data frame containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata, or NULL
#' @param group A string that is the column name for group in the provided metadata, or NULL
#' @param memory.gb A number of gigabytes, Inf for no limit, or NULL for the computer's memory
#' @param infmat.method A string naming the method of generate_infMat
#' @param interactions An integer number of interactions with all genes expressed and at least two distinct genes
#' @returns NULL, invisibly
#' @noRd
check_run_cells <- function(counts, metadata, sample, group, memory.gb, infmat.method, interactions){
  cells <- rownames(counts)[Matrix::rowSums(counts > 0, na.rm=TRUE) > 0]
  factored <- infmat.method == "woodbury" && interactions < length(cells)
  rows <- match(cells, rownames(metadata))
  sample.vals <- label_values(metadata, sample, "all")[rows]
  group.vals <- label_values(metadata, group, "none")[rows]
  check_labels(sample.vals, group.vals)
  kernel_type(cell_pairs(table(sample.vals)), memory.gb, factored)

  return(invisible(NULL))
}

#' Get the counts matrix of a Seurat assay layer
#'
#' `get_seurat_layer` returns the matrix with dimnames for given assay and layer
#'
#' This function takes a Seurat object and returns one layer of one assay, joining split layers of v5 assays. Assays of other classes
#' (v3 assays, SeuratObject 4 and 5) are read from the slot of that name, as GetAssayData does in SeuratObject 5 (in SeuratObject 4
#' GetAssayData has no layer argument and returns the data slot).
#' @param seurat A seurat object
#' @param assay A string naming the assay to get data from
#' @param layer A string naming the layer to get data from
#' @returns A matrix with row names being genes and column names being cells
#' @noRd
get_seurat_layer <- function(seurat, assay="RNA", layer="counts"){
  if(!requireNamespace("SeuratObject", quietly=TRUE)){
    stop("reading a Seurat object needs the SeuratObject package")
  }
  a <- seurat[[assay]]

  if(!inherits(a, "Assay5")){
    if(!(layer %in% c("counts", "data", "scale.data"))){
      stop("layer \"", layer, "\" is not a layer of assay \"", assay, "\": use counts, data or scale.data")
    }
    m <- methods::slot(a, layer)
    if(!length(m)){
      stop("layer \"", layer, "\" is empty or absent in assay \"", assay, "\"")
    }
    return(m)
  }

  lyrs <- suppressWarnings(SeuratObject::Layers(a, search=layer))
  if(!length(lyrs)){
    stop("no layer matching \"", layer, "\" in assay \"", assay, "\"; available: ", paste(SeuratObject::Layers(a), collapse=", "))
  }
  if(length(lyrs) == 1){
    return(SeuratObject::LayerData(a, layer=lyrs))
  }

  return(SeuratObject::LayerData(SeuratObject::JoinLayers(a), layer=layer))
}

#' Run the HyperCom pipeline on a Seurat object
#'
#' `HyperCom_seurat` returns the results of `HyperCom` on a Seurat object
#'
#' This function takes a Seurat object (v5 or older assays) and runs `HyperCom` on one layer of one assay (raw counts by default)
#' and on the cell metadata (`seurat[[]]`). Every other argument is passed on to `HyperCom`; a misspelled one stops with "unused
#' argument". Needs the SeuratObject package.
#' @param seurat A Seurat object
#' @param sample A string that is the column name for sample (library or replicate) in the cell metadata of the Seurat object, or
#' NULL (the default) for one sample; cells are compared only within their sample
#' @param assay A string naming the assay to read counts from, "RNA" by default
#' @param layer A string naming the layer of the assay to read counts from: "counts" (the default) or "data", not the scaled
#' "scale.data"
#' @param group A string that is the column name for group (e.g. condition) in the cell metadata of the Seurat object, or NULL (the
#' default) for one group; groups feed the purity term of the priority (see `prioritize_lr`)
#' @param keep A vector of strings that are additional column names to keep from the cell metadata of the Seurat object, or NULL
#' @inheritParams HyperCom
#' @inheritDotParams HyperCom -lrdb -counts -metadata -sample -group -keep
#' @returns Invisibly, the list of `HyperCom`: priority (the ranked pairs, one data frame per weight), metadata, infMat, adj and lrs,
#' also written to outdir
#' @examplesIf requireNamespace("SeuratObject", quietly=TRUE)
#' seurat <- SeuratObject::CreateSeuratObject(alexanian2021.counts, meta.data=alexanian2021.metadata)
#' lrdb <- cellchat.mm.lrdb[cellchat.mm.lrdb$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' res <- HyperCom_seurat(seurat, lrdb, sample="sample", group="condition",
#'                        outdir=file.path(tempdir(), "HyperCom_seurat"))
#' head(res$priority[["1"]])
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample=NULL, assay="RNA", layer="counts", group=NULL, keep=NULL, ...){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat[[]]

  res <- HyperCom(lrdb, counts, metadata, sample=sample, group=group, keep=keep, ...)

  return(invisible(res))
}

#' Get the counts matrix of a SingleCellExperiment assay
#'
#' `get_sce_assay` returns the matrix with dimnames for given assay
#'
#' This function takes a SingleCellExperiment (or any SummarizedExperiment) and returns one assay with genes as rows and cells as
#' columns. Assays that are neither base nor Matrix matrices (e.g. a DelayedMatrix of HDF5-backed data) are read into a sparse
#' matrix with only the genes of lrdb, so the rest of the assay is never loaded.
#' @param sce A SingleCellExperiment object
#' @param lrdb A data frame of ligand-receptor interactions
#' @param assay A string naming the assay to get data from
#' @returns A matrix with row names being genes and column names being cells
#' @noRd
get_sce_assay <- function(sce, lrdb, assay="counts"){
  if(!requireNamespace("SummarizedExperiment", quietly=TRUE)){
    stop("reading a SingleCellExperiment needs the SummarizedExperiment package")
  }
  if(!methods::is(sce, "SummarizedExperiment")){
    stop("sce must be a SingleCellExperiment (or SummarizedExperiment) object")
  }
  if(is.null(colnames(sce))){
    stop("sce has no cell names: set colnames(sce) to the cell barcodes")
  }
  if(is.character(assay) && !(assay %in% SummarizedExperiment::assayNames(sce))){
    stop("no assay \"", assay, "\" in sce; available: ", paste(SummarizedExperiment::assayNames(sce), collapse=", "))
  }
  genes <- lrdb_genes(lrdb)
  if(!any(rownames(sce) %in% genes)){
    stop("no row name of sce is a gene of lrdb: the row names must be gene symbols of the species of lrdb (cellchat.mm.lrdb mouse, ",
         "cellchat.hs.lrdb human). With Ensembl IDs as row names and symbols in rowData: rownames(sce) <- rowData(sce)$Symbol")
  }
  counts <- read_counts(SummarizedExperiment::assay(sce, assay, withDimnames=TRUE), genes)

  return(counts)
}

#' Run the HyperCom pipeline on a SingleCellExperiment object
#'
#' `HyperCom_sce` returns the results of `HyperCom` on a SingleCellExperiment object
#'
#' This function takes a SingleCellExperiment, or any SummarizedExperiment with genes as rows and cells as columns, and runs
#' `HyperCom` on one assay and on the cell metadata (`colData`). The row names must be gene symbols and the column names cell names.
#' Assays held on disk (DelayedArray, HDF5Array) are read for the genes of lrdb only. Every other argument is passed on to `HyperCom`;
#' a misspelled one stops with "unused argument". Needs the SummarizedExperiment package (installed with SingleCellExperiment).
#' @param sce A SingleCellExperiment object
#' @param sample A string that is the column name for sample (library or replicate) in the cell metadata (`colData`) of the
#' SingleCellExperiment object, or NULL (the default) for one sample; cells are compared only within their sample
#' @param assay A string naming the assay to read counts from: "counts" (the default) or "logcounts", not a scaled assay
#' @param group A string that is the column name for group (e.g. condition) in the cell metadata of the SingleCellExperiment object,
#' or NULL (the default) for one group; groups feed the purity term of the priority (see `prioritize_lr`)
#' @param keep A vector of strings that are additional column names to keep from the cell metadata of the SingleCellExperiment object,
#' or NULL
#' @inheritParams HyperCom
#' @inheritDotParams HyperCom -lrdb -counts -metadata -sample -group -keep
#' @returns Invisibly, the list of `HyperCom`: priority (the ranked pairs, one data frame per weight), metadata, infMat, adj and lrs,
#' also written to outdir
#' @examplesIf requireNamespace("SingleCellExperiment", quietly=TRUE)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=alexanian2021.counts),
#'                                                   colData=alexanian2021.metadata)
#' lrdb <- cellchat.mm.lrdb[cellchat.mm.lrdb$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' res <- HyperCom_sce(sce, lrdb, sample="sample", group="condition",
#'                     outdir=file.path(tempdir(), "HyperCom_sce"))
#' head(res$priority[["1"]])
#' @export
HyperCom_sce <- function(sce, lrdb, sample=NULL, assay="counts", group=NULL, keep=NULL, ...){
  lrdb <- clean_lrdb(lrdb)
  counts <- get_sce_assay(sce, lrdb, assay)
  metadata <- as.data.frame(SummarizedExperiment::colData(sce), optional=TRUE)

  res <- HyperCom(lrdb, counts, metadata, sample=sample, group=group, keep=keep, ...)

  return(invisible(res))
}
