# The whole pipeline: HyperCom() on a counts matrix and cell metadata, HyperCom_seurat() on a Seurat object

#' Run the HyperCom pipeline
#'
#' `HyperCom` runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`,
#' `generate_metadata`, and `prioritize_lr`, writing every step to outdir. By default no matrix with a row and a column for every
#' cell and gene is made, so memory grows with the number of cells rather than its square (see `generate_infMat`).
#' The arguments are checked before any step runs, and the run stops when no interaction of lrdb has all its genes expressed (e.g.
#' lrdb of another species). Cells of counts that are not row names of metadata are left out, with a warning. Priority files of other weights already in outdir are kept, with a warning.
#' If a run stops while pairs are scored, the same call run again repeats the earlier steps and continues the scoring from the cache
#' in `<outdir>/6.priority/cache` (see `prioritize_lr`). The results are returned as
#' well as written, so they need not be read back: `res <- HyperCom(...)`, then `res$priority[["1"]]` or
#' `score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")`.
#' @inheritParams filter_counts
#' @inheritParams generate_metadata
#' @inheritParams prioritize_lr
#' @param outdir A string naming the output directory of the pipeline run
#' @param verbose A boolean of whether to report progress as messages: what each step found, the pairs scored and the top pair of
#' each weight
#' @param infmat.method A string naming the method of `generate_infMat`: "woodbury" (the default; fastest and least memory), "chol" or
#' "lu" (a full matrix of all cells and genes; "lu" gives exactly the numbers of HyperCom 0.0.0.9000). The methods agree to about 1e-15
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency) and lrs (the interactions with all genes expressed)
#' @examples
#' # a few ligands keep the example quick; use the whole database for real data
#' lrdb <- cellchat.mm.lrdb[cellchat.mm.lrdb$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' res <- HyperCom(lrdb, alexanian2021.counts, alexanian2021.metadata, sample="sample",
#'                 group="condition", outdir=file.path(tempdir(), "HyperCom"))
#' head(res$priority[["1"]])
#' list.files(file.path(tempdir(), "HyperCom"), recursive=TRUE)
#' scores <- score_hypercom(res$infMat, res$metadata, ligand="Thbs1", receptors="Sdc4")
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE,
                     outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="woodbury", cache=TRUE){
  check_lrdb(lrdb)
  check_columns(metadata, sample, group, keep)
  check_weights(weights)
  check_outdir(outdir)
  check_cache(cache)
  infmat.method <- match_method(infmat.method)
  threads <- get_threads(parallel, threads)
  in.metadata <- colnames(counts) %in% rownames(metadata)
  if(!any(in.metadata)){
    if(any(rownames(counts) %in% rownames(metadata))){
      stop("counts has cells as row names: give a gene x cell matrix, e.g. t(counts)")
    }
    stop("no cell of counts is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  stale <- setdiff(list.files(file.path(outdir, "6.priority"), pattern="^priority\\..*\\.csv$"), priority_filename(weights))
  if(length(stale) > 0){
    warning("outdir already holds priorities of other weights, which are kept: ", paste(stale, collapse=", "))
  }
  if(!all(in.metadata)){
    warning(sum(!in.metadata), " cells of counts are not in metadata and are left out")
    counts <- counts[, in.metadata, drop=FALSE]
  }
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

  counts <- filter_counts(counts, lrdb, outfile=file.path(outdir, "0.counts.rds"))
  lrs <- generate_lrs_table(counts, lrdb, outfile=file.path(outdir, "1.lrs.csv"))
  if(nrow(lrs) == 0){
    stop("no interaction of lrdb has all its genes expressed in counts (", ncol(counts), " of its genes are): check that lrdb is of the ",
         "species of counts (cellchat.mm.lrdb mouse, cellchat.hs.lrdb human) and that counts is gene x cell")
  }
  if(verbose){
    message(ncol(counts), " ligand and receptor genes expressed in ", nrow(counts), " cells; ", nrow(lrs), " interactions of lrdb have all ",
            "their genes expressed")
  }
  hyperedges <- generate_hyperedges(counts, lrs, outfile=file.path(outdir, "2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=file.path(outdir, "3.adj.rds"))
  rm(counts, hyperedges)
  infMat <- generate_infMat(adj, outfile=file.path(outdir, "4.infMat.rds"), method=infmat.method, threads=threads)
  if(verbose){
    message("Hypergraph of ", nrow(adj), " cells and genes; influence matrix made (method \"", infmat.method, "\")")
  }
  metadata <- generate_metadata(infMat, metadata, sample, group, keep, outfile=file.path(outdir, "5.metadata.csv"))
  priority <- prioritize_lr(lrs, adj, infMat, metadata, weights=weights, significance=significance, outdir=file.path(outdir, "6.priority"),
                            threads=threads, seed=seed, verbose=verbose, cache=cache)
  if(verbose){
    message("Done: results in ", outdir)
  }

  return(invisible(list(priority=priority, metadata=metadata, infMat=infMat, adj=adj, lrs=lrs)))
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
      stop("layer '", layer, "' is not a layer of assay '", assay, "': use counts, data or scale.data")
    }
    m <- methods::slot(a, layer)
    if(!length(m)){
      stop("layer '", layer, "' is empty or absent in assay '", assay, "'")
    }
    return(m)
  }

  lyrs <- SeuratObject::Layers(a, search=layer)
  if(!length(lyrs)){
    stop("no layer matching '", layer, "' in assay '", assay, "'; available: ", paste(SeuratObject::Layers(a), collapse=", "))
  }
  if(length(lyrs) == 1){
    return(SeuratObject::LayerData(a, layer=lyrs))
  }

  return(SeuratObject::LayerData(SeuratObject::JoinLayers(a), layer=layer))
}

#' Run the HyperCom pipeline on a seurat object
#'
#' `HyperCom_seurat` runs the HyperCom pipeline on a Seurat object
#'
#' This function runs `HyperCom` on one assay layer and the cell metadata of a Seurat object (v5 or older assays). Needs the
#' SeuratObject package.
#' @param seurat A Seurat object
#' @param sample A string that is the column name for sample in the cell metadata of the Seurat object
#' @param assay A string naming the assay to read counts from
#' @param layer A string naming the layer of the assay to read counts from
#' @param group A string that is the column name for group in the cell metadata of the Seurat object
#' @param keep A vector of strings that are additional column names to keep from the cell metadata of the Seurat object
#' @inheritParams HyperCom
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency) and lrs (the interactions with all genes expressed)
#' @examplesIf requireNamespace("SeuratObject", quietly=TRUE)
#' seurat <- SeuratObject::CreateSeuratObject(alexanian2021.counts, meta.data=alexanian2021.metadata)
#' lrdb <- cellchat.mm.lrdb[cellchat.mm.lrdb$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' res <- HyperCom_seurat(seurat, lrdb, sample="sample", group="condition",
#'                        outdir=file.path(tempdir(), "HyperCom_seurat"))
#' head(res$priority[["1"]])
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, weights=c(1), significance=FALSE,
                            parallel=FALSE, outdir="HyperCom", threads=NULL, seed=NULL, verbose=TRUE, infmat.method="woodbury", cache=TRUE){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat[[]]

  res <- HyperCom(lrdb, counts, metadata, sample=sample, group=group, keep=keep, weights=weights, significance=significance,
                  parallel=parallel, outdir=outdir, threads=threads, seed=seed, verbose=verbose, infmat.method=infmat.method, cache=cache)

  return(invisible(res))
}
