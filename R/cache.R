# Cache of prioritize_lr: the priority terms of the pairs scored so far, so that a stopped run continues where it stopped

#' Open the cache of a prioritize_lr run
#'
#' `open_cache` returns what the cache folder holds for these inputs
#'
#' This function takes the cache folder and the inputs of prioritize_lr and returns the cache: the folder, the fingerprint of the
#' inputs with the token of the run that made the cache, and the priority terms and background pairs saved by an earlier call on the
#' same inputs. A cache made from other inputs (or damaged) is not used: it is moved to `<dir>.old`, replacing an earlier one, with a
#' warning, so that a call with mistaken inputs does not destroy it. Temporary files of a write that was stopped are deleted. Nothing
#' else is written.
#' @param dir A string naming the cache folder
#' @inheritParams prioritize_lr
#' @returns A list with dir, fingerprint, terms (a matrix of priority terms with one row per pair named by pair_keys, or NULL) and
#' background (a dataframe of background pairs, or NULL)
#' @noRd
open_cache <- function(dir, adj, infMat, metadata){
  fingerprint <- cache_fingerprint(adj, infMat, metadata)
  unlink(file.path(dir, paste0(cache_files(), ".tmp")))
  old <- read_cache_file(dir, "fingerprint.rds")
  terms <- NULL
  background <- NULL

  if(same_fingerprint(old, fingerprint)){
    fingerprint$token <- old$token
    saved <- read_cache_file(dir, "terms.rds")
    if(is.list(saved) && identical(saved$token, old$token)){
      terms <- saved$terms
    }
    background <- read_cache_file(dir, "background.rds")
  } else {
    fingerprint$token <- paste0(Sys.getpid(), "-", format(as.numeric(Sys.time()), digits=17))
    if(any(file.exists(file.path(dir, cache_files())))){
      aside <- paste0(dir, ".old")
      unlink(aside, recursive=TRUE)
      if(!file.rename(dir, aside)){
        unlink(file.path(dir, cache_files()))
        aside <- NULL
      }
      warning("the cache in ", dir, " was made from other inputs, or is damaged, and is not used",
              if(!is.null(aside)) paste0("; it is moved to ", aside, " (rename it back to continue the run it belongs to)"))
    }
  }

  return(list(dir=dir, fingerprint=fingerprint, terms=terms, background=background))
}

#' Fingerprint the inputs of a prioritize_lr run
#'
#' `cache_fingerprint` returns what ties cached priority terms to the inputs they were calculated from
#'
#' This function takes the inputs of prioritize_lr and returns the package version, the class of the influence matrix, the node
#' names, the cell, sample and group of every metadata row as strings (their type does not change the priorities, e.g. metadata read
#' back from 5.metadata.csv), and the entries of the adjacency and the influence matrix between 16 nodes spread over them. Other
#' data, another restart probability or another influence matrix method change these.
#' @inheritParams prioritize_lr
#' @returns A list
#' @noRd
cache_fingerprint <- function(adj, infMat, metadata){
  nodes <- rownames(infMat)
  probe <- nodes[unique(round(seq(1, length(nodes), length.out=16)))]

  return(list(format=1L, version=as.character(utils::packageVersion("HyperCom")), class=class(infMat), nodes=nodes,
              cell=as.character(metadata$cell), sample=as.character(metadata$sample), group=as.character(metadata$group),
              adj=as.matrix(adj[probe, probe, drop=FALSE]), infMat=as.matrix(infMat[probe, probe, drop=FALSE])))
}

#' Compare two cache fingerprints
#'
#' `same_fingerprint` returns whether cached priority terms fit the current inputs
#'
#' This function takes a fingerprint read from the cache (NULL if missing or damaged) and the current one and compares everything but
#' the token exactly, except the entries of a factored influence matrix (method "woodbury"), which are compared within 1e-10: the BLAS
#' (another cluster node, another thread count) can change the last bits of its hyperedge inverse, and no other method makes that
#' class, while other data move the entries far more. A dense influence matrix ("lu" or "chol", which differ by about 1e-15) must
#' match exactly.
#' @param old A list from cache_fingerprint with a token, or NULL
#' @param new A list from cache_fingerprint
#' @returns A boolean
#' @noRd
same_fingerprint <- function(old, new){
  if(!is.list(old) || !identical(setdiff(names(old), "token"), names(new))){
    return(FALSE)
  }
  exact <- setdiff(names(new), "infMat")
  tolerance <- if(identical(new$class, "hypercom_influence")) 1e-10 else 0

  return(identical(old[exact], new[exact]) && isTRUE(all.equal(old$infMat, new$infMat, tolerance=tolerance)))
}

#' Save the inputs of a cached run
#'
#' `save_cache` writes the fingerprint and the background pairs of a run into its cache folder
#'
#' This function takes the cache from open_cache, creates its folder and writes the fingerprint and, if any, the background pairs.
#' It is called when the first batch of pairs is scored, so a call that stops before leaves no folder.
#' @param cached A list from open_cache
#' @returns A boolean of whether the files were written
#' @noRd
save_cache <- function(cached){
  dir.create(cached$dir, showWarnings=FALSE, recursive=TRUE)
  done <- write_cache_file(cached$fingerprint, cached$dir, "fingerprint.rds")
  if(done && !is.null(cached$background)){
    done <- write_cache_file(cached$background, cached$dir, "background.rds")
  }

  return(done)
}

#' Remove the cache of a finished run
#'
#' `close_cache` deletes the cache files and then the cache folder if nothing else is in it
#'
#' This function takes the cache folder and deletes the files prioritize_lr writes there, then the folder if it is empty.
#' @param dir A string naming the cache folder
#' @returns NULL, invisibly
#' @noRd
close_cache <- function(dir){
  unlink(file.path(dir, c(cache_files(), paste0(cache_files(), ".tmp"))))
  if(dir.exists(dir) && length(list.files(dir, all.files=TRUE, no..=TRUE)) == 0){
    unlink(dir, recursive=TRUE)
  }

  return(invisible(NULL))
}

#' Name the files of a cache
#'
#' `cache_files` returns the names of the files prioritize_lr writes into its cache folder
#'
#' This function returns the file names of the fingerprint, the background pairs and the priority terms.
#' @returns A vector of strings of file names
#' @noRd
cache_files <- function(){
  return(c("fingerprint.rds", "background.rds", "terms.rds"))
}

#' Read a cache file
#'
#' `read_cache_file` returns the object saved in a cache file, or NULL
#'
#' This function takes the cache folder and a file name and returns the file's object, or NULL when the file is missing or cannot be
#' read.
#' @param dir A string naming the cache folder
#' @param name A string naming the file
#' @returns The saved object, or NULL
#' @noRd
read_cache_file <- function(dir, name){
  f <- file.path(dir, name)
  if(!file.exists(f)){
    return(NULL)
  }

  return(tryCatch(readRDS(f), error=function(e) NULL))
}

#' Write a cache file
#'
#' `write_cache_file` saves an object into the cache folder so that a stopped write never leaves a damaged file
#'
#' This function takes an object, the cache folder and a file name and saves the object uncompressed under a temporary name, then
#' renames it, replacing the earlier file. A write that fails (full disk, locked file, folder removed) warns and returns FALSE: the
#' cache only saves time, so the run goes on without it.
#' @param x An object
#' @param dir A string naming the cache folder
#' @param name A string naming the file
#' @returns A boolean of whether the file was written
#' @noRd
write_cache_file <- function(x, dir, name){
  f <- file.path(dir, name)
  tmp <- paste0(f, ".tmp")
  done <- tryCatch({
    saveRDS(x, tmp, compress=FALSE)
    file.rename(tmp, f)
  }, error=function(e) FALSE, warning=function(w) FALSE)
  if(!done){
    unlink(tmp)
    warning("could not write the cache file ", f, "; scoring goes on without the cache")
  }

  return(done)
}

#' Name ligand-receptor pairs by their genes
#'
#' `pair_keys` returns one string per pair naming its ligand and receptor subunits
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and returns the ligand and the subunits, in order and
#' without NA, separated by tabs. The priority terms of a pair depend only on these genes, so pairs with the same key share them.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @returns A vector of strings, one per pair
#' @noRd
pair_keys <- function(ligands, receptors){
  subunits <- apply(receptors, 1, function(x) paste(x[!is.na(x)], collapse="\t"))

  return(paste(ligands, subunits, sep="\t"))
}
