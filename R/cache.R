# Cache of prioritize_lr: the priority terms of the pairs scored so far and the term bounds of the pairs screened so far, so that a
# stopped run continues where it stopped

#' Open the cache of a prioritize_lr run
#'
#' `open_cache` returns what the cache folder holds for these inputs
#'
#' This function takes the cache folder and the inputs of prioritize_lr and returns the cache: the folder, the fingerprint of the
#' inputs with the token of the run that made the cache, and the priority terms, term bounds and background pairs saved by an earlier
#' call on the same inputs. A cache made from other inputs (or damaged) is not used: it is moved to `<dir>.old`, replacing an earlier
#' one, with a warning, so that a call with mistaken inputs does not destroy it. Temporary files of a write that was stopped are
#' deleted. Nothing else is written. The cache is an environment, so the steps of a run share what each writes and a failed write
#' turns the cache off for all of them.
#' @param dir A string naming the cache folder
#' @inheritParams prioritize_lr
#' @inheritParams cache_fingerprint
#' @returns An environment with dir, fingerprint, terms (a matrix of priority terms with one row per pair named by pair_keys, or NULL),
#' bounds (a list of two such matrices, lo and hi, or NULL), background (a data frame of background pairs, or NULL), saved (whether
#' this run wrote the fingerprint) and failed (whether a write failed)
#' @noRd
open_cache <- function(dir, adj, infMat, metadata, kernel="exact"){
  fingerprint <- cache_fingerprint(adj, infMat, metadata, kernel)
  unlink(file.path(dir, paste0(cache_files(), ".tmp")))
  old <- read_cache_file(dir, "fingerprint.rds")
  terms <- NULL
  bounds <- NULL
  background <- NULL

  if(same_fingerprint(old, fingerprint)){
    fingerprint$token <- old$token
    saved <- read_cache_file(dir, "terms.rds")
    if(is.list(saved) && identical(saved$token, old$token)){
      terms <- saved$terms
    }
    saved <- read_cache_file(dir, "bounds.rds")
    if(is.list(saved) && identical(saved$token, old$token) && is.matrix(saved$lo) && is.matrix(saved$hi)){
      bounds <- saved[c("lo", "hi")]
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
      warning("the cache in ", dir, " ", cache_mismatch(old, fingerprint), " and is not used",
              if(!is.null(aside)) paste0("; it is moved to ", aside, " (rename it back to continue the run it belongs to)"))
    }
  }

  return(list2env(list(dir=dir, fingerprint=fingerprint, terms=terms, bounds=bounds, background=background, saved=FALSE,
                       failed=FALSE)))
}

#' Say why a cache does not fit a run
#'
#' `cache_mismatch` returns why a cache is not used, for its warning
#'
#' This function takes the fingerprint read from the cache and the current one and returns that the cache was made with the other
#' kernels (exact or memory mode, which follows the memory.gb argument and the computer's memory) when that is the only difference, else
#' that it was made from other inputs or is damaged.
#' @param old A list from cache_fingerprint with a token, or NULL
#' @param new A list from cache_fingerprint
#' @returns A string
#' @noRd
cache_mismatch <- function(old, new){
  kernel_words <- function(type){
    return(if(identical(type, "exact")) "the exact kernels" else "memory mode")
  }
  if(is.list(old) && is.character(old$kernel) && xor(identical(old$kernel, "exact"), identical(new$kernel, "exact")) &&
     same_fingerprint(utils::modifyList(old, list(kernel=new$kernel)), new)){
    return(paste0("was made with ", kernel_words(old$kernel), " but this run uses ", kernel_words(new$kernel),
                  " (it follows memory.gb and the computer's memory),"))
  }

  return("was made from other inputs, or is damaged,")
}

#' Check whether a cache is in use
#'
#' `cache_on` returns whether pairs are read from and written to a cache
#'
#' This function takes the cache from open_cache, or NULL, and returns FALSE for no cache or a cache whose writing failed.
#' @param cached An environment from open_cache, or NULL
#' @returns A boolean
#' @noRd
cache_on <- function(cached){
  return(!is.null(cached) && !cached$failed)
}

#' Write the terms or term bounds of a cache
#'
#' `write_cache_terms` saves the priority terms or the term bounds of all pairs done so far into the cache folder
#'
#' This function takes the cache, the kind of rows and their new contents and writes them into terms.rds or bounds.rds under the run's
#' token, writing the fingerprint and background pairs first if this run has not yet. On success the cache keeps the new rows, so later
#' steps write them again with theirs; a failed write turns the cache off for the rest of the run.
#' @param cached An environment from open_cache
#' @param what A string, "terms" or "bounds"
#' @param x A matrix of priority terms (terms), or a list of two matrices lo and hi (bounds), rows named by pair_keys
#' @returns A boolean of whether the file was written
#' @noRd
write_cache_terms <- function(cached, what, x){
  done <- cached$saved || save_cache(cached)
  cached$saved <- done
  if(done){
    content <- if(what == "terms") list(token=cached$fingerprint$token, terms=x) else c(list(token=cached$fingerprint$token), x)
    done <- write_cache_file(content, cached$dir, paste0(what, ".rds"))
  }
  if(done){
    cached[[what]] <- x
  } else {
    cached$failed <- TRUE
  }

  return(done)
}

#' Fingerprint the inputs of a prioritize_lr run
#'
#' `cache_fingerprint` returns what ties cached priority terms to the inputs they were calculated from
#'
#' This function takes the inputs of prioritize_lr and returns the package version, the class of the influence matrix, the node
#' names, the cell, sample and group of every metadata row as strings (their type does not change the priorities, e.g. metadata read
#' back from 5.metadata.csv), and the entries of the adjacency and the influence matrix between 16 nodes spread over them. Other
#' data, another restart probability or another influence matrix method change these. It also holds the grid kernels: memory mode's
#' terms and bounds differ from the exact kernels' by ~1e-11, so a cache never mixes them.
#' @inheritParams prioritize_lr
#' @param kernel A string naming the grid kernels of the run, "exact" or "memory" (memory mode)
#' @returns A list
#' @noRd
cache_fingerprint <- function(adj, infMat, metadata, kernel="exact"){
  nodes <- rownames(infMat)
  probe <- nodes[unique(round(seq(1, length(nodes), length.out=16)))]

  return(list(format=2L, version=as.character(utils::packageVersion("HyperCom")), kernel=kernel, class=class(infMat), nodes=nodes,
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
  if(!is.list(old) || !identical(setdiff(names(old), "token"), setdiff(names(new), "token"))){
    return(FALSE)
  }
  exact <- setdiff(names(new), c("infMat", "token"))
  tolerance <- if(identical(new$class, "hypercom_influence")) 1e-10 else 0

  return(identical(old[exact], new[exact]) && isTRUE(all.equal(old$infMat, new$infMat, tolerance=tolerance)))
}

#' Save the inputs of a cached run
#'
#' `save_cache` writes the fingerprint and the background pairs of a run into its cache folder
#'
#' This function takes the cache from open_cache, creates its folder and writes the fingerprint and, if any, the background pairs.
#' It is called when the first batch of pairs is scored or screened, so a call that stops before leaves no folder.
#' @param cached An environment from open_cache
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
#' This function returns the file names of the fingerprint, the background pairs, the priority terms and the term bounds.
#' @returns A vector of strings of file names
#' @noRd
cache_files <- function(){
  return(c("fingerprint.rds", "background.rds", "terms.rds", "bounds.rds"))
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
