# Priority terms of ligand-receptor pairs: the pairs are scored (or screened, for bounds) in batches of whole ligands, each grid of a
# ligand and a receptor subunit once, with a cache of the pairs done so far; then each pair's terms from its cells and cell scores

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns the terms of the priority score of each ligand-receptor pair
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata
#' and returns the five terms of each pair's priority score, from the exact cell scores of its grids (pair_terms).
#' @inheritParams pair_terms
#' @returns A matrix with one row per pair and columns frac, pct.min, max, balance and purity
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, threads=1, verbose=FALSE, cached=NULL, kernel=NULL){
  return(pair_terms(lrs, adj, infMat, metadata, threads, verbose, cached, kernel)$terms)
}

#' Score or screen ligand-receptor pairs in batches
#'
#' `pair_terms` returns the priority terms of each ligand-receptor pair, or lower and upper bounds of them
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata and
#' returns the five terms of each pair's priority score: exact from the cell scores of its grids (score_grids, exact_terms), or, with
#' bounds, bounds of them from one screening pass (screen_grids, term_bounds). Pairs are done in batches of whole ligands
#' (batch_pairs): each ligand and receptor subunit grid once, whatever the number of pairs using it, and only one batch of cell scores
#' (about 64 MB at most, at least one grid per thread; batch_size) in memory at a time. With a cache, pairs whose genes are in it
#' take its terms or bounds and are not done again, and those of all pairs done so far (and of cached pairs of other calls) are
#' written to it after each batch; the first write of a run also saves the inputs (save_cache). After a failed write the run goes on
#' without the cache. Screening reuses one set of bucket ids across batches, freed on return; with blocks of 2^28 entries or more,
#' garbage is collected after each batch, as R would otherwise let the score bounds of many batches pile up beside the blocks.
#' @param lrs A data frame of ligand-receptor pairs in the format of lrdb
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param infMat A named influence matrix, dense or a hypercom_influence
#' @param metadata A data frame from generate_metadata
#' @param threads An integer number of threads
#' @param verbose A boolean of whether to report the pairs done, about every 5 s, with an estimate of the time left
#' @param cached An environment from open_cache, or NULL for no cache
#' @param kernel An environment from choose_kernel (exact or memory mode, with the run's tiles), or NULL for the exact kernels
#' @param bounds A boolean of whether to screen the pairs for bounds of their terms instead of scoring them
#' @param label A string naming the pairs in progress messages when screening ("background" or "real")
#' @returns A list with terms (a matrix with one row per pair and columns frac, pct.min, max, balance and purity), or with bounds two
#' such matrices, lo and hi
#' @noRd
pair_terms <- function(lrs, adj, infMat, metadata, threads=1, verbose=FALSE, cached=NULL, kernel=NULL, bounds=FALSE, label="real"){
  cells.idx <- metadata$sample != "gene"
  cells <- metadata$cell[cells.idx]
  groups <- metadata$group[cells.idx]
  sizes <- tabulate(as.integer(factor(metadata$sample[cells.idx])))
  gini.max <- 1-(1/length(unique(groups)))
  eps <- 1e-12 + max(sizes) * 2^-52

  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, -(1:2), drop=FALSE])
  genes <- unique(c(ligands, receptors[!is.na(receptors)]))
  check_genes(genes, rownames(infMat))
  ge <- gene_expression(adj, cells, genes)

  what <- if(bounds) "bounds" else "terms"
  parts <- if(bounds) c("lo", "hi") else "terms"
  pairs.done <- if(bounds) paste(label, "pairs screened") else "pairs scored"
  keys <- pair_keys(ligands, receptors)
  out <- stats::setNames(replicate(length(parts), empty_terms(nrow(lrs)), simplify=FALSE), parts)
  done <- logical(nrow(lrs))
  other <- NULL
  store <- cached_terms(cached, what)
  if(!is.null(store)){
    hit <- match(keys, rownames(store[[1]]))
    done <- !is.na(hit)
    for(p in parts){
      out[[p]][done, ] <- store[[p]][hit[done], , drop=FALSE]
    }
    other <- lapply(store, function(m) m[!(rownames(m) %in% keys), , drop=FALSE])
    if(verbose && any(done)){
      message("  Continuing from the cache in ", cached$dir, ": ", sum(done), " of ", nrow(lrs), " ",
              sub("pairs ", "pairs already ", pairs.done, fixed=TRUE))
    }
  }
  todo <- which(!done)

  blocks <- NULL
  if(length(todo) > 0){
    blocks <- grid_blocks(infMat, grid_cells(rownames(infMat), metadata), threads, kernel)
    if(bounds){
      on.exit(screen_release_cpp(), add=TRUE)
    }
  }
  n.grids <- nrow(pair_grids(ligands[todo], receptors[todo, , drop=FALSE]))
  size <- batch_size(threads, length(cells), sizes, n.grids, if(bounds) 32 else 16)
  progress <- start_progress(nrow(lrs), sum(done), n.grids, pairs.done)
  for(rows in batch_pairs(ligands[todo], receptors[todo, , drop=FALSE], size)){
    rows <- todo[rows]
    pairs <- pair_grids(ligands[rows], receptors[rows, , drop=FALSE])
    grid.keys <- paste(pairs$ligand, pairs$receptor, sep="\t")
    if(bounds){
      scores <- screen_grids(infMat, metadata, pairs$ligand, pairs$receptor, threads, blocks, keep.ids=TRUE)
    } else {
      scores <- score_grids(infMat, metadata, pairs$ligand, pairs$receptor, threads, blocks)
    }

    for(i in rows){
      rs <- stats::na.omit(receptors[i, ])
      cols <- match(paste(ligands[i], rs, sep="\t"), grid.keys)
      pc <- pair_cells(c(ligands[i], rs), ge, groups, gini.max)
      same <- all(ligands[i] == rs)
      if(bounds){
        tb <- term_bounds(scores, cols, same, pc, eps)
        out$lo[i, ] <- tb$lo
        out$hi[i, ] <- tb$hi
      } else {
        out$terms[i, ] <- exact_terms(scores, cols, same, pc)
      }
    }
    rm(scores)
    if(bounds && block_entries(blocks) >= 2^28){
      gc()
    }
    done[rows] <- TRUE
    if(cache_on(cached)){
      write_cache_rows(cached, what, out, keys, done, other)
    }
    report_progress(progress, verbose, length(rows), nrow(pairs))
  }

  return(out)
}

#' Calculate the priority terms of a pair from its cell scores
#'
#' `exact_terms` returns the five priority terms of a ligand-receptor pair
#'
#' This function takes the cell scores of grids, the columns of the pair's subunit grids, whether its ligand is its receptor and its
#' cells from pair_cells, and returns frac, pct.min, the mean of each cell's larger score as sender or receiver (max), how evenly the
#' cells split into senders and receivers (balance: the Gini impurity of receiver or not, over 0.5; a pair of one gene has only
#' senders) and purity.
#' @param sr A list of two matrices, sender and receiver, with one row per cell and one column per grid
#' @param cols An integer vector of the columns of the pair's subunit grids
#' @param same A boolean of whether the ligand is the pair's only receptor subunit
#' @param pc A list from pair_cells
#' @returns A numeric vector of five terms
#' @noRd
exact_terms <- function(sr, cols, same, pc){
  pair.sr <- average_scores(sr, cols)
  ss <- pair.sr$sender
  rr <- pair.sr$receiver
  if(same){
    type <- rep(FALSE, length(ss))
  } else {
    type <- (rr - ss) > 0
  }
  max.score <- mean(pmax(ss, rr)[pc$present])
  balance <- gini_impurity(type[pc$present]) / 0.5

  return(c(pc$frac, pc$pct.min, max.score, balance, pc$purity))
}

#' Name the priority terms
#'
#' `term_names` returns the names of the five terms of a priority score, in their order
#' @returns A vector of strings
#' @noRd
term_names <- function(){
  return(c("frac", "pct.min", "max", "balance", "purity"))
}

#' Make an empty matrix of priority terms
#'
#' `empty_terms` returns a matrix of NA with one row per pair and one column per priority term
#' @param n An integer number of pairs
#' @returns A numeric matrix with columns frac, pct.min, max, balance and purity
#' @noRd
empty_terms <- function(n){
  return(matrix(NA_real_, n, 5, dimnames=list(NULL, term_names())))
}

#' Read the cached terms or bounds of a run
#'
#' `cached_terms` returns the priority terms or term bounds a cache holds, as a list of matrices, or NULL
#' @param cached An environment from open_cache, or NULL
#' @param what A string, "terms" or "bounds"
#' @returns A list with terms, or with lo and hi (matrices with rows named by pair_keys), or NULL when the cache is off or empty
#' @noRd
cached_terms <- function(cached, what){
  if(!cache_on(cached) || is.null(cached[[what]])){
    return(NULL)
  }
  if(what == "terms"){
    return(list(terms=cached$terms))
  }

  return(cached$bounds)
}

#' Find pairs in the cached terms
#'
#' `cached_rows` returns the row of each pair in the cached priority terms, NA for pairs not in them or without a cache
#' @param cached An environment from open_cache, or NULL
#' @param keys A vector of strings of pair keys from pair_keys
#' @returns An integer vector, one per pair
#' @noRd
cached_rows <- function(cached, keys){
  store <- cached_terms(cached, "terms")
  if(is.null(store)){
    return(rep(NA_integer_, length(keys)))
  }

  return(match(keys, rownames(store$terms)))
}

#' Write the terms or bounds of the pairs done so far to a cache
#'
#' `write_cache_rows` saves the terms or bounds of the pairs done so far, after the cached ones of other calls
#'
#' This function takes the cache, what is cached, the terms or bounds of the pairs, their keys, which pairs are done and the cached
#' rows of other pairs, and writes one row per distinct key done (write_cache_terms).
#' @param cached An environment from open_cache
#' @param what A string, "terms" or "bounds"
#' @param out A list of matrices of terms, or lo and hi
#' @param keys A vector of strings of pair keys from pair_keys
#' @param done A logical vector of the pairs done
#' @param other A list of matrices of the cached rows of other pairs, or NULL
#' @returns A boolean of whether the file was written
#' @noRd
write_cache_rows <- function(cached, what, out, keys, done, other){
  keep <- done & !duplicated(keys)
  rows <- lapply(out, function(m){
    m <- m[keep, , drop=FALSE]
    rownames(m) <- keys[keep]
    return(m)
  })
  if(!is.null(other)){
    rows <- Map(rbind, other[names(rows)], rows)
  }

  return(write_cache_terms(cached, what, if(what == "terms") rows$terms else rows))
}

#' Start reporting the progress of a run
#'
#' `start_progress` returns the state of a progress report
#' @param total An integer number of pairs
#' @param done An integer number of pairs already done
#' @param grids An integer number of grids to do
#' @param what A string naming what is done, e.g. "pairs scored"
#' @returns An environment
#' @noRd
start_progress <- function(total, done, grids, what){
  return(list2env(list(total=total, done=done, grids=grids, grids.done=0, what=what, t0=Sys.time(), shown=Sys.time())))
}

#' Report the progress of a run
#'
#' `report_progress` adds a batch to a progress report and reports it about every 5 s and at the end
#'
#' This function takes the progress state, verbose and the pairs and grids of a batch and, every 5 s or so and after the last batch,
#' reports the pairs done with an estimate of the time left from the grids done so far.
#' @param progress An environment from start_progress
#' @param verbose A boolean of whether to report
#' @param pairs An integer number of pairs in the batch
#' @param grids An integer number of grids in the batch
#' @returns NULL, invisibly
#' @noRd
report_progress <- function(progress, verbose, pairs, grids){
  progress$done <- progress$done + pairs
  progress$grids.done <- progress$grids.done + grids
  last <- progress$done == progress$total
  if(verbose && (as.numeric(difftime(Sys.time(), progress$shown, units="secs")) >= 5 || last)){
    left <- as.numeric(difftime(Sys.time(), progress$t0, units="secs")) / progress$grids.done * (progress$grids - progress$grids.done)
    message("  ", progress$done, " of ", progress$total, " ", progress$what, if(!last) paste0(", about ", format_duration(left), " left"))
    progress$shown <- Sys.time()
  }

  return(invisible(NULL))
}

#' Format a duration for a progress message
#'
#' `format_duration` returns a duration in seconds as a short string
#'
#' This function takes a number of seconds and returns it rounded in seconds below 90 s, in minutes below 90 min and in hours above.
#' @param seconds A double number of seconds
#' @returns A string, e.g. "40 s", "12 min" or "2.5 h"
#' @noRd
format_duration <- function(seconds){
  if(seconds < 90){
    return(paste0(max(1, round(seconds)), " s"))
  }
  if(seconds < 5400){
    return(paste0(round(seconds / 60), " min"))
  }

  return(paste0(round(seconds / 3600, 1), " h"))
}

#' Find the cells that express each gene
#'
#' `gene_expression` returns the fraction of cells expressing each gene and the cells that do
#'
#' This function takes the adjacency, the cells and the genes and returns, from the gene x cell counts of the adjacency, the fraction
#' of the cells expressing each gene and their positions. The counts of a symmetric sparse adjacency are read from its stored triangle,
#' without expanding the whole matrix.
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param cells A vector of strings of cell names
#' @param genes A vector of strings of gene names
#' @returns A list with frac (a numeric vector named by gene) and expressed (a list named by gene of integer vectors of cell positions)
#' @noRd
gene_expression <- function(adj, cells, genes){
  if(inherits(adj, "hypercom_adjacency")){
    gene.counts <- adj[cells, genes, drop=FALSE]
  } else if(methods::is(adj, "dsCMatrix") && !any(genes %in% cells)){
    U <- methods::new("dgCMatrix", i=adj@i, p=adj@p, x=adj@x, Dim=adj@Dim, Dimnames=dimnames(adj))
    gene.counts <- as.matrix(Matrix::t(U[genes, cells, drop=FALSE]) + U[cells, genes, drop=FALSE])
    rm(U)
  } else {
    gene.counts <- as.matrix(Matrix::t(adj[genes, cells, drop=FALSE]))
  }
  expressed <- lapply(seq_len(ncol(gene.counts)), function(k) which(gene.counts[, k] > 0))
  names(expressed) <- colnames(gene.counts)

  return(list(frac=colSums(gene.counts) / length(cells), expressed=expressed))
}

#' Find the cells of a ligand-receptor pair
#'
#' `pair_cells` returns the cells expressing any gene of a pair and the priority terms that depend on them alone
#'
#' This function takes the genes of a pair, the expression from gene_expression, the group of each cell and the largest possible Gini
#' impurity of the groups and returns which cells express any of the genes, their fraction of all cells (frac), the smallest fraction of
#' cells expressing one of the genes (pct.min) and the group purity of the cells.
#' @param genes A vector of strings, the ligand and the receptor subunits of the pair
#' @param ge A list from gene_expression
#' @param groups A vector of the group of each cell
#' @param gini.max A double of the largest possible Gini impurity, 1 - 1 / the number of groups
#' @returns A list with present (a logical vector over the cells), frac, pct.min and purity
#' @noRd
pair_cells <- function(genes, ge, groups, gini.max){
  present <- logical(length(groups))
  for(g in genes){
    present[ge$expressed[[g]]] <- TRUE
  }

  return(list(present=present, frac=sum(present) / length(groups), pct.min=min(ge$frac[genes]),
              purity=group_purity(groups[present], gini.max)))
}

#' Calculate the group purity of the cells of a pair
#'
#' `group_purity` returns how far the groups of a pair's cells are from evenly mixed
#'
#' This function takes the groups of the cells expressing a pair's genes and the largest possible Gini impurity and returns their
#' difference, at least 0: cells perfectly balanced over 5, 10 or 11 groups round the difference to -1.1e-16, whose log would make the
#' priority NaN. One group (largest impurity 0) gives 1.
#' @param groups A factor of the groups of the cells
#' @param gini.max A double of the largest possible Gini impurity, 1 - 1 / the number of groups
#' @returns A numeric purity
#' @noRd
group_purity <- function(groups, gini.max){
  if(gini.max <= 0){
    return(1)
  }

  return(max(0, gini.max - gini_impurity(groups)))
}

#' Calculate the Gini impurity of a set of labels
#'
#' `gini_impurity` returns the Gini impurity of a vector of labels
#'
#' This function takes a vector of labels and returns the probability that a label drawn at random is mislabelled
#' when labelled according to the distribution of labels (as mltools::gini_impurity does for labels, not counts).
#' @param vals A vector of labels
#' @returns A numeric Gini impurity
#' @noRd
gini_impurity <- function(vals){
  if(is.factor(vals)){
    counts <- tabulate(vals, nbins=nlevels(vals))
  } else if(is.logical(vals)){
    counts <- c(sum(!vals, na.rm=TRUE), sum(vals, na.rm=TRUE))
    counts <- counts[counts > 0]
  } else {
    counts <- table(vals)
  }
  total <- sum(counts)

  return(sum((counts/total)*(1-counts/total)))
}

#' Average the cell scores of the subunits of a receptor complex
#'
#' `average_scores` returns the sender and receiver scores of a ligand and receptor complex
#'
#' This function takes the scores of grids and the columns of a complex's subunit grids and returns the sum of the columns,
#' added in order, divided by their number (one column is returned as is).
#' @param sr A list of two matrices, sender and receiver, with one row per cell and one column per grid
#' @param cols An integer vector of the columns of the complex's subunits
#' @returns A list of two vectors, sender and receiver
#' @noRd
average_scores <- function(sr, cols){
  ss <- sr$sender[, cols[1]]
  rr <- sr$receiver[, cols[1]]
  if(length(cols) > 1){
    for(k in cols[-1]){
      ss <- ss + sr$sender[, k]
      rr <- rr + sr$receiver[, k]
    }
    ss <- ss / length(cols)
    rr <- rr / length(cols)
  }

  return(list(sender=ss, receiver=rr))
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

#' Name the ligand-receptor pairs of a table by their genes
#'
#' `lrs_keys` returns the pair_keys of the pairs of a table of ligand-receptor pairs
#' @param lrs A data frame of ligand-receptor pairs in the format of lrdb
#' @returns A vector of strings, one per pair
#' @noRd
lrs_keys <- function(lrs){
  return(pair_keys(as.character(lrs$ligand), as.matrix(lrs[, -(1:2), drop=FALSE])))
}

#' List the distinct grids of ligand-receptor pairs
#'
#' `pair_grids` returns the distinct ligand and receptor subunit pairs, one per cell grid
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and returns each distinct (ligand, receptor
#' subunit) pair once, in order of first appearance by receptor column: a complex is scored through one grid per subunit.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @returns A data frame with columns ligand and receptor, one row per grid
#' @noRd
pair_grids <- function(ligands, receptors){
  grids <- unique(data.frame(ligand=rep(ligands, ncol(receptors)), receptor=as.vector(receptors)))

  return(grids[!is.na(grids$receptor), ])
}

#' Group ligand-receptor pairs into batches of whole ligands
#'
#' `batch_pairs` returns the rows of each batch of ligand-receptor pairs
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and splits the pairs into batches of whole ligands
#' with about size distinct ligand and receptor subunit grids each. A grid belongs to one ligand, so no grid is needed by two batches.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @param size An integer number of grids per batch
#' @returns A list of integer vectors of rows, one per batch, rows in their original order
#' @noRd
batch_pairs <- function(ligands, receptors, size){
  grids <- pair_grids(ligands, receptors)
  ligand.order <- unique(ligands)
  per.ligand <- tabulate(match(grids$ligand, ligand.order), nbins=length(ligand.order))

  batch <- integer(length(ligand.order))
  b <- 1
  filled <- 0
  for(k in seq_along(ligand.order)){
    if(filled > 0 && filled + per.ligand[k] > size){
      b <- b + 1
      filled <- 0
    }
    batch[k] <- b
    filled <- filled + per.ligand[k]
  }

  return(split(seq_along(ligands), batch[match(ligands, ligand.order)]))
}

#' Choose the number of grids per batch
#'
#' `batch_size` returns how many grids pair_terms scores or screens per kernel call
#'
#' This function takes the thread count, the number of cells, the cells of each sample, the number of grids and the bytes of cell
#' scores per cell and grid (16 for scores, 32 for bounds) and returns the batch size: at least one grid per thread and at most about
#' 64 MB of cell scores; when grids are big (their kernel buffers above the kernel's cache budget, kernel_constants_cpp), all threads
#' work on one grid at a time, so smaller batches keep every thread busy and there are about 20 batches or more, for progress messages.
#' @param threads An integer number of threads
#' @param cells An integer number of cells
#' @param sample.sizes An integer vector of the number of cells of each sample
#' @param n.grids An integer number of grids
#' @param bytes A number of bytes of cell scores per cell and grid
#' @returns An integer number of grids per batch
#' @noRd
batch_size <- function(threads, cells, sample.sizes, n.grids, bytes=16){
  size <- max(threads, min(max(256, 32 * threads), floor(2^26 / (bytes * cells))))
  kc <- kernel_constants_cpp()
  if(kc$grid.entry.bytes * cell_pairs(sample.sizes) > kc$cache.budget){
    size <- max(threads, min(size, ceiling(n.grids / 20)))
  }

  return(size)
}
