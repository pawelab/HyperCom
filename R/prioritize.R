# Pipeline step 6: priority scores of ligand-receptor pairs and their p-values

#' Generate a background set of random gene interactions
#'
#' `generate_background_lrs` returns a dataframe of random gene interactions.
#'
#' This generates a table of random gene interactions where the number of genes in an interaction is equal to the complex size + 1.
#' @inheritParams prioritize_lr
#' @param count An integer of the number of interactions to generate
#' @returns A dataframe of random gene interactions
#' @noRd
generate_background_lrs <- function(lrs, metadata, count=1000){
  background <- data.frame(interaction_name=rep("background", count), check.names=FALSE)
  nodes <- metadata$cell[metadata$sample == "gene"]
  for(i in seq_len(ncol(lrs))[-1]){
    background[[i]] <- sample(nodes, count, replace=TRUE)
  }
  colnames(background) <- colnames(lrs)

  return(background)
}

#' Perform a permutation test
#'
#' `permutation_test` returns a numeric p-value
#'
#' This function takes a priority and the priorities of the background pairs and returns the fraction of background priorities at
#' least as high in absolute value. Background pairs without a priority (NA or NaN) are left out; NA if x has no priority or no
#' background pair has one.
#' @param x A numeric priority
#' @param background A numeric vector of background priorities
#' @returns A numeric p-value
#' @noRd
permutation_test <- function(x, background){
  background <- background[!is.na(background)]
  if(length(background) == 0){
    return(NA_real_)
  }

  return(sum(abs(background) >= abs(x)) / length(background))
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

#' List the distinct grids of ligand-receptor pairs
#'
#' `pair_grids` returns the distinct ligand and receptor subunit pairs, one per cell grid
#'
#' This function takes the ligands and receptor subunits of ligand-receptor pairs and returns each distinct (ligand, receptor
#' subunit) pair once, in order of first appearance by receptor column: a complex is scored through one grid per subunit.
#' @param ligands A vector of strings of ligand names, one per pair
#' @param receptors A matrix of strings of receptor subunit names, one row per pair, NA where a complex has fewer subunits
#' @returns A dataframe with columns ligand and receptor, one row per grid
#' @noRd
pair_grids <- function(ligands, receptors){
  grids <- unique(data.frame(ligand=rep(ligands, ncol(receptors)), receptor=as.vector(receptors)))

  return(grids[!is.na(grids$receptor), ])
}

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns the terms of the priority score of each ligand-receptor pair
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata
#' and returns the five terms of each pair's priority score. Pairs are scored in batches of whole ligands: each ligand and receptor
#' subunit grid is scored once, whatever the number of pairs using it, and only one batch of cell scores (about 64 MB at most, at least
#' one grid per thread) is in memory at a time. The gene x cell counts of a symmetric sparse adjacency are read from its stored
#' triangle, without expanding the whole matrix. With a cache, pairs whose genes are in it take its terms and are not scored again,
#' and the terms of all pairs scored so far (and of cached pairs of other calls) are written to it after each batch; the first write
#' of a run also saves the inputs (save_cache). After a failed write the run goes on without the cache.
#' @inheritParams prioritize_lr
#' @param threads An integer number of threads
#' @param verbose A boolean of whether to report the pairs scored, about every 5 s, with an estimate of the time left (pairs are
#' scored in about 20 batches or more when grids are bigger than the kernel's 24 MB cache budget; see batch_size)
#' @param cached An environment from open_cache, or NULL for no cache
#' @param kernel An environment from choose_kernel (exact or memory mode, with the run's tiles), or NULL for the exact kernels
#' @returns A matrix with one row per pair and columns frac, pct.min, max, balance and purity
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, threads=1, verbose=FALSE, cached=NULL, kernel=NULL){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells.idx <- metadata$sample != "gene"
  cells <- metadata$cell[cells.idx]
  groups <- metadata$group[cells.idx]

  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, -(1:2), drop=FALSE])
  genes <- unique(c(ligands, receptors[!is.na(receptors)]))
  check_genes(genes, rownames(infMat))
  ge <- gene_expression(adj, cells, genes)

  terms <- matrix(NA_real_, nrow(lrs), 5, dimnames=list(NULL, term_names()))
  pair.keys <- pair_keys(ligands, receptors)
  scored <- logical(nrow(lrs))
  other.terms <- NULL
  if(cache_on(cached) && !is.null(cached$terms)){
    hit <- match(pair.keys, rownames(cached$terms))
    scored <- !is.na(hit)
    terms[scored, ] <- cached$terms[hit[scored], , drop=FALSE]
    other.terms <- cached$terms[!(rownames(cached$terms) %in% pair.keys), , drop=FALSE]
    if(verbose && any(scored)){
      message("  Continuing from the cache in ", cached$dir, ": ", sum(scored), " of ", nrow(lrs), " pairs already scored")
    }
  }
  todo <- which(!scored)

  blocks <- NULL
  if(inherits(infMat, "hypercom_influence") && length(todo) > 0){
    blocks <- grid_blocks(infMat, grid_cells(rownames(infMat), metadata), threads, kernel)
  }

  n.grids <- nrow(pair_grids(ligands[todo], receptors[todo, , drop=FALSE]))
  size <- batch_size(threads, length(cells), tabulate(as.integer(factor(metadata$sample[cells.idx]))), n.grids)
  done <- sum(scored)
  grids.done <- 0
  t0 <- Sys.time()
  shown <- t0
  for(rows in batch_pairs(ligands[todo], receptors[todo, , drop=FALSE], size)){
    rows <- todo[rows]
    pairs <- pair_grids(ligands[rows], receptors[rows, , drop=FALSE])
    keys <-paste(pairs$ligand, pairs$receptor, sep="\t")
    sr <- score_grids(infMat, metadata, pairs$ligand, pairs$receptor, threads, blocks)

    for(i in rows){
      l <- ligands[i]
      rs <- stats::na.omit(receptors[i, ])
      cols <- match(paste(l, rs, sep="\t"), keys)

      pair.sr <- average_scores(sr, cols)
      ss <- pair.sr$sender
      rr <- pair.sr$receiver
      if(all(l == rs)){
        type <- rep(FALSE, length(ss))
      } else {
        type <- (rr - ss) > 0
      }

      pc <- pair_cells(c(l, rs), ge, groups, gini.max)
      max.score <- mean(pmax(ss, rr)[pc$present])
      balance <- gini_impurity(type[pc$present]) / 0.5

      terms[i, ] <- c(pc$frac, pc$pct.min, max.score, balance, pc$purity)
    }
    rm(sr)
    scored[rows] <- TRUE
    if(cache_on(cached)){
      keep <- scored & !duplicated(pair.keys)
      cache.terms <- terms[keep, , drop=FALSE]
      rownames(cache.terms) <- pair.keys[keep]
      write_cache_terms(cached, "terms", rbind(other.terms, cache.terms))
    }
    done <- done + length(rows)
    grids.done <- grids.done + nrow(pairs)
    if(verbose && (as.numeric(difftime(Sys.time(), shown, units="secs")) >= 5 | done == nrow(lrs))){
      left <- as.numeric(difftime(Sys.time(), t0, units="secs")) / grids.done * (n.grids - grids.done)
      message("  ", done, " of ", nrow(lrs), " pairs scored", if(done < nrow(lrs)) paste0(", about ", format_duration(left), " left"))
      shown <- Sys.time()
    }
  }

  return(terms)
}

#' Name the priority terms
#'
#' `term_names` returns the names of the five terms of a priority score, in their order
#' @returns A vector of strings
#' @noRd
term_names <- function(){
  return(c("frac", "pct.min", "max", "balance", "purity"))
}

#' Find the cells that express each gene
#'
#' `gene_expression` returns the fraction of cells expressing each gene and the cells that do
#'
#' This function takes the adjacency, the cells and the genes and returns, from the gene x cell counts of the adjacency, the fraction
#' of the cells expressing each gene and their positions. The counts of a symmetric sparse adjacency are read from its stored triangle,
#' without expanding the whole matrix.
#' @inheritParams prioritize_lr
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

#' Choose the number of grids per batch
#'
#' `batch_size` returns how many grids calculate_priority scores per kernel call
#'
#' This function takes the thread count, the number of cells, the cells of each sample and the number of grids and returns the batch
#' size: at least one grid per thread and at most about 64 MB of cell scores; when a grid is bigger than the kernel's 24 MB cache
#' budget (10 bytes per within-sample cell pair; CACHE_BUDGET in src/score_grids.cpp), all threads score one grid at a time, so smaller batches keep every thread busy and
#' there are about 20 batches or more, for progress messages. Each batch sets up the kernel's buffers again: on 20,000 cells, 104
#' pairs in 20 batches instead of 2 took 5% longer (90 instead of 86 s) and peaked 0.16 GB higher; the full CellChat database already
#' gives about 20 batches there.
#' @param threads An integer number of threads
#' @param cells An integer number of cells
#' @param sample.sizes An integer vector of the number of cells of each sample
#' @param n.grids An integer number of grids to score
#' @returns An integer number of grids per batch
#' @noRd
batch_size <- function(threads, cells, sample.sizes, n.grids){
  size <- max(threads, min(max(256, 32 * threads), floor(2^26 / (16 * cells))))
  if(10 * sum(as.numeric(sample.sizes)^2) > 24 * 2^20){
    size <- max(threads, min(size, ceiling(n.grids / 20)))
  }

  return(size)
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

#' Combine priority terms into priority scores
#'
#' `score_priority` returns the priority score of each ligand-receptor pair for one weight
#'
#' This function takes the priority terms of calculate_priority and returns their weighted geometric mean.
#' @param terms A matrix of priority terms from calculate_priority, one row per pair
#' @param weight A weight of at least 0 for the fraction of cells expressing a pair's genes (frac), the first of the five terms
#' @returns A numeric vector of priority scores
#' @noRd
score_priority <- function(terms, weight=1){
  w <- c(weight, rep(1, 4))
  priority <- numeric(nrow(terms))
  for(i in seq_len(nrow(terms))){
    priority[i] <- exp(stats::weighted.mean(log(terms[i, ]), w))
  }

  return(priority)
}

#' Check the weights of the priority scores
#'
#' `check_weights` stops unless the weights are distinct numbers of at least 0
#'
#' This function takes the weights argument of prioritize_lr and stops unless every weight is a finite number of at least 0 and no two
#' weights print alike (each names a result and an output file).
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_weights <- function(weights){
  if(!is.numeric(weights) || length(weights) == 0 || !all(is.finite(weights)) || any(weights < 0)){
    stop("weights must be numbers of at least 0, e.g. c(1, 0.5)")
  }
  if(anyDuplicated(as.character(weights)) > 0){
    stop("weights must differ: each names a result and an output file")
  }

  return(invisible(NULL))
}

#' Check the cache argument
#'
#' `check_cache` stops unless cache is TRUE or FALSE
#'
#' This function takes the cache argument of prioritize_lr and stops unless it is one TRUE or FALSE.
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_cache <- function(cache){
  if(!(isTRUE(cache) || isFALSE(cache))){
    stop("cache must be TRUE or FALSE")
  }

  return(invisible(NULL))
}

#' Check the memory argument
#'
#' `check_memory` stops unless memory is NULL or one positive number
#'
#' This function takes the memory argument of prioritize_lr and stops unless it is NULL or one positive number of gigabytes (Inf allowed).
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_memory <- function(memory){
  if(!is.null(memory) && !(is.numeric(memory) && length(memory) == 1 && !is.na(memory) && memory > 0)){
    stop("memory must be NULL or one positive number of gigabytes")
  }

  return(invisible(NULL))
}

#' Check the exact.priority argument
#'
#' `check_exact_priority` stops unless exact.priority is TRUE or FALSE
#'
#' This function takes the exact.priority argument of prioritize_lr and stops unless it is one TRUE or one FALSE.
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_exact_priority <- function(exact.priority){
  if(!(isTRUE(exact.priority) || isFALSE(exact.priority))){
    stop("exact.priority must be TRUE or FALSE")
  }

  return(invisible(NULL))
}

#' Name the priority file of a weight
#'
#' `priority_filename` returns the name of the file prioritize_lr writes the priorities of a weight to
#'
#' This function takes weights and returns the file name of each, "priority.<weight>.csv".
#' @param weights A vector of weights
#' @returns A vector of strings of file names
#' @noRd
priority_filename <- function(weights){
  return(paste0("priority.", weights, ".csv"))
}

#' Order ligand-receptor pairs by signaling strength and specificity
#'
#' `prioritize_lr` returns an ordered list of ligand-receptor pairs
#'
#' This function takes the ligand-receptor pairs, the adjacency matrix, a named influence matrix and the corresponding metadata
#' and ranks the pairs by a priority score: the weighted geometric mean of five terms of the cells that express a pair's genes (their
#' fraction of all cells, the smallest fraction of cells expressing one of its genes, the mean of each cell's larger score as sender
#' or receiver, how evenly they split into senders and receivers, and how much they concentrate in some groups). Cells are scored once
#' for all weights; parallel or threads use several CPU cores.
#' With significance, the p-value of a pair is the fraction of 1000 background pairs of random genes whose priority is at least as
#' high, so it moves in steps of 0.001 and is 0 when no background pair scored as high; p.adj is its Benjamini-Hochberg adjustment.
#' Background pairs without a priority (NaN, e.g. when some of their cells have no score) are left out, with a warning. Background
#' pairs are first screened: one pass over their cell grids bounds their priorities without ranking the cells, and only the few whose
#' bounds hold a real pair's priority are scored in full, so the p-values are those of scoring all of them, in less time. Real pairs
#' are screened too unless exact.priority is TRUE, and data whose grids do not fit in memory are scored in memory mode (memory).
#' With outdir set and cache=TRUE, the pairs scored so far are kept in `<outdir>/cache` while pairs are scored, and the folder is
#' removed once the priority files are written. If the run stops (time limit of a cluster job, crash, Esc), the same call run again
#' continues where it stopped, losing at most the last batch of pairs, with the same results as a run that never stopped on the same
#' type of computer (on another type the BLAS can change the last bits of a "woodbury" influence matrix, and priorities by about
#' 1e-13); without seed it also continues with the background pairs drawn before. Weights and threads may differ. A cache made from other inputs (or
#' damaged) is not used: it is moved to `<outdir>/cache.old`, with a warning. If the cache cannot be written, scoring goes on without
#' it, with a warning.
#' @param lrs A dataframe of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param infMat An influence matrix from generate_infMat (a hypercom_influence or a named dense matrix)
#' @param metadata A dataframe from generate_metadata, containing columns named cell, sample and group
#' @param weights A vector of weights of at least 0, one result each, for the fraction of cells expressing a pair's genes in the priority
#' score: 1 weighs it like the other four terms, less than 1 favors pairs of fewer cells, 0 leaves it out
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel: one thread per physical core (one per logical core on Linux, where R does
#' not tell them apart)
#' @param outdir A string naming the output directory, or NULL to write no files
#' @param threads An integer number of threads of at least 1, overriding parallel; more than the available cores warns and uses
#' them all
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state; the session's
#' random state is restored afterwards
#' @param verbose A boolean of whether to report progress as messages: the pairs scored about every 5 s with
#' the time left, and the top pair of each weight
#' @param cache A boolean of whether to keep the pairs scored so far in `<outdir>/cache`, so that a stopped run continues where it
#' stopped when run again; ignored when outdir is NULL
#' @param memory A number of gigabytes the run may use, NULL for 90% of the computer's memory, or Inf for no limit. When the exact
#' grid kernels (about 17.6 bytes per within-sample cell pair, the sum of squared sample sizes) do not fit, a "woodbury" influence matrix
#' is scored in memory mode (log-domain tiles, about 5.5 bytes per pair, priorities within ~1e-11 of exact); if that does not
#' fit either, the run stops with the estimate. On a cluster give the memory of the job: the scheduler may not stop a job that uses more
#' @param exact.priority A boolean of whether to score every real pair exactly. FALSE (with significance) screens the real pairs too and
#' scores exactly only those whose p-value or rank, or whose balance term, the bounds leave undecided: the other priorities are the
#' midpoints of their bounds (measured within ~1e-8 of exact, median ~1e-9), and p-values, p.adj and order are those of exact scoring, in the same
#' or less time. Without significance every pair is scored exactly
#' @returns A list of dataframes named by weight, each the ligand-receptor pairs ordered by their priority column, with columns p and
#' p.adj when significance is TRUE
#' @examples
#' counts <- filter_counts(alexanian2021.counts, cellchat.mm.lrdb, outfile=NULL)
#' lrs <- generate_lrs_table(counts, cellchat.mm.lrdb, outfile=NULL)
#' hyperedges <- generate_hyperedges(counts, lrs, outfile=NULL)
#' adj <- generate_adjacency_matrix(hyperedges, outfile=NULL)
#' infMat <- generate_infMat(adj, outfile=NULL)
#' metadata <- generate_metadata(infMat, alexanian2021.metadata, sample="sample", group="condition",
#'                               outfile=NULL)
#' # a few ligands keep the example quick; significance=TRUE adds p-values from 1000 background pairs
#' pairs <- lrs[lrs$ligand %in% c("Thbs1", "Comp", "Tgfb1"), ]
#' priority <- prioritize_lr(pairs, adj, infMat, metadata, weights=c(1, 0.5), outdir=NULL)
#' head(priority[["1"]])
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, outdir="6.priority",
                          threads=NULL, seed=NULL, verbose=TRUE, cache=TRUE, memory=NULL, exact.priority=FALSE){
  check_weights(weights)
  check_outdir(outdir, null.ok=TRUE)
  check_cache(cache)
  check_memory(memory)
  check_exact_priority(exact.priority)
  kernel <- choose_kernel(infMat, metadata, memory, verbose)
  cached <- NULL
  if(cache && !is.null(outdir)){
    cached <- open_cache(file.path(outdir, "cache"), adj, infMat, metadata, kernel$type)
  }
  if(significance){
    if(is.null(seed) && identical(colnames(cached$background), colnames(lrs))){
      background.lrs <- cached$background
    } else {
      if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir=globalenv(), inherits=FALSE)
        on.exit(restore_seed(old.seed), add=TRUE)
        set.seed(seed)
      }
      background.lrs <- generate_background_lrs(lrs, metadata)
    }
  }
  threads <- get_threads(parallel, threads)
  if(verbose){
    message("Scoring ", nrow(lrs), " ligand-receptor pairs", if(significance) paste0(" and screening ", nrow(background.lrs), " background pairs"),
            " on ", threads, " thread", if(threads > 1) "s")
  }
  if(!is.null(cached) && significance){
    cached$background <- background.lrs
  }

  screened <- significance && !exact.priority
  if(screened){
    sp <- screened_priority(lrs, background.lrs, weights, adj, infMat, metadata, threads, verbose, cached, kernel)
    bp <- sp$background
    if(verbose){
      message("  ", sp$refined.real, " of ", nrow(lrs), " real pairs scored exactly; the priorities of the others are the midpoints of ",
              "their bounds (p-values and order as exact scoring)")
    }
  } else {
    priority.terms <- calculate_priority(lrs, adj, infMat, metadata, threads, verbose, cached, kernel)
    if(significance){
      bp <- background_priority(background.lrs, lrs, priority.terms, weights, adj, infMat, metadata, threads, verbose, cached, kernel)
    }
  }
  if(significance){
    if(bp$no.priority > 0){
      warning(bp$no.priority, " of ", nrow(background.lrs), " background pairs have no priority (NaN) and are left out of the p-values, ",
              "which use the other ", nrow(background.lrs) - bp$no.priority)
    }
  }

  if(!is.null(outdir)){
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  }
  lrs.list <- list()

  for(i in seq_along(weights)){
    weight <- weights[i]
    priority <- if(screened) sp$priority[, i] else score_priority(priority.terms, weight)
    current.lrs <- lrs
    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      current.lrs$p <- vapply(current.lrs$priority, permutation_test, numeric(1), background=bp$priority[, i])
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method="BH")
    }

    outfile <- NULL
    if(!is.null(outdir)){
      outfile <- file.path(outdir, priority_filename(weight))
      utils::write.csv(current.lrs, outfile, row.names=FALSE)
    }
    if(verbose){
      message("Weight ", weight, ": top pair ", current.lrs[[1]][1], if(!is.null(outfile)) paste0(", written to ", outfile))
    }

    lrs.list[[i]] <- current.lrs
  }
  if(!is.null(cached)){
    close_cache(cached$dir)
  }

  names(lrs.list) <- weights

  return(lrs.list)
}
