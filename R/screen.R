# Background screening for the p-values of prioritize_lr: bounds of the priorities of background pairs from one screening pass over their
# grids (screen_grids in R/grids.R), so that only the background pairs whose priority could fall on either side of a real pair's priority
# are scored exactly. The p-values are the same as with every background pair scored.

#' Bound the priority terms of ligand-receptor pairs
#'
#' `bound_terms` returns lower and upper bounds of the five priority terms of each pair
#'
#' This function takes ligand-receptor pairs, the adjacency, a named influence matrix and the corresponding metadata and returns, for
#' each pair, bounds of the terms calculate_priority would return. frac, pct.min and purity are exact. The cell score bounds of
#' screen_grids, averaged over a complex's subunits as the scores are, are widened for rounding by 1e-12 plus the largest sample
#' times 2^-52, more than the rounding of the exact scores (means of up to that many percent ranks, summed in long double, or double on
#' some computers) and of the bounds (exact integer sums divided twice); the mean of each cell's larger score then lies between the
#' means of the larger lower and of the larger upper bounds, widened by the pair's cells times 2^-52 (R's mean), and the balance between
#' the smallest and largest Gini impurity of the possible numbers of receivers (cells whose receiver bounds lie above their sender
#' bounds are receivers, cells whose bounds overlap may be). A pair with a cell without a score bound gets NA. Pairs are screened in
#' batches of whole ligands, in rounds of one batch per foreach worker as calculate_priority scores them; with a cache, pairs whose
#' bounds are in it are not screened again, and the bounds of all pairs screened so far are written to it after each round. With
#' blocks of 2^28 entries or more, garbage is collected after each round: R holds the blocks, so it would otherwise let the score
#' bounds of many batches pile up before collecting.
#' @inheritParams calculate_priority
#' @returns A list of two matrices, lo and hi, with one row per pair and columns frac, pct.min, max, balance and purity
#' @noRd
bound_terms <- function(lrs, adj, infMat, metadata, verbose=FALSE, cached=NULL){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells.idx <- metadata$sample != "gene"
  cells <- metadata$cell[cells.idx]
  groups <- metadata$group[cells.idx]
  unit <- 2^-52
  eps <- 1e-12 + max(tabulate(as.integer(factor(metadata$sample[cells.idx])))) * unit

  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, -(1:2), drop=FALSE])
  genes <- unique(c(ligands, receptors[!is.na(receptors)]))
  check_genes(genes, rownames(infMat))
  ge <- gene_expression(adj, cells, genes)

  lo <- hi <- matrix(NA_real_, nrow(lrs), 5, dimnames=list(NULL, term_names()))
  pair.keys <- pair_keys(ligands, receptors)
  screened <- logical(nrow(lrs))
  other <- NULL
  if(cache_on(cached) && !is.null(cached$bounds)){
    hit <- match(pair.keys, rownames(cached$bounds$lo))
    screened <- !is.na(hit)
    lo[screened, ] <- cached$bounds$lo[hit[screened], , drop=FALSE]
    hi[screened, ] <- cached$bounds$hi[hit[screened], , drop=FALSE]
    rest <- !(rownames(cached$bounds$lo) %in% pair.keys)
    other <- list(lo=cached$bounds$lo[rest, , drop=FALSE], hi=cached$bounds$hi[rest, , drop=FALSE])
    if(verbose && any(screened)){
      message("  Continuing from the cache in ", cached$dir, ": ", sum(screened), " of ", nrow(lrs), " background pairs already screened")
    }
  }
  todo <- which(!screened)

  workers <- foreach::getDoParWorkers()
  sample.sizes <- tabulate(as.integer(factor(metadata$sample[cells.idx])))
  blocks <- NULL
  if(workers == 1 && length(todo) > 0){
    blocks <- grid_blocks(infMat, grid_cells(rownames(infMat), metadata))
  }
  n.grids <- nrow(pair_grids(ligands[todo], receptors[todo, , drop=FALSE]))
  size <- screen_batch_size(workers, length(cells), sample.sizes, n.grids)
  batches <- lapply(batch_pairs(ligands[todo], receptors[todo, , drop=FALSE], size), function(rows) todo[rows])
  done <- sum(screened)
  grids.done <- 0
  t0 <- Sys.time()
  shown <- t0
  for(round in split(seq_along(batches), ceiling(seq_along(batches) / workers))){
    grids <- lapply(batches[round], function(rows) pair_grids(ligands[rows], receptors[rows, , drop=FALSE]))
    sbs <- run_batches(infMat, metadata, grids, blocks, screen=TRUE)
    rows <- unlist(batches[round])
    keys <- unlist(lapply(seq_along(grids), function(b) paste(b, grids[[b]]$ligand, grids[[b]]$receptor, sep="\t")))
    sr.lo <- list(sender=do.call(cbind, lapply(sbs, function(sb) sb$sender.lo)), receiver=do.call(cbind, lapply(sbs, function(sb) sb$receiver.lo)))
    sr.hi <- list(sender=do.call(cbind, lapply(sbs, function(sb) sb$sender.hi)), receiver=do.call(cbind, lapply(sbs, function(sb) sb$receiver.hi)))
    pairs <- do.call(rbind, grids)
    batch.of <- rep(seq_along(round), lengths(batches[round]))
    rm(sbs)

    for(k in seq_along(rows)){
      i <- rows[k]
      l <- ligands[i]
      rs <- stats::na.omit(receptors[i, ])
      cols <- match(paste(batch.of[k], l, rs, sep="\t"), keys)
      pair.lo <- average_scores(sr.lo, cols)
      pair.hi <- average_scores(sr.hi, cols)
      pc <- pair_cells(c(l, rs), ge, groups, gini.max)

      sl <- pmax(0, pair.lo$sender[pc$present] - eps)
      rl <- pmax(0, pair.lo$receiver[pc$present] - eps)
      sh <- pair.hi$sender[pc$present] + eps
      rh <- pair.hi$receiver[pc$present] + eps
      if(anyNA(c(sl, rl, sh, rh))){
        max.score <- balance <- c(NA_real_, NA_real_)
      } else {
        max.score <- c(max(0, mean(pmax(sl, rl)) - length(sl) * unit), mean(pmax(sh, rh)) + length(sl) * unit)
        if(all(l == rs)){
          balance <- rep(gini_impurity(rep(FALSE, sum(pc$present))) / 0.5, 2)
        } else {
          balance <- bound_balance(sum(pc$present), sum(rl > sh), sum(rh > sl))
        }
      }

      lo[i, ] <- c(pc$frac, pc$pct.min, max.score[1], balance[1], pc$purity)
      hi[i, ] <- c(pc$frac, pc$pct.min, max.score[2], balance[2], pc$purity)
    }
    rm(sr.lo, sr.hi)
    if(sum(as.numeric(sample.sizes)^2) >= 2^28){
      gc()
    }
    screened[rows] <- TRUE
    if(cache_on(cached)){
      keep <- screened & !duplicated(pair.keys)
      cache.lo <- lo[keep, , drop=FALSE]
      cache.hi <- hi[keep, , drop=FALSE]
      rownames(cache.lo) <- rownames(cache.hi) <- pair.keys[keep]
      write_cache_terms(cached, "bounds", list(lo=rbind(other$lo, cache.lo), hi=rbind(other$hi, cache.hi)))
    }
    done <- done + length(rows)
    grids.done <- grids.done + nrow(pairs)
    if(verbose && (as.numeric(difftime(Sys.time(), shown, units="secs")) >= 5 | done == nrow(lrs))){
      left <- as.numeric(difftime(Sys.time(), t0, units="secs")) / grids.done * (n.grids - grids.done)
      message("  ", done, " of ", nrow(lrs), " background pairs screened", if(done < nrow(lrs)) paste0(", about ", format_duration(left), " left"))
      shown <- Sys.time()
    }
  }

  return(list(lo=lo, hi=hi))
}

#' Bound the balance term of a pair
#'
#' `bound_balance` returns the smallest and largest balance of n cells of which between t0 and t1 are receivers
#'
#' This function takes the number of cells and the range of the number of receivers among them and returns the range of the balance
#' (Gini impurity of senders and receivers / 0.5) over that range, each value made by gini_impurity as calculate_priority makes it. The
#' impurity is concave in the number of receivers, so the smallest is at an end of the range and the largest at an end or at the number
#' nearest n / 2.
#' @param n An integer number of cells
#' @param t0 An integer, the fewest receivers
#' @param t1 An integer, the most receivers
#' @returns A numeric vector of two, the smallest and the largest balance
#' @noRd
bound_balance <- function(n, t0, t1){
  balance <- function(t) gini_impurity(c(rep(FALSE, n - t), rep(TRUE, t))) / 0.5
  ends <- c(balance(t0), balance(t1))
  middle <- unique(pmin(t1, pmax(t0, c(floor(n / 2), ceiling(n / 2)))))

  return(c(min(ends), max(ends, vapply(middle, balance, numeric(1)))))
}

#' Choose the number of grids per screening batch
#'
#' `screen_batch_size` returns how many grids bound_terms screens per batch
#'
#' This function takes the number of workers, the number of cells, the cells of each sample and the number of grids and returns the
#' batch size as batch_size does, with at most about 64 MB of score bounds (four doubles per cell and grid).
#' @inheritParams batch_size
#' @returns An integer number of grids per batch
#' @noRd
screen_batch_size <- function(workers, cells, sample.sizes, n.grids){
  size <- min(256, floor(2^26 / (32 * cells)), ceiling(n.grids / workers))
  if(8 * sum(as.numeric(sample.sizes)^2) > 24 * 2^20){
    size <- min(size, ceiling(n.grids / (20 * workers)))
  }

  return(max(1, size))
}

#' Find the priority intervals that hold a real priority
#'
#' `ambiguous` returns which intervals (L, U] contain at least one real priority, or have an end that is not a number
#'
#' This function takes the lower and upper bounds of background priorities and the real priorities and returns, for each background
#' pair, whether its exact priority is needed: a pair whose interval holds no real priority x is above or below every x whatever its
#' exact priority. Real priorities that are NA are left out (their p-value is NA).
#' @param L A numeric vector of lower bounds
#' @param U A numeric vector of upper bounds
#' @param x A numeric vector of real priorities
#' @returns A logical vector
#' @noRd
ambiguous <- function(L, U, x){
  x <- sort(unique(x[!is.na(x)]))

  return(!is.finite(L) | !is.finite(U) | findInterval(U, x) > findInterval(L, x))
}

#' Calculate the priorities of background pairs for the p-values
#'
#' `background_priority` returns, for each weight, the priority of each background pair, or a lower bound of it that gives the same
#' p-values
#'
#' This function takes the background pairs, the real pairs and their priority terms, the weights and the inputs of prioritize_lr.
#' Background pairs whose terms are known (the same genes as a real pair, or in the cache) take them. The others are screened
#' (bound_terms), their priority bounds widened by 1e-12 relative for rounding, and those whose interval (L, U] holds a real priority
#' for some weight, or that have no bound, are scored exactly (calculate_priority). For every other pair each real priority x is at
#' most L or above U, so L is on the same side of x as the exact priority and permutation_test counts it the same.
#' @param background.lrs A dataframe of background pairs from generate_background_lrs
#' @inheritParams prioritize_lr
#' @param priority.terms A matrix of the priority terms of lrs from calculate_priority
#' @param verbose A boolean of whether to report progress
#' @param cached An environment from open_cache, or NULL for no cache
#' @returns A list with priority (a matrix with one row per background pair and one column per weight), no.priority (the number of
#' background pairs without terms, as calculate_priority finds them) and refined (the number of screened pairs scored exactly)
#' @noRd
background_priority <- function(background.lrs, lrs, priority.terms, weights, adj, infMat, metadata, verbose=FALSE, cached=NULL){
  margin <- 1e-12
  keys <- pair_keys(as.character(background.lrs$ligand), as.matrix(background.lrs[, -(1:2), drop=FALSE]))
  terms <- matrix(NA_real_, nrow(background.lrs), 5, dimnames=list(NULL, term_names()))
  hit <- match(keys, pair_keys(as.character(lrs$ligand), as.matrix(lrs[, -(1:2), drop=FALSE])))
  known <- !is.na(hit)
  terms[known, ] <- priority.terms[hit[known], , drop=FALSE]
  if(cache_on(cached) && !is.null(cached$terms)){
    hit <- match(keys, rownames(cached$terms))
    from.cache <- !known & !is.na(hit)
    terms[from.cache, ] <- cached$terms[hit[from.cache], , drop=FALSE]
    known <- known | from.cache
  }

  L <- matrix(NA_real_, nrow(background.lrs), length(weights))
  refine <- logical(nrow(background.lrs))
  screen <- which(!known)
  if(length(screen) > 0){
    b <- bound_terms(background.lrs[screen, , drop=FALSE], adj, infMat, metadata, verbose, cached)
    for(i in seq_along(weights)){
      x <- score_priority(priority.terms, weights[i])
      L[screen, i] <- score_priority(b$lo, weights[i]) * (1 - margin)
      U <- score_priority(b$hi, weights[i]) * (1 + margin)
      refine[screen] <- refine[screen] | ambiguous(L[screen, i], U, x)
    }
  }
  if(any(refine)){
    if(verbose){
      message("  Scoring the ", sum(refine), " background pairs whose priority could be on either side of a real pair's")
    }
    terms[refine, ] <- calculate_priority(background.lrs[refine, , drop=FALSE], adj, infMat, metadata, verbose, cached)
  }

  exact <- known | refine
  priority <- L
  for(i in seq_along(weights)){
    priority[exact, i] <- score_priority(terms[exact, , drop=FALSE], weights[i])
  }

  return(list(priority=priority, no.priority=sum(rowSums(is.na(terms[exact, , drop=FALSE])) > 0), refined=sum(refine)))
}
