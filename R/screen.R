# Screening for prioritize_lr with p-values: bounds of the priorities of background and real pairs from one screening pass over their
# grids (src/screen_grids.cpp), so that only the pairs whose priority could fall on either side of another pair's are scored exactly.
# The p-values and the order are those of scoring every pair exactly.

#' Bound the sender and receiver scores of cell grids
#'
#' `screen_grids` returns lower and upper bounds of the sender and receiver scores of cells for ligand and receptor subunit pairs
#'
#' This function takes a named influence matrix, the corresponding metadata, vectors of ligands and receptor subunits and the
#' within-sample blocks, and returns for each grid bounds of the scores score_grids would return, from one screening pass without
#' ranking: each entry is only placed in a value bucket, which bounds its rank. Every grid needs a ligand and a receptor (the network
#' score is not screened). Teams of threads screen one grid each at a time, as small teams as the bucket ids of all teams (2 bytes per
#' within-sample cell pair each) fit in the bytes of the blocks (8 per within-sample cell pair; 2 in memory mode, beside the tiles'
#' 3), and at least 128 MB, so that screening never needs more memory than scoring; the bounds are the same for any team size.
#' @inheritParams score_grids
#' @param ligands A vector of strings of ligand names
#' @param receptors A vector of strings of receptor subunit names
#' @param blocks A numeric vector of the within-sample blocks, or a hypercom_lambda (memory mode), from grid_blocks
#' @param keep.ids A boolean of whether to keep the bucket ids for the next call (freed by screen_release_cpp), as new memory is slow
#' to first touch
#' @returns A list of four matrices, sender.lo, sender.hi, receiver.lo and receiver.hi, with one row per cell and one column per grid
#' @noRd
screen_grids <- function(infMat, metadata, ligands, receptors, threads=1, blocks, keep.ids=FALSE){
  if(anyNA(ligands) || anyNA(receptors)){
    stop("screening needs a ligand and a receptor for every grid: the network score (no ligand or receptor) is scored, not screened")
  }
  gi <- grid_inputs(infMat, metadata, ligands, receptors, threads, blocks)
  tiles <- inherits(blocks, "hypercom_lambda")
  id.budget <- max(2^27, (if(tiles) 2 else 8) * block_entries(blocks))
  if(tiles){
    sb <- screen_lambda_cpp(blocks$bytes, blocks$base, blocks$toff, blocks$D, gi$A, gi$B, gi$g$sample, gi$lc, gi$rc, threads, id.budget,
                            keep.ids, blocks$p)
  } else {
    sb <- screen_blocks_cpp(blocks, gi$A, gi$B, gi$g$sample, gi$lc, gi$rc, threads, id.budget, keep.ids)
  }

  return(metadata_order(list(sender.lo=sb$sender_lo, sender.hi=sb$sender_hi, receiver.lo=sb$receiver_lo, receiver.hi=sb$receiver_hi),
                        gi$g))
}

#' Bound the priority terms of ligand-receptor pairs
#'
#' `bound_terms` returns lower and upper bounds of the five priority terms of each pair
#'
#' This function takes ligand-receptor pairs, the adjacency, a named influence matrix and the corresponding metadata and returns, for
#' each pair, bounds of the terms calculate_priority would return, from screening (pair_terms, term_bounds). A pair with a cell
#' without a score bound gets NA.
#' @inheritParams pair_terms
#' @returns A list of two matrices, lo and hi, with one row per pair and columns frac, pct.min, max, balance and purity
#' @noRd
bound_terms <- function(lrs, adj, infMat, metadata, threads=1, verbose=FALSE, cached=NULL, kernel=NULL, label="background"){
  return(pair_terms(lrs, adj, infMat, metadata, threads, verbose, cached, kernel, bounds=TRUE, label=label))
}

#' Bound the priority terms of a pair from its cells' score bounds
#'
#' `term_bounds` returns lower and upper bounds of the five priority terms of a ligand-receptor pair
#'
#' This function takes the score bounds of grids, the columns of the pair's subunit grids, whether its ligand is its receptor, its
#' cells from pair_cells and the rounding widening, and returns bounds of exact_terms: frac, pct.min and purity are exact. The score
#' bounds, averaged over a complex's subunits as the scores are, are widened by eps (1e-12 plus the largest sample times 2^-52, more
#' than the rounding of the kernels' divisions); the mean of each cell's larger score then lies between the means of the larger lower
#' and of the larger upper bounds, widened by the pair's cells times 2^-52 (R's mean), and the balance between the smallest and
#' largest Gini impurity of the possible numbers of receivers (cells whose receiver bounds lie above their sender bounds are
#' receivers, cells whose bounds overlap may be). A cell without a score bound gives NA for max and balance.
#' @param sb A list from screen_grids
#' @param cols An integer vector of the columns of the pair's subunit grids
#' @param same A boolean of whether the ligand is the pair's only receptor subunit
#' @param pc A list from pair_cells
#' @param eps A double, the widening of each score bound
#' @returns A list of two numeric vectors of five terms, lo and hi
#' @noRd
term_bounds <- function(sb, cols, same, pc, eps){
  unit <- 2^-52
  pair.lo <- average_scores(list(sender=sb$sender.lo, receiver=sb$receiver.lo), cols)
  pair.hi <- average_scores(list(sender=sb$sender.hi, receiver=sb$receiver.hi), cols)
  sl <- pmax(0, pair.lo$sender[pc$present] - eps)
  rl <- pmax(0, pair.lo$receiver[pc$present] - eps)
  sh <- pair.hi$sender[pc$present] + eps
  rh <- pair.hi$receiver[pc$present] + eps
  if(anyNA(c(sl, rl, sh, rh))){
    max.score <- balance <- c(NA_real_, NA_real_)
  } else {
    max.score <- c(max(0, mean(pmax(sl, rl)) - length(sl) * unit), mean(pmax(sh, rh)) + length(sl) * unit)
    if(same){
      balance <- rep(gini_impurity(rep(FALSE, sum(pc$present))) / 0.5, 2)
    } else {
      balance <- bound_balance(sum(pc$present), sum(rl > sh), sum(rh > sl))
    }
  }

  return(list(lo=c(pc$frac, pc$pct.min, max.score[1], balance[1], pc$purity), hi=c(pc$frac, pc$pct.min, max.score[2], balance[2],
                                                                                    pc$purity)))
}

#' Bound the balance term of a pair
#'
#' `bound_balance` returns the smallest and largest balance of n cells of which between t0 and t1 are receivers
#'
#' This function takes the number of cells and the range of the number of receivers among them and returns the range of the balance
#' (Gini impurity of senders and receivers / 0.5) over that range, each value made by gini_impurity as exact_terms makes it. The
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

#' Widen priority bounds for rounding
#'
#' `widen_bounds` returns priority bounds widened by 1e-12 relative, unless they are equal (an exact priority)
#' @param lower A numeric vector of lower bounds
#' @param upper A numeric vector of upper bounds
#' @returns A list with lower, upper and point (whether the bounds are equal and finite)
#' @noRd
widen_bounds <- function(lower, upper){
  margin <- 1e-12
  point <- is.finite(lower) & lower == upper

  return(list(lower=ifelse(point, lower, lower * (1 - margin)), upper=ifelse(point, upper, upper * (1 + margin)), point=point))
}

#' Find the priority intervals that hold a real priority
#'
#' `ambiguous` returns which intervals (lower, upper] contain at least one real priority, or have an end that is not a number
#'
#' This function takes the lower and upper bounds of background priorities and the real priorities and returns, for each background
#' pair, whether its exact priority is needed: a pair whose interval holds no real priority x is above or below every x whatever its
#' exact priority. Real priorities that are NA are left out (their p-value is NA).
#' @param lower A numeric vector of lower bounds
#' @param upper A numeric vector of upper bounds
#' @param x A numeric vector of real priorities
#' @returns A logical vector
#' @noRd
ambiguous <- function(lower, upper, x){
  x <- sort(unique(x[!is.na(x)]))

  return(!is.finite(lower) | !is.finite(upper) | findInterval(upper, x) > findInterval(lower, x))
}

#' Find the intervals that meet another set of intervals
#'
#' `meets` returns, for each interval [lower1, upper1], whether it meets at least one interval [lower2, upper2] (closed ends)
#'
#' This function takes two sets of intervals and returns, for each of the first, whether any of the second overlaps or touches it,
#' from the second set sorted by lower end and the running largest upper end.
#' @param lower1 A numeric vector of lower ends
#' @param upper1 A numeric vector of upper ends
#' @param lower2 A numeric vector of lower ends
#' @param upper2 A numeric vector of upper ends
#' @returns A logical vector, one per interval of the first set
#' @noRd
meets <- function(lower1, upper1, lower2, upper2){
  if(length(lower2) == 0){
    return(logical(length(lower1)))
  }
  o <- order(lower2)
  lower2 <- lower2[o]
  upper2 <- cummax(upper2[o])
  k <- findInterval(upper1, lower2)

  return(k > 0 & upper2[pmax(k, 1)] >= lower1)
}

#' Find the intervals that meet another interval of the same set
#'
#' `meets_other` returns, for each interval [lower, upper], whether another interval of the set overlaps or touches it (closed ends)
#' @param lower A numeric vector of lower ends
#' @param upper A numeric vector of upper ends
#' @returns A logical vector
#' @noRd
meets_other <- function(lower, upper){
  hit <- logical(length(lower))
  if(length(lower) < 2){
    return(hit)
  }
  o <- order(lower)
  before <- c(-Inf, cummax(upper[o])[-length(o)])
  after <- c(lower[o][-1], Inf)
  hit[o] <- before >= lower[o] | after <= upper[o]

  return(hit)
}

#' Calculate the priorities of real and background pairs from screening
#'
#' `screened_priority` returns, for each weight, the priority of each real pair (exact, or the midpoint of its bounds) and of each
#' background pair (exact, or a bound on the same side of every real priority), with the p-values and the order of exact scoring
#'
#' This function takes the real pairs, the background pairs, the weights and the inputs of prioritize_lr. Pairs whose terms are in the
#' cache take them. Without exact.real the real pairs are screened (bound_terms), with it they are scored exactly (calculate_priority).
#' The background pairs take the bounds or terms of a real pair with the same genes, and the others are screened; those without
#' finite bounds are scored at once, since their priority must be known before real pairs are decided. Priority bounds are widened by
#' 1e-12 relative for rounding (bounds that are equal are the exact priority). A screened real pair is scored exactly when its bounds
#' are not finite, its interval meets a background pair's (its p-value could change) or another real pair's (the order could change),
#' or its balance is undecided: a cell whose sender and receiver bounds overlap may count either way, and the midpoint can then be far
#' from the exact priority. Every other real pair's interval meets no other, so its midpoint gives the exact p-value and order; it is
#' typically within ~1e-9 relative of the exact priority (bounds are ~1e-6 wide; the midpoint is far closer). Last,
#' background pairs whose interval (L, U] holds an exact real priority are scored exactly; for every other one each real priority x is
#' at most L or above U, so L is on the same side of x as the exact priority and permutation_test counts it the same.
#' @inheritParams prioritize_lr
#' @param background.lrs A data frame of background pairs from generate_background_lrs
#' @param threads An integer number of threads
#' @param cached An environment from open_cache, or NULL for no cache
#' @param kernel An environment from choose_kernel, or NULL for the exact kernels
#' @param exact.real A boolean of whether to score every real pair exactly
#' @returns A list with priority (a matrix with one row per real pair and one column per weight), refined.real (the number of
#' screened real pairs scored exactly), background (a matrix with one row per background pair and one column per weight) and
#' no.priority (the number of background pairs without a priority, as calculate_priority finds them)
#' @noRd
screened_priority <- function(lrs, background.lrs, weights, adj, infMat, metadata, threads=1, verbose=FALSE, cached=NULL, kernel=NULL,
                              exact.real=FALSE){
  real.keys <- lrs_keys(lrs)
  if(exact.real){
    terms <- calculate_priority(lrs, adj, infMat, metadata, threads, verbose, cached, kernel)
    known <- rep(TRUE, nrow(lrs))
  } else {
    terms <- empty_terms(nrow(lrs))
    hit <- cached_rows(cached, real.keys)
    known <- !is.na(hit)
    if(any(known)){
      terms[known, ] <- cached$terms[hit[known], , drop=FALSE]
    }
  }
  lo <- hi <- terms
  if(any(!known)){
    rb <- bound_terms(lrs[!known, , drop=FALSE], adj, infMat, metadata, threads, verbose, cached, kernel, label="real")
    lo[!known, ] <- rb$lo
    hi[!known, ] <- rb$hi
  }

  keys <- lrs_keys(background.lrs)
  twin <- match(keys, real.keys)
  bg.known <- !is.na(twin)
  bg.lo <- bg.hi <- bg.terms <- empty_terms(nrow(background.lrs))
  bg.lo[bg.known, ] <- lo[twin[bg.known], , drop=FALSE]
  bg.hi[bg.known, ] <- hi[twin[bg.known], , drop=FALSE]
  hit <- cached_rows(cached, keys)
  bg.exact <- !bg.known & !is.na(hit)
  if(any(bg.exact)){
    bg.terms[bg.exact, ] <- cached$terms[hit[bg.exact], , drop=FALSE]
    bg.lo[bg.exact, ] <- bg.hi[bg.exact, ] <- bg.terms[bg.exact, ]
  }
  bg.known <- bg.known | bg.exact
  screen <- which(!bg.known)
  if(length(screen) > 0){
    bb <- bound_terms(background.lrs[screen, , drop=FALSE], adj, infMat, metadata, threads, verbose, cached, kernel)
    bg.lo[screen, ] <- bb$lo
    bg.hi[screen, ] <- bb$hi
  }
  odd <- !bg.known & (rowSums(!is.finite(bg.lo)) > 0 | rowSums(!is.finite(bg.hi)) > 0)
  if(any(odd)){
    bg.terms[odd, ] <- calculate_priority(background.lrs[odd, , drop=FALSE], adj, infMat, metadata, threads, verbose, cached, kernel)
    bg.lo[odd, ] <- bg.hi[odd, ] <- bg.terms[odd, ]
  }
  bg.exact <- bg.exact | odd

  need <- !known & !(hi[, "balance"] <= lo[, "balance"])
  if(!all(known)){
    for(w in weights){
      r <- widen_bounds(score_priority(lo, w), score_priority(hi, w))
      b <- widen_bounds(score_priority(bg.lo, w), score_priority(bg.hi, w))
      fin <- is.finite(r$lower) & is.finite(r$upper)
      bfin <- is.finite(b$lower) & is.finite(b$upper)
      other <- logical(nrow(lrs))
      other[fin] <- meets_other(r$lower[fin], r$upper[fin])
      need <- need | (!known & (!fin | (!r$point & (other | meets(r$lower, r$upper, b$lower[bfin], b$upper[bfin])))))
    }
  }
  need[is.na(need)] <- TRUE
  if(any(need)){
    if(verbose){
      message("  Fully scoring the ", sum(need), " real pairs whose rank or p-value screening could not settle")
    }
    terms[need, ] <- calculate_priority(lrs[need, , drop=FALSE], adj, infMat, metadata, threads, verbose, cached, kernel)
    lo[need, ] <- hi[need, ] <- terms[need, ]
  }
  exact <- known | need
  from.real <- !is.na(twin) & exact[pmax(twin, 1)]
  bg.terms[from.real, ] <- terms[twin[from.real], , drop=FALSE]
  bg.exact <- bg.exact | from.real

  refine <- logical(nrow(background.lrs))
  priority <- matrix(NA_real_, nrow(lrs), length(weights))
  bg.priority <- matrix(NA_real_, nrow(background.lrs), length(weights))
  for(i in seq_along(weights)){
    plo <- score_priority(lo, weights[i])
    phi <- score_priority(hi, weights[i])
    point <- is.finite(plo) & plo == phi
    priority[, i] <- (plo + phi) / 2
    priority[point, i] <- plo[point]
    priority[exact, i] <- score_priority(terms[exact, , drop=FALSE], weights[i])
    b <- widen_bounds(score_priority(bg.lo, weights[i]), score_priority(bg.hi, weights[i]))
    bg.priority[, i] <- b$lower
    refine <- refine | (!bg.exact & !b$point & ambiguous(b$lower, b$upper, priority[exact | point, i]))
  }
  if(any(refine)){
    if(verbose){
      message("  Fully scoring the ", sum(refine), " background pairs that screening could not place above or below a real pair")
    }
    bg.terms[refine, ] <- calculate_priority(background.lrs[refine, , drop=FALSE], adj, infMat, metadata, threads, verbose, cached, kernel)
  }
  bg.done <- bg.exact | refine
  for(i in seq_along(weights)){
    bg.priority[bg.done, i] <- score_priority(bg.terms[bg.done, , drop=FALSE], weights[i])
  }

  return(list(priority=priority, refined.real=sum(need), background=bg.priority,
              no.priority=sum(rowSums(is.na(bg.terms[bg.done, , drop=FALSE])) > 0)))
}
