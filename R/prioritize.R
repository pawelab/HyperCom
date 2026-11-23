# Pipeline step 6: priority scores of ligand-receptor pairs and their p-values

#' Generate a background set of random gene interactions
#'
#' `generate_background_lrs` returns a data frame of random gene interactions.
#'
#' This generates a table of random gene interactions where the number of genes in an interaction is equal to the complex size + 1.
#' @inheritParams prioritize_lr
#' @param count An integer of the number of interactions to generate
#' @returns A data frame of random gene interactions
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

#' Check a TRUE or FALSE argument
#'
#' `check_flag` stops unless an argument is one TRUE or one FALSE
#'
#' This function takes the value and the name of an argument such as significance, exact.priority, resume or verbose and stops,
#' naming it, unless it is one TRUE or one FALSE.
#' @param x The value of the argument
#' @param name A string naming the argument
#' @returns NULL, invisibly
#' @noRd
check_flag <- function(x, name){
  if(!(isTRUE(x) || isFALSE(x))){
    stop(name, " must be TRUE or FALSE")
  }

  return(invisible(NULL))
}

#' Check the memory.gb argument
#'
#' `check_memory` stops unless memory.gb is NULL or one positive number
#'
#' This function takes the memory.gb argument of prioritize_lr and stops unless it is NULL or one positive number of gigabytes (Inf
#' allowed).
#' @inheritParams prioritize_lr
#' @returns NULL, invisibly
#' @noRd
check_memory <- function(memory.gb){
  if(!is.null(memory.gb) && !(is.numeric(memory.gb) && length(memory.gb) == 1 && !is.na(memory.gb) && memory.gb > 0)){
    stop("memory.gb must be NULL or one positive number of gigabytes")
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

#' Rank ligand-receptor pairs by signaling strength and specificity
#'
#' `prioritize_lr` returns the ligand-receptor pairs ranked by a priority score, one table per weight
#'
#' This function takes the interactions, the adjacency, the influence matrix and the metadata made by the earlier steps (as
#' `HyperCom` does) and ranks the ligand-receptor pairs. Every cell gets a sender and a receiver score for each pair (see
#' `score_hypercom`); cells are scored once for all weights.
#' @section Priority score:
#' The priority of a pair is the weighted geometric mean of five terms of the cells that express any of its genes:
#'
#' - frac: their fraction of all cells (its weight is `weights`; the others have weight 1);
#' - pct.min: the smallest fraction of cells expressing one of the pair's genes;
#' - max: the mean of each cell's larger score, as sender or as receiver;
#' - balance: how evenly the cells split into senders and receivers;
#' - purity: how concentrated the cells are in some groups (the group column of metadata). Cells spread exactly evenly over the groups
#'   give purity 0, and so priority 0; without groups purity is 1.
#'
#' There is one ranking over all cells: groups enter only through purity, not as a ranking per group.
#' @section P-values:
#' With significance=TRUE, 1000 background pairs of random genes of the hypergraph are scored like the real pairs. The p-value of a
#' pair (column p) is the fraction of background pairs whose priority is at least as high, so it moves in steps of 0.001 and 0 means
#' none (report p < 0.001); p.adj is p with the Benjamini-Hochberg adjustment. Background pairs without a priority (NaN) are left
#' out, with a warning.
#'
#' To save time, pairs are first screened: a quick pass gives lower and upper bounds of each priority, and only the pairs whose rank
#' or p-value the bounds cannot settle are scored fully. The priorities of the other real pairs are estimates (the middle of their
#' bounds): typically within ~1e-9 of exact, up to ~1e-7 on the bundled data (1.5e-8 measured on the full bundled run). P-values,
#' p.adj and order are exact. exact.priority=TRUE gives exact values.
#' @section Memory:
#' Memory grows with the number of pairs of cells in the same sample: the sum over samples of the squared number of cells (one
#' sample of 30,000 cells: 9e8 pairs). Exact scoring needs about 17.6 bytes per pair (16 bytes plus 10%), plus about 1.1 GB for the
#' rest of the run: about 17 GB for one sample of 30,000 cells. When that is more than memory.gb, the pairs are scored in memory
#' mode, which needs about 5.5 bytes per pair (plus the same 1.1 GB) and gives priorities within ~1e-11 of exact scoring; memory
#' mode needs an influence matrix made with the default method "woodbury". When neither fits, the run stops at the start, giving
#' the estimate. Splitting the cells into more samples saves memory but changes which cells are compared, and so the results.
#'
#' The estimate is an upper bound with a floor of about 1.1 GB (1 GiB for the rest of the run): small runs report about 1.08 GB
#' while using less, and memory.gb below that always stops.
#' @section Stopped runs:
#' With outdir set and resume=TRUE, the pairs scored so far are kept in `<outdir>/cache`. If the run stops (time limit of a cluster
#' job, crash, Esc), run the same call again: it continues where it stopped, losing at most the last batch of pairs. The results
#' are those of a run that never stopped on the same type of computer (another type of CPU can change the last digits, about
#' 1e-13). Without seed, the run keeps the background pairs drawn before it stopped. weights and threads may differ between the
#' calls. The cache folder is removed once the priority files are written.
#'
#' A cache is not used when the run differs in its data (adj, infMat, or the cells' sample or group in metadata; through HyperCom
#' this means counts, lrdb, metadata, r or infmat.method), in its kind of scoring (exact or memory mode, which follows memory.gb and
#' the computer's memory) or in its version of HyperCom, or when the cache is damaged. It is then moved to `<outdir>/cache.old`,
#' with a warning that says why: rename it back to continue the run it belongs to. If the cache cannot be written (e.g. a full
#' disk), scoring goes on without it, with a warning.
#' @param lrs A data frame of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @param adj An adjacency from generate_adjacency_matrix, or a named adjacency matrix
#' @param infMat An influence matrix from generate_infMat (a hypercom_influence or a named dense matrix)
#' @param metadata A data frame from generate_metadata, with columns cell, sample and group and one row per node of infMat (cells and
#' genes; the p-values draw background pairs from the gene rows)
#' @param weights A vector of numbers of at least 0, one result table each: the weight of frac, the fraction of cells expressing a
#' pair's genes, in the priority score. 1 (the default) weighs it like the other four terms, less than 1 favors pairs of fewer cells,
#' 0 leaves it out
#' @param significance A boolean of whether to add p-values (columns p and p.adj) from 1000 background pairs; FALSE by default
#' @param seed An integer seed for drawing the background pairs, so that p-values are reproducible, or NULL (the default) to use the
#' session's random numbers; the session's random state is restored afterwards. Ignored without significance
#' @param exact.priority A boolean of whether to score every real pair exactly when significance is TRUE. FALSE (the default) screens
#' the real pairs, in the same or less time: priorities typically within ~1e-9 of exact, up to ~1e-7 on the bundled data (1.5e-8
#' measured on the full bundled run; more with a very small lrdb, e.g. 3e-5 with 15 pairs); p-values, p.adj and order exact. TRUE
#' gives exact values. Without significance every pair is scored exactly
#' @param threads An integer number of threads, NULL (the default) for one per physical core, or 1 to run on one thread; a number above
#' the cores available to R warns and uses them all. The threads run at a lower priority than other programs, so the computer stays
#' responsive while HyperCom uses all cores. A cluster job or container gets only its allocated cores, and `options(mc.cores=4)` (or
#' the environment variable MC_CORES) caps the default. Results do not depend on it
#' @param memory.gb A number of gigabytes (10^9 bytes) the run may use, NULL (the default) for 90% of the computer's memory (on Linux,
#' of a container's or cluster job's memory limit when that is smaller), or Inf for no limit. It decides between exact scoring and
#' memory mode (see the Memory section of `prioritize_lr`). On a cluster give the job's memory request, as the computer's memory may
#' be more than the job may use
#' @param resume A boolean of whether to keep the pairs scored so far in `<outdir>/cache`, so that a stopped run continues where it
#' stopped when the same call is run again (see Stopped runs); TRUE by default, ignored when outdir is NULL
#' @param verbose A boolean of whether to report progress as messages: the memory needed, the pairs scored about every 5 s with the
#' time left, and the top pair of each weight
#' @param outdir A string naming the output directory for the files priority.<weight>.csv (and the cache), or NULL to write no files
#' @returns A list of data frames named by weight (e.g. "1"), each with the interactions of lrs ordered by decreasing priority: the
#' columns of lrs, priority, and with significance p and p.adj
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
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=1, significance=FALSE, seed=NULL, exact.priority=FALSE, threads=NULL,
                          memory.gb=NULL, resume=TRUE, verbose=TRUE, outdir="6.priority"){
  check_weights(weights)
  check_flag(significance, "significance")
  check_flag(exact.priority, "exact.priority")
  check_flag(resume, "resume")
  check_flag(verbose, "verbose")
  check_memory(memory.gb)
  check_outdir(outdir)
  check_metadata(metadata, significance)
  lrs <- clean_lrdb(lrs, "lrs", empty.ok=TRUE)
  in.graph <- lrs_in_graph(lrs, rownames(infMat))
  if(nrow(lrs) > 0 && !any(in.graph)){
    stop("no interaction of lrs has all its genes in infMat: make lrs, adj and infMat from the same counts")
  }
  if(!all(in.graph)){
    warning(sum(!in.graph), " interactions of lrs have genes that are not nodes of infMat (not connected to the rest of the hypergraph) ",
            "and are left out")
    lrs <- lrs[in.graph, , drop=FALSE]
  }
  kernel <- choose_kernel(infMat, metadata, memory.gb, verbose)
  cached <- NULL
  if(resume && !is.null(outdir)){
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
  threads <- get_threads(threads)
  if(verbose){
    message("Scoring ", nrow(lrs), " ligand-receptor pairs", if(significance) paste0(" and screening ", nrow(background.lrs), " background pairs"),
            " on ", threads, " thread", if(threads > 1) "s")
  }
  if(!is.null(cached) && significance){
    cached$background <- background.lrs
  }

  if(significance){
    sp <- screened_priority(lrs, background.lrs, weights, adj, infMat, metadata, threads, verbose, cached, kernel, exact.priority)
    if(verbose && !exact.priority){
      message("  ", sp$refined.real, " of ", nrow(lrs), " real pairs needed full scoring; the others' priorities are estimated, ",
              "typically within ~1e-9, with exact p-values and order (exact.priority=TRUE scores all)")
    }
    if(sp$no.priority > 0){
      warning(sp$no.priority, " of ", nrow(background.lrs), " background pairs have no priority (NaN) and are left out of the p-values, ",
              "which use the other ", nrow(background.lrs) - sp$no.priority)
    }
  } else {
    terms <- calculate_priority(lrs, adj, infMat, metadata, threads, verbose, cached, kernel)
  }

  if(!is.null(outdir)){
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  }
  lrs.list <- list()

  for(i in seq_along(weights)){
    weight <- weights[i]
    current.lrs <- lrs
    current.lrs$priority <- if(significance) sp$priority[, i] else score_priority(terms, weight)
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(.data$priority))

    if(significance){
      current.lrs$p <- vapply(current.lrs$priority, permutation_test, numeric(1), background=sp$background[, i])
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
