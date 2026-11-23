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
#' subunit grid is scored once, whatever the number of pairs using it, and only one batch of cell scores per worker (about 64 MB at
#' most) is in memory at a time. With a foreach backend of several workers registered (prioritize_lr registers cl), the workers score
#' one batch each at a time, in rounds. The gene x cell counts of a symmetric sparse adjacency are read from its stored triangle,
#' without expanding the whole matrix. With a cache, pairs whose genes are in it take its terms and are not scored again, and the
#' terms of all pairs scored so far (and of cached pairs of other calls) are written to it after each round; the first write of a
#' run also saves the inputs (save_cache). After a failed write the run goes on without the cache.
#' @inheritParams prioritize_lr
#' @param verbose A boolean of whether to report the pairs scored, about every 5 s, with an estimate of the time left (pairs are
#' scored in about 20 rounds or more when grids are big; see batch_size)
#' @param cached An environment from open_cache, or NULL for no cache
#' @returns A matrix with one row per pair and columns frac, pct.min, max, balance and purity
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, verbose=FALSE, cached=NULL){
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

  workers <- foreach::getDoParWorkers()
  blocks <- NULL
  if(workers == 1 && length(todo) > 0){
    blocks <- grid_blocks(infMat, grid_cells(rownames(infMat), metadata))
  }

  n.grids <- nrow(pair_grids(ligands[todo], receptors[todo, , drop=FALSE]))
  size <- batch_size(workers, length(cells), tabulate(as.integer(factor(metadata$sample[cells.idx]))), n.grids)
  batches <- lapply(batch_pairs(ligands[todo], receptors[todo, , drop=FALSE], size), function(rows) todo[rows])
  done <- sum(scored)
  grids.done <- 0
  t0 <- Sys.time()
  shown <- t0
  for(round in split(seq_along(batches), ceiling(seq_along(batches) / workers))){
    grids <- lapply(batches[round], function(rows) pair_grids(ligands[rows], receptors[rows, , drop=FALSE]))
    srs <- run_batches(infMat, metadata, grids, blocks)

    for(b in seq_along(round)){
      rows <- batches[[round[b]]]
      sr <- srs[[b]]
      keys <- paste(grids[[b]]$ligand, grids[[b]]$receptor, sep="\t")
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
    }
    rm(srs, sr)
    rows <- unlist(batches[round])
    pairs <- do.call(rbind, grids)
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
#' `batch_size` returns how many grids calculate_priority scores per batch
#'
#' This function takes the number of workers, the number of cells, the cells of each sample and the number of grids and returns the
#' batch size: at most 256 grids and about 64 MB of cell scores, few enough grids that every worker gets a batch, and, when grids are
#' big (more than about 3 million within-sample cell pairs, 24 MB at 8 bytes each), about 20 rounds of batches or more, for progress
#' messages and the cache.
#' @param workers An integer number of workers
#' @param cells An integer number of cells
#' @param sample.sizes An integer vector of the number of cells of each sample
#' @param n.grids An integer number of grids to score
#' @returns An integer number of grids per batch
#' @noRd
batch_size <- function(workers, cells, sample.sizes, n.grids){
  size <- min(256, floor(2^26 / (16 * cells)), ceiling(n.grids / workers))
  if(8 * sum(as.numeric(sample.sizes)^2) > 24 * 2^20){
    size <- min(size, ceiling(n.grids / (20 * workers)))
  }

  return(max(1, size))
}

#' Score or screen batches of grids
#'
#' `run_batches` returns the cell scores (score_grids) or their bounds (screen_grids) of each batch of grids
#'
#' This function takes an influence matrix, the metadata, batches of grids and the within-sample blocks of the influence matrix and
#' scores or screens each batch. With one foreach worker registered the batches are done in turn on the blocks given; with several,
#' foreach gives each worker a batch, which makes its own blocks, so only the influence matrix and the metadata are sent to the
#' workers.
#' @inheritParams prioritize_lr
#' @param batches A list of dataframes from pair_grids, with columns ligand and receptor
#' @param blocks A list from grid_blocks, or NULL to make it for each batch
#' @param screen A boolean of whether to screen (screen_grids) instead of score (score_grids)
#' @returns A list of the results of score_grids or screen_grids, one per batch
#' @noRd
run_batches <- function(infMat, metadata, batches, blocks=NULL, screen=FALSE){
  grid_fun <- if(screen) screen_grids else score_grids
  if(foreach::getDoParWorkers() == 1){
    return(lapply(batches, function(b) grid_fun(infMat, metadata, b$ligand, b$receptor, blocks)))
  }

  return(foreach::`%dopar%`(foreach::foreach(b=batches, .packages="HyperCom"), {
    grid_fun(infMat, metadata, b$ligand, b$receptor)
  }))
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
#' for all weights; parallel scores batches of grids on a cluster of R processes (cl).
#' With significance, the p-value of a pair is the fraction of 1000 background pairs of random genes whose priority is at least as
#' high, so it moves in steps of 0.001 and is 0 when no background pair scored as high; p.adj is its Benjamini-Hochberg adjustment.
#' Background pairs without a priority (NaN, e.g. when some of their cells have no score) are left out, with a warning. Background
#' pairs are first screened: one pass over their cell grids bounds their priorities without ranking the cells, and only the few whose
#' bounds hold a real pair's priority are scored in full, so the p-values are those of scoring all of them, in less time.
#' With outdir set and cache=TRUE, the pairs scored so far are kept in `<outdir>/cache` while pairs are scored, and the folder is
#' removed once the priority files are written. If the run stops (time limit of a cluster job, crash, Esc), the same call run again
#' continues where it stopped, losing at most the last batch of pairs, with the same results as a run that never stopped on the same
#' type of computer (on another type the BLAS can change the last bits of a "woodbury" influence matrix, and priorities by about
#' 1e-13); without seed it also continues with the background pairs drawn before. Weights and the cluster may differ. A cache made from other inputs (or
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
#' @param parallel A boolean of whether to run in parallel on the cluster cl
#' @param cl A cluster object, which is left running. Default if parallel: inside an SGE parallel job (PE_HOSTFILE set)
#' `parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd = "qrsh", rshopts = c("-inherit", "-nostdin", "-V"), outfile = "")`,
#' elsewhere `parallelly::makeClusterPSOCK(parallelly::availableCores(), outfile = "")`; it is stopped at the end. Each worker holds
#' the influence matrix and the within-sample blocks of the batch it scores (8 bytes per within-sample cell pair)
#' @param outdir A string naming the output directory, or NULL to write no files
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state; the session's
#' random state is restored afterwards
#' @param verbose A boolean of whether to report progress as messages: the pairs scored about every 5 s with
#' the time left, and the top pair of each weight
#' @param cache A boolean of whether to keep the pairs scored so far in `<outdir>/cache`, so that a stopped run continues where it
#' stopped when run again; ignored when outdir is NULL
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
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="6.priority",
                          seed=NULL, verbose=TRUE, cache=TRUE){
  check_weights(weights)
  check_outdir(outdir, null.ok=TRUE)
  check_cache(cache)
  cached <- NULL
  if(cache && !is.null(outdir)){
    cached <- open_cache(file.path(outdir, "cache"), adj, infMat, metadata)
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
  if(parallel){
    if(is.null(cl)){
      cl <- do.call(parallelly::makeClusterPSOCK, cluster_args())
      on.exit(parallel::stopCluster(cl), add=TRUE)
    }
    doParallel::registerDoParallel(cl)
    on.exit(foreach::registerDoSEQ(), add=TRUE)
  } else {
    foreach::registerDoSEQ()
  }
  workers <- foreach::getDoParWorkers()
  if(verbose){
    message("Scoring ", nrow(lrs), " ligand-receptor pairs", if(significance) paste0(" and screening ", nrow(background.lrs), " background pairs"),
            if(workers > 1) paste0(" on ", workers, " workers"))
  }
  if(!is.null(cached) && significance){
    cached$background <- background.lrs
  }

  priority.terms <- calculate_priority(lrs, adj, infMat, metadata, verbose, cached)
  if(significance){
    bp <- background_priority(background.lrs, lrs, priority.terms, weights, adj, infMat, metadata, verbose, cached)
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
    priority <- score_priority(priority.terms, weight)
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
