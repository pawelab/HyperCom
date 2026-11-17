#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression. Missing (NA) counts are ignored.
#' A gene name used by several rows keeps its first row, with a warning. Counts with negative values (scaled data) and counts with
#' genes as column names stop.
#' @param counts A matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column names the interaction, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions, e.g. `cellchat.mm.lrdb`
#' @param outfile A string naming the output file, or NULL to write none
#' @returns A matrix with row names being cells and column names being genes of only ligand and receptor genes with greater than 0 expression,
#' sparse if counts is sparse
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  check_lrdb(lrdb)
  counts <- Matrix::t(counts)
  genes <- unique(unlist(lapply(lrdb[, -1, drop=FALSE], as.character), use.names=FALSE))
  features <- intersect(colnames(counts), genes)
  if(length(features) == 0 && any(rownames(counts) %in% genes)){
    stop("counts has genes as column names: give a gene x cell matrix, e.g. t(counts)")
  }
  dup.genes <- unique(colnames(counts)[duplicated(colnames(counts)) & colnames(counts) %in% features])
  if(length(dup.genes) > 0){
    warning("gene names used by several rows of counts, the first row is kept: ", paste(dup.genes, collapse=", "))
  }

  exp <- counts[, features, drop=FALSE]
  if(any(!is.na(exp)) && min(exp, na.rm=TRUE) < 0){
    stop("counts has negative values (scaled data?): give raw or normalised counts; HyperCom uses whether each gene is above 0 in each cell")
  }
  genes.exp <- Matrix::colMeans(exp, na.rm=TRUE)
  genes.filter <- features[which(genes.exp > 0)]
  cc.counts <- counts[, genes.filter, drop=FALSE]

  if(!is.null(outfile)){
    saveRDS(cc.counts, outfile)
  }
  return(cc.counts)
}

#' Check a table of ligand-receptor interactions
#'
#' `check_lrdb` stops unless lrdb looks like a table of ligand-receptor interactions
#'
#' This function takes the lrdb argument and stops, naming the expected format, unless it is a dataframe with a column named ligand
#' and at least one receptor column, so that swapped arguments fail at once.
#' @param lrdb A dataframe of ligand-receptor interactions
#' @returns NULL, invisibly
#' @noRd
check_lrdb <- function(lrdb){
  if(!is.data.frame(lrdb) || !("ligand" %in% colnames(lrdb)) || ncol(lrdb) < 3){
    stop("lrdb must be a dataframe of ligand-receptor interactions (interaction, ligand and receptor columns) like cellchat.mm.lrdb")
  }

  return(invisible(NULL))
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns dataframe of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptors interactions to produce
#' a table of the ligand-receptor interactions whose genes are all in the counts matrix.
#' @param counts A filtered matrix with row names being cells and column names being genes, from filter_counts
#' @inheritParams filter_counts
#' @returns A dataframe of ligand-receptor interactions in the format of lrdb with only interactions present in the counts matrix, and
#' without receptor columns that are empty for all of them
#' @export
generate_lrs_table <- function(counts, lrdb, outfile="1.lrs.csv"){
  genes <- as.matrix(lrdb[, -1, drop=FALSE])
  found <- matrix(genes %in% colnames(counts), nrow=nrow(genes), ncol=ncol(genes)) | is.na(genes)
  keep <- rowSums(!found) == 0
  lrdb <- lrdb[keep, , drop=FALSE]
  if(nrow(lrdb) > 0){
    lrdb <- Filter(function(x)!all(is.na(x)), lrdb)
  }

  if(!is.null(outfile)){
    utils::write.csv(lrdb, outfile, row.names = FALSE)
  }
  return(lrdb)
}

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns dataframe of hyperedges
#'
#' This function takes a filtered cell x gene counts matrix and table of ligand-receptors interactions to produce
#' a table of hyperedges: one per ligand-receptor interaction (its genes) and one per gene (the gene and the cells that express it).
#' Cells are nodes named by the row names of counts, so they must be unique (prefix barcodes with their sample when merging samples).
#' @inheritParams generate_lrs_table
#' @param lrs A dataframe of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @returns A dataframe of hyperedges where the edge column are the hyperedges and the node column are the node belonging to the hyperedge
#' @export
generate_hyperedges <- function(counts, lrs, outfile="2.hyperedges.csv"){
  dup.cells <- unique(rownames(counts)[duplicated(rownames(counts))])
  if(length(dup.cells) > 0){
    stop("cell names used by several cells of counts (prefix them with their sample): ", paste(utils::head(dup.cells, 10), collapse=", "))
  }

  lrs <- data.frame(lapply(lrs, as.character), check.names=FALSE)
  lr.edges <- unique(tidyr::pivot_longer(lrs, -1, values_to = "node", names_to=NULL, values_drop_na = TRUE))
  lr.edges$weight <- 1
  colnames(lr.edges) <- c("edge", "node", "weight")

  gene.edges <- data.frame(edge=colnames(counts), node=colnames(counts), weight=rep(1, ncol(counts)), check.names=FALSE)

  counts[counts == 0] <- NA
  counts <- data.frame(as.matrix(counts), check.names=FALSE)
  counts$node <- row.names(counts)
  exp.edges <- tidyr::pivot_longer(counts, !node, names_to="edge", values_to="weight", values_drop_na=TRUE)
  exp.edges$weight <- 1
  exp.edges <- exp.edges[,c(2, 1, 3)]

  edges <- rbind(lr.edges, gene.edges, exp.edges)

  if(!is.null(outfile)){
    utils::write.csv(edges, outfile, row.names=FALSE)
  }
  return(edges)
}

#' Generate adjacency matrix of hypergraph
#'
#' `generate_adjacency_matrix` returns the adjacency matrix of the hypergraph
#'
#' This function takes a dataframe of hyperedges and returns its adjacency matrix. Nodes are ordered by their names in the C locale,
#' so the order does not depend on the session's collation. Nodes that share no hyperedge with another node are dropped.
#' @param hyperedges A dataframe of hyperedges where the edge column are the hyperedges and the node column are the node belonging to the hyperedge
#' @inheritParams filter_counts
#' @returns A sparse adjacency matrix where each entry is the number of shared hyperedges between nodes
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds"){
  node <- as.character(hyperedges$node)
  edge <- as.character(hyperedges$edge)
  he <- data.frame(node=factor(node, levels=sort(unique(node), method="radix")), edge=factor(edge, levels=sort(unique(edge), method="radix")))
  incidence.matrix <- stats::xtabs(~node + edge, he, sparse=TRUE)
  adj <- Matrix::tcrossprod(incidence.matrix)
  diag(adj) <- 0

  rm.indices <- which(Matrix::rowSums(adj) == 0)
  if(length(rm.indices) > 0) {
    adj <- adj[-rm.indices, -rm.indices]
  }

  if(!is.null(outfile)){
    saveRDS(adj, outfile)
  }
  return(adj)
}

#' Calculate influence matrix of hypergraph
#'
#' `generate_infMat` returns the influence matrix of the hypergraph
#'
#' This function takes an adjacency matrix and returns its influence matrix using heat diffusion.
#' @param adj A named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @inheritParams filter_counts
#' @returns A dense influence matrix with row and column names being the nodes of adj
#' @export
generate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds"){
  len <- dim(adj)[1]

  D <- as.matrix(Matrix::Diagonal(len))
  W <- adj/Matrix::rowSums(adj)
  W <- (1 - r) * W
  W <- as.matrix(W)
  M <- D - W

  infMat <- r * Matrix::solve(M)
  rownames(infMat) <- rownames(adj)
  colnames(infMat) <- rownames(adj)

  if(!is.null(outfile)){
    saveRDS(infMat, outfile)
  }
  return(infMat)
}

#' Generate cleaned metadata table
#'
#' `generate_metadata` returns a minimal metadata dataframe
#'
#' This function takes an influence matrix and metadata table and creates a metadata data table containing the
#' node name, sample, and group. If group is unspecified, it is assumed all nodes are part of the same group.
#' Nodes of the influence matrix that are not cells of the metadata are genes: their sample, group and kept columns are "gene", so a
#' cell missing from the metadata is treated as a gene (`HyperCom` leaves such cells out, with a warning). Kept columns named cell, sample or group are renamed with make.unique (e.g. sample.1); kept dates and other classed columns
#' other than factors become strings. Cells are compared only within their sample, so every sample needs at least 2 cells in the
#' graph (cells that express a ligand or receptor gene).
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing cell metadata with row names being cells
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @inheritParams filter_counts
#' @returns A dataframe with columns cell, sample, group and the kept columns, one row per node of the influence matrix
#' @export
generate_metadata <- function(infMat, metadata, sample, group=NULL, keep=NULL, outfile="5.metadata.csv"){
  check_columns(metadata, sample, group, keep)
  index <- data.frame(cell = row.names(infMat), check.names=FALSE)
  cells <- rownames(metadata)

  if(is.null(group)){
    group.vals <- "none"
  } else {
    group.vals <- as.character(metadata[[group]])
  }
  sample.vals <- as.character(metadata[[sample]])

  in.graph <- cells %in% index$cell
  if(!any(in.graph)){
    stop("no node of infMat is a row name of metadata: the row names of metadata must be the cell names of the counts")
  }
  missing.sample <- sum(is.na(sample.vals[in.graph]))
  missing.group <- sum(is.na(rep(group.vals, length.out=nrow(metadata))[in.graph]))
  if(missing.sample > 0 || missing.group > 0){
    stop("sample and group must not be missing for cells in the graph: ", missing.sample, " cells miss a sample and ", missing.group,
         " a group")
  }
  if(any(sample.vals == "gene", na.rm=TRUE) || any(group.vals == "gene", na.rm=TRUE)){
    stop("\"gene\" labels gene nodes and cannot name a sample or group")
  }
  sizes <- table(sample.vals[in.graph])
  if(any(sizes < 2)){
    stop("every sample needs at least 2 cells in the graph (cells that express a ligand or receptor gene), as cells are compared only ",
         "within their sample; remove the cells of: ", paste(names(sizes)[sizes < 2], collapse=", "))
  }

  clean.df <- data.frame(cell=cells, sample=sample.vals, group=group.vals, check.names=FALSE)
  clean.df$sample <- factor(clean.df$sample, levels=c(unique(clean.df$sample), "gene"))
  clean.df$group <- factor(clean.df$group, levels=c(unique(clean.df$group), "gene"))

  if(!is.null(keep)){
    additional <- metadata[, keep, drop=FALSE]
    colnames(additional) <- make.unique(c(colnames(clean.df), colnames(additional)))[-seq_len(ncol(clean.df))]
    for(col in colnames(additional)){
      if(is.factor(additional[[col]])){
        additional[[col]] <- factor(additional[[col]], levels=unique(c(levels(additional[[col]]), "gene")))
      } else if(is.object(additional[[col]])){
        additional[[col]] <- as.character(additional[[col]])
      }
    }
    clean.df <- cbind(clean.df, additional)
  }

  clean.df <- dplyr::left_join(index, clean.df, by=dplyr::join_by(cell))
  genes <- !(clean.df$cell %in% cells)
  for(col in colnames(clean.df)[-1]){
    clean.df[genes, col] <- "gene"
  }

  if(!is.null(outfile)){
    utils::write.csv(clean.df, outfile, row.names = FALSE)
  }

  return(clean.df)
}

#' Check the metadata columns named by the arguments
#'
#' `check_columns` stops unless sample, group and keep name columns of the metadata
#'
#' This function takes the metadata and the column arguments of generate_metadata and stops, naming the problem, unless metadata is a
#' dataframe, sample and group are each one string and every named column exists.
#' @inheritParams generate_metadata
#' @returns NULL, invisibly
#' @noRd
check_columns <- function(metadata, sample, group=NULL, keep=NULL){
  if(!is.data.frame(metadata)){
    stop("metadata must be a dataframe with row names being cells")
  }
  for(arg in list(sample, group)){
    if(!is.null(arg) && !(is.character(arg) && length(arg) == 1)){
      stop("sample and group must each be one column name of metadata, e.g. sample=\"sample\"")
    }
  }
  absent <- setdiff(c(sample, group, keep[is.character(keep)]), colnames(metadata))
  if(length(absent) > 0){
    stop("not columns of metadata: ", paste(absent, collapse=", "), " (its columns: ", paste(utils::head(colnames(metadata), 20), collapse=", "),
         ")")
  }

  return(invisible(NULL))
}

#' Perform element wise multiplication of two matrices
#'
#' `fast_mult` returns the element wise product of two matrices
#'
#' This function takes two matrices and multiplies them together element-wise.
#' @param A,B Matrices
#' @returns A matrix
#' @noRd
fast_mult <- function(A, B){
  setup <- {AA <- A*1}
  return(collapse::setop(AA, "*", B))
}

#' Perform element wise division of two matrices
#'
#' `fast_div` returns the element wise quotient of two matrices
#'
#' This function takes two matrices and divides them together element-wise.
#' @param A,B Matrices
#' @returns A matrix
#' @noRd
fast_div <- function(A, B){
  setup <- {AA <- A/1}
  return(collapse::setop(AA, "/", B))
}

#' Check that genes are nodes of the influence matrix
#'
#' `check_genes` stops if a gene is not a node of the influence matrix
#'
#' This function takes gene names and the node names of an influence matrix and stops naming the genes that are not among them.
#' @param genes A vector of strings of gene names, NA allowed
#' @param nodes A vector of strings of the node names of the influence matrix
#' @returns NULL, invisibly
#' @noRd
check_genes <- function(genes, nodes){
  missing <- setdiff(stats::na.omit(as.character(genes)), nodes)
  if(length(missing) > 0){
    stop("not nodes of infMat: ", paste(missing, collapse=", "))
  }

  return(invisible(NULL))
}

#' Calculate vectors used in HyperCom scoring
#'
#' `score_sr` returns a dataframe of sender and receiver scores use by HyperCom scoring
#'
#' This function takes a named matrix and a 2 column dataframe with the sender and receiver score. If ligand and receptor args are NULL calculates network score.
#' Cells and genes are found in the influence matrix by name and scored in its node order, so the metadata rows may be in any order.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligand A string of the ligand name
#' @param receptor A string of the receptor name
#' @returns A dataframe with the sender and receiver score, one row per cell in the order of the metadata
#' @noRd
score_sr <- function(infMat, metadata, ligand=NULL, receptor=NULL){
  cells.idx <- metadata$group != "gene"
  nodes <- rownames(infMat)
  cell.pos <- match(metadata$cell[cells.idx], nodes)
  if(anyNA(cell.pos) || anyDuplicated(cell.pos) > 0){
    stop("metadata cells must be distinct nodes of infMat: make it with generate_metadata")
  }
  o <- order(cell.pos)
  cell.pos <- cell.pos[o]
  mat.names <- metadata$cell[cells.idx][o]
  samples <- as.numeric(factor(metadata$sample[cells.idx][o]))
  length <- length(mat.names)

  infMat.cells <- infMat[cell.pos, cell.pos]

  if(!is.null(ligand) & !is.null(receptor)){
    idx.1 <- match(ligand, nodes)
    idx.2 <- match(receptor, nodes)

    m <- matrix(infMat[,idx.1])
    n <- matrix(infMat[idx.2,])
    trans <- Rfast::Tcrossprod(m, n)

    trans <- trans[cell.pos, cell.pos]
    trans <- fast_div(trans, infMat.cells)
  } else {
    trans <- infMat.cells
  }

  same.sample <- matrix(Rfast::Outer(samples, samples,"/"), nrow=length, ncol=length) == 1
  trans <- fast_mult(trans, same.sample)
  diag(trans) <- 0
  trans[trans == 0] <- NA

  trans <- dplyr::percent_rank(as.vector(trans))
  trans <- matrix(trans, nrow=length, ncol=length, byrow=TRUE)

  sender <- colMeans(trans, na.rm = TRUE)
  receiver <- rowMeans(trans, na.rm = TRUE)

  df <- data.frame(cell=mat.names, sender=sender, receiver=receiver)
  if(is.unsorted(o)){
    df <- df[order(o), ]
    rownames(df) <- NULL
  }

  return(df)
}

#' Calculate the HyperCom score for each cell
#'
#' `score_hypercom` returns a dataframe containing the HyperCom score for each cell
#'
#' This function takes a named influence matrix, the corresponding metadata, the ligand, and receptor
#' to give each cell a HyperCom score. If both the ligand and receptor are specified this is the LR Pair Score;
#' if neither is, the score is the Network Score. If a receptor complex is specified, the average of the score vectors is used.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligand A string of the ligand name, NULL for the network score
#' @param receptors A vector containing string elements of the receptor names, NULL for the network score
#' @returns A dataframe with one row per cell, sorted by hypercom (highest first), with columns cell, sender, receiver, dif
#' (receiver - sender), avg (their mean), max (the larger of the two) and hypercom (the HyperCom score)
#' @export
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL){
  if(is.null(ligand) != is.null(receptors)){
    stop("give both ligand and receptors for a pair score, or neither for the network score")
  }
  if(!is.null(ligand) && length(ligand) != 1){
    stop("ligand must be one gene")
  }
  if(!is.null(receptors)){
    receptors <- receptors[!is.na(receptors)]
  }
  if(!is.null(receptors) && length(receptors) == 0){
    stop("receptors must name at least one gene")
  }
  check_genes(c(ligand, receptors), rownames(infMat))

  complex.size <- length(receptors)
  if(complex.size > 1){
    ss <- c()
    rs <- c()
    for(receptor in receptors){
      sr <- score_sr(infMat, metadata, ligand, receptor)
      if(length(ss) == 0){
        ss <- sr$sender
        rs <- sr$receiver
      } else {
        ss <- ss + sr$sender
        rs <- rs + sr$receiver
      }
    }
    sr$sender <- ss / complex.size
    sr$receiver <- rs / complex.size
  } else {
    sr <- score_sr(infMat, metadata, ligand, receptors)
  }
  if(!is.null(receptors) && all(ligand == receptors)){
    sr$dif <- 0
  } else {
    sr$dif <- sr$receiver - sr$sender
  }
  sr$avg <- (sr$receiver + sr$sender) / 2
  sr$max <- pmax(sr$sender, sr$receiver)
  sr$hypercom <- dplyr::percent_rank(abs(sr$dif))
  sr$hypercom <- sign(sr$dif) * sr$hypercom * sr$avg
  sr <- dplyr::arrange(sr, dplyr::desc(hypercom))

  return(sr)
}

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
  nodes <- dplyr::filter(metadata, sample=="gene")$cell
  for (i in 2:dim(lrs)[2]) {
    col <- sample(nodes, count, replace = TRUE)
    background <- cbind(background, col)
  }

  colnames(background) <- colnames(lrs)

  return(background)
}

#' Perform a permutation test
#'
#' `permutation_test` returns a numeric p-value
#'
#' This function takes a priority and the priorities of the background pairs and returns the fraction of background priorities at
#' least as high in absolute value. Background pairs without a priority (NA or NaN) are left out; NA if no background pair has one.
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

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to score ligand-receptor pairs. Genes are read by name, also from factor columns. The group purity is at least 0: cells perfectly
#' balanced over the groups round it to about -1e-16, whose log would make the priority NaN.
#' @inheritParams prioritize_lr
#' @param weight A weight of at least 0 for the fraction of cells expressing a pair's genes
#' @returns A numeric vector of priority scores
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, weight=1){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells <- dplyr::filter(metadata, sample != "gene")$cell
  ligands <- as.character(lrs$ligand)
  receptors <- as.matrix(lrs[, 3:dim(lrs)[2], drop=FALSE])
  check_genes(c(ligands, receptors), rownames(infMat))

  priority <- foreach::`%dopar%`(foreach::foreach(i=seq_along(ligands), .packages=c("HyperCom", "Matrix"), .combine=c), {
    l <- ligands[i]
    rs <- stats::na.omit(receptors[i, ])
    cc.score <- score_hypercom(infMat, metadata, l, rs)
    cc.score <- dplyr::left_join(cc.score, metadata, by=dplyr::join_by(cell))

    cc.score <- dplyr::mutate(cc.score, max=pmax(sender, receiver))
    cc.score$type <- cc.score$dif > 0
    rownames(cc.score) <- cc.score$cell
    cc.score <- cc.score[cells,]

    if(all(l == rs)){
      pair.rows <- data.frame(as.matrix(adj[l, cells]), check.names=FALSE)
      colnames(pair.rows) <- l
    }else{
      pair.rows <- data.frame(as.matrix(Matrix::t(adj[c(l, rs), cells])), check.names=FALSE)
    }

    cc.score$present <- as.logical(do.call(pmax, pair.rows))
    pct.min <- min(colSums(pair.rows) / length(cells))

    cc.expressed <- dplyr::filter(cc.score, present > 0)

    frac <- length(cc.expressed$cell) / length(cells)
    max.score <- mean(cc.expressed$max)

    if(gini.max > 0){
      purity <- max(0, gini.max - mltools::gini_impurity(cc.expressed$group))
    }else{purity <- 1}

    balance <- mltools::gini_impurity(cc.expressed$type) / 0.5

    score <- exp(stats::weighted.mean(log(c(frac, pct.min, max.score, balance, purity)), c(weight, rep(1, 4))))
    score
  })

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

#' Check the output directory argument
#'
#' `check_outdir` stops unless outdir is one string, or NULL where it means writing no files
#'
#' This function takes the outdir argument of prioritize_lr or HyperCom and stops unless it is one string (or NULL when null.ok).
#' @param outdir The outdir argument
#' @param null.ok A boolean of whether NULL (write no files) is allowed
#' @returns NULL, invisibly
#' @noRd
check_outdir <- function(outdir, null.ok=FALSE){
  if(is.null(outdir) && null.ok){
    return(invisible(NULL))
  }
  if(!(is.character(outdir) && length(outdir) == 1 && !is.na(outdir))){
    stop("outdir must be one string naming the output directory", if(null.ok) " (or NULL to write no files)")
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

#' Choose the arguments of the default cluster
#'
#' `cluster_args` returns the arguments of parallelly::makeClusterPSOCK for the cluster prioritize_lr makes when none is given
#'
#' This function returns, inside an SGE parallel job (PE_HOSTFILE set), the job's workers started with qrsh, which can span several
#' nodes; elsewhere, one local worker per available core.
#' @returns A list of arguments of parallelly::makeClusterPSOCK
#' @noRd
cluster_args <- function(){
  if(nzchar(Sys.getenv("PE_HOSTFILE"))){
    return(list(parallelly::availableWorkers(), rshcmd="qrsh", rshopts=c("-inherit", "-nostdin", "-V"), autoStop=TRUE, outfile=""))
  }

  return(list(parallelly::availableCores(), autoStop=TRUE, outfile=""))
}

#' Restore the state of the random number generator
#'
#' `restore_seed` puts back the random number state saved before a set.seed
#'
#' This function takes the value of .Random.seed saved before set.seed (NULL if there was none) and restores it, so that a seed
#' argument does not change the random numbers of the calling session.
#' @param old.seed An integer vector saved from .Random.seed, or NULL
#' @returns NULL, invisibly
#' @noRd
restore_seed <- function(old.seed){
  if(is.null(old.seed)){
    if(exists(".Random.seed", envir=globalenv(), inherits=FALSE)){
      rm(".Random.seed", envir=globalenv())
    }
  } else {
    assign(".Random.seed", old.seed, envir=globalenv())
  }

  return(invisible(NULL))
}

#' Order ligand receptor pairs based on signalling strength and specificity
#'
#' `prioritize_lr` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to prioritize ligand-receptor pairs.
#' With significance, the p-value of a pair is the fraction of 1000 background pairs of random genes whose priority is at least as
#' high, so it moves in steps of 0.001 and is 0 when no background pair scored as high; p.adj is its Benjamini-Hochberg adjustment.
#' Background pairs without a priority (NaN) are left out, with a warning.
#' @param lrs A dataframe of ligand-receptor interactions from generate_lrs_table, in the format of lrdb with only interactions present
#' in the counts matrix
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param weights A vector of distinct weights of at least 0, one result each, for the importance of number of cells involved in an
#' interaction: 1 weighs it like the other four terms, less than 1 favors pairs of fewer cells, 0 leaves it out
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel
#' @param cl A cluster object, which is left running. Default if parallel: inside an SGE parallel job (PE_HOSTFILE set)
#' `parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd = "qrsh", rshopts = c("-inherit", "-nostdin", "-V"), outfile = "")`,
#' elsewhere `parallelly::makeClusterPSOCK(parallelly::availableCores(), outfile = "")`; it is stopped at the end
#' @param outdir A string naming the output directory, or NULL to write no files
#' @param seed An integer seed for the background interactions of the p-values, NULL to use the current random state; the session's
#' random state is restored afterwards
#' @param verbose A boolean of whether to report progress as messages
#' @returns A list of dataframes named by weight, each the ligand-receptor pairs ordered by their priority column, with columns p and
#' p.adj when significance is TRUE
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="6.priority",
                          seed=NULL, verbose=TRUE){
  check_weights(weights)
  check_outdir(outdir, null.ok=TRUE)
  check_genes(c(as.character(lrs$ligand), as.matrix(lrs[, -(1:2), drop=FALSE])), rownames(infMat))
  if(significance){
    if(!is.null(seed)){
      old.seed <- get0(".Random.seed", envir=globalenv(), inherits=FALSE)
      on.exit(restore_seed(old.seed), add=TRUE)
      set.seed(seed)
    }
    background.lrs <- generate_background_lrs(lrs, metadata)
  }

  if(parallel){
    if(is.null(cl)){
      cl <- do.call(parallelly::makeClusterPSOCK, cluster_args())
      on.exit(parallel::stopCluster(cl), add=TRUE)
    }
    doParallel::registerDoParallel(cl)
    on.exit(foreach::registerDoSEQ(), add=TRUE)
  } else{
    foreach::registerDoSEQ()
  }

  if(!is.null(outdir)){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }
  lrs.list <- list()
  warned <- FALSE

  for(i in seq_along(weights)){
    weight <- weights[i]
    if(verbose){
      message("Weight ", weight, ": scoring ", nrow(lrs), " ligand-receptor pairs", if(significance) " and 1000 background pairs")
    }
    priority <- calculate_priority(lrs, adj, infMat, metadata, weight)

    current.lrs <- lrs

    if(significance){
      current.background <- background.lrs
      background.priority <- calculate_priority(current.background, adj, infMat, metadata, weight)
      no.priority <- sum(is.na(background.priority))
      if(no.priority > 0 && !warned){
        warning(no.priority, " of ", length(background.priority), " background pairs have no priority (NaN) and are left out of the p-values")
        warned <- TRUE
      }
    }

    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      current.lrs$p <- vapply(current.lrs$priority, permutation_test, numeric(1), background=background.priority)
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method = "BH")
    }

    if(!is.null(outdir)){
      outfile <- file.path(outdir, priority_filename(weight))
      utils::write.csv(current.lrs, outfile, row.names = FALSE)
    }
    if(verbose){
      message("Weight ", weight, ": top pair ", current.lrs[[1]][1])
    }

    lrs.list[[i]] <- current.lrs
  }

  names(lrs.list) <- weights

  return(lrs.list)
}

#' Run the HyperCom pipeline
#'
#' `HyperCom` Runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency_matrix`, `generate_infMat`,
#' `generate_metadata`, and `prioritize_lr`, writing every step to outdir (the working directory is not changed). The arguments are
#' checked before any step runs, and the run stops when no interaction of lrdb has all its genes expressed (e.g. lrdb of another
#' species). Cells of counts that are not row names of metadata are left out, with a warning. Priority files of other weights already
#' in outdir are kept, with a warning.
#' @inheritParams filter_counts
#' @inheritParams generate_metadata
#' @inheritParams prioritize_lr
#' @param outdir A string naming the output directory of the pipeline run
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency matrix) and lrs (the interactions with all genes expressed)
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL,
                     outdir="HyperCom", seed=NULL, verbose=TRUE){
  check_lrdb(lrdb)
  check_columns(metadata, sample, group, keep)
  check_weights(weights)
  check_outdir(outdir)
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
  hyperedges <- generate_hyperedges(counts, lrs, outfile=file.path(outdir, "2.hyperedges.csv"))
  adj <- generate_adjacency_matrix(hyperedges, outfile=file.path(outdir, "3.adj.rds"))
  infMat <- generate_infMat(adj, outfile=file.path(outdir, "4.infMat.rds"))
  metadata <- generate_metadata(infMat, metadata, sample, group, keep, outfile=file.path(outdir, "5.metadata.csv"))
  priority <- prioritize_lr(lrs, adj, infMat, metadata, weights=weights, significance=significance, parallel=parallel, cl=cl,
                            outdir=file.path(outdir, "6.priority"), seed=seed, verbose=verbose)

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
#' `HyperCom_seurat` Runs the HyperCom pipeline on a seurat object
#'
#' This function runs `HyperCom` on one assay layer and the cell metadata of a Seurat object (v5 or older assays).
#' @param seurat A seurat object
#' @param assay A string naming the assay to read counts from
#' @param layer A string naming the layer of the assay to read counts from
#' @inheritParams HyperCom
#' @returns Invisibly, a list of the results also written to outdir: priority (from prioritize_lr, one dataframe per weight), metadata
#' (from generate_metadata), infMat (the influence matrix), adj (the adjacency matrix) and lrs (the interactions with all genes expressed)
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, weights=c(1), significance=FALSE,
                            parallel=FALSE, cl=NULL, outdir="HyperCom", seed=NULL, verbose=TRUE){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat@meta.data

  res <- HyperCom(lrdb, counts, metadata, sample=sample, group=group, keep=keep, weights=weights, significance=significance,
                  parallel=parallel, cl=cl, outdir=outdir, seed=seed, verbose=verbose)

  return(invisible(res))
}
