#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression.
#' @param counts A matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A dense matrix with row names being genes and column names being cells of only ligand and receptor genes with greater than 0 expression
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  counts <- Matrix::t(counts)
  genes <- unique(c(unlist(lrdb[,2:length(lrdb)])))
  features <- intersect(colnames(counts), genes)
  exp <- counts[, features]
  genes.exp <- Matrix::colMeans(exp)
  genes.filter <- names(genes.exp[genes.exp > 0])
  cc.counts <- counts[, genes.filter]
  if(!is.null(outfile)){
    saveRDS(cc.counts, outfile)
  }
  return(cc.counts)
}

#' Load Seurat object and filter for HyperCom analysis
#'
#' `load_seurat` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a Seurat object and filters it to
#' just ligand and receptor genes with greater than 0 expression.
#' @param seurat A seurat v5 object with a counts matrix
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' @param assay assay to get data from
#' @param layer layer to get data from
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A dense matrix with row names being genes and column names being cells of only ligand and receptor genes with greater than 0 expression
#' @export
load_seurat <- function(seurat, assay="RNA", layer="counts", lrdb, outfile="0.counts.rds"){
  counts <- seurat@assays[[assay]]@layers[[layer]]
  rownames(counts) <- rownames(seurat)
  colnames(counts) <- colnames(seurat)
  cc.counts <- filter_counts(counts, lrdb, outfile)

  return(cc.counts)
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns dataframe of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered gene x cell counts matrix and table of ligand-receptors interactions to produce
#' a table of ligand-receptors interactions in the counts matrix.
#' @param counts A filtered matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' with only interactions present in the counts matrix
#' @export
generate_lrs_table <- function(counts, lrdb, outfile="1.lrs.csv"){
  length <- dim(lrdb)[2]
  keep <- c()
  for (i in 1:dim(lrdb)[1]) {
    lrs <- stats::na.omit(unlist(lrdb[i, 2:length]))
    keep <- c(keep, all(lrs %in% colnames(counts)))
  }
  lrdb <- lrdb[keep,]
  lrdb <- Filter(function(x)!all(is.na(x)), lrdb)

  if(!is.null(outfile)){
    utils::write.csv(lrdb, outfile, row.names = FALSE)
  }
  return(lrdb)
}

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns dataframe of hyperedges
#'
#' This function takes a filtered gene x cell counts matrix and table of ligand-receptors interactions to produce
#' a table of ligand-receptors interactions in the counts matrix.
#' @param counts A filtered matrix with row names being genes and column names being cells
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param outfile A string naming the output file
#' @returns A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @export
generate_hyperedges <- function(counts, lrs, outfile="2.hyperedges.csv"){
  lr.edges <- unique(tidyr::pivot_longer(lrs, !interaction, values_to = "node", names_to=NULL, values_drop_na = TRUE))
  lr.edges$weight <- 1
  colnames(lr.edges) <- c("edge", "node", "weight")

  gene.edges <- data.frame(edge=colnames(counts), node=colnames(counts), weight=1, check.names=FALSE)

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
#' This function takes a dataframe of hyperedges and returns its adjacency matrix.
#' @param hyperedges A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @param outfile A string naming the output file
#' @returns An adjacency matrix where each entry is the number of shared hyperedges between nodes
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds"){
  incidence.matrix <- stats::xtabs(~node + edge, hyperedges, sparse=TRUE)
  adj <- Matrix::tcrossprod(incidence.matrix)
  diag(adj) <- 0

  rm.indices <- which(Matrix::rowSums(adj) == 0)
  if(length(rm.indices > 0)) {
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
#' @param outfile A string naming the output file
#' @returns An influence matrix
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
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing cell metadata
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param outfile A string naming the output file
#' @returns An influence matrix
#' @export
generate_metadata <- function(infMat, metadata, sample, group=NULL, keep=NULL, outfile="5.metadata.csv"){
  index <- data.frame(cell = row.names(infMat), check.names=FALSE)
  metadata$cell <- rownames(metadata)

  clean.df <- data.frame(cell=metadata$cell, sample=as.character(metadata[, sample]), group=as.character(metadata[, group]), check.names=FALSE)
  clean.df$sample <- factor(clean.df$sample, levels=c(unique(clean.df$sample), "gene"))
  if(is.null(group)){
    clean.df$group <- "none"
  }
  else{
    clean.df$group <- factor(clean.df$group, levels=c(unique(clean.df$group), "gene"))
  }

  if(!is.null(keep)){
    additional <- metadata[, keep]
    clean.df <- cbind(clean.df, additional)
  }

  clean.df <- dplyr::left_join(index, clean.df, by=dplyr::join_by(cell))
  clean.df[is.na(clean.df)] <- "gene"

  if(!is.null(outfile)){
    utils::write.csv(clean.df, outfile, row.names = FALSE)
  }
  return(clean.df)
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

#' Calculate vectors used in HyperCom scoring
#'
#' `score_sr` returns a dataframe of sender and receiver scores use by HyperCom scoring
#'
#' This function takes a named matrix and a 2 column dataframe with the sender and receiver score. If ligand and receptor args are NULL calculates network score.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named sample and group
#' @param ligand A string of the ligand name
#' @param receptor A string of the ligand name
#' @returns A dataframe with the sender and receiver score
#' @noRd
score_sr <- function(infMat, metadata, ligand=NULL, receptor=NULL){
  cells.idx <- metadata$group != "gene"
  mat.names <- metadata$cell[cells.idx]
  samples <- as.numeric(factor(metadata$sample[cells.idx]))
  length <- length(mat.names)

  infMat.cells <- infMat[cells.idx, cells.idx]

  if(!is.null(ligand) & !is.null(receptor)){
    idx.1 <- which(metadata$cell == ligand)
    idx.2 <- which(metadata$cell == receptor)

    m <- matrix(infMat[,idx.1])
    n <- matrix(infMat[idx.2,])
    trans <- Rfast::Tcrossprod(m, n)

    trans <- trans[cells.idx, cells.idx]
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

  return(df)
}

#' Calculate the HyperCom score for each cell
#'
#' `score_hypercom` returns a dataframe containing the HyperCom score for each cell
#'
#' This function takes a named influence matrix, the corresponding metadata, the ligand, and receptor
#' to give each cell a HyperCom score. If both the ligand and receptor are specified this is the LR Pair Score
#' otherwise the score is the Network Score. If a receptor complex is specified, the average of the score vectors is used.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligand A string of the ligand name
#' @param receptors A vector containing string elements of the receptor names
#' @returns A dataframe containing the sender score, receiver score, avg, dif, and HyperCom score
#' @export
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL){
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
  if(all(ligand == receptors) & !is.null(receptors)){
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
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param count An integer of the number of interactions to generate
#' @return A dataframe of random gene interactions
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
#' This function takes a numeric and list of numeric vector and returns how many values in the vector are less than the numeric.
#' @param x A numeric
#' @param background A numeric vector
#' @return A numeric p-value
#' @noRd
permutation_test <- function(x, background){
  sum(abs(background) >= abs(x)) / length(background)
}

#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to score ligand-receptor pairs.
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param weight a value between 0 to 1 for the importance of number of cells involved in an interaction
#' @return A numeric vector of priority scores
#' @noRd
calculate_priority <- function(lrs, adj, infMat, metadata, weight=1){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells <- dplyr::filter(metadata, sample != "gene")$cell

  priority <- foreach::`%dopar%`(foreach::foreach(i=1:length(lrs$ligand), .packages=c("HyperCom", "Matrix"), .combine=c), {
    l <- lrs$ligand[i]
    rs <- stats::na.omit(unlist(lrs[i, 3:dim(lrs)[2]]))
    title <- paste(c(l, rs), collapse = "_")
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
    max <- mean(cc.expressed$max)

    if(gini.max > 0){
      purity <- gini.max - mltools::gini_impurity(cc.expressed$group)
    }else{purity <- 1}

    balance <- mltools::gini_impurity(cc.expressed$type) / 0.5

    score <- exp(weighted.mean(log(c(frac, pct.min, max, balance, purity)), c(weight, rep(1, 4))))
    print(paste0(title, ",", score))
    score
  })

  return(priority)
}

#' Order ligand receptor pairs based on signalling strength and specificity
#'
#' `prioritize_lr` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the adjacency matrix, a named influence matrix, and the corresponding metadata
#' to prioritize ligand-receptor pairs.
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix.
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param weight a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel
#' @param cl A cluster object. Default if parallel: `parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd = "qrsh", rshopts = c("-inherit", "-nostdin", "-V"), outfile = "")`
#' @param outdir A string naming the output directory
#' @returns A dataframe of ligand-receptor pairs ordered by priority score for each weight
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="6.priority"){
  if(significance){
    background.lrs <- HyperCom:::generate_background_lrs(lrs, metadata)
  }

  if(parallel){
    if(is.null(cl)){
      cl <- parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd="qrsh", rshopts=c("-inherit", "-nostdin", "-V"), autoStop=TRUE, outfile="")
    }
    doParallel::registerDoParallel(cl)
  } else{
    foreach::registerDoSEQ()
  }

  lrs.list <- list()

  for(i in seq(1, length(weights))){
    weight <- weights[i]
    print(weight)
    priority <- calculate_priority(lrs, adj, infMat, metadata, weight)

    current.lrs <- lrs
    current.background <- background.lrs

    if(significance){
      background.priority <- calculate_priority(current.background, adj, infMat, metadata, weight)
    }

    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      current.lrs$p <- as.numeric(lapply(current.lrs$priority, FUN=HyperCom:::permutation_test, background=background.priority))
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method = "BH")
    }

    if(!is.null(outdir)){
      dir.create(file.path(outdir), showWarnings = FALSE)
      outfile <- paste0(outdir, "/priority.", weight, ".csv")
      utils::write.csv(current.lrs, outfile, row.names = FALSE)
    }

    lrs.list[[i]] <- current.lrs
  }

  if(parallel){
    parallel::stopCluster(cl)
    gc()
  }

  names(lrs.list) <- weights

  return(lrs.list)
}

#' Run the HyperCom pipeline
#'
#' `HyperCom` Runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency-matrix`, `generate_infMat`, `generate-metadata`, and `prioritize_lr`.
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions.
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel
#' @param cl A cluster object. Default if parallel: `parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd = "qrsh", rshopts = c("-inherit", "-nostdin", "-V"), outfile = "")`
#' @returns 0 if run successfully
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=1, significance=FALSE, parallel=FALSE, cl=NULL){
  counts <- filter_counts(counts, lrdb)
  lrs <- generate_lrs_table(counts, lrdb)
  hyperedges <- generate_hyperedges(counts, lrs)
  adj <- generate_adjacency_matrix(hyperedges)
  infMat <- generate_infMat(adj)
  metadata <- generate_metadata(infMat, metadata, sample, group, keep)
  priority <- prioritize_lr(lrs, adj, infMat, metadata, weights, significance, parallel, cl)
  return(0)
}
