#' Filter counts matrix to ligand and receptor genes
#'
#' `filter_counts` returns the count matrix with only ligand and receptor genes
#'
#' This function takes a gene x cell counts matrix and filters it to
#' just ligand and receptor genes with greater than 0 expression.
#' @param counts A dense matrix with row names being genes and column names being cells
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param outfile A string naming the output file
#' @returns A dense matrix with row names being genes and column names being cells of only ligand and receptor genes with greater than 0 expression
#' @export
filter_counts <- function(counts, lrdb, outfile="0.counts.rds"){
  counts <- t(counts)
  genes <- unique(c(unlist(lrdb[,2:length(lrdb)])))
  features <- intersect(colnames(counts), genes)
  exp <- counts[, features]
  genes.exp <- collapse::fmean(exp)   
  genes.filter <- names(genes.exp[genes.exp > 0])
  cc.counts <- counts[, genes.filter]
  saveRDS(cc.counts, outfile)
  return(cc.counts)
}

#' Generate a table of ligand-receptor interactions present in a counts matrix
#'
#' `generate_lrs_table` returns dataframe of ligand-receptor interactions present in a counts matrix
#'
#' This function takes a filtered gene x cell counts matrix and table of ligand-receptors interactions to produce
#' a table of ligand-receptors interactions in the counts matrix.
#' @param counts A filtered dense matrix with row names being genes and column names being cells
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
    lrs <- na.omit(unlist(lrdb[i, 2:length]))
    keep <- c(keep, all(lrs %in% colnames(counts)))
  }
  lrdb <- lrdb[keep,]
  lrdb <- Filter(function(x)!all(is.na(x)), lrdb)
  write.csv(lrdb, outfile, row.names = FALSE)
  return(lrdb)
}

#' Generate a table of hyperedges
#'
#' `generate_hyperedges` returns dataframe of hyperedges
#'
#' This function takes a filtered gene x cell counts matrix and table of ligand-receptors interactions to produce
#' a table of ligand-receptors interactions in the counts matrix.
#' @param counts A filtered dense matrix with row names being genes and column names being cells
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param outfile A string naming the output file
#' @returns A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @export
generate_hyperedges <- function(counts, lrs, outfile="2.hyperedges.csv"){
  lr.edges <- unique(tidyr::pivot_longer(lrs, !edges, values_to = "nodes", names_to=NULL, values_drop_na = TRUE))
  lr.edges$weight <- 1
  
  gene.edges <- data.frame(edges=colnames(counts), nodes=colnames(counts), weight=1)
  
  complexes <- unique(dplyr::select(lrs, 3:length(lrs)))
  complexes <- complexes <- dplyr::filter(complexes, !is.na(complexes$r2))
  
  for (i in 1:dim(complexes)[1]) {
    genes <- na.omit(unlist(complexes[i, 1:dim(complexes)[2]]))
    complex.size <- length(genes)
    edge.name <- paste0(genes, collapse="_")
    complex.nodes <- counts[, genes] > 0
    complex.df <- data.frame(edges=edge.name, nodes=c(genes, rownames(counts)), weight=c(rep(complex.size, complex.size), rowSums(complex.nodes)))
    complex.df <- dplyr::filter(complex.df, weight == complex.size)
    complex.df$weight <- 1
    lr.edges <- dplyr::bind_rows(lr.edges, complex.df)
  }
  
  counts[counts == 0] <- NA
  counts <- as.data.frame(counts)
  counts$nodes <- row.names(counts)
  exp.edges <- tidyr::pivot_longer(counts, !nodes, names_to="edges", values_to="weight", values_drop_na=TRUE)
  exp.edges$weight <- 1
  exp.edges <- exp.edges[,c(2, 1, 3)]
  
  edges <- rbind(lr.edges, gene.edges, exp.edges)
  
  write.csv(edges, outfile, row.names=FALSE)
  return(edges)
}

#' Generate adjacency matrix of hypergraph
#'
#' `generate_adjacency_matrix` returns the adjacency matrix of the hypergraph
#'
#' This function takes a dataframe of hyperedges and returns its adjacency matrix.
#' @param hyperedgelist A dataframe of hyperedges where the edges column are the hyperedges and the nodes column are the node belonging to the hyperedge
#' @param outfile A string naming the output file
#' @returns An adjacency matrix where each entry is the number of shared hyperedges between nodes
#' @export
generate_adjacency_matrix <- function(hyperedgelist, outfile="3.adj.rds"){
  incidence.matrix <- xtabs(~nodes + edges, hyperedgelist, sparse=TRUE)
  adj <- Matrix::tcrossprod(incidence.matrix)
  diag(adj) <- 0
  adj <- as.matrix(adj)
  
  saveRDS(adj, outfile)
  return(adj)
}

#' Calculate influence matrix of hypergraph
#'
#' `calculate_infMat` returns the influence matrix of the hypergraph
#'
#' This function takes an adjacency matrix and returns its influence matrix using heat diffusion.
#' @param adj A named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @param outfile A string naming the output file
#' @returns An influence matrix
#' @export
calculate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds"){
  rm.indices <- which(rowSums(adj) == 0)
  if(length(rm.indices > 0)) {
    adj <- adj[-rm.indices, -rm.indices]
  }
  len <- dim(adj)[1]
  
  D <- as.matrix(Matrix::Diagonal(len))
  W <- adj/Matrix::rowSums(adj)
  W <- (1 - r) * W
  W <- as.matrix(W)
  M <- D - W
  
  Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends="RcppArmadillo")
  infMat <- r*armaInv(M)
  rownames(infMat) <- rownames(adj)
  colnames(infMat) <- rownames(adj)
  
  saveRDS(infMat, outfile)
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
#' @param outfile A string naming the output file
#' @returns An influence matrix
#' @export
generate_metadata <- function(infMat, metadata, sample, group=NULL, outfile="5.metadata.csv"){
  index <- data.frame(cell = row.names(infMat))
  metadata$cell <- rownames(metadata)
  
  clean.df <- data.frame(cell=metadata$cell, sample=as.character(metadata[, sample]), group=as.character(metadata[, group]))
  clean.df$sample <- factor(clean.df$sample, levels=c(unique(clean.df$sample), "gene"))
  if(is.null(group)){
    clean.df$group <- "none"
  }
  else{
    clean.df$group <- factor(clean.df$group, levels=c(unique(clean.df$group), "gene"))
  }
  
  clean.df <- dplyr::left_join(index, clean.df, by=dplyr::join_by(cell))
  clean.df[is.na(clean.df)] <- "gene"
  
  write.csv(clean.df, outfile, row.names = FALSE)
  return(clean.df)
}

#' Perform element wise multiplication of two matrices
#'
#' `fast_mult` returns the element wise product of two matrices
#'
#' This function takes two matrices and multiplies them together element-wise.
#' @param A,B Matrices
#' @returns A matrix
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
fast_div <- function(A, B){
  setup <- {AA <- A/1}
  return(collapse::setop(AA, "/", B))
}

#' Calculate vectors used in HyperCom scoring
#'
#' `process_transition` returns a list of four vectors used by HyperCom scoring
#'
#' This function takes a named matrix and returns four vectors, a forward vector, a backward vector,
#' an element wise average of those vectors, and an element wise difference of those two vectors.
#' @param trans A named matrix
#' @param metadata A dataframe containing columns named sample and group
#' @returns A list containing four elements: forward, backward, average, and difference vectors
process_transition <- function(trans, metadata){
  cells.idx <- metadata$group != "gene"
  mat.names <- metadata$cell[cells.idx]
  samples <- as.numeric(factor(metadata$sample[cells.idx]))
  length <- length(mat.names)
  
  same.sample <- matrix(Rfast::Outer(samples, samples,"/"), nrow=length, ncol=length) == 1
  trans <- fast_mult(trans, same.sample)
  
  tri <- upper.tri(trans, diag=TRUE)
  forward.trans <- fast_mult(trans, tri)
  diag(tri) <- FALSE
  backward.trans <- fast_mult(trans, !tri)
  backward.trans <- Rfast::transpose(backward.trans)
  
  forward.trans[forward.trans == 0] <- NA
  backward.trans[backward.trans == 0] <- NA
  
  forward.vec <- dplyr::percent_rank(c(forward.trans))
  backward.vec <- dplyr::percent_rank(c(backward.trans))
  avg.vec <- (forward.vec + backward.vec) / 2
  dif.vec <- backward.vec - forward.vec
  
  vecs <- list(forward=forward.vec, backward=backward.vec, avg=avg.vec, dif=dif.vec)
  return(vecs)
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
  cells.idx <- metadata$group != "gene"
  mat.names <- metadata$cell[cells.idx]
  length <- length(mat.names)
  
  complex.forward.vec <- c()
  complex.backward.vec <- c()
  complex.avg.vec <- c()
  complex.dif.vec <- c()
  
  if(is.null(ligand) & is.null(receptors)){
    trans <- infMat[cells.idx, cells.idx]
    
    vecs <- process_transition(trans, metadata)
    
    complex.forward.vec <- vecs$forward
    complex.backward.vec <- vecs$backward
    
    complex.avg.vec <- vecs$avg
    complex.dif.vec <- vecs$dif
  }
  else{
    complex.size <- length(receptors)
    idx.1 <- which(metadata$cell == ligand)
    
    for (receptor in receptors){
      idx.2 <- which(metadata$cell == receptor)
      
      m <- matrix(infMat[,idx.1])
      n <- matrix(infMat[idx.2,])
      trans <- Rfast::Tcrossprod(m, n)
      
      trans <- trans[cells.idx, cells.idx]
      infMat.cells <- infMat[cells.idx, cells.idx]
      
      trans <- fast_div(trans, infMat.cells)
      
      vecs <- process_transition(trans, metadata)
      
      forward.vec <- vecs$forward
      backward.vec <- vecs$backward
      avg.vec <- vecs$avg
      dif.vec <- vecs$dif
      
      if(length(complex.avg.vec) == 0){
        complex.forward.vec <- forward.vec
        complex.backward.vec <- backward.vec
        complex.avg.vec <- avg.vec
        complex.dif.vec <- dif.vec
      }
      else{
        complex.forward.vec <- complex.forward.vec + forward.vec
        complex.backward.vec <- complex.backward.vec + backward.vec
        complex.avg.vec <- complex.avg.vec + avg.vec
        complex.dif.vec <- complex.dif.vec + dif.vec
      }
    }
    
    complex.forward.vec <- complex.forward.vec / complex.size
    complex.backward.vec <- complex.backward.vec / complex.size
    complex.avg.vec <- complex.avg.vec / complex.size
    complex.dif.vec <- complex.dif.vec / complex.size
  }
  
  f.mat <- matrix(complex.forward.vec, nrow=length, ncol=length, byrow=TRUE)
  b.mat <- matrix(complex.backward.vec, nrow=length, ncol=length, byrow=TRUE)
  avg.mat <- matrix(complex.avg.vec, nrow=length, ncol=length, byrow=TRUE)
  dif.mat <- matrix(complex.dif.vec, nrow=length, ncol=length, byrow=TRUE)
  
  scores.df <- data.frame(cell=mat.names, sender=collapse::fmean(f.mat, na.rm=TRUE), receiver=collapse::fmean(b.mat, na.rm=TRUE), dif=collapse::fmean(dif.mat, na.rm=TRUE), avg=collapse::fmean(avg.mat, na.rm=TRUE))
  scores.df$hypercom <- dplyr::percent_rank(abs(scores.df$dif))
  scores.df$hypercom <- sign(scores.df$dif) * scores.df$hypercom * scores.df$avg
  scores.df <- dplyr::arrange(scores.df, desc(hypercom))
  return(scores.df)
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
#' interactions with only interactions present in the counts matrix
#' @param adj A named adjacency matrix
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param parallel A boolean of whether to run in parallel
#' @param outfile A string naming the output file
#' @returns An dataframe of ligand-receptor pairs ordered by priority score
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, parallel=TRUE, outfile="6.priority.csv"){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells <- dplyr::filter(metadata, sample != "gene")$cell
  
  if(parallel){
    cl <- parallelly::makeClusterPSOCK(parallelly::availableWorkers(), rshcmd = "qrsh", rshopts = c("-inherit", "-nostdin", "-V"), outfile = "")
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("process_transition", "score_hypercom", "fast_div", "fast_mult"), envir=environment())
  } else{
    foreach::registerDoSEQ()
  }
  
  statistics <- foreach::`%dopar%`(foreach::foreach(i=1:length(lrs$ligand), .combine=c),  {
    l <- lrs$ligand[i]
    rs <- na.omit(unlist(lrs[i, 3:dim(lrs)[2]]))
    
    percent.done <- i / length(lrs$ligand)
    print(paste(percent.done, "-", paste(c(l, rs), collapse = " ")))
    
    cc.score <- score_hypercom(infMat, metadata, l, rs)
    cc.score <- dplyr::left_join(cc.score, metadata, by=dplyr::join_by(cell))
    
    cc.score <- dplyr::mutate(cc.score, max=pmax(sender, receiver))
    cc.score$type <- cc.score$dif > 0
    rownames(cc.score) <- cc.score$cell
    cc.score <- cc.score[cells,]
    
    if(all(l == rs)){
      pair.rows <- as.data.frame(adj[l, cells])
      colnames(pair.rows) <- l
    }else{
      pair.rows <- as.data.frame(t(adj[c(l, rs), cells]))
    }
    
    cc.score$present <- do.call(pmax, pair.rows)
    pct.min <- min(colSums(pair.rows) / length(cells))
    
    cc.expressed <- dplyr::filter(cc.score, present > 0)
    
    frac <- length(cc.expressed$cell) / length(cells)
    max <- mean(cc.expressed$max)
    
    if(gini.max > 0){
      purity <- gini.max - mltools::gini_impurity(cc.expressed$group)
    }else{purity <- 1}
    
    balance <- mltools::gini_impurity(cc.expressed$type) / 0.5
    
    score <- (pct.min * frac * max * balance * purity) ^ (1/5)
    score
  })
  
  if(parallel){
    parallel::stopCluster(cl)
    gc()
  }

lrs$priority <- statistics
lrs <- dplyr::arrange(lrs, dplyr::desc(priority))
write.csv(lrs, outfile, row.names = FALSE)
return(lrs)
}

#' Load outputs from running HyperCom
#'
#' `read_vars` loads outputs from HyperCom
#'
#' This function loads the influence matrix, metadata, and ligand-receptor priority
#' from a directory and appends a supplied prefix to the variable names.
#' @param prefix A string to append to the variable name
#' @param dir A string of the directory to load files from
read_vars <- function(prefix="", dir="."){
  if(prefix != ""){
    prefix <- paste0(prefix, ".")
  }
  dir <- paste0(dir, "/")
  assign(paste0(prefix, "infMat"), readRDS(paste0(dir, "4.infMat.rds")), envir = parent.frame())
  assign(paste0(prefix, "metadata"), read.csv(paste0(dir, "5.metadata.csv")), envir = parent.frame())
  assign(paste0(prefix, "priority"), read.csv(paste0(dir, "6.priority.csv")), envir = parent.frame())
}