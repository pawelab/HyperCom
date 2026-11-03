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
#' @param sort A boolean; if TRUE, reorders adj so all genes come first, then all cells
#' @param counts A matrix with row names being genes and column names being cells; required only if sort=TRUE
#' @returns An adjacency matrix where each entry is the number of shared hyperedges between nodes
#' @export
generate_adjacency_matrix <- function(hyperedges, outfile="3.adj.rds", sort=FALSE, counts=NULL){
  
  incidence.matrix <- stats::xtabs(~node + edge, hyperedges, sparse=TRUE)
  adj <- Matrix::tcrossprod(incidence.matrix)
  diag(adj) <- 0
  rm.indices <- which(Matrix::rowSums(adj) == 0)
  
  if(length(rm.indices) > 0) {
    adj <- adj[-rm.indices, -rm.indices]
  }

  if(sort){
    if(is.null(counts)){
      stop("counts must be provided (rownames=genes, colnames=cells) when sort=TRUE")
    }
    cell.ids <- colnames(counts)
    node.ids <- rownames(adj)
    is.cell <- node.ids %in% cell.ids
    
    # sanity check: gene/cell namespaces should never collide
    if(any(node.ids[!is.cell] %in% cell.ids)){
      warning("some node names classified as genes also match cell names; check for ID collisions")
    }
    
    # genes first (FALSE), cells last (TRUE); alphabetical within each group for determinism
    ord <- order(is.cell, node.ids)
    adj <- adj[ord, ord]
  } 

  if(!is.null(outfile)){
    saveRDS(adj, outfile)
  }

  return(adj)
}

#' Apply spatial decay to the cell-cell block of a sorted adjacency matrix
#'
#' `score_spatial_decay` reweights the cell-cell block of a sorted adjacency matrix by physical distance.
#' @param adj A sorted adjacency matrix (genes first, then cells; see generate_adjacency_matrix(sort=TRUE))
#' @param coords A matrix/data.frame of cell coordinates, rownames = cell IDs, two columns (x, y), same units as radius/nugget (e.g. microns)
#' @param mode 1 = binary threshold, 2 = exponential with a small-value floor, 3 = exponential with a flat nugget region
#' @param radius Distance (same units as coords) at which decay reaches `target`. Also the hard cutoff for mode 1.
#' @param nugget Distance below which weight is held at 1 (mode 3 only)
#' @param target Weight value decay should reach at `radius` (modes 2, 3)
#' @param floor Weights below this value are snapped to 0 (mode 2 only)
#' @param outfile A string naming the output file
#' @returns The adjacency matrix A with cell-cell block (A4) spatially reweighted
#' @export
score_spatial_decay <- function(adj, coords, mode=1, radius=100, nugget=10, target=0.01, floor=1e-3, outfile=NULL){
  node.ids <- rownames(adj)
  stopifnot(identical(node.ids, colnames(adj)))
  
  cell.ids <- rownames(coords)
  is.cell <- node.ids %in% cell.ids
  n.cell <- sum(is.cell)
  n.total <- length(node.ids)
  
  # must be sorted: all cell nodes occupy the trailing block, in this exact order
  expected.pos <- (n.total - n.cell + 1):n.total
  is.sorted <- identical(which(is.cell), expected.pos)
  
  if(!is.sorted){
    stop("adj is not sorted gene-then-cell; run generate_adjacency_matrix(..., sort=TRUE) first")
  }

  cell.idx <- expected.pos
  cell.names <- node.ids[cell.idx]
  coords <- coords[cell.names, , drop=FALSE]  # align coords to adj's cell order
  
  # A4 refers to cell-cell block of adjacency matrix A
  A4 <- adj[cell.idx, cell.idx]
  
  # force full (non-symmetric-class) storage so both (i,j) and (j,i) are explicit triplets
  A4 <- as(as(A4, "generalMatrix"), "CsparseMatrix")
  A4.t <- as(A4, "TsparseMatrix")
  
  i <- A4.t@i + 1
  j <- A4.t@j + 1
  x <- A4.t@x
  
  keep <- (x > 0) & (i != j)
  i <- i[keep]; j <- j[keep]; x <- x[keep]
  
  d <- sqrt((coords[i, 1] - coords[j, 1])^2 + (coords[i, 2] - coords[j, 2])^2)
  
  # FLAG:TBD collapse 2&3 as one, and add gaussian decay. and also consider other spatial decay functions
  w <- switch(mode,
    `1` = {
      # spatial mode 1: binary classification with threshold parameter r
      wt <- as.numeric(d <= radius)
      wt
    },
    `2` = {
      # spatial mode 2: exponential decay with threshold parameter r
      lambda <- -log(target) / radius
      wt <- exp(-lambda * d)
      wt[d > radius] <- 0
      wt
    },
    `3` = {
      # spatial mode 3: exponential decay with threshold parameter r 
      lambda <- -log(target) / (radius - nugget)
      wt <- ifelse(d <= nugget, 1, exp(-lambda * (d - nugget)))
      wt[d > radius] <- 0
      wt
    },
    stop("mode must be 1, 2, or 3")
  )

  A4.new <- Matrix::sparseMatrix(i=i, j=j, x=x*w, dims=dim(A4), dimnames=dimnames(A4))
  
  adj[cell.idx, cell.idx] <- A4.new
  
  if(!is.null(outfile)){
    saveRDS(adj, outfile)  
  }

  return(adj)
}

# FLAG:MODIFIED_20261025 (v1-fast)
# ----------------------------------------------------------------------------------------------
# CHANGES 20261025 (generate_infMat)
#   The influence matrix is the same F = r (I - (1 - r) W)^-1, W = D^-1 A. What changed is how the
#   inverse is computed:
#   1. Fewer dense L x L copies. The old code held a dense identity D, the dense W, M = D - W,
#      solve()'s identity right-hand side, its LU copy and the result (~6 L^2 doubles, ~20 GB at
#      L = 20k). Both methods below build M straight from the sparse adjacency and scale in place
#      (in place, C++: scale_inplace_cpp / scale_rows_cols_inplace_cpp, which unlike collapse::setop
#      accept matrices of more than 2^31 - 1 elements), and the Cholesky inverse is taken in place
#      (chol_inverse_inplace_cpp): peak ~3 L^2 for "lu", ~1 L^2 for "chol".
#   2. method = "chol" (default when adj is symmetric, which it always is for HyperCom: adj is
#      B B^T, and the spatial kernel is symmetric). M is similar to the symmetric positive definite
#      S = D^1/2 M D^-1/2 = I - (1 - r) D^-1/2 A D^-1/2 (eigenvalues in [r, 2 - r]), so
#      F = r D^-1/2 S^-1 D^1/2, with S^-1 from a Cholesky factorisation: ~L^3 flops instead of
#      the ~8/3 L^3 of an LU-based full inverse. S has condition number <= (2 - r) / r (3 at
#      r = 0.5), so the result matches "lu" to ~1e-15 relative, but not bit for bit.
#   3. method = "lu" is the old computation (same M, same LAPACK call), only without the extra
#      copies: with the same BLAS it reproduces the old infMat bit for bit.
#   4. `threads`: BLAS threads for the inversion (default: the cores this session was granted).
#      On the SCC, R uses FlexiBLAS; its default NETLIB backend is single-threaded reference BLAS
#      (~1 GFLOP/s), which is what made the 20k inversion take hours. Start R with
#      FLEXIBLAS=OpenBLASOpenMP (or MklOpenMP) to get an optimised, multithreaded backend.
# ----------------------------------------------------------------------------------------------
#' Calculate influence matrix of hypergraph
#'
#' `generate_infMat` returns the influence matrix of the hypergraph
#'
#' This function takes an adjacency matrix and returns its influence matrix using heat diffusion.
#' @param adj A named adjacency matrix
#' @param r A double between 0 and 1 representing the restart probability
#' @param outfile A string naming the output file
#' @param coords,spatial.mode,spatial.radius,spatial.nugget,spatial.target,spatial.floor Optional spatial
#' decay applied to adj first; see `score_spatial_decay`
#' @param method "auto" (Cholesky if adj is symmetric, else LU), "chol", or "lu" (the original computation)
#' @param threads BLAS threads used for the inversion; default all cores granted to this session
#' @returns An influence matrix
#' @export
generate_infMat <- function(adj, r=0.5, outfile="4.infMat.rds",
                            coords=NULL, spatial.mode=1, spatial.radius=100,
                            spatial.nugget=10, spatial.target=0.01, spatial.floor=1e-3,
                            method=c("auto", "chol", "lu"),
                            threads=as.integer(parallelly::availableCores())){
  method <- match.arg(method)

  if(!is.null(coords)){
    adj <- score_spatial_decay(adj, coords, mode=spatial.mode, radius=spatial.radius,
                               nugget=spatial.nugget, target=spatial.target,
                               floor=spatial.floor, outfile=NULL)
  }
  len <- dim(adj)[1]
  node.names <- rownames(adj)

  if(method == "auto"){
    method <- if(Matrix::isSymmetric(adj)) "chol" else "lu"
  }

  old.threads <- blas_threads_cpp(as.integer(threads))
  if(old.threads > 0){
    on.exit(blas_threads_cpp(old.threads), add = TRUE)
  }
  diag.idx <- seq.int(1, as.numeric(len) * len, by = len + 1)
  d <- Matrix::rowSums(adj)

  if(method == "lu"){
    # identical operations to the old D - (1 - r) * W, without the dense D
    M <- as.matrix(-((1 - r) * (adj / d)))
    M[diag.idx] <- M[diag.idx] + 1
    infMat <- solve(M)
    rm(M)
    scale_inplace_cpp(infMat, r, as.integer(threads))          # in place, any size (collapse::setop refuses > 2^31 elements)
  } else {
    # S = I - (1 - r) D^-1/2 A D^-1/2, symmetric positive definite
    s <- 1 / sqrt(d)
    S <- as.matrix(Matrix::Diagonal(x = s) %*% adj %*% Matrix::Diagonal(x = s))
    scale_inplace_cpp(S, -(1 - r), as.integer(threads))
    S[diag.idx] <- S[diag.idx] + 1
    # S^-1 in place: LAPACK dpotrf + dpotri on the one buffer (src/infmat_lapack.cpp), the same
    # routines chol() + chol2inv() call on copies -- one dense L x L copy alive instead of two,
    # bit-identical result
    chol_inverse_inplace_cpp(S, as.integer(threads))
    infMat <- S
    rm(S)
    # F = r D^-1/2 S^-1 D^1/2: scale row i by r s_i, column j by 1 / s_j, in place
    scale_rows_cols_inplace_cpp(infMat, r * s, 1 / s, as.integer(threads))
  }

  dimnames(infMat) <- list(node.names, node.names)

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

  if(is.null(group)){
    group.vals <- "none"
  }
  else{
    group.vals <- as.character(metadata[, group])
  }

  clean.df <- data.frame(cell=metadata$cell, sample=as.character(metadata[, sample]), group=group.vals, check.names=FALSE)
  clean.df$sample <- factor(clean.df$sample, levels=c(unique(clean.df$sample), "gene"))
  clean.df$group <- factor(clean.df$group, levels=c(unique(clean.df$group), "gene"))

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

#' Run the HyperCom pipeline
#'
#' `HyperCom` Runs the HyperCom pipeline
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency-matrix`, `generate_infMat`, `generate-metadata`, and `prioritize_lr`.
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param counts A matrix with row names being genes and column names being cells
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param sample A string that is the column name for sample in the provided metadata
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel (forked workers; see `prioritize_lr`)
#' @param cl No longer used (kept so existing calls do not break); parallelism is set by `parallel`
#' @param outdir output directory of pipeline run
#' @param coords Optional matrix/data.frame of cell coordinates: rownames = cell IDs (colnames of counts),
#' two columns (x, y), same units as spatial.radius / spatial.nugget (e.g. microns). NULL = non-spatial run.
#' @param spatial.mode Decay kernel (see `score_spatial_decay`): 1 = binary cutoff, 2 = exponential with floor, 3 = nugget + exponential
#' @param spatial.radius Distance at which the decay reaches spatial.target (also the cutoff for mode 1)
#' @param spatial.nugget Distance below which the weight is held at 1 (mode 3 only)
#' @param spatial.target Weight the decay reaches at spatial.radius (modes 2, 3)
#' @param spatial.floor Weights below this are set to 0 (mode 2 only)
#' @param infmat.method Inversion method for `generate_infMat`: "auto", "chol" or "lu"
#' @param cores Cores used for the inversion (BLAS threads) and the LR scoring (see `prioritize_lr`)
#' @param workers Forked workers for the LR scoring; 1 = one grid at a time on `cores` threads
#' @param grid.mem.gb Memory budget in GB per grid for the LR scoring (see `prioritize_lr`)
#' @returns 0 if run successfully
#' @export
HyperCom <- function(lrdb, counts, metadata, sample, group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="HyperCom",
                     coords=NULL, spatial.mode=3, spatial.radius=100, spatial.nugget=10,
                     spatial.target=0.01, spatial.floor=1e-3,
                     infmat.method="auto",
                     cores=if(parallel) as.integer(parallelly::availableCores()) else 1L,
                     workers=1L, grid.mem.gb=1){
  spatial <- !is.null(coords)

  # [20261016] NEW: check coords before any heavy work
  if(spatial){
    if(is.null(rownames(coords))){
      stop("coords must have rownames = cell IDs (the colnames of counts)")
    }
    missing.cells <- setdiff(colnames(counts), rownames(coords))
    if(length(missing.cells) > 0){
      stop(length(missing.cells), " of ", ncol(counts), " cells have no row in coords (e.g. '",
           missing.cells[1], "'); coords rownames must be the colnames of counts")
    }
  }

  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  setwd(outdir)

  # [20261016] keep the original genes x cells matrix: filter_counts() returns CELLS x genes,
  # and the sort below needs the cell IDs, i.e. colnames of the ORIGINAL orientation
  counts.raw <- counts
  counts <- filter_counts(counts, lrdb)
  lrs <- generate_lrs_table(counts, lrdb)
  hyperedges <- generate_hyperedges(counts, lrs)

  # [20261016] CHANGED: sort only when a spatial kernel will be applied (score_spatial_decay needs
  # the cell nodes as one trailing block); with coords = NULL this is the original non-spatial
  # call, so the node order, and therefore every number, matches the original HyperCom.
  adj <- generate_adjacency_matrix(hyperedges, sort = spatial,
                                   counts = if(spatial) counts.raw else NULL)
  rm(counts.raw)

  # [20261016] CHANGED: coords and kernel settings come from the arguments (coords = NULL skips
  # the decay inside generate_infMat)
  # [20261025] infmat.method / cores forwarded (Cholesky inverse, multithreaded BLAS)
  infMat <- generate_infMat(adj, coords=coords, spatial.mode=spatial.mode, spatial.radius=spatial.radius,
                            spatial.nugget=spatial.nugget, spatial.target=spatial.target,
                            spatial.floor=spatial.floor, method=infmat.method, threads=cores)

  metadata <- generate_metadata(infMat, metadata, sample, group, keep)

  # [20261016] NEW: everything that depends only on metadata, computed once per run.
  # [20261025] infMat.cell is no longer built: the scoring kernel reads the cell block of infMat
  #            in place, which saves a dense n x n copy (3.1 GB at 20k cells).
  cell.info <- generate_cell_info(metadata)

  # [20261016] CHANGED: arguments by name; cell.info passed down.
  # NOTE: adj here is the UNDECAYED adjacency (generate_infMat decays its own copy). That is
  # numerically irrelevant: calculate_priority() reads only adj[genes, cells], which the kernel
  # never touches.
  priority <- prioritize_lr(lrs=lrs, adj=adj, infMat=infMat, metadata=metadata,
                            weights=weights, significance=significance,
                            parallel=parallel, cl=cl, cell.info=cell.info,
                            cores=cores, workers=workers, grid.mem.gb=grid.mem.gb)

  return(0)
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

#' Pre-compute the cell-level information used by every LR pair
#'
#' `generate_cell_info` returns everything in HyperCom scoring that depends only on metadata.
#'
#' Computed once per run instead of once per receptor subunit. `sample.id` is NULL when every cell
#' comes from one sample: the same-sample mask would be all TRUE, so it is skipped.
#' @param metadata A dataframe containing columns named cell, sample and group (from `generate_metadata`)
#' @returns A list: cells.idx (logical, which infMat rows/cols are cells), mat.names (cell names in
#' infMat order), n (number of cells), sample.id (integer sample per cell, or NULL with one sample)
#' @export
generate_cell_info <- function(metadata){
  cells.idx <- metadata$group != "gene"
  mat.names <- metadata$cell[cells.idx]
  samples <- as.integer(factor(metadata$sample[cells.idx]))
  n <- length(mat.names)

  # [20261025] one id per cell; the kernel tests sample_i == sample_j (was an n x n logical mask)
  sample.id <- if(length(unique(samples)) > 1) samples else NULL

  list(cells.idx=cells.idx, mat.names=mat.names, n=n, sample.id=sample.id)
}

#' Receptor subunits of one row of the LR table
#'
#' `pair_receptors` returns the receptor subunits of row i, read exactly as `calculate_priority`
#' always has (columns 3 onward, NA dropped), so every function agrees on what a pair's subunits are.
#' @param lrs A dataframe of ligand-receptor interactions (see `prioritize_lr`)
#' @param i Row index
#' @returns A character vector of receptor subunit names
#' @noRd
pair_receptors <- function(lrs, i){
  stats::na.omit(unlist(lrs[i, 3:dim(lrs)[2]]))
}

# FLAG:NEW_20261016
#' Key naming one (ligand, receptor subunit) grid
#'
#' @param ligand A string of the ligand name
#' @param receptors A vector of receptor subunit names
#' @returns A character vector "ligand|receptor", one per subunit
#' @noRd
subunit_key <- function(ligand, receptors){
  paste(ligand, as.character(receptors), sep = "|")
}

# FLAG:NEW_20261016
#' Every unique (ligand, receptor subunit) grid a table of LR pairs needs
#'
#' `all_subunit_combos` lists each (ligand, subunit) combination once. Many LR pairs share
#' subunits (e.g. the same ligand with several complexes containing ITGB1), and the old code
#' recomputed the full n x n grid for every pair that used it. Here each grid is computed once and
#' reused, which is where part of the speed-up comes from.
#' @param lrs A dataframe of ligand-receptor interactions (see `prioritize_lr`)
#' @returns A dataframe with columns ligand and receptor, one row per unique grid
#' @noRd
all_subunit_combos <- function(lrs){
  combos <- do.call(rbind, lapply(seq_len(nrow(lrs)), function(i){
    data.frame(ligand = lrs$ligand[i], receptor = as.character(pair_receptors(lrs, i)),
               stringsAsFactors = FALSE)
  }))
  combos[!duplicated(subunit_key(combos$ligand, combos$receptor)), , drop = FALSE]
}

# FLAG:NEW_20261025
#' Split the cores between forked workers and kernel threads
#'
#' `grid_schedule` returns how many grids run at once (workers) and how many OpenMP threads each
#' uses. GNU OpenMP is not fork-safe, so forked workers always run the kernel on one thread.
#' @param cores Total cores to use
#' @param workers Requested forked workers (1 = one grid at a time on all cores)
#' @param n.grids Number of grids to score
#' @returns A list with workers and threads
#' @noRd
grid_schedule <- function(cores, workers, n.grids){
  cores <- max(1L, as.integer(cores))
  workers <- max(1L, min(as.integer(workers), cores, n.grids))
  if(workers > 1 && .Platform$OS.type == "windows"){
    warning("forked workers are not available on Windows; using workers = 1")
    workers <- 1L
  }
  list(workers = workers, threads = if(workers > 1) 1L else cores)
}

# FLAG:MODIFIED_20261025
#' Sender and receiver scores for every unique (ligand, subunit) grid
#'
#' `score_subunits` scores Q[i, j] = infMat[i, ligand] * infMat[subunit, j] / infMat[i, j] over
#' cells with the C++ kernel `rank_grid_cpp`, which reads infMat in place and keeps at most
#' `grid.mem.gb` per grid.
#'
#' Every ligand column and every receptor row is read out of infMat once, restricted to cells.
#' With workers > 1 the grids are spread over forked workers (parallel::mclapply), which share
#' infMat through copy-on-write memory; otherwise grids run one at a time on `cores` threads.
#' @param combos A dataframe from `all_subunit_combos`
#' @param infMat A named influence matrix
#' @param cell.info A list from `generate_cell_info`
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param cores Total cores
#' @param workers Forked workers (1 = serial over grids, multithreaded within each)
#' @param grid.mem.gb Memory budget per grid, in GB
#' @returns A named list (names from `subunit_key`) of lists with sender and receiver
#' @noRd
score_subunits <- function(combos, infMat, cell.info, metadata, cores=1L, workers=1L, grid.mem.gb=1){
  ligs <- unique(combos$ligand)
  recs <- unique(combos$receptor)
  idx.l <- match(ligs, metadata$cell)
  idx.r <- match(recs, metadata$cell)
  if(anyNA(idx.l) || anyNA(idx.r)){
    stop("not nodes of infMat: ", paste(c(ligs[is.na(idx.l)], recs[is.na(idx.r)]), collapse = ", "))
  }
  if(!identical(rownames(infMat), metadata$cell)){
    stop("metadata rows must be in infMat order (use generate_metadata(infMat, ...))")
  }
  lig.cols <- infMat[cell.info$cells.idx, idx.l, drop = FALSE]
  rec.rows <- infMat[idx.r, cell.info$cells.idx, drop = FALSE]
  cell.pos <- which(cell.info$cells.idx) - 1L          # 0-based positions of the cells in infMat
  budget <- grid.mem.gb * 1024^3

  sched <- grid_schedule(cores, workers, nrow(combos))
  one_grid <- function(k){
    s <- rank_grid_cpp(infMat, cell.pos,
                       a = lig.cols[, match(combos$ligand[k], ligs)],
                       b = rec.rows[match(combos$receptor[k], recs), ],
                       sample_id = cell.info$sample.id,
                       nthreads = sched$threads, budget_bytes = budget)
    list(sender = s$sender, receiver = s$receiver)
  }

  idx <- seq_len(nrow(combos))
  if(sched$workers > 1){
    # mc.preschedule = FALSE: one fork per grid, handed out as workers free up, so a slow grid
    # does not hold back a whole pre-assigned batch
    res <- parallel::mclapply(idx, one_grid, mc.cores = sched$workers, mc.preschedule = FALSE)
    # a grid that raised an error comes back as a "try-error"; a worker killed by the OS (e.g. out
    # of memory) comes back as NULL. Either way, stop instead of scoring pairs with missing grids.
    failed <- vapply(res, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
    if(any(failed)){
      first <- res[[which(failed)[1]]]
      stop(sum(failed), " grid(s) failed in workers; first: ",
           if(is.null(first)) "worker returned NULL (killed? check memory per worker)" else as.character(first))
    }
  } else {
    res <- lapply(idx, one_grid)
  }
  names(res) <- subunit_key(combos$ligand, combos$receptor)
  res
}

# FLAG:NEW_20261016
#' Combine one pair's subunit grids into its sender and receiver scores
#'
#' `combine_subunits` sums the subunits' sender/receiver in subunit order and divides by the number
#' of subunits, the same arithmetic as the old `score_hypercom` loop; a single subunit is passed
#' through untouched, as before.
#' @param subunit.scores A named list from `score_subunits`
#' @param ligand A string of the ligand name
#' @param receptors A vector of receptor subunit names
#' @returns A list with sender and receiver
#' @noRd
combine_subunits <- function(subunit.scores, ligand, receptors){
  keys <- subunit_key(ligand, receptors)
  missing <- setdiff(keys, names(subunit.scores))
  if(length(missing) > 0) stop("no precomputed grid for: ", paste(missing, collapse = ", "))
  parts <- subunit.scores[keys]
  if(length(parts) == 1) return(parts[[1]])
  k <- length(parts)
  list(sender   = Reduce(`+`, lapply(parts, `[[`, "sender"))   / k,
       receiver = Reduce(`+`, lapply(parts, `[[`, "receiver")) / k)
}

# FLAG:MODIFIED_20261016
# ----------------------------------------------------------------------------------------------
# CHANGES 20261016 (score_sr)
#   - Now scores a whole receptor COMPLEX in one call (`receptors` may have several subunits),
#     averaging over subunits as score_hypercom used to. A single subunit gives the old result.
#   - The work moved into generate_cell_info / score_subunits / rank_grid (see the header above).
#   - Signature: the old arguments keep their positions (infMat, metadata, ligand, receptor[s]), so
#     old positional calls still work, and `receptor =` still matches `receptors` by partial
#     matching. New optional arguments infMat.cell and cell.info are appended; when NULL they are
#     built here, so standalone use (e.g. plotting one pair) needs nothing new.
# ----------------------------------------------------------------------------------------------
#' Calculate vectors used in HyperCom scoring
#'
#' `score_sr` returns a dataframe of sender and receiver scores use by HyperCom scoring
#'
#' This function takes a named matrix and a 2 column dataframe with the sender and receiver score. If ligand and receptor args are NULL calculates network score.
#' @param infMat A named influence matrix
#' @param metadata A dataframe containing columns named sample and group
#' @param ligand A string of the ligand name
#' @param receptors A vector of receptor subunit names (one or more)
#' @param infMat.cell Ignored (kept so existing calls work); the kernel reads infMat in place
#' @param cell.info Optional list from `generate_cell_info`; built here if NULL
#' @param cores Kernel threads
#' @returns A dataframe with the sender and receiver score
#' @noRd
score_sr <- function(infMat, metadata, ligand=NULL, receptors=NULL, infMat.cell=NULL, cell.info=NULL, cores=1L){
  if(is.null(cell.info)){
    cell.info <- generate_cell_info(metadata)
  }

  if(is.null(ligand) | is.null(receptors)){
    # network score: rank the cell block of infMat itself (a = b = NULL)
    # [20261025] kernel instead of a private copy of infMat.cell
    s <- rank_grid_cpp(infMat, which(cell.info$cells.idx) - 1L, sample_id = cell.info$sample.id,
                       nthreads = cores)
  } else {
    combos <- data.frame(ligand = ligand, receptor = as.character(receptors), stringsAsFactors = FALSE)
    combos <- combos[!duplicated(combos$receptor), , drop = FALSE]
    s <- combine_subunits(score_subunits(combos, infMat, cell.info, metadata, cores = cores),
                          ligand, receptors)
  }

  data.frame(cell=cell.info$mat.names, sender=s$sender, receiver=s$receiver)
}

# FLAG:MODIFIED_20261016
# ----------------------------------------------------------------------------------------------
# CHANGES 20261016 (score_hypercom)
#   - The per-subunit score_sr() loop is gone: inside prioritize_lr the subunit grids are already
#     computed (subunit.scores) and only combined here; called on its own, it gets sender and
#     receiver from one score_sr() call, which now handles complexes itself.
#   - Everything after sender/receiver (dif, avg, max, hypercom, ordering) is unchanged.
#   - Signature: old arguments keep their positions; infMat.cell, cell.info, subunit.scores are
#     appended and optional, so existing calls (e.g. build_ridge_plot in section0.R) still work.
# ----------------------------------------------------------------------------------------------
#' Calculate the HyperCom score for each cell
#'
#' `score_hypercom` returns a dataframe containing the HyperCom score for each cell
#'
#' This function takes a named influence matrix, the corresponding metadata, the ligand, and receptor
#' to give each cell a HyperCom score. If both the ligand and receptor are specified this is the LR Pair Score
#' otherwise the score is the Network Score. If a receptor complex is specified, the average of the score vectors is used.
#' @param infMat A named influence matrix (may be NULL when subunit.scores is given)
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param ligand A string of the ligand name
#' @param receptors A vector containing string elements of the receptor names
#' @param infMat.cell Ignored (kept so existing calls work)
#' @param cell.info Optional list from `generate_cell_info`; built if NULL
#' @param subunit.scores Optional named list from the internal `score_subunits`; when given, no grid is computed here
#' @returns A dataframe containing the sender score, receiver score, avg, dif, and HyperCom score
#' @export
score_hypercom <- function(infMat, metadata, ligand=NULL, receptors=NULL, infMat.cell=NULL, cell.info=NULL, subunit.scores=NULL){
  if(is.null(subunit.scores)){
    # [20261016] one call for the whole complex (was: one score_sr() call per subunit)
    sr <- score_sr(infMat, metadata, ligand, receptors, infMat.cell=infMat.cell, cell.info=cell.info)
  } else {
    # [20261016] precomputed grids, just combined
    if(is.null(cell.info)){
      cell.info <- generate_cell_info(metadata)
    }
    s <- combine_subunits(subunit.scores, ligand, receptors)
    sr <- data.frame(cell=cell.info$mat.names, sender=s$sender, receiver=s$receiver)
  }

  # ---- unchanged from here ----
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

# FLAG:MODIFIED_20261016
# ----------------------------------------------------------------------------------------------
# CHANGES 20261016 (calculate_priority)
#   - Only the cheap O(n) per-pair part is left here. All n x n work (the grids) is done once,
#     before this is called, by score_subunits() in prioritize_lr(), and read from subunit.scores.
#   - It therefore no longer needs infMat; signature is now
#       calculate_priority(lrs, subunit.scores, cell.info, adj, metadata, weight)
#     (internal, @noRd; prioritize_lr is the only caller).
#   - foreach %dopar% replaced by a serial vapply: parallelism now lives in score_subunits(),
#     where the time is. The per-pair body is otherwise the original line for line.
# ----------------------------------------------------------------------------------------------
#' Score ligand receptor pairs based on signalling strength and specificity
#'
#' `calculate_priority` returns an ordered list of ligand-receptor pairs
#'
#' This function takes list of all ligand-receptor pairs in the influence matrix,
#' the precomputed subunit grids, the adjacency matrix, and the corresponding metadata
#' to score ligand-receptor pairs.
#' @param lrs A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor
#' interactions with only interactions present in the counts matrix
#' @param subunit.scores A named list from `score_subunits` covering every (ligand, subunit) in lrs
#' @param cell.info A list from `generate_cell_info`
#' @param adj A named adjacency matrix
#' @param metadata A dataframe containing columns named cell, sample and group
#' @param weight a value between 0 to 1 for the importance of number of cells involved in an interaction
#' @return A numeric vector of priority scores
#' @noRd
calculate_priority <- function(lrs, subunit.scores, cell.info, adj, metadata, weight=1){
  num.groups <- length(unique(metadata$group)) - 1
  gini.max <- 1-(1/num.groups)
  cells <- dplyr::filter(metadata, sample != "gene")$cell

  # [20261016] vapply instead of foreach %dopar%
  priority <- vapply(seq_len(nrow(lrs)), function(i){
    l <- lrs$ligand[i]
    rs <- pair_receptors(lrs, i)                     # [20261016] same read as before, shared helper
    title <- paste(c(l, rs), collapse = "_")
    # [20261016] reads the precomputed grids; no score_sr() call here
    cc.score <- score_hypercom(infMat=NULL, metadata=metadata, ligand=l, receptors=rs,
                               cell.info=cell.info, subunit.scores=subunit.scores)
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
  }, numeric(1))

  return(priority)
}

# FLAG:MODIFIED_20261016
# ----------------------------------------------------------------------------------------------
# CHANGES 20261016 (prioritize_lr)
#   1. Two phases instead of one loop of full score_sr() calls per pair:
#        (a) score_subunits(): every UNIQUE (ligand, subunit) grid in lrs (and in the background
#            set when significance = TRUE) is computed once, serially or over forked workers;
#        (b) calculate_priority(): the cheap per-pair scoring, reading those grids.
#      The grids do not depend on `weight`, so with several weights they are computed once, not
#      once per weight as before.
#   2. Parallelism: forked workers (parallel::mclapply) over grids, `cores` of them. The old PSOCK
#      cluster (started over qrsh) copied infMat to every worker and, through
#      foreach(.packages = "HyperCom"), could run the INSTALLED package's functions instead of these.
#      `cl` is no longer used; a supplied `cl` triggers a warning and is not stopped here (whoever
#      created it owns it). On Windows, where forking is unavailable, cores falls back to 1.
#   3. data.table's thread count is set to the cores this session was granted while this function
#      runs (frankv sorts with it in the serial path) and restored on exit.
#   4. New optional arguments, appended so old positional calls still work:
#        infMat.cell  dense cell x cell block of infMat (built here if NULL)
#        cell.info    from generate_cell_info(metadata) (built here if NULL)
#        cores        workers; default = parallelly::availableCores() (NSLOTS on SGE, cgroups,
#                     mc.cores) when parallel = TRUE, else 1
#   Output (columns, ordering, priority values, p / p.adj, CSVs in outdir) is unchanged.
#
# CHANGES 20261025 (prioritize_lr)
#   - infMat.cell is ignored (the kernel reads infMat in place); cell.info carries sample.id
#     instead of the n x n mask.
#   - `cores` is now the total core count. New `workers` (default 1: one grid at a time on `cores`
#     threads, memory independent of cores) and `grid.mem.gb` (memory budget per grid, default 1).
#   - data.table is no longer used.
# ----------------------------------------------------------------------------------------------
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
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel (sets the default of `cores`)
#' @param cl No longer used; kept so existing calls do not break (a warning is given if supplied)
#' @param outdir A string naming the output directory
#' @param infMat.cell Ignored (kept so existing calls work)
#' @param cell.info Optional list from `generate_cell_info(metadata)`; built here if NULL
#' @param cores Total cores for the per-grid work. Defaults to the cores this R session was granted
#' (`parallelly::availableCores()`) when parallel = TRUE, else 1.
#' @param workers Forked workers. 1 (default) scores one grid at a time on `cores` threads; k > 1
#' scores k grids at once on one thread each (needs k * grid.mem.gb).
#' @param grid.mem.gb Memory budget per grid in GB. Smaller budgets mean more passes over infMat
#' per grid, not different results.
#' @returns A dataframe of ligand-receptor pairs ordered by priority score for each weight
#' @export
prioritize_lr <- function(lrs, adj, infMat, metadata, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="6.priority",
                          infMat.cell=NULL, cell.info=NULL,
                          cores=if(parallel) as.integer(parallelly::availableCores()) else 1L,
                          workers=1L, grid.mem.gb=1){

  if(is.null(cell.info)){
    cell.info <- generate_cell_info(metadata)
  }
  if(!is.list(cell.info) || !all(c("cells.idx", "mat.names", "n", "sample.id") %in% names(cell.info))){
    stop("cell.info must come from generate_cell_info(metadata)")
  }
  if(!is.null(infMat.cell)){
    message("prioritize_lr: infMat.cell is no longer used; the cell block is read from infMat in place")
  }

  if(!is.null(cl)){
    warning("prioritize_lr: `cl` is no longer used; parallelism is set with `cores` / `workers`")
  }

  if(significance){
    background.lrs <- generate_background_lrs(lrs, metadata)
  }

  # ---- [20261016] phase (a): every unique (ligand, subunit) grid, once --------------------------
  combos <- all_subunit_combos(lrs)
  if(significance){
    combos <- rbind(combos, all_subunit_combos(background.lrs))
    combos <- combos[!duplicated(subunit_key(combos$ligand, combos$receptor)), , drop = FALSE]
  }
  n.calls <- sum(vapply(seq_len(nrow(lrs)), function(i) length(pair_receptors(lrs, i)), integer(1)))
  sched <- grid_schedule(cores, workers, nrow(combos))
  message(sprintf("<prioritize_lr> %d pairs, %d subunit calls -> %d unique grids, %d worker(s) x %d thread(s), %.2f GB per grid, %d cells",
                  nrow(lrs), n.calls, nrow(combos), sched$workers, sched$threads, grid.mem.gb, cell.info$n))

  subunit.scores <- score_subunits(combos, infMat, cell.info, metadata,
                                   cores = cores, workers = workers, grid.mem.gb = grid.mem.gb)

  # ---- phase (b): per-weight scoring, unchanged apart from the calculate_priority() arguments ----
  lrs.list <- list()

  for(i in seq(1, length(weights))){
    weight <- weights[i]
    print(weight)
    priority <- calculate_priority(lrs, subunit.scores, cell.info, adj, metadata, weight)

    current.lrs <- lrs

    if(significance){
      current.background <- background.lrs
      background.priority <- calculate_priority(current.background, subunit.scores, cell.info, adj, metadata, weight)
    }

    current.lrs$priority <- priority
    current.lrs <- dplyr::arrange(current.lrs, dplyr::desc(priority))

    if(significance){
      current.lrs$p <- as.numeric(lapply(current.lrs$priority, FUN=permutation_test, background=background.priority))
      current.lrs$p.adj <- stats::p.adjust(current.lrs$p, method = "BH")
    }

    if(!is.null(outdir)){
      dir.create(file.path(outdir), showWarnings = FALSE)
      outfile <- paste0(outdir, "/priority.", weight, ".csv")
      utils::write.csv(current.lrs, outfile, row.names = FALSE)
    }

    lrs.list[[i]] <- current.lrs
  }

  # [20261016] removed: parallel::stopCluster(cl) -- there is no cluster to stop; forked workers
  # exit when mclapply returns.

  names(lrs.list) <- weights

  return(lrs.list)
}

#' Get matrix with dimnames for given layer and slot
#'
#' `get_seurat_layer` returns the matrix with dimnames for given layer and slot
#'
#' This function takes a Seurat object and filters it to
#' just ligand and receptor genes with greater than 0 expression.
#' @param seurat A seurat object
#' @param assay assay to get data from
#' @param layer layer to get data from
#' @returns A matrix with dimnames for given layer and slot
get_seurat_layer <- function(seurat, assay = "RNA", layer = "counts") {
  a <- seurat[[assay]]

  # legacy v3/v4: slot-based storage, use the long-standing accessor
  if (!inherits(a, "Assay5")) {
    m <- SeuratObject::GetAssayData(a, layer = layer)
    if (!length(m))
      stop("layer '", layer, "' is empty or absent in assay '", assay, "'")
    return(m)
  }

  lyrs <- SeuratObject::Layers(a, search = layer)
  if (!length(lyrs))
    stop("no layer matching '", layer, "' in assay '", assay,
         "'; available: ", paste(SeuratObject::Layers(a), collapse = ", "))

  if (length(lyrs) == 1L)
    return(SeuratObject::LayerData(a, layer = lyrs))

  return(SeuratObject::LayerData(SeuratObject::JoinLayers(a), layer = layer))
}

#' Run the HyperCom pipeline on a seurat object
#'
#' `HyperCom_seurat` Runs the HyperCom pipeline on a seurat object
#'
#' This function runs `filter_counts`, `generate_lrs_table`, `generate_hyperedges`, `generate_adjacency-matrix`, `generate_infMat`, `generate-metadata`, and `prioritize_lr`.
#' @param seurat A seurat object
#' @param lrdb A dataframe of ligand-receptor interactions where the first column is named edges, the second column is named ligand,
#' and subsequent columns are named r1, ..., rn where n is the maximum number of genes in a receptor complex in the list of ligand-receptor interactions
#' @param sample A string that is the column name for sample in the provided metadata
#' @param assay assay to get data from
#' @param layer layer to get data from
#' @param group A string that is the column name for group in the provided metadata
#' @param keep A vector of strings that are additional column names to keep in the provided metadata
#' @param weights a vector of values between 0 to 1 for the importance of number of cells involved in an interaction
#' @param significance A boolean of whether to generate p-values for ligand-receptor interactions
#' @param parallel A boolean of whether to run in parallel (forked workers; see `prioritize_lr`)
#' @param cl No longer used (kept so existing calls do not break)
#' @param outdir output directory of pipeline run
#' @param coords Optional cell coordinates, rownames = cell IDs (colnames of the Seurat object); NULL = non-spatial. See `HyperCom`.
#' @param spatial.mode,spatial.radius,spatial.nugget,spatial.target,spatial.floor Decay kernel settings; see `HyperCom`
#' @returns 0 if run successfully
#' @export
HyperCom_seurat <- function(seurat, lrdb, sample, assay="RNA", layer="counts", group=NULL, keep=NULL, weights=c(1), significance=FALSE, parallel=FALSE, cl=NULL, outdir="HyperCom",
                            coords=NULL, spatial.mode=3, spatial.radius=100, spatial.nugget=10,
                            spatial.target=0.01, spatial.floor=1e-3){
  counts <- get_seurat_layer(seurat, assay, layer)
  metadata <- seurat@meta.data

  # [20261016] CHANGED: spatial arguments forwarded by name
  status <- HyperCom(lrdb, counts, metadata, sample, group, keep, weights, significance, parallel, cl, outdir,
                     coords=coords, spatial.mode=spatial.mode, spatial.radius=spatial.radius,
                     spatial.nugget=spatial.nugget, spatial.target=spatial.target,
                     spatial.floor=spatial.floor)

  return(status)
}
