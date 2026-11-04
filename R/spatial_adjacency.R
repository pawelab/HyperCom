# FLAG:NEW_20261104 (v1-fast)
# ----------------------------------------------------------------------------------------------
# The spatial adjacency, built without the dense pre-decay cell-cell block.
#
# generate_adjacency_matrix() forms B B' over every node, and score_spatial_decay() afterwards
# zeroes the cell pairs beyond the kernel's reach. Before the decay two cells are adjacent whenever
# they share any expressed LR gene, so that block is nearly dense: 76.5% at 20k cells (3.0e8
# nonzeros), ~1.9e9 at 50k, which is the 2^31 - 1 limit of R's sparse matrices; past ~40k cells
# the pipeline cannot run that way. Every kernel of score_spatial_decay() is 0 beyond `radius`, so
# generate_spatial_adjacency_matrix() forms the cell-cell block only for pairs that can be within
# `radius`: the cells are binned into square tiles of side `radius`, and each tile is multiplied
# against itself and its 8 neighbours. The shared-hyperedge counts, distances and weights are
# computed with the same arithmetic as before, so the result equals
# score_spatial_decay(generate_adjacency_matrix(sort = TRUE)) entry for entry.
# ----------------------------------------------------------------------------------------------

#' Spatial decay weights
#'
#' `spatial_weights` returns the weight of a cell pair at distance `d` for one of the kernels of
#' `score_spatial_decay`. Every mode is 0 beyond `radius`.
#' @param d Distances (same units as radius / nugget)
#' @param mode 1 = binary cutoff at radius; 2 = exponential reaching `target` at radius;
#' 3 = 1 up to `nugget`, then exponential reaching `target` at radius
#' @param radius,nugget,target Kernel parameters (see `score_spatial_decay`)
#' @param floor Not used by the current kernels (kept so existing calls run unchanged)
#' @returns A numeric vector of weights, one per distance
#' @export
spatial_weights <- function(d, mode=1, radius=100, nugget=10, target=0.01, floor=1e-3){
  if(!(length(mode) == 1 && mode %in% c(1, 2, 3))){
    stop("mode must be 1, 2, or 3")
  }
  # FLAG:TBD collapse 2&3 as one, and add gaussian decay. and also consider other spatial decay functions
  if(mode == 1){
    # spatial mode 1: binary classification with threshold parameter r
    wt <- as.numeric(d <= radius)
  } else if(mode == 2){
    # spatial mode 2: exponential decay with threshold parameter r
    lambda <- -log(target) / radius
    wt <- exp(-lambda * d)
    wt[d > radius] <- 0
  } else {
    # spatial mode 3: exponential decay with threshold parameter r
    lambda <- -log(target) / (radius - nugget)
    wt <- ifelse(d <= nugget, 1, exp(-lambda * (d - nugget)))
    wt[d > radius] <- 0
  }
  wt
}

#' Spatially decayed adjacency matrix of the hypergraph, built tile by tile
#'
#' `generate_spatial_adjacency_matrix` returns the same matrix as
#' `score_spatial_decay(generate_adjacency_matrix(hyperedges, sort = TRUE, counts = counts), coords, ...)`
#' (genes first, then cells, alphabetical within each), but never forms the nearly dense cell-cell
#' block of B B' that those two steps go through. Cell pairs are formed only within neighbouring
#' tiles of side `radius`, beyond which every kernel is 0. Memory and time grow with the number of
#' cell pairs within `radius`, not with the square of the number of cells.
#'
#' Normally called through `generate_adjacency_matrix(..., coords = coords)` or `HyperCom(coords = )`.
#' @param hyperedges A dataframe of hyperedges (see `generate_hyperedges`)
#' @param cell.ids The cell IDs (colnames of the genes x cells counts matrix)
#' @param coords Cell coordinates: rownames = cell IDs, first two columns x, y
#' @param mode,radius,nugget,target,floor Kernel settings (see `spatial_weights`)
#' @param verbose Print a one-line summary
#' @returns A sparse adjacency matrix (dgCMatrix)
#' @export
generate_spatial_adjacency_matrix <- function(hyperedges, cell.ids, coords, mode=1, radius=100, nugget=10,
                                              target=0.01, floor=1e-3, verbose=TRUE){
  if(!(length(mode) == 1 && mode %in% c(1, 2, 3))){
    stop("mode must be 1, 2, or 3")
  }
  if(!(is.numeric(radius) && length(radius) == 1 && is.finite(radius) && radius > 0)){
    stop("radius must be a positive number")
  }
  if(is.null(rownames(coords))){
    stop("coords must have rownames = cell IDs")
  }
  t0 <- Sys.time()

  # nodes x hyperedges, exactly as generate_adjacency_matrix()
  inc <- stats::xtabs(~node + edge, hyperedges, sparse=TRUE)
  node.ids <- rownames(inc)
  is.cell <- node.ids %in% cell.ids

  # generate_adjacency_matrix() drops nodes of degree 0 (diagonal excluded). A cell always shares
  # the hyperedge of each gene it expresses with that gene's node, so only gene nodes can drop;
  # their degrees are row sums of the gene rows of B B'.
  G <- Matrix::tcrossprod(inc[!is.cell, , drop=FALSE], inc)
  gene.pos <- which(!is.cell)
  self <- G[cbind(seq_len(nrow(G)), gene.pos)]
  drop <- node.ids[gene.pos][(Matrix::rowSums(G) - self) == 0]
  keep <- !(node.ids %in% drop)
  node.ids <- node.ids[keep]; is.cell <- is.cell[keep]; inc <- inc[keep, , drop=FALSE]

  # the order of generate_adjacency_matrix(sort = TRUE)
  ord <- order(is.cell, node.ids)
  node.ids <- node.ids[ord]; is.cell <- is.cell[ord]; inc <- inc[ord, , drop=FALSE]
  L <- length(node.ids)
  gene.idx <- which(!is.cell); cell.idx <- which(is.cell)
  nG <- length(gene.idx); nC <- length(cell.idx)

  # gene rows (gene-gene and gene-cell blocks), unchanged by the decay; mirrored into cell -> gene
  G <- as(Matrix::tcrossprod(inc[gene.idx, , drop=FALSE], inc), "TsparseMatrix")
  gi <- gene.idx[G@i + 1]; gj <- G@j + 1; gx <- G@x
  off <- gi != gj
  gi <- gi[off]; gj <- gj[off]; gx <- gx[off]
  to.cell <- gj > nG
  ii <- c(gi, gj[to.cell]); jj <- c(gj, gi[to.cell]); xx <- c(gx, gx[to.cell])
  rm(G)

  # cell-cell block: pairs within neighbouring tiles only
  cells <- node.ids[cell.idx]
  missing.cells <- setdiff(cells, rownames(coords))
  if(length(missing.cells) > 0){
    stop(length(missing.cells), " cell(s) have no row in coords, e.g. '", missing.cells[1], "'")
  }
  xy <- coords[cells, , drop=FALSE]
  cx1 <- as.numeric(xy[, 1]); cy1 <- as.numeric(xy[, 2])
  if(anyNA(cx1) || anyNA(cy1)){
    stop("coords has NA for some cells")
  }
  BT <- Matrix::t(inc[cell.idx, , drop=FALSE])          # hyperedges x cells: a tile is a column subset
  side <- radius * (1 + 1e-9)                           # a hair wider than radius, against rounding at tile edges
  tx <- floor(cx1 / side); ty <- floor(cy1 / side)
  tiles <- split(seq_len(nC), paste(tx, ty))
  txy <- do.call(rbind, lapply(strsplit(names(tiles), " ", fixed=TRUE), as.numeric))
  ci <- vector("list", length(tiles)); cj <- ci; cv <- ci
  for(k in seq_along(tiles)){
    cT <- tiles[[k]]
    nb <- unlist(lapply(-1:1, function(dx) lapply(-1:1, function(dy) tiles[[paste(txy[k, 1] + dx, txy[k, 2] + dy)]])),
                 use.names=FALSE)
    nb <- sort(nb)
    P <- as(Matrix::crossprod(BT[, cT, drop=FALSE], BT[, nb, drop=FALSE]), "TsparseMatrix")
    i <- cT[P@i + 1]; j <- nb[P@j + 1]; x <- P@x
    ok <- (x > 0) & (i != j)
    i <- i[ok]; j <- j[ok]; x <- x[ok]
    # as score_spatial_decay() computes them
    d <- sqrt((cx1[i] - cx1[j])^2 + (cy1[i] - cy1[j])^2)
    w <- spatial_weights(d, mode=mode, radius=radius, nugget=nugget, target=target, floor=floor)
    ok <- w != 0
    ci[[k]] <- cell.idx[i[ok]]; cj[[k]] <- cell.idx[j[ok]]; cv[[k]] <- x[ok] * w[ok]
  }
  ci <- unlist(ci); cj <- unlist(cj); cv <- unlist(cv)

  # dimnames named "node", as the xtabs-based matrix of generate_adjacency_matrix()
  adj <- Matrix::sparseMatrix(i=c(ii, ci), j=c(jj, cj), x=c(xx, cv), dims=c(L, L),
                              dimnames=list(node=node.ids, node=node.ids))
  if(verbose){
    message(sprintf("<generate_spatial_adjacency_matrix> %d nodes (%d genes, %d cells), %s cell pairs within %g, %d tiles, %.1f s",
                    L, nG, nC, format(length(ci), big.mark=","), radius, length(tiles),
                    as.numeric(difftime(Sys.time(), t0, units="secs"))))
  }
  adj
}
