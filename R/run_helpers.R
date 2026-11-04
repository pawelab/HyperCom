# FLAG:NEW_20261104 (v1-fast)
# ----------------------------------------------------------------------------------------------
# What a long run needs besides the computation: a memory budget for the LR scoring chosen from
# the machine, and checkpoints so that a run stopped part way resumes instead of starting over.
# ----------------------------------------------------------------------------------------------

#' Memory available to this R session, in GB
#'
#' `available_memory_gb` returns the smaller of the memory the operating system reports as available
#' (`MemAvailable` in /proc/meminfo) and what is left under a cgroup memory limit, if one is set
#' (containers, some batch schedulers). NA where neither can be read (e.g. not Linux).
#' @returns A number of GB, or NA
#' @export
available_memory_gb <- function(){
  read_num <- function(f){
    if(!file.exists(f)) return(NA_real_)
    v <- suppressWarnings(as.numeric(readLines(f, n = 1, warn = FALSE)))
    if(length(v) == 0) NA_real_ else v
  }
  out <- NA_real_
  if(file.exists("/proc/meminfo")){
    l <- grep("^MemAvailable:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
    if(length(l) == 1) out <- as.numeric(gsub("[^0-9]", "", l)) / 1024^2       # kB -> GB
  }
  # cgroup v2, then v1; "max" or a huge number means no limit
  lim <- read_num("/sys/fs/cgroup/memory.max"); use <- read_num("/sys/fs/cgroup/memory.current")
  if(is.na(lim)){
    lim <- read_num("/sys/fs/cgroup/memory/memory.limit_in_bytes"); use <- read_num("/sys/fs/cgroup/memory/memory.usage_in_bytes")
  }
  if(!is.na(lim) && !is.na(use) && lim < 2^60){
    out <- min(out, (lim - use) / 1024^3, na.rm = TRUE)
  }
  out
}

#' Memory budget per grid for the LR scoring
#'
#' `auto_grid_mem_gb` returns the budget `prioritize_lr(grid.mem.gb = "auto")` uses: the memory one
#' grid needs to be ranked in a single pass, 16 n^2 bytes for n cells, if that fits in 40% of the
#' memory still available divided by the number of workers; otherwise that 40% share (the kernel
#' then ranks the grid in several passes, which is slower but gives the same numbers). The 40%
#' leaves room for the rest of the session and for memory other processes take later. Pass a
#' number to `grid.mem.gb` instead to set the budget yourself, e.g. when a batch scheduler limits
#' the job to less memory than the machine reports as available.
#' @param n Number of cells
#' @param workers Grids scored at once (see `prioritize_lr`)
#' @param share Fraction of the available memory to use
#' @returns A number of GB
#' @export
auto_grid_mem_gb <- function(n, workers = 1L, share = 0.4){
  need <- 16 * as.numeric(n)^2 / 1024^3 * 1.001 + 0.01
  avail <- available_memory_gb()
  if(is.na(avail)){
    cap <- 8
    message(sprintf("<prioritize_lr> available memory unknown: grid.mem.gb = min(%.2f, %g); pass grid.mem.gb to choose", need, cap))
  } else {
    cap <- max(0.5, share * avail / max(1L, workers))
  }
  min(need, cap)
}

# file name of one grid's checkpoint: ligand and receptor, characters outside [A-Za-z0-9._-] replaced
grid_checkpoint_name <- function(ligand, receptor){
  clean <- function(x) gsub("[^A-Za-z0-9._-]", "_", as.character(x))
  paste0(clean(ligand), "__", clean(receptor), ".rds")
}

# what ties checkpoints to one influence matrix and one set of cells: the cell names, the
# per-cell sample ids and 64 entries of infMat spread over the matrix (any change to the data,
# the kernel or r changes them)
grid_checkpoint_fingerprint <- function(infMat, cell.info){
  L <- nrow(infMat)
  i <- unique(round(seq(1, L, length.out = 64)))
  list(format = 1L, L = L, cells = cell.info$mat.names, sample.id = cell.info$sample.id,
       probe = infMat[cbind(i, rev(i))])
}

# create the checkpoint folder, or check that the grids already in it come from the same data;
# if not, they are discarded (with a warning) and the folder starts empty
prepare_grid_checkpoints <- function(dir, infMat, cell.info){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fp <- grid_checkpoint_fingerprint(infMat, cell.info)
  f <- file.path(dir, "fingerprint.rds")
  if(file.exists(f)){
    old <- tryCatch(readRDS(f), error = function(e) NULL)
    if(!identical(old, fp)){
      stale <- list.files(dir, pattern = "\\.rds(\\.tmp[0-9]+)?$", full.names = TRUE)
      warning(sprintf("prioritize_lr: the %d checkpoint(s) in %s were made from different data and are discarded",
                      length(stale) - 1L, dir))
      unlink(stale)
    }
  }
  unlink(list.files(dir, pattern = "\\.tmp[0-9]+$", full.names = TRUE))     # half-written grids of a killed run
  saveRDS(fp, f)
  invisible(dir)
}

# remove the checkpoint files this package wrote, then the folder if it is empty
remove_grid_checkpoints <- function(dir){
  unlink(list.files(dir, pattern = "\\.rds(\\.tmp[0-9]+)?$", full.names = TRUE))
  if(length(list.files(dir, all.files = TRUE, no.. = TRUE)) == 0) unlink(dir, recursive = TRUE)
  invisible(NULL)
}
