// rank_grid.cpp -- sender / receiver scores of one cell x cell grid, without materialising it.
//
// For one (ligand l, receptor subunit r) the old R code built the dense n x n grid
//
//   Q[i, j] = a[i] * b[j] / C[i, j] * 1[sample_i == sample_j]      (a = F[cells, l], b = F[r, cells])
//   Q[i, j] = NA   if Q == 0, Q is NaN, or i == j
//
// percent-ranked all non-NA entries globally (min rank on ties, (rank - 1) / (K - 1)), and took
// row means (sender) and column means (receiver) of the ranks. That needs >= 12 n^2 bytes per grid
// (the grid, the sort keys, the integer ranks, the double ranks): ~4.8 GB at 20k cells, per worker.
//
// This kernel returns the same sender / receiver vectors with memory that does not grow with n^2:
//
//   pass 1   scan the grid once (computing Q on the fly from C, a, b) and histogram the values by
//            the top HIST_BITS bits of an order-preserving integer key. Count K, per-row and
//            per-column numbers of non-NA entries.
//   chunks   group consecutive histogram bins so that each group fits in `budget_bytes`. Because
//            bins are value ranges, equal values always land in the same bin, so ties never span a
//            chunk and the rank of everything below a chunk is known from the histogram alone.
//   pass 2   for each chunk, rescan the grid, scatter the entries that fall in the chunk into their
//            bins (exact offsets from the per-block histograms, so no atomics and a deterministic
//            layout), sort each bin, and add (rank - 1) to the row and column of every entry.
//
// Rank sums are accumulated as 64-bit integers, so they are exact and independent of the thread
// count or the order entries are visited in. The final means are computed in long double:
//   sender_i = sum_j (rank_ij - 1) / ((K - 1) * count_i)
// The old code rounded every (rank - 1) / (K - 1) to double before summing, so values agree with
// it to ~1e-16 relative, not bit for bit.
//
// C is read in place from the full influence matrix through `cell_idx`, so the n x n cell block is
// never copied. Do not compile with -ffast-math: the order of a * b / c must match the R code.
//
// Threads: the grid's columns are split into fixed blocks processed by `nthreads` OpenMP threads.
// Inside forked R workers, call with nthreads = 1 (GNU OpenMP is not fork-safe).

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>
#include <cstdio>
#ifndef _WIN32
#include <dlfcn.h>
#include <fcntl.h>
#include <unistd.h>
#endif

using namespace Rcpp;

namespace {

const int HIST_BITS = 18;                      // 2^18 bins: 64 per octave of |Q|
const uint64_t N_BINS = (uint64_t) 1 << HIST_BITS;

struct Entry {
  uint64_t key;                                // order-preserving image of the double
  uint32_t i, j;                               // row, column in the cell grid
};

// IEEE-754 doubles ordered as unsigned integers: flip all bits of negatives, set the sign bit of
// positives. Equal doubles give equal keys (0 and NaN never get here).
inline uint64_t sortable_key(double v) {
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | ((uint64_t) 1 << 63));
}

inline uint64_t bin_of(uint64_t key) { return key >> (64 - HIST_BITS); }

inline bool key_less(const Entry &x, const Entry &y) { return x.key < y.key; }

// Sort one bin by key. The values of a dataset occupy only a few thousand of the 2^HIST_BITS bins,
// so a bin can hold 10^4 - 10^6 entries; a counting sort on the next key bits (enough of them for
// ~1 entry per sub-bucket) leaves std::sort only tiny groups. `tmp` and `count` are scratch space.
void sort_bin(Entry *e, uint64_t m, std::vector<Entry> &tmp, std::vector<uint64_t> &count) {
  if (m < 256) {
    std::sort(e, e + m, key_less);
    return;
  }
  int sub_bits = 0;
  while (((uint64_t) 1 << sub_bits) < m && sub_bits < 20) ++sub_bits;
  const int shift = 64 - HIST_BITS - sub_bits;
  const uint64_t n_sub = (uint64_t) 1 << sub_bits, mask = n_sub - 1;
  if (tmp.size() < m) tmp.resize(m);
  count.assign(n_sub + 1, 0);
  for (uint64_t p = 0; p < m; ++p) ++count[((e[p].key >> shift) & mask) + 1];
  for (uint64_t s = 0; s < n_sub; ++s) count[s + 1] += count[s];
  for (uint64_t p = 0; p < m; ++p) tmp[count[(e[p].key >> shift) & mask]++] = e[p];
  // count[s] is now the end of sub-bucket s
  uint64_t start = 0;
  for (uint64_t s = 0; s < n_sub; ++s) {
    const uint64_t end = count[s];
    if (end - start > 1) std::sort(tmp.data() + start, tmp.data() + end, key_less);
    start = end;
  }
  std::memcpy(e, tmp.data(), m * sizeof(Entry));
}

struct Grid {
  const double *F;                             // full L x L influence matrix
  R_xlen_t L;
  const int *idx;                              // 0-based positions of the cells in F
  R_xlen_t n;
  const double *a, *b;                         // ligand column, receptor row (both over cells), or null
  const int *sid;                              // sample id per cell, or null (one sample)

  // the value of Q[i, j], or false if it is NA (zero, NaN or diagonal)
  inline bool value(R_xlen_t i, R_xlen_t j, double &v) const {
    if (i == j) return false;
    const double c = F[(R_xlen_t) idx[i] + (R_xlen_t) idx[j] * L];
    v = a ? (a[i] * b[j]) / c : c;
    if (sid) v = v * (double) (sid[i] == sid[j]);
    return !(v == 0 || ISNAN(v));
  }
};

int resolve_threads(int nthreads) {
#ifdef _OPENMP
  return nthreads < 1 ? 1 : nthreads;
#else
  (void) nthreads;
  return 1;
#endif
}

}  // namespace

// [[Rcpp::export]]
List rank_grid_cpp(NumericMatrix F, IntegerVector cell_idx,
                   Nullable<NumericVector> a = R_NilValue,
                   Nullable<NumericVector> b = R_NilValue,
                   Nullable<IntegerVector> sample_id = R_NilValue,
                   int nthreads = 1, double budget_bytes = 1e9) {
  const R_xlen_t L = F.nrow();
  if (F.ncol() != L) stop("rank_grid_cpp: F must be square");
  const R_xlen_t n = cell_idx.size();
  if (n >= 4294967295.0) stop("rank_grid_cpp: too many cells");
  for (R_xlen_t k = 0; k < n; ++k) {
    if (cell_idx[k] < 0 || cell_idx[k] >= L) stop("rank_grid_cpp: cell_idx out of range (must be 0-based)");
  }
  if (a.isNotNull() != b.isNotNull()) stop("rank_grid_cpp: give both a and b, or neither");

  NumericVector av, bv;
  IntegerVector sv;
  Grid g{REAL(F), L, INTEGER(cell_idx), n, nullptr, nullptr, nullptr};
  if (a.isNotNull()) {
    av = as<NumericVector>(a.get());
    bv = as<NumericVector>(b.get());
    if (av.size() != n || bv.size() != n) stop("rank_grid_cpp: a and b must have one value per cell");
    g.a = REAL(av);
    g.b = REAL(bv);
  }
  if (sample_id.isNotNull()) {
    sv = as<IntegerVector>(sample_id.get());
    if (sv.size() != n) stop("rank_grid_cpp: sample_id must have one value per cell");
    g.sid = INTEGER(sv);
  }

  const int nt = resolve_threads(nthreads);
  // fixed column blocks, independent of how many threads OpenMP actually provides
  const int n_blocks = (int) std::min<R_xlen_t>(std::max(nt, 1), std::max<R_xlen_t>(n, 1));
  std::vector<R_xlen_t> col_start(n_blocks + 1);
  for (int k = 0; k <= n_blocks; ++k) col_start[k] = (R_xlen_t) ((double) n * k / n_blocks);

  // ---- pass 1: histograms and counts ---------------------------------------------------------
  std::vector<std::vector<uint64_t>> hist(n_blocks, std::vector<uint64_t>(N_BINS, 0));
  std::vector<std::vector<int64_t>> row_cnt_blk(n_blocks, std::vector<int64_t>(n, 0));
  std::vector<int64_t> col_cnt(n, 0);

  #pragma omp parallel for num_threads(nt) if(nt > 1) schedule(dynamic, 1)
  for (int blk = 0; blk < n_blocks; ++blk) {
    uint64_t *h = hist[blk].data();
    int64_t *rc = row_cnt_blk[blk].data();
    for (R_xlen_t j = col_start[blk]; j < col_start[blk + 1]; ++j) {
      int64_t cc = 0;
      for (R_xlen_t i = 0; i < n; ++i) {
        double v;
        if (!g.value(i, j, v)) continue;
        ++h[bin_of(sortable_key(v))];
        ++rc[i];
        ++cc;
      }
      col_cnt[j] = cc;
    }
  }

  std::vector<int64_t> row_cnt(n, 0);
  for (int blk = 0; blk < n_blocks; ++blk) {
    for (R_xlen_t i = 0; i < n; ++i) row_cnt[i] += row_cnt_blk[blk][i];
    std::vector<int64_t>().swap(row_cnt_blk[blk]);
  }

  std::vector<uint64_t> total(N_BINS, 0), below(N_BINS + 1, 0);
  for (uint64_t bin = 0; bin < N_BINS; ++bin) {
    for (int blk = 0; blk < n_blocks; ++blk) total[bin] += hist[blk][bin];
    below[bin + 1] = below[bin] + total[bin];
  }
  const uint64_t K = below[N_BINS];

  // ---- chunks of consecutive bins, each within the memory budget -----------------------------
  uint64_t max_entries = (uint64_t) std::max(1.0, budget_bytes / (double) sizeof(Entry));
  std::vector<uint64_t> chunk_lo, chunk_hi;    // bin ranges [lo, hi)
  uint64_t largest = 0;
  for (uint64_t bin = 0; bin < N_BINS;) {
    while (bin < N_BINS && total[bin] == 0) ++bin;
    if (bin == N_BINS) break;
    const uint64_t lo = bin;
    uint64_t size = 0;
    // a single bin larger than the budget becomes its own (oversized) chunk
    while (bin < N_BINS && (size == 0 || size + total[bin] <= max_entries)) size += total[bin++];
    chunk_lo.push_back(lo);
    chunk_hi.push_back(bin);
    largest = std::max(largest, size);
  }

  std::unique_ptr<Entry[]> buf(largest > 0 ? new Entry[largest] : nullptr);

  // per-thread exact rank sums
  std::vector<std::vector<int64_t>> row_sum_t(nt, std::vector<int64_t>(n, 0));
  std::vector<std::vector<int64_t>> col_sum_t(nt, std::vector<int64_t>(n, 0));

  // ---- pass 2: one chunk at a time ---------------------------------------------------------
  for (size_t c = 0; c < chunk_lo.size(); ++c) {
    const uint64_t lo = chunk_lo[c], hi = chunk_hi[c], nb = hi - lo;
    const uint64_t base = below[lo];           // entries below this chunk

    // write position of (block, bin): bins in order, blocks in order within a bin
    #pragma omp parallel for num_threads(nt) if(nt > 1) schedule(dynamic, 1)
    for (int blk = 0; blk < n_blocks; ++blk) {
      std::vector<uint64_t> cursor(nb);
      for (uint64_t k = 0; k < nb; ++k) {
        uint64_t pos = below[lo + k] - base;
        for (int other = 0; other < blk; ++other) pos += hist[other][lo + k];
        cursor[k] = pos;
      }
      for (R_xlen_t j = col_start[blk]; j < col_start[blk + 1]; ++j) {
        for (R_xlen_t i = 0; i < n; ++i) {
          double v;
          if (!g.value(i, j, v)) continue;
          const uint64_t key = sortable_key(v), bin = bin_of(key);
          if (bin < lo || bin >= hi) continue;
          buf[cursor[bin - lo]++] = Entry{key, (uint32_t) i, (uint32_t) j};
        }
      }
    }

    // sort each bin and add the min rank (minus 1) of every entry to its row and column
    #pragma omp parallel num_threads(nt) if(nt > 1)
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      int64_t *rs = row_sum_t[tid].data();
      int64_t *cs = col_sum_t[tid].data();
      std::vector<Entry> tmp;
      std::vector<uint64_t> count;
      #pragma omp for schedule(dynamic, 1)
      for (uint64_t k = 0; k < nb; ++k) {
        const uint64_t bin = lo + k, cnt = total[bin];
        if (cnt == 0) continue;
        Entry *e = buf.get() + (below[bin] - base);
        sort_bin(e, cnt, tmp, count);
        uint64_t first = 0;
        for (uint64_t p = 0; p < cnt; ++p) {
          if (p > 0 && e[p].key != e[p - 1].key) first = p;
          const int64_t r = (int64_t) (below[bin] + first);     // min rank - 1
          rs[e[p].i] += r;
          cs[e[p].j] += r;
        }
      }
    }
  }

  // ---- means -------------------------------------------------------------------------------
  NumericVector sender(n), receiver(n);
  const long double km1 = (long double) K - 1.0L;
  for (R_xlen_t i = 0; i < n; ++i) {
    int64_t rs = 0, cs = 0;
    for (int t = 0; t < nt; ++t) { rs += row_sum_t[t][i]; cs += col_sum_t[t][i]; }
    sender[i]   = row_cnt[i] > 0 ? (double) ((long double) rs / (km1 * (long double) row_cnt[i])) : R_NaN;
    receiver[i] = col_cnt[i] > 0 ? (double) ((long double) cs / (km1 * (long double) col_cnt[i])) : R_NaN;
  }

  return List::create(_["sender"] = sender, _["receiver"] = receiver,
                      _["K"] = (double) K, _["chunks"] = (int) chunk_lo.size(),
                      _["chunk_bytes"] = (double) largest * sizeof(Entry));
}

// Thread count of the BLAS R is linked to, set at run time. R on the SCC links FlexiBLAS, whose
// backend (OpenBLAS, MKL, ...) reads OMP_NUM_THREADS once at start-up; this lets generate_infMat()
// use all cores for the inversion even when OMP_NUM_THREADS = 1 was exported for forked workers.
// Returns the previous thread count, or -1 if no known BLAS control function is found.

// [[Rcpp::export]]
int blas_threads_cpp(int n = -1) {
#ifndef _WIN32
  typedef int (*getter_t)(void);
  typedef void (*setter_t)(int);
  const char *get_names[] = {"flexiblas_get_num_threads", "openblas_get_num_threads", "MKL_Get_Max_Threads"};
  const char *set_names[] = {"flexiblas_set_num_threads", "openblas_set_num_threads", "MKL_Set_Num_Threads"};
  for (int k = 0; k < 3; ++k) {
    void *gs = dlsym(RTLD_DEFAULT, get_names[k]);
    void *ss = dlsym(RTLD_DEFAULT, set_names[k]);
    if (gs == nullptr || ss == nullptr) continue;
    const int old = ((getter_t) gs)();
    if (n > 0) {
      // FlexiBLAS 3.3 prints a stray "BLAS" line on stdout here; keep it out of the logs
      std::fflush(stdout);
      const int saved = dup(1), devnull = open("/dev/null", O_WRONLY);
      if (saved >= 0 && devnull >= 0) dup2(devnull, 1);
      ((setter_t) ss)(n);
      std::fflush(stdout);
      if (saved >= 0 && devnull >= 0) dup2(saved, 1);
      if (devnull >= 0) close(devnull);
      if (saved >= 0) close(saved);
    }
    return old;
  }
#endif
  (void) n;
  return -1;
}
