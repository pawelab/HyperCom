// score_grids.cpp: sender and receiver scores of cell x cell grids, within one ulp of those of HyperCom 0.0.0.9000
//
// For a ligand l and one receptor subunit r, HyperCom 0.0.0.9000 builds the cell x cell grid
//   Q[i, j] = (F[i, l] * F[r, j]) / F[i, j]      (F = influence matrix; Q[i, j] = F[i, j] for the network score)
// sets cross-sample entries, the diagonal, zeros and NaN to NA, percent-ranks the remaining K entries together
// ((min rank - 1) / (K - 1)) and averages the ranks over j (sender of i, colMeans) and over i (receiver of j, rowMeans).
//
// Here only the within-sample blocks are visited (every cross-sample entry is NA). Entries are ranked by an
// order-preserving integer image of their value: each entry falls in one of the grid's value buckets in key order, or in one
// of two catch-alls below and above the sampled key range. The bounds come from a sample of the entries: 2^b1 buckets of
// equal width (b1 = 11-15, growing with the grid; fewer when the sampled key range is narrow), or, for a grid above the
// cache budget scored by several threads, buckets of equal counts from a histogram of the sample, so every bucket sorts in
// L2. Each entry's key bits go to its bucket, then each bucket is sorted in cache as packed (key bits, index) words and its
// entries get their min ranks in place. Any partition in key order gives the same ranks. Each score is its exact integer rank
// sum divided once, sum / (count x (K - 1)), correctly rounded: the same on every computer, and within one ulp of the long double
// means of HyperCom 0.0.0.9000. A team of threads scores one grid at a time; several teams run at once when their grids fit the
// cache together, else one team takes all threads and splits the row ranges of a big block among them. Threads run as threads.h
// runs them (errors and interrupts stop the run at the next barrier). Memory per within-sample entry per team: its bucket id (2 B)
// and its key bits, 8 B, or 6 B for a grid above the cache budget on several threads (a bucket then keeps at most 48 key bits; the
// few entries of the catch-alls keep 8 B), plus each thread's sort buffers (16 B per entry of the biggest bucket it ranks).
//
// score_blocks_cpp is given every within-sample block of F (from the factored influence matrix, influence_blocks_cpp, or copied from a
// dense one, dense_blocks_cpp) and F[cell, ligand] / F[receptor, cell] for the genes of the pairs.
//
// Memory mode (score_lambda_cpp): no F block either. Entries are ranked by integer keys ln(a b / F) = A[i] + B'[j] + Λ[i, j] (f2.h),
// Λ = -ln F of the within-sample blocks as 6-byte offsets in tiles of the upper triangle (influence_lambda_cpp, about 3 B per entry),
// A, B' the logs of F[cell, ligand] and F[receptor, cell] / D, all made from integers alone; the network grid is ranked by
// ln D[j] - Λ[i, j], the order of F. Only the diagonal is skipped. Scores are within ~1e-11 of the exact kernel's.
//
// Passes (score_blocks_cpp, score_lambda_cpp): a big grid on several threads can place, rank and sum its buckets in P ranges, one
// after the other, recomputing each entry's bucket at each pass instead of storing it: no bucket ids, and a key pool of ~6.75 / P B per
// entry instead of 6 B. Same integer sums, so the same scores at any P.
//
// The stages are in headers included here in order: grid_team.h (limits, knobs, Cells, Buckets, Team, Work), grid_sort.h (ranking a
// bucket), grid_input.h (the entries: Q or keys), grid_buckets.h (sampling, the buckets), grid_place.h (counting, placing, ranking the
// buckets), grid_sums.h (rank sums, scores), grid_passes.h (passes > 1). This file runs a team (run_team) and the teams (run_grids).

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include "f2.h"
#include "grid_common.h"
#include "r_alloc.h"
#include "threads.h"
// The kernel's stages, in the order they build on each other; one translation unit, so that the hot loops inline across them
#include "grid_team.h"
#include "grid_sort.h"
#include "grid_input.h"
#include "grid_buckets.h"
#include "grid_place.h"
#include "grid_sums.h"
#include "grid_passes.h"

namespace {

// After each grid: frees per-thread sort buffers above KEEP_SORT_BUFFER entries
void trim_work(Work &w){
  if(w.x.size() > KEEP_SORT_BUFFER){
    std::vector<uint64_t>().swap(w.x);
    std::vector<uint64_t>().swap(w.y);
    std::vector<uint32_t>().swap(w.z);
  }
  if(w.pb.size() > KEEP_SORT_BUFFER) std::vector<uint64_t>().swap(w.pb);
  if(w.sb.size() > KEEP_SORT_BUFFER) std::vector<uint64_t>().swap(w.sb);
  if(w.group.capacity() > KEEP_SORT_BUFFER) std::vector<std::pair<uint64_t, uint64_t>>().swap(w.group);
}

// One team member: score grids until none are left; me = 0 leads
void run_team(const Cells &c, Team &tm, size_t me, const int *lig, const int *rec, size_t G, std::atomic<size_t> &next,
              const std::atomic<bool> &failed, double *snd, double *rcv){
  const size_t n = c.out.size(), nchunks = tm.chunks.size();
  Work w(c, tm);
  for(;;){
    tm.bar.sync(failed);
    if(me == 0){
      tm.g = failed.load() ? G : next.fetch_add(1);
      if(tm.g < G){
        tm.l = lig[tm.g];
        tm.r = rec[tm.g];
      }
      tm.next_count = 0;
      tm.next_split = 0;
      tm.next_place = 0;
      tm.next_bucket = 0;
      tm.next_sum = 0;
    }
    tm.bar.sync(failed);
    const size_t g = tm.g;
    if(g >= G) break;
    prepare_ab(c, tm, me);
    tm.bar.sync(failed);
    sample_keys(c, tm, me);
    tm.bar.sync(failed);
    uint64_t lo, hi;
    int range;
    const bool sampled = key_range(tm, lo, hi, range);
    const int tc = range > TABLE_BITS ? range - TABLE_BITS : 0;  // big grids: histogram cells of 2^tc keys from lo
    const size_t ncell = (size_t) ((hi - lo) >> tc) + 1;
    if(tm.big && sampled){
      sample_cells(c, tm, me, lo, tc, ncell);
      tm.bar.sync(failed);
    }
    if(me == 0) make_buckets(tm, c.cap, sampled, lo, range, tc, ncell);
    tm.bar.sync(failed);
    const Buckets &bk = tm.bk;
    const size_t NB = bk.nb + 2;
    if(w.start.size() < NB + 1){
      w.start.resize(NB + 1);
      w.at.resize(NB);
      w.rp.resize(NB);
    }
    w.start[0] = 0;
    if(tm.slabbed){
      if(tm.next_count.fetch_add(1) == 0) fill_slabs(c, tm, bk, w);
      tm.bar.sync(failed);
      for(size_t d = 0; d < NB; ++d) w.start[d + 1] = w.start[d] + tm.full[d] * SLAB + tm.slab_pos[d] % SLAB;
    } else {
      for(size_t ch; (ch = tm.next_count.fetch_add(1)) < nchunks; ){
        if(bk.map) count_chunk<true>(c, tm, bk, ch, w);
        else count_chunk<false>(c, tm, bk, ch, w);
      }
      tm.bar.sync(failed);
      split_buckets(tm, NB);
      tm.bar.sync(failed);
      for(size_t d = 0; d < NB; ++d) w.start[d + 1] = w.start[d] + tm.tot[d];
      if(tm.multi){
        run_passes(c, tm, me, bk, w, failed, snd + g * n, rcv + g * n);
        trim_work(w);
        continue;
      }
      if(tm.big){
        if(me == 0) six_storage(tm, bk, 0, NB);
        tm.bar.sync(failed);
      }
      for(size_t ch; (ch = tm.next_place.fetch_add(1)) < nchunks; ){
        if(!tm.big) place_chunk_ids(c, tm, bk, ch, w);
        else if(bk.map) place_chunk_six_nt<true>(c, tm, bk, ch, w);
        else place_chunk_six_nt<false>(c, tm, bk, ch, w);
      }
      tm.bar.sync(failed);
    }
    rank_buckets(tm, bk, w, 0, NB);
    tm.bar.sync(failed);
    const double denom = (double) w.start[NB] - 1.0;
    for(size_t ch; (ch = tm.next_sum.fetch_add(1)) < nchunks; ){
      if(tm.slabbed) sum_chunk<Store::SLABS>(c, tm, bk, ch, denom, w, snd + g * n, rcv + g * n);
      else sum_chunk<Store::RUNS>(c, tm, bk, ch, denom, w, snd + g * n, rcv + g * n);
    }
    if(tm.has_split){
      tm.bar.sync(failed);
      finish_columns(c, tm, me, denom, rcv + g * n);
    }
    trim_work(w);
  }
}


// The cells grouped into the within-sample blocks of grid_common.h (out: the cell at each position). The rank sums must fit uint64
// (RANK_SUM_MAX) and the ranks, counts and their products stay exact as doubles (GRID_MAX).
void group_cells(Cells &c, const int *sample, size_t n){
  c.block_start.push_back(0);
  for(const std::vector<int> &members : sample_blocks(sample, n)){
    for(int k : members) c.out.push_back(k);
    const size_t size = members.size();
    c.block_off.push_back(c.cap);
    c.block_start.push_back(c.out.size());
    c.cap += size * size;
    c.max_block = std::max(c.max_block, size);
  }
  if(c.cap >= GRID_MAX) Rcpp::stop("grid scoring: too many cells per sample (the sum of squared sample sizes must stay below 2^52)");
  if((double) c.max_block * (double) c.cap >= RANK_SUM_MAX){
    Rcpp::stop("grid scoring: a sample too big for 64-bit rank sums (its cells times the sum of squared sample sizes must stay below "
               "9.2e18)");
  }
}

// Sender and receiver scores (cell x pair) of the grids of ligands lp and receptor subunits rp (negative: the network score); a big grid
// on several threads in `passes` ranges of buckets when passes > 1
Rcpp::List run_grids(const Cells &c, const int *lp, const int *rp, size_t G, int threads, size_t passes){
  const size_t n = c.out.size();
  Rcpp::NumericMatrix sender = alloc_numeric_matrix(n, G), receiver = alloc_numeric_matrix(n, G);
  double *snd = REAL(sender), *rcv = REAL(receiver);

  const size_t cores = std::thread::hardware_concurrency();
  size_t nt = G == 0 || threads < 1 ? 1 : threads;
  if(cores > 0) nt = std::min(nt, cores);
  const size_t bytes = std::max<size_t>(1, grid_bytes(c));
  const size_t fit = std::max<size_t>(1, CACHE_BUDGET / bytes);
  const size_t nteams = std::max<size_t>(1, std::min(std::min(nt, fit), G));
  std::vector<std::unique_ptr<Team>> teams;
  std::vector<std::pair<Team *, size_t>> member;  // each thread's team and place in it
  try {
    for(size_t k = 0; k < nteams; ++k){
      teams.emplace_back(new Team(c, nt / nteams + (k < nt % nteams ? 1 : 0), nt == 1, passes));
      for(size_t me = 0; me < teams[k]->size; ++me) member.emplace_back(teams[k].get(), me);
    }
  } catch(std::bad_alloc &){
    Rcpp::stop("grid scoring: not enough memory for the grid buffers (" + std::to_string(nteams * bytes >> 20) + " MB)");
  }

  const unsigned mxcsr = get_mxcsr();
  std::atomic<size_t> next(0);
  hc::Stop stop;
  hc::run_threads(member.size(), stop, [&](size_t k){
    set_mxcsr(mxcsr);
    run_team(c, *member[k].first, member[k].second, lp, rp, G, next, stop.failed, snd, rcv);
  });
  stop.raise("grid scoring");

  return Rcpp::List::create(Rcpp::Named("sender") = sender, Rcpp::Named("receiver") = receiver);
}

// The checks of the arguments of score_blocks_cpp and score_lambda_cpp: one row of A and one column of B per cell, one receptor per
// ligand (both negative for the network grid), passes at least 1, sample ids positive
void check_grids(const Rcpp::NumericMatrix &A, const Rcpp::NumericMatrix &B, const Rcpp::IntegerVector &sample,
                 const Rcpp::IntegerVector &lig, const Rcpp::IntegerVector &rec, int passes){
  const size_t n = sample.size(), G = lig.size();
  if((size_t) A.nrow() != n || (size_t) B.ncol() != n) Rcpp::stop("grid scoring: A needs one row and B one column per cell");
  if((size_t) rec.size() != G) Rcpp::stop("grid scoring: one receptor per ligand");
  if(passes < 1) Rcpp::stop("grid scoring: passes must be at least 1");
  check_samples(sample, "grid scoring");
  for(size_t g = 0; g < G; ++g){
    if(lig[g] >= A.ncol() || rec[g] >= B.nrow() || (lig[g] < 0) != (rec[g] < 0)) Rcpp::stop("grid scoring: bad ligand or receptor index");
  }
}

// The cells of the grids by position (A's rows, B's columns) in the blocks of their samples, each block given whole (in Fb, or as
// Λ's tiles)
void block_cells(Cells &c, const Rcpp::NumericMatrix &A, const Rcpp::NumericMatrix &B, const Rcpp::IntegerVector &sample){
  c.A = REAL(A);
  c.B = REAL(B);
  c.nrec = B.nrow();
  group_cells(c, INTEGER(sample), sample.size());
}

}  // namespace

// The exact kernel on the within-sample blocks of F (influence_blocks_cpp or dense_blocks_cpp), A = F[cell, ligand genes] and
// B = F[receptor genes, cell]; lig and rec: each grid's column of A and row of B, both -1 for the network grid. passes > 1 ranks a big
// grid on several threads in that many ranges of buckets (the multi-pass path memory mode takes; the same scores at any passes, which
// the tests check)
// [[Rcpp::export]]
Rcpp::List score_blocks_cpp(Rcpp::NumericVector blocks, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample,
                            Rcpp::IntegerVector lig, Rcpp::IntegerVector rec, int threads, int passes = 1){
  check_grids(A, B, sample, lig, rec, passes);
  Cells c;
  block_cells(c, A, B, sample);
  if((size_t) blocks.size() != c.cap) Rcpp::stop("grid scoring: the blocks must hold every within-sample entry");
  c.Fb = REAL(blocks);

  return run_grids(c, INTEGER(lig), INTEGER(rec), lig.size(), threads, (size_t) passes);
}

// Memory mode: score_blocks_cpp ranking by the keys A + B' + Λ (f2.h). bytes, base, toff: Λ from influence_lambda_cpp for the same
// cells and samples at the same p; D: the adjacency's row sums of the cells, in the order of sample; A and B as score_blocks_cpp's,
// every entry positive and finite (their logs are taken)
// [[Rcpp::export]]
Rcpp::List score_lambda_cpp(Rcpp::RawVector bytes, Rcpp::NumericVector base, Rcpp::NumericVector toff, Rcpp::NumericVector D,
                            Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::IntegerVector sample, Rcpp::IntegerVector lig,
                            Rcpp::IntegerVector rec, int threads, int p, int passes){
  check_grids(A, B, sample, lig, rec, passes);
  const size_t n = sample.size();
  Cells c;
  block_cells(c, A, B, sample);
  const size_t nblocks = c.block_start.size() - 1;
  std::vector<size_t> sizes(nblocks);
  for(size_t blk = 0; blk < nblocks; ++blk) sizes[blk] = c.block_start[blk + 1] - c.block_start[blk];
  check_lambda(bytes, base, toff, D, A, B, sizes, n, p, "grid scoring");
  c.LamT = reinterpret_cast<const char *>(RAW(bytes));
  c.lbase.resize(nblocks);
  if(nblocks > 0) std::memcpy(c.lbase.data(), REAL(base), nblocks * sizeof(int64_t));  // int64 bits in doubles
  c.toff.resize(nblocks);
  for(size_t blk = 0; blk < nblocks; ++blk) c.toff[blk] = (size_t) toff[blk];
  c.p = p;
  c.lD52.resize(n);
  c.kD.resize(n);
  for(size_t k = 0; k < n; ++k){
    c.lD52[k] = f2::fixlog52(D[c.out[k]]);
    c.kD[k] = f2::round_shift(c.lD52[k], 52 - p);
  }

  return run_grids(c, INTEGER(lig), INTEGER(rec), lig.size(), threads, (size_t) passes);
}

// The kernel constants R's memory estimates and batch sizes use, so R reads them here instead of keeping copies: cache.budget (bytes
// of team buffers in flight at once), grid.entry.bytes (team buffer bytes per within-sample entry), pass.slack.div and pass.slack.max
// (passes > 1: a pass's key pool holds cap / passes entries plus cap / (passes x pass.slack.div), at most pass.slack.max, 6 B each),
// key.bits (memory mode: fractional bits of the keys at which R makes Λ)
// [[Rcpp::export]]
Rcpp::List kernel_constants_cpp(){
  return Rcpp::List::create(Rcpp::Named("cache.budget") = (double) CACHE_BUDGET, Rcpp::Named("grid.entry.bytes") = (double) GRID_ENTRY_BYTES,
                            Rcpp::Named("pass.slack.div") = (double) PASS_SLACK_DIV, Rcpp::Named("pass.slack.max") = (double) PASS_SLACK_MAX,
                            Rcpp::Named("key.bits") = f2::KEY_BITS);
}
