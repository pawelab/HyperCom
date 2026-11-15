// adjacency.cpp: adjacency matrix of a hypergraph, the same dsCMatrix as generate_adjacency_matrix of HyperCom 0.0.0.9000
//
// That version counts node x edge incidences with xtabs, takes tcrossprod (CHOLMOD builds the whole product, then its upper
// triangle), sets the diagonal to 0 (a copy) and drops nodes without neighbours (another copy). Here each column j of the
// upper triangle is accumulated directly from the edges of node j (sum over shared edges of the product of incidence counts,
// exact integers in any order), nodes without neighbours are left out from the start, and the diagonal is written as the
// explicit zeros that version keeps. Two passes (count, then fill) so the output is allocated once at its final size.

#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <vector>
#include "r_alloc.h"

const size_t SCAN_RATIO = 16;          // column j: its touched rows are sorted when fewer than j / 16, else rows 0..j-1 are scanned

// [[Rcpp::export]]
Rcpp::List adjacency_cpp(Rcpp::IntegerVector node, Rcpp::IntegerVector edge, int n_nodes, int n_edges){
  const size_t R = node.size(), N = n_nodes, E = n_edges;
  if((size_t) edge.size() != R) Rcpp::stop("adjacency_cpp: one edge per node entry");
  for(size_t k = 0; k < R; ++k){
    if(node[k] < 0 || node[k] >= n_nodes || edge[k] < 0 || edge[k] >= n_edges) Rcpp::stop("adjacency_cpp: index out of range");
  }

  // incidence by edge: nodes in ascending order with their counts
  std::vector<size_t> eptr(E + 1, 0);
  for(size_t k = 0; k < R; ++k) ++eptr[edge[k] + 1];
  for(size_t e = 0; e < E; ++e) eptr[e + 1] += eptr[e];
  std::vector<int> enode(R);
  std::vector<size_t> fill(eptr.begin(), eptr.end() - 1);
  for(size_t k = 0; k < R; ++k) enode[fill[edge[k]]++] = node[k];
  std::vector<size_t> cptr(E + 1, 0);
  std::vector<int> cnode;
  std::vector<double> ccount;
  cnode.reserve(R);
  ccount.reserve(R);
  for(size_t e = 0; e < E; ++e){
    std::sort(enode.begin() + eptr[e], enode.begin() + eptr[e + 1]);
    for(size_t k = eptr[e]; k < eptr[e + 1]; ++k){
      if(k > eptr[e] && enode[k] == enode[k - 1]){
        ccount.back() += 1;
      } else {
        cnode.push_back(enode[k]);
        ccount.push_back(1);
      }
    }
    cptr[e + 1] = cnode.size();
  }
  std::vector<int>().swap(enode);

  // incidence by node: edges in ascending order with their counts; a node is kept when one of its edges has another node
  std::vector<size_t> nptr(N + 1, 0);
  for(size_t k = 0; k < cnode.size(); ++k) ++nptr[cnode[k] + 1];
  for(size_t j = 0; j < N; ++j) nptr[j + 1] += nptr[j];
  std::vector<int> nedge(cnode.size());
  std::vector<double> ncount(cnode.size());
  fill.assign(nptr.begin(), nptr.end() - 1);
  std::vector<int> keep(N, 0);
  for(size_t e = 0; e < E; ++e){
    const bool shared = cptr[e + 1] - cptr[e] >= 2;
    for(size_t k = cptr[e]; k < cptr[e + 1]; ++k){
      const size_t j = cnode[k];
      nedge[fill[j]] = e;
      ncount[fill[j]++] = ccount[k];
      if(shared) keep[j] = 1;
    }
  }
  std::vector<int> new_idx(N, -1);
  int kept = 0;
  for(size_t j = 0; j < N; ++j) if(keep[j]) new_idx[j] = kept++;

  // column j: rows i < j sharing an edge with j, then the diagonal
  std::vector<size_t> mark(N, 0);
  std::vector<double> acc(N, 0);
  std::vector<int> touched;
  Rcpp::IntegerVector p = alloc_vector<INTSXP>(kept + 1);
  size_t nnz = 0;
  for(size_t j = 0; j < N; ++j){
    if(!keep[j]) continue;
    size_t m = 0;
    for(size_t k = nptr[j]; k < nptr[j + 1]; ++k){
      const size_t e = nedge[k];
      for(size_t q = cptr[e]; q < cptr[e + 1] && (size_t) cnode[q] < j; ++q){
        const size_t i = cnode[q];
        if(mark[i] != j + 1){
          mark[i] = j + 1;
          ++m;
        }
      }
    }
    nnz += m + 1;
    if(nnz > (size_t) INT_MAX) Rcpp::stop("adjacency_cpp: too many node pairs for a sparse matrix");
    p[new_idx[j] + 1] = (int) nnz;
  }

  Rcpp::IntegerVector out_i = alloc_vector<INTSXP>(nnz);
  Rcpp::NumericVector out_x = alloc_vector<REALSXP>(nnz);
  std::fill(mark.begin(), mark.end(), 0);
  size_t pos = 0;
  for(size_t j = 0; j < N; ++j){
    if(!keep[j]) continue;
    touched.clear();
    for(size_t k = nptr[j]; k < nptr[j + 1]; ++k){
      const size_t e = nedge[k];
      const double cj = ncount[k];
      for(size_t q = cptr[e]; q < cptr[e + 1] && (size_t) cnode[q] < j; ++q){
        const size_t i = cnode[q];
        if(mark[i] != j + 1){
          mark[i] = j + 1;
          acc[i] = 0;
          touched.push_back((int) i);
        }
        acc[i] += cj * ccount[q];
      }
    }
    if(touched.size() * SCAN_RATIO < j){
      std::sort(touched.begin(), touched.end());
    } else {
      touched.clear();
      for(size_t i = 0; i < j; ++i) if(mark[i] == j + 1) touched.push_back((int) i);
    }
    for(int i : touched){
      out_i[pos] = new_idx[i];
      out_x[pos++] = acc[i];
    }
    out_i[pos] = new_idx[j];
    out_x[pos++] = 0;
  }

  return Rcpp::List::create(Rcpp::Named("i") = out_i, Rcpp::Named("p") = p, Rcpp::Named("x") = out_x,
                            Rcpp::Named("keep") = Rcpp::LogicalVector(keep.begin(), keep.end()));
}
