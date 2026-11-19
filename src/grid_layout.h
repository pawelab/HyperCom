// grid_layout.h: the within-sample blocks of the cell x cell grids, shared by the grid kernel (score_grids.cpp) and the blocks
// of the factored influence matrix made for it (influence.cpp)
#ifndef HYPERCOM_GRID_LAYOUT_H
#define HYPERCOM_GRID_LAYOUT_H

#include <algorithm>
#include <cstddef>
#include <vector>

// The cells of each within-sample block, by their position 0..n-1: samples ascending (ids >= 1, checked by the caller), empty
// samples skipped, the cells of a sample in the given order. The blocks F[cells, cells] follow one another column-major in
// this order.
inline std::vector<std::vector<int>> sample_blocks(const int *sample, size_t n){
  int n_samples = 0;
  for(size_t k = 0; k < n; ++k) n_samples = std::max(n_samples, sample[k]);
  std::vector<std::vector<int>> members(n_samples);
  for(size_t k = 0; k < n; ++k) members[sample[k] - 1].push_back((int) k);
  members.erase(std::remove_if(members.begin(), members.end(), [](const std::vector<int> &m){ return m.empty(); }), members.end());

  return members;
}

#endif
