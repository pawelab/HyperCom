// grid_sort.h: ranking one value bucket. Its key bits are sorted as packed (key bits, index) words by radix passes and an insertion
// sort, and each entry gets its min rank. One thread ranks a bucket, in its own buffers (Work). Included by score_grids.cpp.
#ifndef HYPERCOM_GRID_SORT_H
#define HYPERCOM_GRID_SORT_H

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>
#include "grid_team.h"

namespace {

// Insertion sort with x[0] made the minimum first (searched in x[0, e0), where it is), so the inner loop needs no bounds test
void insertion_sort_sentinel(uint64_t *x, size_t m, size_t e0){
  size_t k = 0;
  for(size_t i = 1; i < e0; ++i){
    if(x[i] < x[k]) k = i;
  }
  std::swap(x[0], x[k]);
  for(size_t i = 1; i < m; ++i){
    const uint64_t v = x[i];
    if(x[i - 1] <= v) continue;
    size_t j = i;
    do {
      x[j] = x[j - 1];
      --j;
    } while(x[j - 1] > v);
    x[j] = v;
  }
}

// After a radix pass into y (cnt[coff + b] = end of sub-bucket b of nb, the biggest holds most): deeper passes for crowded
// sub-buckets, then one insertion sort over the nearly sorted result
void sort_rest(uint64_t *x, uint64_t *y, size_t m, int P, size_t nb, uint32_t most, size_t coff, std::vector<uint32_t> &cnt);

// One radix pass of the m words of x into y by digit(word) in [0, nb), with the nb + 1 counts c: afterwards c[b] is the end of
// sub-bucket b; returns the size of the biggest
template <class Digit> uint32_t radix_pass(const uint64_t *x, uint64_t *y, size_t m, size_t nb, uint32_t *c, Digit digit){
  std::fill(c, c + nb + 1, 0u);
  for(size_t p = 0; p < m; ++p) ++c[digit(x[p]) + 1];
  uint32_t most = 0;
  for(size_t b = 0; b < nb; ++b){
    most = std::max(most, c[b + 1]);
    c[b + 1] += c[b];
  }
  for(size_t p = 0; p < m; ++p) y[c[digit(x[p])]++] = x[p];
  return most;
}

// Sort packed words (key bits << P | index) of x into y, x is scratch: one radix pass over the key range of the words, then sort_rest
void sort_words(uint64_t *x, uint64_t *y, size_t m, int P, std::vector<uint32_t> &cnt){
  uint64_t lo = x[0] >> P, hi = lo;
  for(size_t p = 1; p < m; ++p){
    const uint64_t k = x[p] >> P;
    lo = std::min(lo, k);
    hi = std::max(hi, k);
  }
  if(lo == hi){
    std::memcpy(y, x, m * sizeof(uint64_t));
    return;
  }
  const int bits = bit_width(hi - lo), want = std::min(RADIX_BITS, bit_width(m)), shift = bits > want ? bits - want : 0;
  const size_t nb = (size_t) ((hi - lo) >> shift) + 1, base = cnt.size();
  cnt.resize(base + nb + 1);
  const uint32_t most = radix_pass(x, y, m, nb, cnt.data() + base, [&](uint64_t v){ return ((v >> P) - lo) >> shift; });
  if(shift > 0) sort_rest(x, y, m, P, nb, most, base, cnt);
  cnt.resize(base);
}

void sort_rest(uint64_t *x, uint64_t *y, size_t m, int P, size_t nb, uint32_t most, size_t coff, std::vector<uint32_t> &cnt){
  if(most > INSERTION_MAX){
    for(size_t b = 0, s = 0; b < nb; ++b){
      const size_t e = cnt[coff + b];
      if(e - s > INSERTION_MAX){
        std::memcpy(x + s, y + s, (e - s) * sizeof(uint64_t));
        sort_words(x + s, y + s, e - s, P, cnt);
      }
      s = e;
    }
  }
  size_t b = 0;
  while(cnt[coff + b] == 0) ++b;
  insertion_sort_sentinel(y, m, cnt[coff + b]);
}

// Byte offset in the pool of the rank of slab entry p: a slab's entries hold key bits until ranked, then the first half of the
// slab holds their ranks as uint32
inline size_t slab_rank_offset(size_t p){
  return (p - p % SLAB) * sizeof(uint64_t) + p % SLAB * sizeof(uint32_t);
}

inline uint32_t get_rank(const uint64_t *pool, size_t p){
  uint32_t rank;
  std::memcpy(&rank, (const char *) pool + slab_rank_offset(p), sizeof rank);
  return rank;
}

// Keys of a bucket held in one run of the pool; their ranks go to the run's first half as uint32
struct RunKeys {
  uint64_t *keys;

  void copy(uint64_t *x, size_t m) const {
    std::memcpy(x, keys, m * sizeof(uint64_t));
  }
  uint64_t hist(uint32_t *c, size_t m, int cut, int shift, uint64_t top) const {
    uint64_t seen = 0;
    for(size_t i = 0; i < m; ++i){
      seen |= keys[i];
      ++c[((keys[i] >> cut >> shift) & top) + 1];
    }
    return seen;
  }
  void scatter(uint64_t *y, uint32_t *c, size_t m, int cut, int I, int shift) const {
    for(size_t i = 0; i < m; ++i){
      const uint64_t key = keys[i] >> cut;
      y[c[key >> shift]++] = (key << I) | i;
    }
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy((char *) keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(){}
};

// A 6-byte pool slot holds the low 48 bits of a key, least significant byte first on every platform (f2::put48; a 6-byte memcpy of a
// uint64_t would keep the high bytes on a big-endian machine), read by one 8-byte load (f2::get48: pool6 has f2::SLACK bytes of
// slack), except the last slot of a bucket's run (SixKeys)
using f2::get48;
using f2::put48;

// Keys of a bucket held in one run of 6-byte slots (at most NARROW_BITS key bits); their ranks go to the run's first 4m bytes. Key i
// is one 8-byte load (get48), except the last (m >= 1): its load would reach 2 bytes into the next bucket's run, where another thread
// may be writing ranks, so it is read without crossing the run's end (get48_exact)
struct SixKeys {
  char *keys;

  uint64_t key(size_t i, size_t m) const {
    return i + 1 < m ? get48(keys + 6 * i) : f2::get48_exact(keys + 6 * i);
  }
  void copy(uint64_t *x, size_t m) const {
    for(size_t i = 0; i + 1 < m; ++i) x[i] = get48(keys + 6 * i);
    x[m - 1] = key(m - 1, m);
  }
  uint64_t hist(uint32_t *c, size_t m, int cut, int shift, uint64_t top) const {
    uint64_t seen = 0;
    auto add = [&](uint64_t k){
      seen |= k;
      ++c[((k >> cut >> shift) & top) + 1];
    };
    for(size_t i = 0; i + 1 < m; ++i) add(get48(keys + 6 * i));
    add(key(m - 1, m));
    return seen;
  }
  void scatter(uint64_t *y, uint32_t *c, size_t m, int cut, int I, int shift) const {
    auto place = [&](size_t i, uint64_t k){
      k >>= cut;
      y[c[k >> shift]++] = (k << I) | i;
    };
    for(size_t i = 0; i + 1 < m; ++i) place(i, get48(keys + 6 * i));
    place(m - 1, key(m - 1, m));
  }
  void put(size_t i, size_t rank){
    const uint32_t r = (uint32_t) rank;
    std::memcpy(keys + i * sizeof(uint32_t), &r, sizeof r);
  }
  void finish(){}
};

// Keys of a bucket held in slab segments (offset, length); their ranks go to z, then to each slab's first half
struct SlabKeys {
  uint64_t *pool;
  const std::vector<std::pair<size_t, size_t>> &segs;
  uint32_t *z;

  void copy(uint64_t *x, size_t) const {
    for(const auto &sg : segs){
      std::memcpy(x, pool + sg.first, sg.second * sizeof(uint64_t));
      x += sg.second;
    }
  }
  uint64_t hist(uint32_t *c, size_t, int cut, int shift, uint64_t top) const {
    uint64_t seen = 0;
    for(const auto &sg : segs){
      const uint64_t *src = pool + sg.first;
      for(size_t j = 0; j < sg.second; ++j){
        seen |= src[j];
        ++c[((src[j] >> cut >> shift) & top) + 1];
      }
    }
    return seen;
  }
  void scatter(uint64_t *y, uint32_t *c, size_t, int cut, int I, int shift) const {
    size_t i = 0;
    for(const auto &sg : segs){
      const uint64_t *src = pool + sg.first;
      for(size_t j = 0; j < sg.second; ++j, ++i){
        const uint64_t key = src[j] >> cut;
        y[c[key >> shift]++] = (key << I) | i;
      }
    }
  }
  void put(size_t i, size_t rank){
    z[i] = (uint32_t) rank;
  }
  void finish(){
    size_t i = 0;
    for(const auto &sg : segs){
      std::memcpy((char *) pool + slab_rank_offset(sg.first), z + i, sg.second * sizeof(uint32_t));
      i += sg.second;
    }
  }
};

// A bucket's keys hold at most its t key bits, by construction. A bit above them (seen: the keys or-ed together) can only come from a
// memory error; the key would index past the counts of the sort, so hist masks the index and this stops the run before the scatter
inline void check_key_width(uint64_t seen, int t){
  if(t < 64 && (seen >> t) != 0) throw std::runtime_error("a key is wider than its value bucket (a memory error? check the computer's RAM)");
}

// The min rank - 1 of each of the n sorted words y (key bits << I | index) to keys: first plus the position of the first word with
// the same key bits
template <class K> void put_ranks(K &keys, const uint64_t *y, size_t n, int I, uint64_t imask, size_t first){
  uint64_t prev = y[0] >> I, at = first;
  for(size_t v = 0; v < n; ++v){
    if((y[v] >> I) != prev){
      prev = y[v] >> I;
      at = first + v;
    }
    keys.put(y[v] & imask, at);
  }
}

// Ranks of one bucket of m entries, held in position order as their low t key bits: sort (key bits, index) words, then
// give each entry its min rank - 1: s0 plus its rank in the bucket (s0 = the bucket's start, or 0 for local ranks). A bucket of
// more than 2^(PIECE_BITS + 2) entries is cut by a top digit into pieces (equal keys share one), then each piece is sorted by a
// second digit, deeper passes and an insertion sort, and ranked, in one reused buffer; a smaller one gets one radix pass of up
// to TOP_BITS straight from its keys, then sort_rest.
template <class K> void rank_bucket(K &keys, size_t m, size_t s0, int t, Work &w){
  const int I = bit_width(m - 1), cut = t + I > 64 ? t + I - 64 : 0, width = t - cut;
  const uint64_t imask = I ? (~(uint64_t) 0 >> (64 - I)) : 0;
  const bool pieces = cut == 0 && bit_width(m) > PIECE_BITS + 2;
  if(w.x.size() < m) w.x.resize(m);
  if(!pieces && w.y.size() < m) w.y.resize(m);
  uint64_t *x = w.x.data(), *y = w.y.data();
  if(pieces){
    const int b1 = std::min(width, bit_width(m) - PIECE_BITS), s1 = width - b1;
    const size_t nb1 = (size_t) 1 << b1;
    w.bnd.assign(nb1 + 1, 0);
    uint32_t *c1 = w.bnd.data();
    check_key_width(keys.hist(c1, m, 0, s1, nb1 - 1), t);
    for(size_t b = 0; b < nb1; ++b) c1[b + 1] += c1[b];
    keys.scatter(x, c1, m, 0, I, s1);
    for(size_t b = 0, s = 0; b < nb1; ++b){
      const size_t e = c1[b], n = e - s;
      if(n == 0) continue;
      if(w.pb.size() < n) w.pb.resize(n);
      uint64_t *pc = w.pb.data(), *xs = x + s;
      if(n == 1 || s1 == 0){
        std::memcpy(pc, xs, n * sizeof(uint64_t));
      } else {
        // the second digit: the piece's next key bits; with none left below it (s2 = 0) the pass leaves the piece sorted
        const int b2 = std::min(s1, bit_width(n)), s2 = s1 - b2, sh = I + s2;
        const uint64_t mask2 = ((uint64_t) 1 << b2) - 1;
        const size_t nb2 = (size_t) mask2 + 1;
        w.cnt.resize(nb2 + 1);
        const uint32_t most = radix_pass(xs, pc, n, nb2, w.cnt.data(), [&](uint64_t v){ return (v >> sh) & mask2; });
        if(s2 > 0) sort_rest(xs, pc, n, I, nb2, most, 0, w.cnt);
      }
      put_ranks(keys, pc, n, I, imask, s0 + s);
      s = e;
    }
    keys.finish();
    return;
  }
  const int want = std::min(TOP_BITS, bit_width(m)), shift = width > want ? width - want : 0;
  const size_t nb = (size_t) 1 << (width - shift);
  w.cnt.assign(nb + 1, 0);
  uint32_t *c = w.cnt.data();
  check_key_width(keys.hist(c, m, cut, shift, nb - 1), t);
  uint32_t most = 0;
  for(size_t b = 0; b < nb; ++b){
    most = std::max(most, c[b + 1]);
    c[b + 1] += c[b];
  }
  keys.scatter(y, c, m, cut, I, shift);
  if(m > 1 && shift > 0) sort_rest(x, y, m, I, nb, most, 0, w.cnt);
  if(cut == 0){
    put_ranks(keys, y, m, I, imask, s0);
  } else {
    uint64_t *key = x;
    keys.copy(key, m);
    for(size_t q = 0; q < m; ){
      size_t u = q + 1;
      while(u < m && (y[u] >> I) == (y[q] >> I)) ++u;
      w.group.clear();
      for(size_t v = q; v < u; ++v) w.group.emplace_back(key[y[v] & imask], y[v] & imask);
      std::sort(w.group.begin(), w.group.end());
      for(size_t g = 0, h = 0; g < w.group.size(); ++g){
        if(w.group[g].first != w.group[h].first) h = g;
        keys.put(w.group[g].second, s0 + q + h);
      }
      q = u;
    }
  }
  keys.finish();
}

}  // namespace

#endif
