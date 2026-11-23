// f2.h: log-domain integer keys (memory mode), shared by influence.cpp (Λ = log F tiles from the factors), score_grids.cpp and
// screen_grids.cpp. A key orders ln(a b / F) as an integer made from integers only (fixlog52), so every computer gives the same bits:
// a grid entry's key is A[i] + B'[j] + Λ[i, j] (Λ = -ln F, rounded, stored as 6-byte offsets from a per-block base in tiles of the
// upper triangle; F is symmetric). Scores are within ~1e-11 of the exact kernel's (J. Moy 2026-11-20).
#ifndef HYPERCOM_F2_H
#define HYPERCOM_F2_H

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace f2 {

// 2^52 ln(x) of a positive finite double, rounded, from integers alone: the mantissa's log2 bit by bit (54 squarings), times ln 2
inline int64_t fixlog52(double x){
  uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  int e = (int) ((u >> 52) & 0x7FF);
  uint64_t frac = u & (((uint64_t) 1 << 52) - 1);
  if(e == 0){
    const int sh = __builtin_clzll(frac) - 11;
    frac = (frac << sh) & (((uint64_t) 1 << 52) - 1);
    e = 1 - sh;
  }
  uint64_t y = (((uint64_t) 1 << 52) | frac) << 10, f = 0;
  for(int k = 0; k < 54; ++k){
    y = (uint64_t) (((unsigned __int128) y * y) >> 62);
    f <<= 1;
    if(y >> 63){
      y >>= 1;
      f |= 1;
    }
  }
  const uint64_t LN2 = 0xB17217F7D1CF79ACULL;
  const __int128 total = (__int128) (e - 1023) * (__int128) LN2 + (__int128) (((unsigned __int128) f * LN2) >> 54);
  return (int64_t) ((total + ((__int128) 1 << 11)) >> 12);
}

// v / 2^sh rounded half up (sh >= 0)
inline int64_t round_shift(int64_t v, int sh){
  return sh > 0 ? (v + ((int64_t) 1 << (sh - 1))) >> sh : v;
}

// 6-byte little-endian slots
inline void put48(char *p, uint64_t v){
  unsigned char *u = reinterpret_cast<unsigned char *>(p);
  u[0] = (unsigned char) v;
  u[1] = (unsigned char) (v >> 8);
  u[2] = (unsigned char) (v >> 16);
  u[3] = (unsigned char) (v >> 24);
  u[4] = (unsigned char) (v >> 32);
  u[5] = (unsigned char) (v >> 40);
}

// One 8-byte load and a mask: every buffer of 6-byte slots read this way (Λ, the key pool) is made 8 bytes longer than its slots
inline uint64_t get48(const char *p){
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  const unsigned char *u = reinterpret_cast<const unsigned char *>(p);
  return (uint64_t) u[0] | (uint64_t) u[1] << 8 | (uint64_t) u[2] << 16 | (uint64_t) u[3] << 24 | (uint64_t) u[4] << 32 |
         (uint64_t) u[5] << 40;
#else
  uint64_t v;
  std::memcpy(&v, p, sizeof v);
  return v & (((uint64_t) 1 << 48) - 1);
#endif
}

// Slack bytes after the last 6-byte slot of a buffer read by get48
const size_t SLACK = 8;

// Λ tiles: TILE x TILE entries, 6 bytes each, tile (I, J) with I <= J of a block at the block's offset + TILE_BYTES (J (J + 1) / 2 + I),
// entry (a, b) of a tile at 6 (TILE b + a)
const size_t TILE = 8;
const size_t TILE_BYTES = 6 * TILE * TILE;

}  // namespace f2

#endif
