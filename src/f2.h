// f2.h: log-domain integer keys (memory mode), shared by influence.cpp (Λ = log F tiles from the factors), score_grids.cpp and
// screen_grids.cpp. A key orders ln(a b / F) as an integer made from integers only (fixlog52), so every computer gives the same bits:
// a grid entry's key is A[i] + B'[j] + Λ[i, j] (Λ = -ln F, rounded, stored as 6-byte offsets from a per-block base in tiles of the
// upper triangle; F is symmetric). Scores are within ~1e-11 of the exact kernel's.
#ifndef HYPERCOM_F2_H
#define HYPERCOM_F2_H

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace f2 {

// 128-bit integers (a GCC and clang extension on 64-bit targets; __extension__ keeps -pedantic quiet)
#ifndef __SIZEOF_INT128__
#error "HyperCom needs 64-bit R (a compiler with __int128): install the 64-bit build of R 4.2 or later"
#endif
__extension__ typedef __int128 i128;
__extension__ typedef unsigned __int128 u128;

// Fractional bits of the keys at which R makes Λ (lambda_blocks): ln F to 2^-44, while a block's span still fits 48 bits
const int KEY_BITS = 44;

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
    y = (uint64_t) (((u128) y * y) >> 62);
    f <<= 1;
    if(y >> 63){
      y >>= 1;
      f |= 1;
    }
  }
  const uint64_t LN2 = 0xB17217F7D1CF79ACULL;
  const i128 total = (i128) (e - 1023) * (i128) LN2 + (i128) (((u128) f * LN2) >> 54);
  return (int64_t) ((total + ((i128) 1 << 11)) >> 12);
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

// get48 reading only the slot's 6 bytes: for a slot whose 8-byte load would reach into memory another thread may be writing
inline uint64_t get48_exact(const char *p){
  const unsigned char *u = reinterpret_cast<const unsigned char *>(p);
  return (uint64_t) u[0] | (uint64_t) u[1] << 8 | (uint64_t) u[2] << 16 | (uint64_t) u[3] << 24 | (uint64_t) u[4] << 32 |
         (uint64_t) u[5] << 40;
}

// Slack bytes after the last 6-byte slot of a buffer read by get48
const size_t SLACK = 8;

// Λ tiles: TILE x TILE entries, 6 bytes each, tile (I, J) with I <= J of a block at the block's offset + TILE_BYTES (J (J + 1) / 2 + I),
// entry (a, b) of a tile at 6 (TILE b + a)
const size_t TILE = 8;
const size_t TILE_BYTES = 6 * TILE * TILE;

// Byte offset in its block of entry (i, j) of a tile at or above the diagonal tiles (i / TILE <= j / TILE)
inline size_t tile_offset(size_t i, size_t j){
  const size_t I = i / TILE, J = j / TILE;
  return TILE_BYTES * (J * (J + 1) / 2 + I) + 6 * (TILE * (j % TILE) + i % TILE);
}

// Entry (i, j) of the block whose tiles start at T: an entry below the diagonal tiles reads its mirror (j, i); diagonal tiles are
// whole and read as stored
inline int64_t tile_value(const char *T, size_t i, size_t j){
  if(i / TILE > j / TILE) return (int64_t) get48(T + tile_offset(j, i));
  return (int64_t) get48(T + tile_offset(i, j));
}

// Rows [i0, i1) of the ncol columns of tile column J of the block whose tiles start at T, column-major into strip (column b at
// strip + b (i1 - i0)), from whole tiles: tile (I, J) at and above the diagonal, tile (J, I) transposed below it, so each tile is
// read once for all its columns
inline void tile_strip(const char *T, size_t J, size_t ncol, size_t i0, size_t i1, int64_t *strip){
  const size_t m = i1 - i0;
  for(size_t I = i0 / TILE; TILE * I < i1; ++I){
    const size_t r0 = i0 > TILE * I ? i0 : TILE * I, r1 = i1 < TILE * I + TILE ? i1 : TILE * I + TILE;
    if(I <= J){
      const char *t = T + TILE_BYTES * (J * (J + 1) / 2 + I);
      for(size_t b = 0; b < ncol; ++b){
        int64_t *col = strip + b * m;
        for(size_t r = r0; r < r1; ++r) col[r - i0] = (int64_t) get48(t + 6 * (TILE * b + r % TILE));
      }
    } else {
      const char *t = T + TILE_BYTES * (I * (I + 1) / 2 + J);
      for(size_t b = 0; b < ncol; ++b){
        int64_t *col = strip + b * m;
        for(size_t r = r0; r < r1; ++r) col[r - i0] = (int64_t) get48(t + 6 * (TILE * (r % TILE) + b));
      }
    }
  }
}

}  // namespace f2

#endif
