// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Bit depth independent read-only data.

#include "f265/bdi.h"

// Hadamard matrix.
const int8_t f265_h8x8[F265_H_DIM][F265_H_DIM] =
{
    {1,  1,  1,  1,  1,  1,  1,  1 },
    {1, -1,  1, -1,  1, -1,  1, -1 },
    {1,  1, -1, -1,  1,  1, -1, -1 },
    {1, -1, -1,  1,  1, -1, -1,  1 },
    {1,  1,  1,  1, -1, -1, -1, -1 },
    {1, -1,  1, -1, -1,  1, -1,  1 },
    {1,  1, -1, -1, -1, -1,  1,  1 },
    {1, -1, -1,  1, -1,  1,  1, -1 }
};

// Luma interpolation filter. The first index is the quarterpel component. The
// second index is the filter coefficient.
const int8_t f265_if_luma[4][8] =
{
    {  0, 0,   0, 64,  0,   0, 0,  0 },
    { -1, 4, -10, 58, 17,  -5, 1,  0 },
    { -1, 4, -11, 40, 40, -11, 4, -1 },
    {  0, 1,  -5, 17, 58, -10, 4, -1 }
};

// Chroma interpolation filter (with 1/8 pixel components).
const int8_t f265_if_chroma[8][4] =
{
    {  0, 64,  0,  0 },
    { -2, 58, 10, -2 },
    { -4, 54, 16, -2 },
    { -6, 46, 28, -4 },
    { -4, 36, 36, -4 },
    { -4, 28, 46, -6 },
    { -2, 16, 54, -4 },
    { -2, 10, 58, -2 }
};

// CABAC context initialization values for I, P, B frames.
const uint8_t f265_cabac_ctx_init_table[3][F265_NB_CABAC_CTX] =
{
    {
        153, 200, 139, 141, 157, 154, 154, 154, 154, 154, 184, 154, 154, 154, 184,  63,
        154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 153, 138, 138, 154, 111, 141,
         94, 138, 182, 154, 154, 154, 154, 154, 139, 139, 110, 110, 124, 125, 140, 153,
        125, 127, 140, 109, 111, 143, 127, 111,  79, 108, 123,  63, 110, 110, 124, 125,
        140, 153, 125, 127, 140, 109, 111, 143, 127, 111,  79, 108, 123,  63,  91, 171,
        134, 141, 111, 111, 125, 110, 110,  94, 124, 108, 124, 107, 125, 141, 179, 153,
        125, 107, 125, 141, 179, 153, 125, 107, 125, 141, 179, 153, 125, 140, 139, 182,
        182, 152, 136, 152, 136, 153, 136, 139, 111, 136, 139, 111, 140,  92, 137, 138,
        140, 152, 138, 139, 153,  74, 149,  92, 139, 107, 122, 152, 140, 179, 166, 182,
        140, 227, 122, 197, 138, 153, 136, 167, 152, 152,
    },
    {
        153, 185, 107, 139, 126, 154, 197, 185, 201, 149, 154, 139, 154, 154, 154, 152,
        110, 122,  95,  79,  63,  31,  31, 153, 153, 168, 124, 138,  94,  79, 153, 111,
        149, 107, 167, 154, 140, 198, 154, 154, 139, 139, 125, 110,  94, 110,  95,  79,
        125, 111, 110,  78, 110, 111, 111,  95,  94, 108, 123, 108, 125, 110,  94, 110,
         95,  79, 125, 111, 110,  78, 110, 111, 111,  95,  94, 108, 123, 108, 121, 140,
         61, 154, 155, 154, 139, 153, 139, 123, 123,  63, 153, 166, 183, 140, 136, 153,
        154, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 170, 153, 123,
        123, 107, 121, 107, 121, 167, 151, 183, 140, 151, 183, 140, 154, 196, 196, 167,
        154, 152, 167, 182, 182, 134, 149, 136, 153, 121, 136, 137, 169, 194, 166, 167,
        154, 167, 137, 182, 107, 167,  91, 122, 107, 167,
    },
    {
        153, 160, 107, 139, 126, 154, 197, 185, 201, 134, 154, 139, 154, 154, 183, 152,
        154, 137,  95,  79,  63,  31,  31, 153, 153, 168, 224, 167, 122,  79, 153, 111,
        149,  92, 167, 154, 169, 198, 154, 154, 139, 139, 125, 110, 124, 110,  95,  94,
        125, 111, 111,  79, 125, 126, 111, 111,  79, 108, 123,  93, 125, 110, 124, 110,
         95,  94, 125, 111, 111,  79, 125, 126, 111, 111,  79, 108, 123,  93, 121, 140,
         61, 154, 170, 154, 139, 153, 139, 123, 123,  63, 124, 166, 183, 140, 136, 153,
        154, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 170, 153, 138,
        138, 122, 121, 122, 121, 167, 151, 183, 140, 151, 183, 140, 154, 196, 167, 167,
        154, 152, 167, 182, 182, 134, 149, 136, 153, 121, 136, 122, 169, 208, 166, 167,
        154, 152, 167, 182, 107, 167,  91, 107, 107, 167,
    },
};

// Each CABAC context has 6 bits of probability state and 1 bit for MPS symbol
// identification (whether bit 0 or 1 is the most probable symbol), laid out as
// follow: [prob|mps].
//
// This table yields the LPS value of the current bin, which is used to update
// the interval range. The first index is the context probability. The second
// index is the range scale (bits 6 and 7 of the range value).
const uint8_t f265_cabac_range_table[64][4] =
{
    { 128, 176, 208, 240 },
    { 128, 167, 197, 227 },
    { 128, 158, 187, 216 },
    { 123, 150, 178, 205 },
    { 116, 142, 169, 195 },
    { 111, 135, 160, 185 },
    { 105, 128, 152, 175 },
    { 100, 122, 144, 166 },
    {  95, 116, 137, 158 },
    {  90, 110, 130, 150 },
    {  85, 104, 123, 142 },
    {  81,  99, 117, 135 },
    {  77,  94, 111, 128 },
    {  73,  89, 105, 122 },
    {  69,  85, 100, 116 },
    {  66,  80,  95, 110 },
    {  62,  76,  90, 104 },
    {  59,  72,  86,  99 },
    {  56,  69,  81,  94 },
    {  53,  65,  77,  89 },
    {  51,  62,  73,  85 },
    {  48,  59,  69,  80 },
    {  46,  56,  66,  76 },
    {  43,  53,  63,  72 },
    {  41,  50,  59,  69 },
    {  39,  48,  56,  65 },
    {  37,  45,  54,  62 },
    {  35,  43,  51,  59 },
    {  33,  41,  48,  56 },
    {  32,  39,  46,  53 },
    {  30,  37,  43,  50 },
    {  29,  35,  41,  48 },
    {  27,  33,  39,  45 },
    {  26,  31,  37,  43 },
    {  24,  30,  35,  41 },
    {  23,  28,  33,  39 },
    {  22,  27,  32,  37 },
    {  21,  26,  30,  35 },
    {  20,  24,  29,  33 },
    {  19,  23,  27,  31 },
    {  18,  22,  26,  30 },
    {  17,  21,  25,  28 },
    {  16,  20,  23,  27 },
    {  15,  19,  22,  25 },
    {  14,  18,  21,  24 },
    {  14,  17,  20,  23 },
    {  13,  16,  19,  22 },
    {  12,  15,  18,  21 },
    {  12,  14,  17,  20 },
    {  11,  14,  16,  19 },
    {  11,  13,  15,  18 },
    {  10,  12,  15,  17 },
    {  10,  12,  14,  16 },
    {   9,  11,  13,  15 },
    {   9,  11,  12,  14 },
    {   8,  10,  12,  14 },
    {   8,   9,  11,  13 },
    {   7,   9,  11,  12 },
    {   7,   9,  10,  12 },
    {   7,   8,  10,  11 },
    {   6,   8,   9,  11 },
    {   6,   7,   9,  10 },
    {   6,   7,   8,   9 },
    {   2,   2,   2,   2 }
};

// This table yields the updated context after encoding the current bin. The
// first index is the context. The second is the bin value.
//
// The entries 0 and 1 correspond to probability 50% with the MPS values 0 and 1
// respectively. The MPS symbol inverts when the probability falls below 50%.
// The entries 124 and 125 correspond to the saturated probabilities. The
// entries 126 and 127 are self-reflecting. The table follows the "exponential
// aging" model.
//
// Mapping of the 2D array f265_cabac_transit_table[128][2] (bdi_ro.c) to Table
// 9-41 - State transition table. In our code, the even indices relate to the case
// where the most probable symbol (MPS) is zero. These are the values presented
// below. The same mapping is observed for the even indices (when MPS is one).
//
// The first column relates to the index used to access f265_cabac_transit_table,
// while the second columns refers directly to the pStateIdx value in Table 9-41.
// The third and fourth columns indicate the transitions. The first value is the
// transition used in f265_cabac_transit_table, while the value between
// parentheses is found in Table 9-41. Note that the least probable symbol (LPS)
// in the table below refers to a bin of 1, while the MPS refers to a bin of 0.
// When comparing with the values found in Table 9-41, the order is inverted.
//
// +-----+----------+--------+-------------+--------+-------------+
// | idx | pSateIdx | LPS(1) | transIdxLps | MPS(0) | transIdxMps |
// +-----+----------+--------+-------------+--------+-------------+
// |   0 |        0 |      1*|           0 |      2 |           1 |
// |   2 |        1 |      0 |           0 |      4 |           2 |
// |   4 |        2 |      2 |           1 |      6 |           3 |
// |   6 |        3 |      4 |           2 |      8 |           4 |
// | ... |      ... |    ... |         ... |    ... |         ... |
// | 124 |       62 |     76 |          38 |    124^|          62 |
// +-----+----------+--------+-------------+--------+-------------+
// * Using LPS shifts to the odd indices (MPS = 1)
// ^ MPS probability saturation (self-reflecting state)
const uint8_t f265_cabac_transit_table[128][2] =
{
    {2,1},    {0,3},    {4,0},    {1,5},    {6,2},    {3,7},    {8,4},    {5,9},
    {10,4},   {5,11},   {12,8},   {9,13},   {14,8},   {9,15},   {16,10},  {11,17},
    {18,12},  {13,19},  {20,14},  {15,21},  {22,16},  {17,23},  {24,18},  {19,25},
    {26,18},  {19,27},  {28,22},  {23,29},  {30,22},  {23,31},  {32,24},  {25,33},
    {34,26},  {27,35},  {36,26},  {27,37},  {38,30},  {31,39},  {40,30},  {31,41},
    {42,32},  {33,43},  {44,32},  {33,45},  {46,36},  {37,47},  {48,36},  {37,49},
    {50,38},  {39,51},  {52,38},  {39,53},  {54,42},  {43,55},  {56,42},  {43,57},
    {58,44},  {45,59},  {60,44},  {45,61},  {62,46},  {47,63},  {64,48},  {49,65},
    {66,48},  {49,67},  {68,50},  {51,69},  {70,52},  {53,71},  {72,52},  {53,73},
    {74,54},  {55,75},  {76,54},  {55,77},  {78,56},  {57,79},  {80,58},  {59,81},
    {82,58},  {59,83},  {84,60},  {61,85},  {86,60},  {61,87},  {88,60},  {61,89},
    {90,62},  {63,91},  {92,64},  {65,93},  {94,64},  {65,95},  {96,66},  {67,97},
    {98,66},  {67,99},  {100,66}, {67,101}, {102,68}, {69,103}, {104,68}, {69,105},
    {106,70}, {71,107}, {108,70}, {71,109}, {110,70}, {71,111}, {112,72}, {73,113},
    {114,72}, {73,115}, {116,72}, {73,117}, {118,74}, {75,119}, {120,74}, {75,121},
    {122,74}, {75,123}, {124,76}, {77,125}, {124,76}, {77,125}, {126,126},{127,127}
};

// Entropy lost when encoding a bit with the context specified. The index is
// context^bin.
//
// FIXME: this is straight from HM. We'll have to convert it to 16 bits
// eventually.
//
// Explanations needed: the first version seems to have been derived
// experimentally, the entropy for the two 50% probability entries are off. The
// second version matches the entropy math directly.
//
// Cabac probability range = [0.01875, 0.5].
// Entropy = -lg(p).
// Last entry, MPS: -lg(1-0.01875) ~= 0x0037f/32768.
// Last entry, LPS: -lg(0.01875) ~= 0x2de55/32768.
const uint32_t f265_cabac_entropy_table[128] =
{
#if 1
  // Corrected table, most notably for the last state.
  0x07b23, 0x085f9, 0x074a0, 0x08cbc, 0x06ee4, 0x09354, 0x067f4, 0x09c1b, 0x060b0, 0x0a62a, 0x05a9c, 0x0af5b, 0x0548d, 0x0b955, 0x04f56, 0x0c2a9,
  0x04a87, 0x0cbf7, 0x045d6, 0x0d5c3, 0x04144, 0x0e01b, 0x03d88, 0x0e937, 0x039e0, 0x0f2cd, 0x03663, 0x0fc9e, 0x03347, 0x10600, 0x03050, 0x10f95,
  0x02d4d, 0x11a02, 0x02ad3, 0x12333, 0x0286e, 0x12cad, 0x02604, 0x136df, 0x02425, 0x13f48, 0x021f4, 0x149c4, 0x0203e, 0x1527b, 0x01e4d, 0x15d00,
  0x01c99, 0x166de, 0x01b18, 0x17017, 0x019a5, 0x17988, 0x01841, 0x18327, 0x016df, 0x18d50, 0x015d9, 0x19547, 0x0147c, 0x1a083, 0x0138e, 0x1a8a3,
  0x01251, 0x1b418, 0x01166, 0x1bd27, 0x01068, 0x1c77b, 0x00f7f, 0x1d18e, 0x00eda, 0x1d91a, 0x00e19, 0x1e254, 0x00d4f, 0x1ec9a, 0x00c90, 0x1f6e0,
  0x00c01, 0x1fef8, 0x00b5f, 0x208b1, 0x00ab6, 0x21362, 0x00a15, 0x21e46, 0x00988, 0x2285d, 0x00934, 0x22ea8, 0x008a8, 0x239b2, 0x0081d, 0x24577,
  0x007c9, 0x24ce6, 0x00763, 0x25663, 0x00710, 0x25e8f, 0x006a0, 0x26a26, 0x00672, 0x26f23, 0x005e8, 0x27ef8, 0x005ba, 0x284b5, 0x0055e, 0x29057,
  0x0050c, 0x29bab, 0x004c1, 0x2a674, 0x004a7, 0x2aa5e, 0x0046f, 0x2b32f, 0x0041f, 0x2c0ad, 0x003e7, 0x2ca8d, 0x003ba, 0x2d323, 0x0010c, 0x3bfbb
#else
  0x08000, 0x08000, 0x076da, 0x089a0, 0x06e92, 0x09340, 0x0670a, 0x09cdf, 0x06029, 0x0a67f, 0x059dd, 0x0b01f, 0x05413, 0x0b9bf, 0x04ebf, 0x0c35f,
  0x049d3, 0x0ccff, 0x04546, 0x0d69e, 0x0410d, 0x0e03e, 0x03d22, 0x0e9de, 0x0397d, 0x0f37e, 0x03619, 0x0fd1e, 0x032ee, 0x106be, 0x02ffa, 0x1105d,
  0x02d37, 0x119fd, 0x02aa2, 0x1239d, 0x02836, 0x12d3d, 0x025f2, 0x136dd, 0x023d1, 0x1407c, 0x021d2, 0x14a1c, 0x01ff2, 0x153bc, 0x01e2f, 0x15d5c,
  0x01c87, 0x166fc, 0x01af7, 0x1709b, 0x0197f, 0x17a3b, 0x0181d, 0x183db, 0x016d0, 0x18d7b, 0x01595, 0x1971b, 0x0146c, 0x1a0bb, 0x01354, 0x1aa5a,
  0x0124c, 0x1b3fa, 0x01153, 0x1bd9a, 0x01067, 0x1c73a, 0x00f89, 0x1d0da, 0x00eb7, 0x1da79, 0x00df0, 0x1e419, 0x00d34, 0x1edb9, 0x00c82, 0x1f759,
  0x00bda, 0x200f9, 0x00b3c, 0x20a99, 0x00aa5, 0x21438, 0x00a17, 0x21dd8, 0x00990, 0x22778, 0x00911, 0x23118, 0x00898, 0x23ab8, 0x00826, 0x24458,
  0x007ba, 0x24df7, 0x00753, 0x25797, 0x006f2, 0x26137, 0x00696, 0x26ad7, 0x0063f, 0x27477, 0x005ed, 0x27e17, 0x0059f, 0x287b6, 0x00554, 0x29156,
  0x0050e, 0x29af6, 0x004cc, 0x2a497, 0x0048d, 0x2ae35, 0x00451, 0x2b7d6, 0x00418, 0x2c176, 0x003e2, 0x2cb15, 0x003af, 0x2d4b5, 0x0037f, 0x2de55
#endif
};

// Unskewed version of the above for testing purposes.
const uint32_t f265_cabac_unskewed_entropy_table[128] =
{
  0x08000, 0x08000, 0x076da, 0x089a0, 0x06e92, 0x09340, 0x0670a, 0x09cdf, 0x06029, 0x0a67f, 0x059dd, 0x0b01f, 0x05413, 0x0b9bf, 0x04ebf, 0x0c35f,
  0x049d3, 0x0ccff, 0x04546, 0x0d69e, 0x0410d, 0x0e03e, 0x03d22, 0x0e9de, 0x0397d, 0x0f37e, 0x03619, 0x0fd1e, 0x032ee, 0x106be, 0x02ffa, 0x1105d,
  0x02d37, 0x119fd, 0x02aa2, 0x1239d, 0x02836, 0x12d3d, 0x025f2, 0x136dd, 0x023d1, 0x1407c, 0x021d2, 0x14a1c, 0x01ff2, 0x153bc, 0x01e2f, 0x15d5c,
  0x01c87, 0x166fc, 0x01af7, 0x1709b, 0x0197f, 0x17a3b, 0x0181d, 0x183db, 0x016d0, 0x18d7b, 0x01595, 0x1971b, 0x0146c, 0x1a0bb, 0x01354, 0x1aa5a,
  0x0124c, 0x1b3fa, 0x01153, 0x1bd9a, 0x01067, 0x1c73a, 0x00f89, 0x1d0da, 0x00eb7, 0x1da79, 0x00df0, 0x1e419, 0x00d34, 0x1edb9, 0x00c82, 0x1f759,
  0x00bda, 0x200f9, 0x00b3c, 0x20a99, 0x00aa5, 0x21438, 0x00a17, 0x21dd8, 0x00990, 0x22778, 0x00911, 0x23118, 0x00898, 0x23ab8, 0x00826, 0x24458,
  0x007ba, 0x24df7, 0x00753, 0x25797, 0x006f2, 0x26137, 0x00696, 0x26ad7, 0x0063f, 0x27477, 0x005ed, 0x27e17, 0x0059f, 0x287b6, 0x00554, 0x29156,
  0x0050e, 0x29af6, 0x004cc, 0x2a497, 0x0048d, 0x2ae35, 0x00451, 0x2b7d6, 0x00418, 0x2c176, 0x003e2, 0x2cb15, 0x003af, 0x2d4b5, 0x0037f, 0x2de55
};

// Table of precomputed VLC values. The first index is the code. The second
// index is 0 for the VLC value, 1 for the VLC length.
const uint8_t f265_vlc_table[64][2] =
{
    {1,1},
    {2,3}, {3,3},
    {4,5}, {5,5}, {6,5}, {7,5},
    {8,7}, {9,7}, {10,7},{11,7},{12,7},{13,7},{14,7},{15,7},
    {16,9},{17,9},{18,9},{19,9},{20,9},{21,9},{22,9},{23,9},
    {24,9},{25,9},{26,9},{27,9},{28,9},{29,9},{30,9},{31,9},
    {32,11},{33,11},{34,11},{35,11},{36,11},{37,11},{38,11},{39,11},
    {40,11},{41,11},{42,11},{43,11},{44,11},{45,11},{46,11},{47,11},
    {48,11},{49,11},{50,11},{51,11},{52,11},{53,11},{54,11},{55,11},
    {56,11},{57,11},{58,11},{59,11},{60,11},{61,11},{62,11},{63,11},
    {64,13}
};

// Default quantization/dequantization multipliers (QP%6).
const int16_t f265_quant_mult[6] = { 26214, 23302, 20560, 18396, 16384, 14564 };
const int16_t f265_dequant_mult[6] = { 40, 45, 51, 57, 64, 72 };

// DST matrix.
const int8_t f265_dst_mat[4][4] =
{
    { 29, 55, 74, 84},
    { 74, 74,  0,-74},
    { 84,-29,-74, 55},
    { 55,-84, 74,-29}
};

// DCT matrices.
const int8_t f265_dct_mat_4[4][4] =
{
    { 64, 64, 64, 64},
    { 83, 36,-36,-83},
    { 64,-64,-64, 64},
    { 36,-83, 83,-36}
};

const int8_t f265_dct_mat_8[8][8] =
{
    { 64, 64, 64, 64, 64, 64, 64, 64},
    { 89, 75, 50, 18,-18,-50,-75,-89},
    { 83, 36,-36,-83,-83,-36, 36, 83},
    { 75,-18,-89,-50, 50, 89, 18,-75},
    { 64,-64,-64, 64, 64,-64,-64, 64},
    { 50,-89, 18, 75,-75,-18, 89,-50},
    { 36,-83, 83,-36,-36, 83,-83, 36},
    { 18,-50, 75,-89, 89,-75, 50,-18}
};

const int8_t f265_dct_mat_16[16][16] =
{
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
    { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
    { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
    { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
    { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
    { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
    { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
    { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
    { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
    { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
    { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
    { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
    { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
    { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
    { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
    {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

const int8_t f265_dct_mat_32[32][32] =
{
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
    { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4, -4,-13,-22,-31,-38,-46,-54,-61,-67,-73,-78,-82,-85,-88,-90,-90},
    { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90,-90,-87,-80,-70,-57,-43,-25, -9,  9, 25, 43, 57, 70, 80, 87, 90},
    { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31,  4,-22,-46,-67,-82,-90},
    { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89, 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
    { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22,-22,-61,-85,-90,-73,-38,  4, 46, 78, 90, 82, 54, 13,-31,-67,-88},
    { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87,-87,-57, -9, 43, 80, 90, 70, 25,-25,-70,-90,-80,-43,  9, 57, 87},
    { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31, 31, 78, 90, 61,  4,-54,-88,-82,-38, 22, 73, 90, 67, 13,-46,-85},
    { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
    { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38,-38,-88,-73, -4, 67, 90, 46,-31,-85,-78,-13, 61, 90, 54,-22,-82},
    { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80,-80, -9, 70, 87, 25,-57,-90,-43, 43, 90, 57,-25,-87,-70,  9, 80},
    { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46, 46, 90, 38,-54,-90,-31, 61, 88, 22,-67,-85,-13, 73, 82,  4,-78},
    { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75, 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
    { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54,-54,-85,  4, 88, 46,-61,-82, 13, 90, 38,-67,-78, 22, 90, 31,-73},
    { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70,-70, 43, 87, -9,-90,-25, 80, 57,-57,-80, 25, 90,  9,-87,-43, 70},
    { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61, 61, 73,-46,-82, 31, 88,-13,-90, -4, 90, 22,-85,-38, 78, 54,-67},
    { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
    { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67,-67,-54, 78, 38,-85,-22, 90,  4,-90, 13, 88,-31,-82, 46, 73,-61},
    { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57,-57, 80, 25,-90,  9, 87,-43,-70, 70, 43,-87, -9, 90,-25,-80, 57},
    { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73, 73, 31,-90, 22, 78,-67,-38, 90,-13,-82, 61, 46,-88,  4, 85,-54},
    { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50, 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
    { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78,-78, -4, 82,-73,-13, 85,-67,-22, 88,-61,-31, 90,-54,-38, 90,-46},
    { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43,-43, 90,-57,-25, 87,-70, -9, 80,-80,  9, 70,-87, 25, 57,-90, 43},
    { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82, 82,-22,-54, 90,-61,-13, 78,-85, 31, 46,-90, 67,  4,-73, 88,-38},
    { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
    { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85,-85, 46, 13,-67, 90,-73, 22, 38,-82, 88,-54, -4, 61,-90, 78,-31},
    { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25,-25, 70,-90, 80,-43, -9, 57,-87, 87,-57,  9, 43,-80, 90,-70, 25},
    { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88, 88,-67, 31, 13,-54, 82,-90, 78,-46,  4, 38,-73, 90,-85, 61,-22},
    { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18, 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
    { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90,-90, 82,-67, 46,-22, -4, 31,-54, 73,-85, 90,-88, 78,-61, 38,-13},
    {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9, -9, 25,-43, 57,-70, 80,-87, 90,-90, 87,-80, 70,-57, 43,-25,  9},
    {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90, 90,-90, 88,-85, 82,-78, 73,-67, 61,-54, 46,-38, 31,-22, 13, -4}
};

// Map to convert from a position in last-to-first encoding order to a position
// in raster scan. The first index is the block size (1x1, 2x2, 4x4, 8x8). The
// second index is the encoding order (up-right, horizontal, vertical). The
// output value is an index in venc_scan_map_data[], as follow.
//
// venc_scan_map_data + index is the address of the position map. This is not a
// symmetrical map, i.e. map[map[X]] != X for up-right.
// Example: raster   up-right/last-to-first
//          [0123]     [FBE7]
//          [4567] <== [AD36]
//          [89AB]     [9C25]
//          [CDEF]     [8140]
//
// venc_scan_map_data + index + 256 is the address of a map yielding the
// position of the first coefficient of a subblock inside a transform block,
// given a subblock position in encoding order. The map value must be multiplied
// by four prior to indexing the coefficient array of the transform block.
const uint8_t f265_scan_map_idx[4][3] =
{
    { 0xfc, 0xfd, 0xfe },
    { 0xf0, 0xf4, 0xf8 },
    { 0xc0, 0xd0, 0xe0 },
    { 0x00, 0x40, 0x80 },
};
const uint8_t f265_scan_map_data[2*256] =
{
    0x3f, 0x37, 0x3e, 0x2f, 0x36, 0x3d, 0x27, 0x2e, 0x35, 0x3c, 0x1f, 0x26, 0x2d, 0x34, 0x3b, 0x17,
    0x1e, 0x25, 0x2c, 0x33, 0x3a, 0x0f, 0x16, 0x1d, 0x24, 0x2b, 0x32, 0x39, 0x07, 0x0e, 0x15, 0x1c,
    0x23, 0x2a, 0x31, 0x38, 0x06, 0x0d, 0x14, 0x1b, 0x22, 0x29, 0x30, 0x05, 0x0c, 0x13, 0x1a, 0x21,
    0x28, 0x04, 0x0b, 0x12, 0x19, 0x20, 0x03, 0x0a, 0x11, 0x18, 0x02, 0x09, 0x10, 0x01, 0x08, 0x00,
    0x3f, 0x3e, 0x3d, 0x3c, 0x3b, 0x3a, 0x39, 0x38, 0x37, 0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30,
    0x2f, 0x2e, 0x2d, 0x2c, 0x2b, 0x2a, 0x29, 0x28, 0x27, 0x26, 0x25, 0x24, 0x23, 0x22, 0x21, 0x20,
    0x1f, 0x1e, 0x1d, 0x1c, 0x1b, 0x1a, 0x19, 0x18, 0x17, 0x16, 0x15, 0x14, 0x13, 0x12, 0x11, 0x10,
    0x0f, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09, 0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
    0x3f, 0x37, 0x2f, 0x27, 0x1f, 0x17, 0x0f, 0x07, 0x3e, 0x36, 0x2e, 0x26, 0x1e, 0x16, 0x0e, 0x06,
    0x3d, 0x35, 0x2d, 0x25, 0x1d, 0x15, 0x0d, 0x05, 0x3c, 0x34, 0x2c, 0x24, 0x1c, 0x14, 0x0c, 0x04,
    0x3b, 0x33, 0x2b, 0x23, 0x1b, 0x13, 0x0b, 0x03, 0x3a, 0x32, 0x2a, 0x22, 0x1a, 0x12, 0x0a, 0x02,
    0x39, 0x31, 0x29, 0x21, 0x19, 0x11, 0x09, 0x01, 0x38, 0x30, 0x28, 0x20, 0x18, 0x10, 0x08, 0x00,
    0x0f, 0x0b, 0x0e, 0x07, 0x0a, 0x0d, 0x03, 0x06, 0x09, 0x0c, 0x02, 0x05, 0x08, 0x01, 0x04, 0x00,
    0x0f, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09, 0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
    0x0f, 0x0b, 0x07, 0x03, 0x0e, 0x0a, 0x06, 0x02, 0x0d, 0x09, 0x05, 0x01, 0x0c, 0x08, 0x04, 0x00,
    0x03, 0x01, 0x02, 0x00, 0x03, 0x02, 0x01, 0x00, 0x03, 0x01, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00,
    0xe7, 0xc7, 0xe6, 0xa7, 0xc6, 0xe5, 0x87, 0xa6, 0xc5, 0xe4, 0x67, 0x86, 0xa5, 0xc4, 0xe3, 0x47,
    0x66, 0x85, 0xa4, 0xc3, 0xe2, 0x27, 0x46, 0x65, 0x84, 0xa3, 0xc2, 0xe1, 0x07, 0x26, 0x45, 0x64,
    0x83, 0xa2, 0xc1, 0xe0, 0x06, 0x25, 0x44, 0x63, 0x82, 0xa1, 0xc0, 0x05, 0x24, 0x43, 0x62, 0x81,
    0xa0, 0x04, 0x23, 0x42, 0x61, 0x80, 0x03, 0x22, 0x41, 0x60, 0x02, 0x21, 0x40, 0x01, 0x20, 0x00,
    0xe7, 0xe6, 0xe5, 0xe4, 0xe3, 0xe2, 0xe1, 0xe0, 0xc7, 0xc6, 0xc5, 0xc4, 0xc3, 0xc2, 0xc1, 0xc0,
    0xa7, 0xa6, 0xa5, 0xa4, 0xa3, 0xa2, 0xa1, 0xa0, 0x87, 0x86, 0x85, 0x84, 0x83, 0x82, 0x81, 0x80,
    0x67, 0x66, 0x65, 0x64, 0x63, 0x62, 0x61, 0x60, 0x47, 0x46, 0x45, 0x44, 0x43, 0x42, 0x41, 0x40,
    0x27, 0x26, 0x25, 0x24, 0x23, 0x22, 0x21, 0x20, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
    0xe7, 0xc7, 0xa7, 0x87, 0x67, 0x47, 0x27, 0x07, 0xe6, 0xc6, 0xa6, 0x86, 0x66, 0x46, 0x26, 0x06,
    0xe5, 0xc5, 0xa5, 0x85, 0x65, 0x45, 0x25, 0x05, 0xe4, 0xc4, 0xa4, 0x84, 0x64, 0x44, 0x24, 0x04,
    0xe3, 0xc3, 0xa3, 0x83, 0x63, 0x43, 0x23, 0x03, 0xe2, 0xc2, 0xa2, 0x82, 0x62, 0x42, 0x22, 0x02,
    0xe1, 0xc1, 0xa1, 0x81, 0x61, 0x41, 0x21, 0x01, 0xe0, 0xc0, 0xa0, 0x80, 0x60, 0x40, 0x20, 0x00,
    0x33, 0x23, 0x32, 0x13, 0x22, 0x31, 0x03, 0x12, 0x21, 0x30, 0x02, 0x11, 0x20, 0x01, 0x10, 0x00,
    0x33, 0x32, 0x31, 0x30, 0x23, 0x22, 0x21, 0x20, 0x13, 0x12, 0x11, 0x10, 0x03, 0x02, 0x01, 0x00,
    0x33, 0x23, 0x13, 0x03, 0x32, 0x22, 0x12, 0x02, 0x31, 0x21, 0x11, 0x01, 0x30, 0x20, 0x10, 0x00,
    0x09, 0x01, 0x08, 0x00, 0x09, 0x08, 0x01, 0x00, 0x09, 0x01, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
};

// Table providing the encoding for the last coefficient position. The index is
// the position. The value is a bitfield where the low 4 bits is the number of
// "1" in the prefix and the high 4 bits is the number of bits in the suffix.
const uint8_t f265_last_coeff_table[32] =
{
    0x00, 0x01, 0x02, 0x03, 0x14, 0x14, 0x15, 0x15,
    0x26, 0x26, 0x26, 0x26, 0x27, 0x27, 0x27, 0x27,
    0x38, 0x38, 0x38, 0x38, 0x38, 0x38, 0x38, 0x38,
    0x39, 0x39, 0x39, 0x39, 0x39, 0x39, 0x39, 0x39,
};

// Coefficient non-zero flag context-derivation table. The first four rows
// correspond to the neighbour index. The last row is specific to the 4x4
// transform. The row is indexed by the coefficient position in raster scan.
const uint8_t f265_coeff_nz_ctx_table[5*16] =
{
    2, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    2, 1, 0, 0, 2, 1, 0, 0, 2, 1, 0, 0, 2, 1, 0, 0,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    0, 1, 4, 5, 2, 3, 4, 5, 6, 6, 8, 8, 7, 7, 8, 8,
};

// Coefficient greater-than-1 context counter update table.
const uint8_t f265_gt1_ctx_counter_table[8] = { 0, 2, 3, 3, 0, 0, 0, 0 };

// Intra MPM binarization table. The first two bits are the binarization string,
// the next two bits are the number of bins.
const uint8_t f265_mpm_bin_table[3] = { 0x4, 0xa, 0xb };

// Common chroma mode table, in order (Planar, Vertical, Horizontal, DC).
const uint8_t f265_chroma_mode_table[4] = { 0, 26, 10, 1 };

// Intra prediction angle and inverse angle tables. See usage.
const int8_t f265_intra_angle_table[17] = { -32, -26, -21, -17, -13, -9, -5, -2, 0, 2, 5, 9, 13, 17, 21, 26, 32 };
const int16_t f265_intra_inv_angle_table[8] = { 4096, 1638, 910, 630, 482, 390, 315, 256 };

// Intra prediction mode distances tweaked for comparison against the threshold
// values. The index is the mode.
const uint8_t f265_intra_mode_dist[35] =
{
    // Modes 2, 10, 18, 26, 34.
    //    v
    8, 0, 8, 7, 6, 5, 4, 3,
    2, 1, 0, 1, 2, 3, 4, 5,
    6, 7, 8, 7, 6, 5, 4, 3,
    2, 1, 0, 1, 2, 3, 4, 5,
    6, 7, 8
};

// Intra prediction mode distance thresholds. The index is the log of the block
// size minus 2.
const uint8_t f265_intra_dist_thresholds[4] = { 8, 7, 1, 0 };

// Reference index bit cost table. The first index is derived from the number of
// indices in the reference list: 0 if 1 index, 1 if 2 indices, 2 otherwise. The
// second index is the reference index. The costs are integer.
const int8_t f265_ref_bits[3][16] =
{
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
    { 1, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9 }
};

// Block partition table. The first index is the partitioning mode. The second
// index is the partition index. The third index is 0 for the partition sizes, 1
// for the partition offsets. The sizes and offsets are encoded as two packed
// 4-bit values (X,Y). X and Y represent 1/4 fractions of the block, e.g. X=2
// means the partition is half the size of the block horizontally.
const uint8_t f265_part_table[8][4][2] =
{
    { { 0x44, 0x00 }, { 0x00, 0x00 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // UN.
    { { 0x24, 0x00 }, { 0x24, 0x20 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // H2.
    { { 0x42, 0x00 }, { 0x42, 0x02 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // V2.
    { { 0x22, 0x00 }, { 0x22, 0x02 }, { 0x22, 0x20 }, { 0x22, 0x22 } }, // HV.
    { { 0x14, 0x00 }, { 0x34, 0x10 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // H1.
    { { 0x34, 0x00 }, { 0x14, 0x30 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // H3.
    { { 0x41, 0x00 }, { 0x43, 0x01 }, { 0x00, 0x00 }, { 0x00, 0x00 } }, // V1.
    { { 0x43, 0x00 }, { 0x41, 0x03 }, { 0x00, 0x00 }, { 0x00, 0x00 } }  // V3.
};

// Partition count table. The index is the partitioning mode.
const uint8_t f265_nb_parts[8] = { 1, 2, 2, 4, 2, 2, 2, 2 };

// Table yielding the AMP width index corresponding to the block width specified
// (0xff if unmapped).
const uint8_t f265_width_to_awi[65] =
{
    0xff, 0xff, 0x00, 0xff, 0x01, 0xff, 0x06, 0xff,
    0x02, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0xff,
    0x03, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x08, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x04, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x09, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x05
};

// Deblocking filter tC table. The index is the QP.
const uint8_t f265_tc_table[54] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3,
    3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 11, 13,
    14, 16, 18, 20, 22, 24
};

// Motion vector cost table. The first index is the QP. The second index is
// ActualMv - PredictedMv + F265_NB_MV_COSTS/2. This is initialized when the
// library is loaded.
uint16_t f265_mv_costs[52][F265_NB_MV_COSTS];

// Distortion lambda for MV cost. The index represent the QP.
const int16_t f265_lambdas[52] =
{
    1 * 16,   1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,  1 * 16,
    1 * 16,   1 * 16,  1 * 16,  2 * 16,  2 * 16,  2 * 16,  2 * 16,  3 * 16,  3 * 16,  3 * 16,  4 * 16,  4 * 16,  4 * 16,
    5 * 16,   6 * 16,  6 * 16,  7 * 16,  8 * 16,  9 * 16, 10 * 16, 11 * 16, 13 * 16, 14 * 16, 16 * 16, 18 * 16, 20 * 16,
    23 * 16, 25 * 16, 29 * 16, 32 * 16, 36 * 16, 40 * 16, 45 * 16, 51 * 16, 57 * 16, 64 * 16, 72 * 16, 81 * 16, 91 * 16,
};

// The following tables are used to determine the sources used in the
// interpolation of a luma block given the subpel components of its motion
// vector. The table index is mvy*4 + mvx (quarterpel component only).
//
// Luma interpolation layout:
// FA HA FB   FA, FB, FC, FD: fullpels A, B, C, D.
// VA DA VB   HA, HC:         horizontal halfpel of A and C.
// FC HC FD   VA, VB:         vertical halfpel of A and B.
//            DA:             diagonal halfpel of A.
//
// Interpolation logic table, by index:
// 0000 0 fullpel.
// 0001 1 qpel FA HA.
// 0010 2 hpel HA.
// 0011 3 qpel HA FB (bump src1).
// 0100 4 qpel FA VA.
// 0101 5 qpel HA VA.
// 0110 6 qpel HA DA.
// 0111 7 qpel HA VB (bump src1).
// 1000 8 hpel VA.
// 1001 9 qpel VA DA.
// 1010 a hpel DA.
// 1011 b qpel DA VB (bump src1).
// 1100 c qpel FC VA (bump src0).
// 1101 d qpel HC VA (bump src0).
// 1110 e hpel HC DA (bump src0).
// 1111 f qpel HC VB (bump src0 and src1).
//
// Notice that the values associated to the fullpels B, C, D are only used when
// a subpel component is equal to 3.
//
//                                  0  1  2  3  4  5  6  7  8  9  a  b  c  d  e  f
const int8_t f265_hpel_src0[16] = { 0, 0, 1, 1, 0, 1, 1, 1, 2, 2, 3, 3, 0, 1, 1, 1 };
const int8_t f265_hpel_src1[16] = { 9, 1, 9, 0, 2, 2, 3, 2, 9, 3, 9, 2, 2, 2, 3, 2 };

