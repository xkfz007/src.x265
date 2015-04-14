/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *          Mandar Gurav <mandar@multicorewareinc.com>
 *          Deepthi Devaki Akkoorath <deepthidevaki@multicorewareinc.com>
 *          Mahesh Pittala <mahesh@multicorewareinc.com>
 *          Rajesh Paulraj <rajesh@multicorewareinc.com>
 *          Praveen Kumar Tiwari <praveen@multicorewareinc.com>
 *          Min Chen <chenm003@163.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com.
 *****************************************************************************/

#ifndef X265_PRIMITIVES_H
#define X265_PRIMITIVES_H

#include "common.h"
#include "cpu.h"

namespace x265 {
// x265 private namespace

enum LumaPartitions
{
    // Square
    LUMA_4x4,   LUMA_8x8,   LUMA_16x16, LUMA_32x32, LUMA_64x64,
    // Rectangular
    LUMA_8x4,   LUMA_4x8,
    LUMA_16x8,  LUMA_8x16,  
    LUMA_32x16, LUMA_16x32,
    LUMA_64x32, LUMA_32x64,
    // Asymmetrical (0.75, 0.25)
    LUMA_16x12, LUMA_12x16, LUMA_16x4,  LUMA_4x16,
    LUMA_32x24, LUMA_24x32, LUMA_32x8,  LUMA_8x32,
    LUMA_64x48, LUMA_48x64, LUMA_64x16, LUMA_16x64,
    NUM_LUMA_PARTITIONS
};

// 4:2:0 chroma partition sizes. These enums are just a convenience for indexing into the
// chroma primitive arrays when instantiating templates. The function tables should always
// be indexed by the luma partition enum
enum Chroma420Partitions
{
    CHROMA_2x2,   CHROMA_4x4,   CHROMA_8x8,   CHROMA_16x16, CHROMA_32x32,
    CHROMA_4x2,   CHROMA_2x4,
    CHROMA_8x4,   CHROMA_4x8,
    CHROMA_16x8,  CHROMA_8x16,
    CHROMA_32x16, CHROMA_16x32,
    CHROMA_8x6,   CHROMA_6x8,   CHROMA_8x2,  CHROMA_2x8,
    CHROMA_16x12, CHROMA_12x16, CHROMA_16x4, CHROMA_4x16,
    CHROMA_32x24, CHROMA_24x32, CHROMA_32x8, CHROMA_8x32,
    NUM_CHROMA_PARTITIONS
};

enum Chroma422Partitions
{
    CHROMA422_2x4,   CHROMA422_4x8,   CHROMA422_8x16,  CHROMA422_16x32, CHROMA422_32x64,
    CHROMA422_4x4,   CHROMA422_2x8,
    CHROMA422_8x8,   CHROMA422_4x16,
    CHROMA422_16x16, CHROMA422_8x32,
    CHROMA422_32x32, CHROMA422_16x64,
    CHROMA422_8x12,  CHROMA422_6x16,  CHROMA422_8x4,   CHROMA422_2x16,
    CHROMA422_16x24, CHROMA422_12x32, CHROMA422_16x8,  CHROMA422_4x32,
    CHROMA422_32x48, CHROMA422_24x64, CHROMA422_32x16, CHROMA422_8x64,
    NUM_CHROMA_PARTITIONS422
};

enum SquareBlocks   // Routines can be indexed using log2n(width)-2
{
    BLOCK_4x4,
    BLOCK_8x8,
    BLOCK_16x16,
    BLOCK_32x32,
    BLOCK_64x64,
    NUM_SQUARE_BLOCKS
};

enum { NUM_TR_SIZE = 4 };

// NOTE: Not all DCT functions support dest stride
enum Dcts
{
    DST_4x4,
    DCT_4x4,
    DCT_8x8,
    DCT_16x16,
    DCT_32x32,
    NUM_DCTS
};

enum IDcts
{
    IDST_4x4,
    IDCT_4x4,
    IDCT_8x8,
    IDCT_16x16,
    IDCT_32x32,
    NUM_IDCTS
};

// Returns a LumaPartitions enum for the given size, always expected to return a valid enum
inline int partitionFromSizes(int width, int height)
{
    X265_CHECK(((width | height) & ~(4 | 8 | 16 | 32 | 64)) == 0, "Invalid block width/height\n");
    extern const uint8_t lumaPartitionMapTable[];
    int w = (width >> 2) - 1;
    int h = (height >> 2) - 1;
    int part = (int)lumaPartitionMapTable[(w << 4) + h];
    X265_CHECK(part != 255, "Invalid block width %d height %d\n", width, height);
    return part;
}

inline int partitionFromLog2Size(int log2Size)
{
    X265_CHECK(2 <= log2Size && log2Size <= 6, "Invalid block size\n");
    return log2Size - 2;
}

typedef int  (*pixelcmp_t)(pixel *fenc, intptr_t fencstride, pixel *fref, intptr_t frefstride); // fenc is aligned
typedef int  (*pixelcmp_ss_t)(int16_t *fenc, intptr_t fencstride, int16_t *fref, intptr_t frefstride);
typedef int  (*pixelcmp_sp_t)(int16_t *fenc, intptr_t fencstride, pixel *fref, intptr_t frefstride);
typedef int  (*pixel_ssd_s_t)(int16_t *fenc, intptr_t fencstride);
typedef void (*pixelcmp_x4_t)(pixel *fenc, pixel *fref0, pixel *fref1, pixel *fref2, pixel *fref3, intptr_t frefstride, int32_t *res);
typedef void (*pixelcmp_x3_t)(pixel *fenc, pixel *fref0, pixel *fref1, pixel *fref2, intptr_t frefstride, int32_t *res);
typedef void (*blockcpy_sp_t)(int bx, int by, int16_t *dst, intptr_t dstride, pixel *src, intptr_t sstride); // dst is aligned
typedef void (*blockcpy_sc_t)(int bx, int by, int16_t *dst, intptr_t dstride, uint8_t *src, intptr_t sstride); // dst is aligned
typedef void (*pixelsub_ps_t)(int bx, int by, int16_t *dst, intptr_t dstride, pixel *src0, pixel *src1, intptr_t sstride0, intptr_t sstride1);
typedef void (*pixelavg_pp_t)(pixel *dst, intptr_t dstride, pixel *src0, intptr_t sstride0, pixel *src1, intptr_t sstride1, int weight);
typedef void (*blockfill_s_t)(int16_t *dst, intptr_t dstride, int16_t val);

typedef void (*intra_pred_t)(pixel* dst, intptr_t dstStride, pixel *refLeft, pixel *refAbove, int dirMode, int bFilter);
typedef void (*intra_allangs_t)(pixel *dst, pixel *above0, pixel *left0, pixel *above1, pixel *left1, int bLuma);

typedef void (*cvt16to32_shl_t)(int32_t *dst, int16_t *src, intptr_t, int, int);
typedef void (*cvt16to32_shr_t)(int32_t *dst, int16_t *src, intptr_t, int, int);
typedef void (*cvt32to16_shr_t)(int16_t *dst, int32_t *src, intptr_t, int, int);
typedef void (*cvt32to16_shl_t)(int16_t *dst, int32_t *src, intptr_t, int);
typedef uint32_t (*copy_cnt_t)(int16_t* coeff, int16_t* residual, intptr_t stride);
typedef void (*copy_shr_t)(int16_t *dst, int16_t *src, intptr_t stride, int shift, int size);
typedef void (*copy_shl_t)(int16_t *dst, int16_t *src, intptr_t stride, int shift);

typedef void (*dct_t)(int16_t *src, int32_t *dst, intptr_t stride);
typedef void (*idct_t)(int32_t *src, int16_t *dst, intptr_t stride);
typedef void (*denoiseDct_t)(int32_t* dctCoef, uint32_t* resSum, uint16_t* offset, int numCoeff);

typedef void (*calcresidual_t)(pixel *fenc, pixel *pred, int16_t *residual, intptr_t stride);
typedef void (*calcrecon_t)(pixel* pred, int16_t* residual, int16_t* reconqt, pixel *reconipred, int stride, int strideqt, int strideipred);
typedef void (*transpose_t)(pixel* dst, pixel* src, intptr_t stride);
typedef uint32_t (*quant_t)(int32_t *coef, int32_t *quantCoeff, int32_t *deltaU, int16_t *qCoef, int qBits, int add, int numCoeff);
typedef uint32_t (*nquant_t)(int32_t *coef, int32_t *quantCoeff, int16_t *qCoef, int qBits, int add, int numCoeff);
typedef void (*dequant_scaling_t)(const int16_t* src, const int32_t *dequantCoef, int32_t* dst, int num, int mcqp_miper, int shift);
typedef void (*dequant_normal_t)(const int16_t* quantCoef, int32_t* coef, int num, int scale, int shift);
typedef int  (*count_nonzero_t)(const int16_t *quantCoeff, int numCoeff);

typedef void (*weightp_pp_t)(pixel *src, pixel *dst, intptr_t stride, int width, int height, int w0, int round, int shift, int offset);
typedef void (*weightp_sp_t)(int16_t *src, pixel *dst, intptr_t srcStride, intptr_t dstStride, int width, int height, int w0, int round, int shift, int offset);
typedef void (*scale_t)(pixel *dst, pixel *src, intptr_t stride);
typedef void (*downscale_t)(pixel *src0, pixel *dstf, pixel *dsth, pixel *dstv, pixel *dstc,
                            intptr_t src_stride, intptr_t dst_stride, int width, int height);
typedef void (*extendCURowBorder_t)(pixel* txt, intptr_t stride, int width, int height, int marginX);
typedef void (*ssim_4x4x2_core_t)(const pixel *pix1, intptr_t stride1, const pixel *pix2, intptr_t stride2, int sums[2][4]);
typedef float (*ssim_end4_t)(int sum0[5][4], int sum1[5][4], int width);
typedef uint64_t (*var_t)(pixel *pix, intptr_t stride);
typedef void (*plane_copy_deinterleave_t)(pixel *dstu, intptr_t dstuStride, pixel *dstv, intptr_t dstvStride, pixel *src, intptr_t srcStride, int w, int h);

typedef void (*filter_pp_t) (pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx);
typedef void (*filter_hps_t) (pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt);
typedef void (*filter_ps_t) (pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx);
typedef void (*filter_sp_t) (int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx);
typedef void (*filter_ss_t) (int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx);
typedef void (*filter_hv_pp_t) (pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int idxX, int idxY);
typedef void (*filter_p2s_t)(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height);

typedef void (*copy_pp_t)(pixel *dst, intptr_t dstride, pixel *src, intptr_t sstride); // dst is aligned
typedef void (*copy_sp_t)(pixel *dst, intptr_t dstStride, int16_t *src, intptr_t srcStride);
typedef void (*copy_ps_t)(int16_t *dst, intptr_t dstStride, pixel *src, intptr_t srcStride);
typedef void (*copy_ss_t)(int16_t *dst, intptr_t dstStride, int16_t *src, intptr_t srcStride);

typedef void (*pixel_sub_ps_t)(int16_t *dst, intptr_t dstride, pixel *src0, pixel *src1, intptr_t sstride0, intptr_t sstride1);
typedef void (*pixel_add_ps_t)(pixel *a, intptr_t dstride, pixel *b0, int16_t *b1, intptr_t sstride0, intptr_t sstride1);
typedef void (*addAvg_t)(int16_t* src0, int16_t* src1, pixel* dst, intptr_t src0Stride, intptr_t src1Stride, intptr_t dstStride);

typedef void (*saoCuOrgE0_t)(pixel * rec, int8_t * offsetEo, int width, int8_t signLeft);
typedef void (*planecopy_cp_t) (uint8_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int width, int height, int shift);
typedef void (*planecopy_sp_t) (uint16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int width, int height, int shift, uint16_t mask);

typedef void (*cutree_propagate_cost) (int *dst, uint16_t *propagateIn, int32_t *intraCosts, uint16_t *interCosts, int32_t *invQscales, double *fpsFactor, int len);

/* Define a structure containing function pointers to optimized encoder
 * primitives.  Each pointer can reference either an assembly routine,
 * a vectorized primitive, or a C function. */
struct EncoderPrimitives
{
    pixelcmp_t      sad[NUM_LUMA_PARTITIONS];        // Sum of Differences for each size
    pixelcmp_x3_t   sad_x3[NUM_LUMA_PARTITIONS];     // Sum of Differences 3x for each size
    pixelcmp_x4_t   sad_x4[NUM_LUMA_PARTITIONS];     // Sum of Differences 4x for each size
    pixelcmp_t      sse_pp[NUM_LUMA_PARTITIONS];     // Sum of Square Error (pixel, pixel) fenc alignment not assumed
    pixelcmp_ss_t   sse_ss[NUM_LUMA_PARTITIONS];     // Sum of Square Error (short, short) fenc alignment not assumed
    pixelcmp_sp_t   sse_sp[NUM_LUMA_PARTITIONS];     // Sum of Square Error (short, pixel) fenc alignment not assumed
    pixel_ssd_s_t   ssd_s[NUM_SQUARE_BLOCKS - 1];    // Sum of Square Error (short) fenc alignment not assumed
    pixelcmp_t      satd[NUM_LUMA_PARTITIONS];       // Sum of Transformed differences (HADAMARD)
    pixelcmp_t      sa8d_inter[NUM_LUMA_PARTITIONS]; // sa8d primitives for motion search partitions
    pixelcmp_t      sa8d[NUM_SQUARE_BLOCKS];         // sa8d primitives for square intra blocks
    pixelcmp_t      psy_cost_pp[NUM_SQUARE_BLOCKS];     // difference in AC energy between two blocks
    pixelcmp_ss_t   psy_cost_ss[NUM_SQUARE_BLOCKS];

    blockfill_s_t   blockfill_s[NUM_SQUARE_BLOCKS];  // block fill with value
    cvt16to32_shl_t cvt16to32_shl;
    cvt16to32_shr_t cvt16to32_shr[NUM_SQUARE_BLOCKS - 1];
    cvt32to16_shr_t cvt32to16_shr;
    cvt32to16_shl_t cvt32to16_shl[NUM_SQUARE_BLOCKS - 1];
    copy_cnt_t      copy_cnt[NUM_SQUARE_BLOCKS - 1];
    copy_shr_t      copy_shr;
    copy_shl_t      copy_shl[NUM_SQUARE_BLOCKS - 1];

    copy_pp_t       luma_copy_pp[NUM_LUMA_PARTITIONS];
    copy_sp_t       luma_copy_sp[NUM_LUMA_PARTITIONS];
    copy_ps_t       luma_copy_ps[NUM_LUMA_PARTITIONS];
    copy_ss_t       luma_copy_ss[NUM_LUMA_PARTITIONS];
    pixel_sub_ps_t  luma_sub_ps[NUM_SQUARE_BLOCKS];
    pixel_add_ps_t  luma_add_ps[NUM_SQUARE_BLOCKS];
    copy_pp_t       square_copy_pp[NUM_SQUARE_BLOCKS];
    copy_sp_t       square_copy_sp[NUM_SQUARE_BLOCKS];
    copy_ps_t       square_copy_ps[NUM_SQUARE_BLOCKS];
    copy_ss_t       square_copy_ss[NUM_SQUARE_BLOCKS];

    filter_pp_t     luma_hpp[NUM_LUMA_PARTITIONS];
    filter_hps_t    luma_hps[NUM_LUMA_PARTITIONS];
    filter_pp_t     luma_vpp[NUM_LUMA_PARTITIONS];
    filter_ps_t     luma_vps[NUM_LUMA_PARTITIONS];
    filter_sp_t     luma_vsp[NUM_LUMA_PARTITIONS];
    filter_ss_t     luma_vss[NUM_LUMA_PARTITIONS];
    filter_hv_pp_t  luma_hvpp[NUM_LUMA_PARTITIONS];
    filter_p2s_t    luma_p2s;
    filter_p2s_t    chroma_p2s[X265_CSP_COUNT];

    weightp_sp_t    weight_sp;
    weightp_pp_t    weight_pp;
    pixelavg_pp_t   pixelavg_pp[NUM_LUMA_PARTITIONS];
    addAvg_t        luma_addAvg[NUM_LUMA_PARTITIONS];

    intra_pred_t    intra_pred[NUM_INTRA_MODE][NUM_TR_SIZE];
    intra_allangs_t intra_pred_allangs[NUM_TR_SIZE];
    scale_t         scale1D_128to64;
    scale_t         scale2D_64to32;

    dct_t           dct[NUM_DCTS];
    idct_t          idct[NUM_IDCTS];
    quant_t         quant;
    nquant_t        nquant;
    dequant_scaling_t dequant_scaling;
    dequant_normal_t dequant_normal;
    count_nonzero_t count_nonzero;
    denoiseDct_t    denoiseDct;

    calcresidual_t  calcresidual[NUM_SQUARE_BLOCKS];
    transpose_t     transpose[NUM_SQUARE_BLOCKS];

    var_t           var[NUM_SQUARE_BLOCKS];
    ssim_4x4x2_core_t ssim_4x4x2_core;
    ssim_end4_t     ssim_end_4;

    downscale_t     frame_init_lowres_core;
    plane_copy_deinterleave_t plane_copy_deinterleave_c;
    extendCURowBorder_t extendRowBorder;
    // sao primitives
    saoCuOrgE0_t      saoCuOrgE0;
    planecopy_cp_t    planecopy_cp;
    planecopy_sp_t    planecopy_sp;

    cutree_propagate_cost    propagateCost;

    struct
    {
        filter_pp_t     filter_vpp[NUM_LUMA_PARTITIONS];
        filter_ps_t     filter_vps[NUM_LUMA_PARTITIONS];
        filter_sp_t     filter_vsp[NUM_LUMA_PARTITIONS];
        filter_ss_t     filter_vss[NUM_LUMA_PARTITIONS];
        filter_pp_t     filter_hpp[NUM_LUMA_PARTITIONS];
        filter_hps_t    filter_hps[NUM_LUMA_PARTITIONS];
        addAvg_t        addAvg[NUM_LUMA_PARTITIONS];
        copy_pp_t       copy_pp[NUM_LUMA_PARTITIONS];
        copy_sp_t       copy_sp[NUM_LUMA_PARTITIONS];
        copy_ps_t       copy_ps[NUM_LUMA_PARTITIONS];
        copy_ss_t       copy_ss[NUM_LUMA_PARTITIONS];
        pixel_sub_ps_t  sub_ps[NUM_SQUARE_BLOCKS];
        pixel_add_ps_t  add_ps[NUM_SQUARE_BLOCKS];
    } chroma[4]; // X265_CSP_COUNT - do not want to include x265.h here
};

void extendPicBorder(pixel* recon, intptr_t stride, int width, int height, int marginX, int marginY);

/* This copy of the table is what gets used by the encoder.
 * It must be initialized before the encoder begins. */
extern EncoderPrimitives primitives;

void Setup_C_Primitives(EncoderPrimitives &p);
void Setup_Instrinsic_Primitives(EncoderPrimitives &p, int cpuMask);
void Setup_Assembly_Primitives(EncoderPrimitives &p, int cpuMask);
void Setup_Alias_Primitives(EncoderPrimitives &p);
}

#endif // ifndef X265_PRIMITIVES_H
