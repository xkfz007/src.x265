/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *          Praveen Kumar Tiwari <praveen@multicorewareinc.com>
 *          Min Chen <chenm003@163.com> <min.chen@multicorewareinc.com>
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

#include "common.h"
#include "primitives.h"
#include "x265.h"
#include "cpu.h"

extern "C" {
#include "pixel.h"
#include "pixel-util.h"
#include "mc.h"
#include "ipfilter8.h"
#include "loopfilter.h"
#include "blockcopy8.h"
#include "intrapred.h"
#include "dct8.h"
}

#define ALL_LUMA_CU_TYPED(prim, fncdef, fname, cpu) \
    p.cu[BLOCK_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.cu[BLOCK_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.cu[BLOCK_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.cu[BLOCK_64x64].prim = fncdef x265_ ## fname ## _64x64_ ## cpu
#define ALL_LUMA_CU_TYPED_S(prim, fncdef, fname, cpu) \
    p.cu[BLOCK_8x8].prim   = fncdef x265_ ## fname ## 8_ ## cpu; \
    p.cu[BLOCK_16x16].prim = fncdef x265_ ## fname ## 16_ ## cpu; \
    p.cu[BLOCK_32x32].prim = fncdef x265_ ## fname ## 32_ ## cpu; \
    p.cu[BLOCK_64x64].prim = fncdef x265_ ## fname ## 64_ ## cpu
#define ALL_LUMA_TU_TYPED(prim, fncdef, fname, cpu) \
    p.cu[BLOCK_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.cu[BLOCK_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.cu[BLOCK_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.cu[BLOCK_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu
#define ALL_LUMA_TU_TYPED_S(prim, fncdef, fname, cpu) \
    p.cu[BLOCK_4x4].prim   = fncdef x265_ ## fname ## 4_ ## cpu; \
    p.cu[BLOCK_8x8].prim   = fncdef x265_ ## fname ## 8_ ## cpu; \
    p.cu[BLOCK_16x16].prim = fncdef x265_ ## fname ## 16_ ## cpu; \
    p.cu[BLOCK_32x32].prim = fncdef x265_ ## fname ## 32_ ## cpu
#define ALL_LUMA_BLOCKS_TYPED(prim, fncdef, fname, cpu) \
    p.cu[BLOCK_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.cu[BLOCK_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.cu[BLOCK_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.cu[BLOCK_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.cu[BLOCK_64x64].prim = fncdef x265_ ## fname ## _64x64_ ## cpu;
#define ALL_LUMA_CU(prim, fname, cpu)      ALL_LUMA_CU_TYPED(prim, , fname, cpu)
#define ALL_LUMA_CU_S(prim, fname, cpu)    ALL_LUMA_CU_TYPED_S(prim, , fname, cpu)
#define ALL_LUMA_TU(prim, fname, cpu)      ALL_LUMA_TU_TYPED(prim, , fname, cpu)
#define ALL_LUMA_BLOCKS(prim, fname, cpu)  ALL_LUMA_BLOCKS_TYPED(prim, , fname, cpu)
#define ALL_LUMA_TU_S(prim, fname, cpu)    ALL_LUMA_TU_TYPED_S(prim, , fname, cpu)

#define ALL_LUMA_PU_TYPED(prim, fncdef, fname, cpu) \
    p.pu[LUMA_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.pu[LUMA_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.pu[LUMA_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.pu[LUMA_64x64].prim = fncdef x265_ ## fname ## _64x64_ ## cpu; \
    p.pu[LUMA_8x4].prim   = fncdef x265_ ## fname ## _8x4_ ## cpu; \
    p.pu[LUMA_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.pu[LUMA_16x8].prim  = fncdef x265_ ## fname ## _16x8_ ## cpu; \
    p.pu[LUMA_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.pu[LUMA_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.pu[LUMA_32x16].prim = fncdef x265_ ## fname ## _32x16_ ## cpu; \
    p.pu[LUMA_64x32].prim = fncdef x265_ ## fname ## _64x32_ ## cpu; \
    p.pu[LUMA_32x64].prim = fncdef x265_ ## fname ## _32x64_ ## cpu; \
    p.pu[LUMA_16x12].prim = fncdef x265_ ## fname ## _16x12_ ## cpu; \
    p.pu[LUMA_12x16].prim = fncdef x265_ ## fname ## _12x16_ ## cpu; \
    p.pu[LUMA_16x4].prim  = fncdef x265_ ## fname ## _16x4_ ## cpu; \
    p.pu[LUMA_4x16].prim  = fncdef x265_ ## fname ## _4x16_ ## cpu; \
    p.pu[LUMA_32x24].prim = fncdef x265_ ## fname ## _32x24_ ## cpu; \
    p.pu[LUMA_24x32].prim = fncdef x265_ ## fname ## _24x32_ ## cpu; \
    p.pu[LUMA_32x8].prim  = fncdef x265_ ## fname ## _32x8_ ## cpu; \
    p.pu[LUMA_8x32].prim  = fncdef x265_ ## fname ## _8x32_ ## cpu; \
    p.pu[LUMA_64x48].prim = fncdef x265_ ## fname ## _64x48_ ## cpu; \
    p.pu[LUMA_48x64].prim = fncdef x265_ ## fname ## _48x64_ ## cpu; \
    p.pu[LUMA_64x16].prim = fncdef x265_ ## fname ## _64x16_ ## cpu; \
    p.pu[LUMA_16x64].prim = fncdef x265_ ## fname ## _16x64_ ## cpu
#define ALL_LUMA_PU(prim, fname, cpu) ALL_LUMA_PU_TYPED(prim, , fname, cpu)

#define ALL_LUMA_PU_T(prim, fname) \
    p.pu[LUMA_8x8].prim   = fname<LUMA_8x8>; \
    p.pu[LUMA_16x16].prim = fname<LUMA_16x16>; \
    p.pu[LUMA_32x32].prim = fname<LUMA_32x32>; \
    p.pu[LUMA_64x64].prim = fname<LUMA_64x64>; \
    p.pu[LUMA_8x4].prim   = fname<LUMA_8x4>; \
    p.pu[LUMA_4x8].prim   = fname<LUMA_4x8>; \
    p.pu[LUMA_16x8].prim  = fname<LUMA_16x8>; \
    p.pu[LUMA_8x16].prim  = fname<LUMA_8x16>; \
    p.pu[LUMA_16x32].prim = fname<LUMA_16x32>; \
    p.pu[LUMA_32x16].prim = fname<LUMA_32x16>; \
    p.pu[LUMA_64x32].prim = fname<LUMA_64x32>; \
    p.pu[LUMA_32x64].prim = fname<LUMA_32x64>; \
    p.pu[LUMA_16x12].prim = fname<LUMA_16x12>; \
    p.pu[LUMA_12x16].prim = fname<LUMA_12x16>; \
    p.pu[LUMA_16x4].prim  = fname<LUMA_16x4>; \
    p.pu[LUMA_4x16].prim  = fname<LUMA_4x16>; \
    p.pu[LUMA_32x24].prim = fname<LUMA_32x24>; \
    p.pu[LUMA_24x32].prim = fname<LUMA_24x32>; \
    p.pu[LUMA_32x8].prim  = fname<LUMA_32x8>; \
    p.pu[LUMA_8x32].prim  = fname<LUMA_8x32>; \
    p.pu[LUMA_64x48].prim = fname<LUMA_64x48>; \
    p.pu[LUMA_48x64].prim = fname<LUMA_48x64>; \
    p.pu[LUMA_64x16].prim = fname<LUMA_64x16>; \
    p.pu[LUMA_16x64].prim = fname<LUMA_16x64>

#define ALL_CHROMA_420_CU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu
#define ALL_CHROMA_420_CU_TYPED_S(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_4x4].prim   = fncdef x265_ ## fname ## _4_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_8x8].prim   = fncdef x265_ ## fname ## _8_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].prim = fncdef x265_ ## fname ## _16_ ## cpu; \
    p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].prim = fncdef x265_ ## fname ## _32_ ## cpu
#define ALL_CHROMA_420_CU(prim, fname, cpu) ALL_CHROMA_420_CU_TYPED(prim, , fname, cpu)
#define ALL_CHROMA_420_CU_S(prim, fname, cpu) ALL_CHROMA_420_CU_TYPED_S(prim, , fname, cpu)

#define ALL_CHROMA_420_PU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].prim   = fncdef x265_ ## fname ## _4x2_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].prim   = fncdef x265_ ## fname ## _2x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].prim   = fncdef x265_ ## fname ## _8x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].prim  = fncdef x265_ ## fname ## _16x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].prim = fncdef x265_ ## fname ## _32x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].prim   = fncdef x265_ ## fname ## _8x6_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].prim   = fncdef x265_ ## fname ## _6x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].prim   = fncdef x265_ ## fname ## _8x2_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].prim   = fncdef x265_ ## fname ## _2x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].prim = fncdef x265_ ## fname ## _16x12_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].prim = fncdef x265_ ## fname ## _12x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].prim  = fncdef x265_ ## fname ## _16x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].prim  = fncdef x265_ ## fname ## _4x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].prim = fncdef x265_ ## fname ## _32x24_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].prim = fncdef x265_ ## fname ## _24x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].prim  = fncdef x265_ ## fname ## _32x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].prim  = fncdef x265_ ## fname ## _8x32_ ## cpu
#define ALL_CHROMA_420_PU(prim, fname, cpu) ALL_CHROMA_420_PU_TYPED(prim, , fname, cpu)

#define ALL_CHROMA_420_4x4_PU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].prim   = fncdef x265_ ## fname ## _8x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].prim  = fncdef x265_ ## fname ## _16x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].prim = fncdef x265_ ## fname ## _32x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].prim = fncdef x265_ ## fname ## _16x12_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].prim = fncdef x265_ ## fname ## _12x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].prim  = fncdef x265_ ## fname ## _16x4_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].prim  = fncdef x265_ ## fname ## _4x16_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].prim = fncdef x265_ ## fname ## _32x24_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].prim = fncdef x265_ ## fname ## _24x32_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].prim  = fncdef x265_ ## fname ## _32x8_ ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].prim  = fncdef x265_ ## fname ## _8x32_ ## cpu
#define ALL_CHROMA_420_4x4_PU(prim, fname, cpu) ALL_CHROMA_420_4x4_PU_TYPED(prim, , fname, cpu)

#define ALL_CHROMA_422_CU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].prim = fncdef x265_ ## fname ## _32x64_ ## cpu
#define ALL_CHROMA_422_CU(prim, fname, cpu) ALL_CHROMA_422_CU_TYPED(prim, , fname, cpu)

#define ALL_CHROMA_422_PU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].prim = fncdef x265_ ## fname ## _32x64_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].prim   = fncdef x265_ ## fname ## _2x8_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].prim  = fncdef x265_ ## fname ## _4x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].prim  = fncdef x265_ ## fname ## _8x32_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].prim = fncdef x265_ ## fname ## _16x64_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].prim  = fncdef x265_ ## fname ## _8x12_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].prim  = fncdef x265_ ## fname ## _6x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].prim   = fncdef x265_ ## fname ## _8x4_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].prim  = fncdef x265_ ## fname ## _2x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].prim = fncdef x265_ ## fname ## _16x24_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].prim = fncdef x265_ ## fname ## _12x32_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].prim  = fncdef x265_ ## fname ## _16x8_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].prim  = fncdef x265_ ## fname ## _4x32_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].prim = fncdef x265_ ## fname ## _32x48_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].prim = fncdef x265_ ## fname ## _24x64_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].prim = fncdef x265_ ## fname ## _32x16_ ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].prim  = fncdef x265_ ## fname ## _8x64_ ## cpu
#define ALL_CHROMA_422_PU(prim, fname, cpu) ALL_CHROMA_422_PU_TYPED(prim, , fname, cpu)

#define ALL_CHROMA_444_PU_TYPED(prim, fncdef, fname, cpu) \
    p.chroma[X265_CSP_I444].pu[LUMA_4x4].prim   = fncdef x265_ ## fname ## _4x4_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_8x8].prim   = fncdef x265_ ## fname ## _8x8_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x16].prim = fncdef x265_ ## fname ## _16x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_32x32].prim = fncdef x265_ ## fname ## _32x32_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_64x64].prim = fncdef x265_ ## fname ## _64x64_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_8x4].prim   = fncdef x265_ ## fname ## _8x4_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_4x8].prim   = fncdef x265_ ## fname ## _4x8_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x8].prim  = fncdef x265_ ## fname ## _16x8_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_8x16].prim  = fncdef x265_ ## fname ## _8x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x32].prim = fncdef x265_ ## fname ## _16x32_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_32x16].prim = fncdef x265_ ## fname ## _32x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_64x32].prim = fncdef x265_ ## fname ## _64x32_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_32x64].prim = fncdef x265_ ## fname ## _32x64_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x12].prim = fncdef x265_ ## fname ## _16x12_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_12x16].prim = fncdef x265_ ## fname ## _12x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x4].prim  = fncdef x265_ ## fname ## _16x4_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_4x16].prim  = fncdef x265_ ## fname ## _4x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_32x24].prim = fncdef x265_ ## fname ## _32x24_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_24x32].prim = fncdef x265_ ## fname ## _24x32_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_32x8].prim  = fncdef x265_ ## fname ## _32x8_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_8x32].prim  = fncdef x265_ ## fname ## _8x32_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_64x48].prim = fncdef x265_ ## fname ## _64x48_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_48x64].prim = fncdef x265_ ## fname ## _48x64_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_64x16].prim = fncdef x265_ ## fname ## _64x16_ ## cpu; \
    p.chroma[X265_CSP_I444].pu[LUMA_16x64].prim = fncdef x265_ ## fname ## _16x64_ ## cpu
#define ALL_CHROMA_444_PU(prim, fname, cpu) ALL_CHROMA_444_PU_TYPED(prim, , fname, cpu)

#define AVC_LUMA_PU(name, cpu) \
    p.pu[LUMA_16x16].name = x265_pixel_ ## name ## _16x16_ ## cpu; \
    p.pu[LUMA_16x8].name  = x265_pixel_ ## name ## _16x8_ ## cpu; \
    p.pu[LUMA_8x16].name  = x265_pixel_ ## name ## _8x16_ ## cpu; \
    p.pu[LUMA_8x8].name   = x265_pixel_ ## name ## _8x8_ ## cpu; \
    p.pu[LUMA_8x4].name   = x265_pixel_ ## name ## _8x4_ ## cpu; \
    p.pu[LUMA_4x8].name   = x265_pixel_ ## name ## _4x8_ ## cpu; \
    p.pu[LUMA_4x4].name   = x265_pixel_ ## name ## _4x4_ ## cpu; \
    p.pu[LUMA_4x16].name  = x265_pixel_ ## name ## _4x16_ ## cpu

#define HEVC_SAD(cpu) \
    p.pu[LUMA_8x32].sad  = x265_pixel_sad_8x32_ ## cpu; \
    p.pu[LUMA_16x4].sad  = x265_pixel_sad_16x4_ ## cpu; \
    p.pu[LUMA_16x12].sad = x265_pixel_sad_16x12_ ## cpu; \
    p.pu[LUMA_16x32].sad = x265_pixel_sad_16x32_ ## cpu; \
    p.pu[LUMA_16x64].sad = x265_pixel_sad_16x64_ ## cpu; \
    p.pu[LUMA_32x8].sad  = x265_pixel_sad_32x8_ ## cpu; \
    p.pu[LUMA_32x16].sad = x265_pixel_sad_32x16_ ## cpu; \
    p.pu[LUMA_32x24].sad = x265_pixel_sad_32x24_ ## cpu; \
    p.pu[LUMA_32x32].sad = x265_pixel_sad_32x32_ ## cpu; \
    p.pu[LUMA_32x64].sad = x265_pixel_sad_32x64_ ## cpu; \
    p.pu[LUMA_64x16].sad = x265_pixel_sad_64x16_ ## cpu; \
    p.pu[LUMA_64x32].sad = x265_pixel_sad_64x32_ ## cpu; \
    p.pu[LUMA_64x48].sad = x265_pixel_sad_64x48_ ## cpu; \
    p.pu[LUMA_64x64].sad = x265_pixel_sad_64x64_ ## cpu; \
    p.pu[LUMA_48x64].sad = x265_pixel_sad_48x64_ ## cpu; \
    p.pu[LUMA_24x32].sad = x265_pixel_sad_24x32_ ## cpu; \
    p.pu[LUMA_12x16].sad = x265_pixel_sad_12x16_ ## cpu

#define HEVC_SAD_X3(cpu) \
    p.pu[LUMA_16x8].sad_x3  = x265_pixel_sad_x3_16x8_ ## cpu; \
    p.pu[LUMA_16x12].sad_x3 = x265_pixel_sad_x3_16x12_ ## cpu; \
    p.pu[LUMA_16x16].sad_x3 = x265_pixel_sad_x3_16x16_ ## cpu; \
    p.pu[LUMA_16x32].sad_x3 = x265_pixel_sad_x3_16x32_ ## cpu; \
    p.pu[LUMA_16x64].sad_x3 = x265_pixel_sad_x3_16x64_ ## cpu; \
    p.pu[LUMA_32x8].sad_x3  = x265_pixel_sad_x3_32x8_ ## cpu; \
    p.pu[LUMA_32x16].sad_x3 = x265_pixel_sad_x3_32x16_ ## cpu; \
    p.pu[LUMA_32x24].sad_x3 = x265_pixel_sad_x3_32x24_ ## cpu; \
    p.pu[LUMA_32x32].sad_x3 = x265_pixel_sad_x3_32x32_ ## cpu; \
    p.pu[LUMA_32x64].sad_x3 = x265_pixel_sad_x3_32x64_ ## cpu; \
    p.pu[LUMA_24x32].sad_x3 = x265_pixel_sad_x3_24x32_ ## cpu; \
    p.pu[LUMA_48x64].sad_x3 = x265_pixel_sad_x3_48x64_ ## cpu; \
    p.pu[LUMA_64x16].sad_x3 = x265_pixel_sad_x3_64x16_ ## cpu; \
    p.pu[LUMA_64x32].sad_x3 = x265_pixel_sad_x3_64x32_ ## cpu; \
    p.pu[LUMA_64x48].sad_x3 = x265_pixel_sad_x3_64x48_ ## cpu; \
    p.pu[LUMA_64x64].sad_x3 = x265_pixel_sad_x3_64x64_ ## cpu

#define HEVC_SAD_X4(cpu) \
    p.pu[LUMA_16x8].sad_x4  = x265_pixel_sad_x4_16x8_ ## cpu; \
    p.pu[LUMA_16x12].sad_x4 = x265_pixel_sad_x4_16x12_ ## cpu; \
    p.pu[LUMA_16x16].sad_x4 = x265_pixel_sad_x4_16x16_ ## cpu; \
    p.pu[LUMA_16x32].sad_x4 = x265_pixel_sad_x4_16x32_ ## cpu; \
    p.pu[LUMA_16x64].sad_x4 = x265_pixel_sad_x4_16x64_ ## cpu; \
    p.pu[LUMA_32x8].sad_x4  = x265_pixel_sad_x4_32x8_ ## cpu; \
    p.pu[LUMA_32x16].sad_x4 = x265_pixel_sad_x4_32x16_ ## cpu; \
    p.pu[LUMA_32x24].sad_x4 = x265_pixel_sad_x4_32x24_ ## cpu; \
    p.pu[LUMA_32x32].sad_x4 = x265_pixel_sad_x4_32x32_ ## cpu; \
    p.pu[LUMA_32x64].sad_x4 = x265_pixel_sad_x4_32x64_ ## cpu; \
    p.pu[LUMA_24x32].sad_x4 = x265_pixel_sad_x4_24x32_ ## cpu; \
    p.pu[LUMA_48x64].sad_x4 = x265_pixel_sad_x4_48x64_ ## cpu; \
    p.pu[LUMA_64x16].sad_x4 = x265_pixel_sad_x4_64x16_ ## cpu; \
    p.pu[LUMA_64x32].sad_x4 = x265_pixel_sad_x4_64x32_ ## cpu; \
    p.pu[LUMA_64x48].sad_x4 = x265_pixel_sad_x4_64x48_ ## cpu; \
    p.pu[LUMA_64x64].sad_x4 = x265_pixel_sad_x4_64x64_ ## cpu

#define ASSIGN_SSE_PP(cpu) \
    p.cu[BLOCK_8x8].sse_pp   = x265_pixel_ssd_8x8_ ## cpu; \
    p.cu[BLOCK_16x16].sse_pp = x265_pixel_ssd_16x16_ ## cpu; \
    p.cu[BLOCK_32x32].sse_pp = x265_pixel_ssd_32x32_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_8x16].sse_pp = x265_pixel_ssd_8x16_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].sse_pp = x265_pixel_ssd_16x32_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].sse_pp = x265_pixel_ssd_32x64_ ## cpu;

#define ASSIGN_SSE_SS(cpu) ALL_LUMA_BLOCKS(sse_ss, pixel_ssd_ss, cpu)

#define ASSIGN_SA8D(cpu) \
    ALL_LUMA_CU(sa8d, pixel_sa8d, cpu); \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_8x16].sa8d = x265_pixel_sa8d_8x16_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].sa8d = x265_pixel_sa8d_16x32_ ## cpu; \
    p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].sa8d = x265_pixel_sa8d_32x64_ ## cpu

#define PIXEL_AVG(cpu) \
    p.pu[LUMA_64x64].pixelavg_pp = x265_pixel_avg_64x64_ ## cpu; \
    p.pu[LUMA_64x48].pixelavg_pp = x265_pixel_avg_64x48_ ## cpu; \
    p.pu[LUMA_64x32].pixelavg_pp = x265_pixel_avg_64x32_ ## cpu; \
    p.pu[LUMA_64x16].pixelavg_pp = x265_pixel_avg_64x16_ ## cpu; \
    p.pu[LUMA_48x64].pixelavg_pp = x265_pixel_avg_48x64_ ## cpu; \
    p.pu[LUMA_32x64].pixelavg_pp = x265_pixel_avg_32x64_ ## cpu; \
    p.pu[LUMA_32x32].pixelavg_pp = x265_pixel_avg_32x32_ ## cpu; \
    p.pu[LUMA_32x24].pixelavg_pp = x265_pixel_avg_32x24_ ## cpu; \
    p.pu[LUMA_32x16].pixelavg_pp = x265_pixel_avg_32x16_ ## cpu; \
    p.pu[LUMA_32x8].pixelavg_pp  = x265_pixel_avg_32x8_ ## cpu; \
    p.pu[LUMA_24x32].pixelavg_pp = x265_pixel_avg_24x32_ ## cpu; \
    p.pu[LUMA_16x64].pixelavg_pp = x265_pixel_avg_16x64_ ## cpu; \
    p.pu[LUMA_16x32].pixelavg_pp = x265_pixel_avg_16x32_ ## cpu; \
    p.pu[LUMA_16x16].pixelavg_pp = x265_pixel_avg_16x16_ ## cpu; \
    p.pu[LUMA_16x12].pixelavg_pp = x265_pixel_avg_16x12_ ## cpu; \
    p.pu[LUMA_16x8].pixelavg_pp  = x265_pixel_avg_16x8_ ## cpu; \
    p.pu[LUMA_16x4].pixelavg_pp  = x265_pixel_avg_16x4_ ## cpu; \
    p.pu[LUMA_12x16].pixelavg_pp = x265_pixel_avg_12x16_ ## cpu; \
    p.pu[LUMA_8x32].pixelavg_pp  = x265_pixel_avg_8x32_ ## cpu; \
    p.pu[LUMA_8x16].pixelavg_pp  = x265_pixel_avg_8x16_ ## cpu; \
    p.pu[LUMA_8x8].pixelavg_pp   = x265_pixel_avg_8x8_ ## cpu; \
    p.pu[LUMA_8x4].pixelavg_pp   = x265_pixel_avg_8x4_ ## cpu;

#define PIXEL_AVG_W4(cpu) \
    p.pu[LUMA_4x4].pixelavg_pp  = x265_pixel_avg_4x4_ ## cpu; \
    p.pu[LUMA_4x8].pixelavg_pp  = x265_pixel_avg_4x8_ ## cpu; \
    p.pu[LUMA_4x16].pixelavg_pp = x265_pixel_avg_4x16_ ## cpu;

#define CHROMA_420_FILTERS(cpu) \
    ALL_CHROMA_420_PU(filter_hpp, interp_4tap_horiz_pp, cpu); \
    ALL_CHROMA_420_PU(filter_hps, interp_4tap_horiz_ps, cpu); \
    ALL_CHROMA_420_PU(filter_vpp, interp_4tap_vert_pp, cpu); \
    ALL_CHROMA_420_PU(filter_vps, interp_4tap_vert_ps, cpu);

#define CHROMA_422_FILTERS(cpu) \
    ALL_CHROMA_422_PU(filter_hpp, interp_4tap_horiz_pp, cpu); \
    ALL_CHROMA_422_PU(filter_hps, interp_4tap_horiz_ps, cpu); \
    ALL_CHROMA_422_PU(filter_vpp, interp_4tap_vert_pp, cpu); \
    ALL_CHROMA_422_PU(filter_vps, interp_4tap_vert_ps, cpu);

#define CHROMA_444_FILTERS(cpu) \
    ALL_CHROMA_444_PU(filter_hpp, interp_4tap_horiz_pp, cpu); \
    ALL_CHROMA_444_PU(filter_hps, interp_4tap_horiz_ps, cpu); \
    ALL_CHROMA_444_PU(filter_vpp, interp_4tap_vert_pp, cpu); \
    ALL_CHROMA_444_PU(filter_vps, interp_4tap_vert_ps, cpu);

#define SETUP_CHROMA_420_VSP_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vsp = x265_interp_4tap_vert_sp_ ## W ## x ## H ## cpu;

#define CHROMA_420_VSP_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_420_VSP_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(4, 2, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(2, 4, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(6, 8, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(16, 12, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(12, 16, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(16, 4, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(32, 24, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(24, 32, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(32, 8, cpu);

#define CHROMA_420_VSP_FILTERS(cpu) \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 2, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 6, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_420_VSP_FUNC_DEF(8, 32, cpu);

#define SETUP_CHROMA_422_VSP_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vsp = x265_interp_4tap_vert_sp_ ## W ## x ## H ## cpu;

#define CHROMA_422_VSP_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_422_VSP_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(6, 16, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(2, 16, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(16, 24, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(12, 32, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(4, 32, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(32, 64, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(16, 64, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(32, 48, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(24, 64, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(32, 16, cpu);

#define CHROMA_422_VSP_FILTERS(cpu) \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 12, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 32, cpu); \
    SETUP_CHROMA_422_VSP_FUNC_DEF(8, 64, cpu);

#define SETUP_CHROMA_444_VSP_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I444].pu[LUMA_ ## W ## x ## H].filter_vsp = x265_interp_4tap_vert_sp_ ## W ## x ## H ## cpu;

#define CHROMA_444_VSP_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_444_VSP_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 12, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(12, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 4, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(32, 24, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(24, 32, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(32, 8, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(64, 64, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(64, 32, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(32, 64, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(64, 48, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(48, 64, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(64, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(16, 64, cpu);

#define CHROMA_444_VSP_FILTERS(cpu) \
    SETUP_CHROMA_444_VSP_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_444_VSP_FUNC_DEF(8, 32, cpu);

#define SETUP_CHROMA_420_VSS_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vss = x265_interp_4tap_vert_ss_ ## W ## x ## H ## cpu;

#define CHROMA_420_VSS_FILTERS(cpu) \
    SETUP_CHROMA_420_VSS_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(4, 2, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 6, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 2, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(16, 12, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(12, 16, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(16, 4, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(32, 24, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(24, 32, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(32, 8, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(8, 32, cpu);

#define CHROMA_420_VSS_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_420_VSS_FUNC_DEF(2, 4, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_420_VSS_FUNC_DEF(6, 8, cpu);

#define SETUP_CHROMA_422_VSS_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vss = x265_interp_4tap_vert_ss_ ## W ## x ## H ## cpu;

#define CHROMA_422_VSS_FILTERS(cpu) \
    SETUP_CHROMA_422_VSS_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 12, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 32, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(16, 24, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(12, 32, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(4, 32, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(32, 64, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(16, 64, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(32, 48, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(24, 64, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(8, 64, cpu);

#define CHROMA_422_VSS_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_422_VSS_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(2, 16, cpu); \
    SETUP_CHROMA_422_VSS_FUNC_DEF(6, 16, cpu);

#define CHROMA_444_VSS_FILTERS(cpu) ALL_CHROMA_444_PU(filter_vss, interp_4tap_vert_ss, cpu)

#define LUMA_FILTERS(cpu) \
    ALL_LUMA_PU(luma_hpp, interp_8tap_horiz_pp, cpu); p.pu[LUMA_4x4].luma_hpp = x265_interp_8tap_horiz_pp_4x4_ ## cpu; \
    ALL_LUMA_PU(luma_hps, interp_8tap_horiz_ps, cpu); p.pu[LUMA_4x4].luma_hps = x265_interp_8tap_horiz_ps_4x4_ ## cpu; \
    ALL_LUMA_PU(luma_vpp, interp_8tap_vert_pp, cpu); p.pu[LUMA_4x4].luma_vpp = x265_interp_8tap_vert_pp_4x4_ ## cpu; \
    ALL_LUMA_PU(luma_vps, interp_8tap_vert_ps, cpu); p.pu[LUMA_4x4].luma_vps = x265_interp_8tap_vert_ps_4x4_ ## cpu; \
    ALL_LUMA_PU(luma_vsp, interp_8tap_vert_sp, cpu); p.pu[LUMA_4x4].luma_vsp = x265_interp_8tap_vert_sp_4x4_ ## cpu; \
    ALL_LUMA_PU_T(luma_hvpp, interp_8tap_hv_pp_cpu); p.pu[LUMA_4x4].luma_hvpp = interp_8tap_hv_pp_cpu<LUMA_4x4>;

#define LUMA_VSS_FILTERS(cpu) ALL_LUMA_PU(luma_vss, interp_8tap_vert_ss, cpu); p.pu[LUMA_4x4].luma_vss = x265_interp_8tap_vert_ss_4x4_ ## cpu

#define LUMA_CU_BLOCKCOPY(type, cpu) \
    p.cu[BLOCK_4x4].copy_ ## type = x265_blockcopy_ ## type ## _4x4_ ## cpu; \
    ALL_LUMA_CU(copy_ ## type, blockcopy_ ## type, cpu);

#define CHROMA_420_CU_BLOCKCOPY(type, cpu) ALL_CHROMA_420_CU(copy_ ## type, blockcopy_ ## type, cpu)
#define CHROMA_422_CU_BLOCKCOPY(type, cpu) ALL_CHROMA_422_CU(copy_ ## type, blockcopy_ ## type, cpu)

#define LUMA_PU_BLOCKCOPY(type, cpu)       ALL_LUMA_PU(copy_ ## type, blockcopy_ ## type, cpu); p.pu[LUMA_4x4].copy_ ## type = x265_blockcopy_ ## type ## _4x4_ ## cpu
#define CHROMA_420_PU_BLOCKCOPY(type, cpu) ALL_CHROMA_420_PU(copy_ ## type, blockcopy_ ## type, cpu)
#define CHROMA_422_PU_BLOCKCOPY(type, cpu) ALL_CHROMA_422_PU(copy_ ## type, blockcopy_ ## type, cpu)

#define LUMA_PIXELSUB(cpu) \
    p.cu[BLOCK_4x4].sub_ps = x265_pixel_sub_ps_4x4_ ## cpu; \
    p.cu[BLOCK_4x4].add_ps = x265_pixel_add_ps_4x4_ ## cpu; \
    ALL_LUMA_CU(sub_ps, pixel_sub_ps, cpu); \
    ALL_LUMA_CU(add_ps, pixel_add_ps, cpu);

#define CHROMA_420_PIXELSUB_PS(cpu) \
    ALL_CHROMA_420_CU(sub_ps, pixel_sub_ps, cpu); \
    ALL_CHROMA_420_CU(add_ps, pixel_add_ps, cpu);

#define CHROMA_422_PIXELSUB_PS(cpu) \
    ALL_CHROMA_422_CU(sub_ps, pixel_sub_ps, cpu); \
    ALL_CHROMA_422_CU(add_ps, pixel_add_ps, cpu);

#define LUMA_VAR(cpu)          ALL_LUMA_CU(var, pixel_var, cpu)

#define LUMA_ADDAVG(cpu)       ALL_LUMA_PU(addAvg, addAvg, cpu); p.pu[LUMA_4x4].addAvg = x265_addAvg_4x4_ ## cpu
#define CHROMA_420_ADDAVG(cpu) ALL_CHROMA_420_PU(addAvg, addAvg, cpu);
#define CHROMA_422_ADDAVG(cpu) ALL_CHROMA_422_PU(addAvg, addAvg, cpu);

#define SETUP_INTRA_ANG_COMMON(mode, fno, cpu) \
    p.cu[BLOCK_4x4].intra_pred[mode] = x265_intra_pred_ang4_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_8x8].intra_pred[mode] = x265_intra_pred_ang8_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_16x16].intra_pred[mode] = x265_intra_pred_ang16_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_32x32].intra_pred[mode] = x265_intra_pred_ang32_ ## fno ## _ ## cpu;

#define SETUP_INTRA_ANG4(mode, fno, cpu) \
    p.cu[BLOCK_4x4].intra_pred[mode] = x265_intra_pred_ang4_ ## fno ## _ ## cpu;

#define SETUP_INTRA_ANG16_32(mode, fno, cpu) \
    p.cu[BLOCK_16x16].intra_pred[mode] = x265_intra_pred_ang16_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_32x32].intra_pred[mode] = x265_intra_pred_ang32_ ## fno ## _ ## cpu;

#define SETUP_INTRA_ANG4_8(mode, fno, cpu) \
    p.cu[BLOCK_4x4].intra_pred[mode] = x265_intra_pred_ang4_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_8x8].intra_pred[mode] = x265_intra_pred_ang8_ ## fno ## _ ## cpu;

#define INTRA_ANG_SSSE3(cpu) \
    SETUP_INTRA_ANG_COMMON(2, 2, cpu); \
    SETUP_INTRA_ANG_COMMON(34, 2, cpu);

#define INTRA_ANG_SSE4_COMMON(cpu) \
    SETUP_INTRA_ANG_COMMON(3,  3,  cpu); \
    SETUP_INTRA_ANG_COMMON(4,  4,  cpu); \
    SETUP_INTRA_ANG_COMMON(5,  5,  cpu); \
    SETUP_INTRA_ANG_COMMON(6,  6,  cpu); \
    SETUP_INTRA_ANG_COMMON(7,  7,  cpu); \
    SETUP_INTRA_ANG_COMMON(8,  8,  cpu); \
    SETUP_INTRA_ANG_COMMON(9,  9,  cpu); \
    SETUP_INTRA_ANG_COMMON(10, 10, cpu); \
    SETUP_INTRA_ANG_COMMON(11, 11, cpu); \
    SETUP_INTRA_ANG_COMMON(12, 12, cpu); \
    SETUP_INTRA_ANG_COMMON(13, 13, cpu); \
    SETUP_INTRA_ANG_COMMON(14, 14, cpu); \
    SETUP_INTRA_ANG_COMMON(15, 15, cpu); \
    SETUP_INTRA_ANG_COMMON(16, 16, cpu); \
    SETUP_INTRA_ANG_COMMON(17, 17, cpu); \
    SETUP_INTRA_ANG_COMMON(18, 18, cpu);

#define SETUP_INTRA_ANG_HIGH(mode, fno, cpu) \
    p.cu[BLOCK_8x8].intra_pred[mode] = x265_intra_pred_ang8_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_16x16].intra_pred[mode] = x265_intra_pred_ang16_ ## fno ## _ ## cpu; \
    p.cu[BLOCK_32x32].intra_pred[mode] = x265_intra_pred_ang32_ ## fno ## _ ## cpu;

#define INTRA_ANG_SSE4_HIGH(cpu) \
    SETUP_INTRA_ANG_HIGH(19, 19, cpu); \
    SETUP_INTRA_ANG_HIGH(20, 20, cpu); \
    SETUP_INTRA_ANG_HIGH(21, 21, cpu); \
    SETUP_INTRA_ANG_HIGH(22, 22, cpu); \
    SETUP_INTRA_ANG_HIGH(23, 23, cpu); \
    SETUP_INTRA_ANG_HIGH(24, 24, cpu); \
    SETUP_INTRA_ANG_HIGH(25, 25, cpu); \
    SETUP_INTRA_ANG_HIGH(26, 26, cpu); \
    SETUP_INTRA_ANG_HIGH(27, 27, cpu); \
    SETUP_INTRA_ANG_HIGH(28, 28, cpu); \
    SETUP_INTRA_ANG_HIGH(29, 29, cpu); \
    SETUP_INTRA_ANG_HIGH(30, 30, cpu); \
    SETUP_INTRA_ANG_HIGH(31, 31, cpu); \
    SETUP_INTRA_ANG_HIGH(32, 32, cpu); \
    SETUP_INTRA_ANG_HIGH(33, 33, cpu); \
    SETUP_INTRA_ANG4(19, 17, cpu); \
    SETUP_INTRA_ANG4(20, 16, cpu); \
    SETUP_INTRA_ANG4(21, 15, cpu); \
    SETUP_INTRA_ANG4(22, 14, cpu); \
    SETUP_INTRA_ANG4(23, 13, cpu); \
    SETUP_INTRA_ANG4(24, 12, cpu); \
    SETUP_INTRA_ANG4(25, 11, cpu); \
    SETUP_INTRA_ANG4(26, 26, cpu); \
    SETUP_INTRA_ANG4(27, 9, cpu); \
    SETUP_INTRA_ANG4(28, 8, cpu); \
    SETUP_INTRA_ANG4(29, 7, cpu); \
    SETUP_INTRA_ANG4(30, 6, cpu); \
    SETUP_INTRA_ANG4(31, 5, cpu); \
    SETUP_INTRA_ANG4(32, 4, cpu); \
    SETUP_INTRA_ANG4(33, 3, cpu);

#define INTRA_ANG_SSE4(cpu) \
    SETUP_INTRA_ANG4_8(19, 17, cpu); \
    SETUP_INTRA_ANG4_8(20, 16, cpu); \
    SETUP_INTRA_ANG4_8(21, 15, cpu); \
    SETUP_INTRA_ANG4_8(22, 14, cpu); \
    SETUP_INTRA_ANG4_8(23, 13, cpu); \
    SETUP_INTRA_ANG4_8(24, 12, cpu); \
    SETUP_INTRA_ANG4_8(25, 11, cpu); \
    SETUP_INTRA_ANG4_8(26, 26, cpu); \
    SETUP_INTRA_ANG4_8(27, 9, cpu); \
    SETUP_INTRA_ANG4_8(28, 8, cpu); \
    SETUP_INTRA_ANG4_8(29, 7, cpu); \
    SETUP_INTRA_ANG4_8(30, 6, cpu); \
    SETUP_INTRA_ANG4_8(31, 5, cpu); \
    SETUP_INTRA_ANG4_8(32, 4, cpu); \
    SETUP_INTRA_ANG4_8(33, 3, cpu); \
    SETUP_INTRA_ANG16_32(19, 19, cpu); \
    SETUP_INTRA_ANG16_32(20, 20, cpu); \
    SETUP_INTRA_ANG16_32(21, 21, cpu); \
    SETUP_INTRA_ANG16_32(22, 22, cpu); \
    SETUP_INTRA_ANG16_32(23, 23, cpu); \
    SETUP_INTRA_ANG16_32(24, 24, cpu); \
    SETUP_INTRA_ANG16_32(25, 25, cpu); \
    SETUP_INTRA_ANG16_32(26, 26, cpu); \
    SETUP_INTRA_ANG16_32(27, 27, cpu); \
    SETUP_INTRA_ANG16_32(28, 28, cpu); \
    SETUP_INTRA_ANG16_32(29, 29, cpu); \
    SETUP_INTRA_ANG16_32(30, 30, cpu); \
    SETUP_INTRA_ANG16_32(31, 31, cpu); \
    SETUP_INTRA_ANG16_32(32, 32, cpu); \
    SETUP_INTRA_ANG16_32(33, 33, cpu);

#define CHROMA_420_VERT_FILTERS(cpu) \
    ALL_CHROMA_420_4x4_PU(filter_vss, interp_4tap_vert_ss, cpu); \
    ALL_CHROMA_420_4x4_PU(filter_vpp, interp_4tap_vert_pp, cpu); \
    ALL_CHROMA_420_4x4_PU(filter_vps, interp_4tap_vert_ps, cpu); \
    ALL_CHROMA_420_4x4_PU(filter_vsp, interp_4tap_vert_sp, cpu)

#define SETUP_CHROMA_420_VERT_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vss = x265_interp_4tap_vert_ss_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vpp = x265_interp_4tap_vert_pp_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vps = x265_interp_4tap_vert_ps_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I420].pu[CHROMA_420_ ## W ## x ## H].filter_vsp = x265_interp_4tap_vert_sp_ ## W ## x ## H ## cpu;

#define CHROMA_420_VERT_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_420_VERT_FUNC_DEF(2, 4, cpu); \
    SETUP_CHROMA_420_VERT_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_420_VERT_FUNC_DEF(4, 2, cpu); \
    SETUP_CHROMA_420_VERT_FUNC_DEF(6, 8, cpu);

#define SETUP_CHROMA_422_VERT_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vss = x265_interp_4tap_vert_ss_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vpp = x265_interp_4tap_vert_pp_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vps = x265_interp_4tap_vert_ps_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_vsp = x265_interp_4tap_vert_sp_ ## W ## x ## H ## cpu;

#define CHROMA_422_VERT_FILTERS(cpu) \
    SETUP_CHROMA_422_VERT_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 12, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 32, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(16, 24, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(12, 32, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(4, 32, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(32, 64, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(16, 64, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(32, 48, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(24, 64, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(8, 64, cpu);

#define CHROMA_422_VERT_FILTERS_SSE4(cpu) \
    SETUP_CHROMA_422_VERT_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(2, 16, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_422_VERT_FUNC_DEF(6, 16, cpu);

#define CHROMA_444_VERT_FILTERS(cpu) \
    ALL_CHROMA_444_PU(filter_vss, interp_4tap_vert_ss, cpu); \
    ALL_CHROMA_444_PU(filter_vpp, interp_4tap_vert_pp, cpu); \
    ALL_CHROMA_444_PU(filter_vps, interp_4tap_vert_ps, cpu); \
    ALL_CHROMA_444_PU(filter_vsp, interp_4tap_vert_sp, cpu)

#define CHROMA_420_HORIZ_FILTERS(cpu) \
    ALL_CHROMA_420_PU(filter_hpp, interp_4tap_horiz_pp, cpu); \
    ALL_CHROMA_420_PU(filter_hps, interp_4tap_horiz_ps, cpu);

#define SETUP_CHROMA_422_HORIZ_FUNC_DEF(W, H, cpu) \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_hpp = x265_interp_4tap_horiz_pp_ ## W ## x ## H ## cpu; \
    p.chroma[X265_CSP_I422].pu[CHROMA_422_ ## W ## x ## H].filter_hps = x265_interp_4tap_horiz_ps_ ## W ## x ## H ## cpu;

#define CHROMA_422_HORIZ_FILTERS(cpu) \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(4, 8, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(4, 4, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(2, 8, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 8, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(4, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 12, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(6, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 4, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(2, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(16, 32, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(16, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 32, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(16, 24, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(12, 32, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(16, 8, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(4, 32, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(32, 64, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(32, 32, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(16, 64, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(32, 48, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(24, 64, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(32, 16, cpu); \
    SETUP_CHROMA_422_HORIZ_FUNC_DEF(8, 64, cpu);

#define CHROMA_444_HORIZ_FILTERS(cpu) \
    ALL_CHROMA_444_PU(filter_hpp, interp_4tap_horiz_pp, cpu); \
    ALL_CHROMA_444_PU(filter_hps, interp_4tap_horiz_ps, cpu);

namespace x265 {
// private x265 namespace

template<int size>
void interp_8tap_hv_pp_cpu(const pixel* src, intptr_t srcStride, pixel* dst, intptr_t dstStride, int idxX, int idxY)
{
    ALIGN_VAR_32(int16_t, immed[MAX_CU_SIZE * (MAX_CU_SIZE + NTAPS_LUMA)]);
    const int filterSize = NTAPS_LUMA;
    const int halfFilterSize = filterSize >> 1;

    x265::primitives.pu[size].luma_hps(src, srcStride, immed, MAX_CU_SIZE, idxX, 1);
    x265::primitives.pu[size].luma_vsp(immed + (halfFilterSize - 1) * MAX_CU_SIZE, MAX_CU_SIZE, dst, dstStride, idxY);
}

#if HIGH_BIT_DEPTH

void setupAssemblyPrimitives(EncoderPrimitives &p, int cpuMask) // 16bpp
{
#if !defined(X86_64)
#error "Unsupported build configuration (32bit x86 and HIGH_BIT_DEPTH), you must configure ENABLE_ASSEMBLY=OFF"
#endif

#if X86_64
    p.scanPosLast = x265_scanPosLast_x64;
#endif

    if (cpuMask & X265_CPU_SSE2)
    {
        /* We do not differentiate CPUs which support MMX and not SSE2. We only check
         * for SSE2 and then use both MMX and SSE2 functions */
        AVC_LUMA_PU(sad, mmx2);

        p.pu[LUMA_16x16].sad = x265_pixel_sad_16x16_sse2;
        p.pu[LUMA_16x8].sad  = x265_pixel_sad_16x8_sse2;
        HEVC_SAD(sse2);

        p.pu[LUMA_4x4].sad_x3   = x265_pixel_sad_x3_4x4_mmx2;
        p.pu[LUMA_4x8].sad_x3   = x265_pixel_sad_x3_4x8_mmx2;
        p.pu[LUMA_4x16].sad_x3  = x265_pixel_sad_x3_4x16_mmx2;
        p.pu[LUMA_8x4].sad_x3   = x265_pixel_sad_x3_8x4_sse2;
        p.pu[LUMA_8x8].sad_x3   = x265_pixel_sad_x3_8x8_sse2;
        p.pu[LUMA_8x16].sad_x3  = x265_pixel_sad_x3_8x16_sse2;
        p.pu[LUMA_8x32].sad_x3  = x265_pixel_sad_x3_8x32_sse2;
        p.pu[LUMA_16x4].sad_x3  = x265_pixel_sad_x3_16x4_sse2;
        p.pu[LUMA_12x16].sad_x3 = x265_pixel_sad_x3_12x16_mmx2;
        HEVC_SAD_X3(sse2);

        p.pu[LUMA_4x4].sad_x4   = x265_pixel_sad_x4_4x4_mmx2;
        p.pu[LUMA_4x8].sad_x4   = x265_pixel_sad_x4_4x8_mmx2;
        p.pu[LUMA_4x16].sad_x4  = x265_pixel_sad_x4_4x16_mmx2;
        p.pu[LUMA_8x4].sad_x4   = x265_pixel_sad_x4_8x4_sse2;
        p.pu[LUMA_8x8].sad_x4   = x265_pixel_sad_x4_8x8_sse2;
        p.pu[LUMA_8x16].sad_x4  = x265_pixel_sad_x4_8x16_sse2;
        p.pu[LUMA_8x32].sad_x4  = x265_pixel_sad_x4_8x32_sse2;
        p.pu[LUMA_16x4].sad_x4  = x265_pixel_sad_x4_16x4_sse2;
        p.pu[LUMA_12x16].sad_x4 = x265_pixel_sad_x4_12x16_mmx2;
        HEVC_SAD_X4(sse2);

        p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_mmx2;
        ALL_LUMA_PU(satd, pixel_satd, sse2);

        ASSIGN_SA8D(sse2);
        LUMA_PIXELSUB(sse2);
        CHROMA_420_PIXELSUB_PS(sse2);
        CHROMA_422_PIXELSUB_PS(sse2);

        LUMA_CU_BLOCKCOPY(ss, sse2);
        CHROMA_420_CU_BLOCKCOPY(ss, sse2);
        CHROMA_422_CU_BLOCKCOPY(ss, sse2);

        p.pu[LUMA_4x4].copy_pp = (copy_pp_t)x265_blockcopy_ss_4x4_sse2;
        ALL_LUMA_PU_TYPED(copy_pp, (copy_pp_t), blockcopy_ss, sse2);
        ALL_CHROMA_420_PU_TYPED(copy_pp, (copy_pp_t), blockcopy_ss, sse2);
        ALL_CHROMA_422_PU_TYPED(copy_pp, (copy_pp_t), blockcopy_ss, sse2);

        CHROMA_420_VERT_FILTERS(sse2);
        CHROMA_422_VERT_FILTERS(_sse2);
        CHROMA_444_VERT_FILTERS(sse2);

        p.ssim_4x4x2_core = x265_pixel_ssim_4x4x2_core_sse2;
        p.ssim_end_4 = x265_pixel_ssim_end4_sse2;
        PIXEL_AVG(sse2);
        PIXEL_AVG_W4(mmx2);
        LUMA_VAR(sse2);


        ALL_LUMA_TU(blockfill_s, blockfill_s, sse2);
        ALL_LUMA_TU_S(cpy1Dto2D_shr, cpy1Dto2D_shr_, sse2);
        ALL_LUMA_TU_S(cpy1Dto2D_shl, cpy1Dto2D_shl_, sse2);
        ALL_LUMA_TU_S(cpy2Dto1D_shr, cpy2Dto1D_shr_, sse2);
        ALL_LUMA_TU_S(cpy2Dto1D_shl, cpy2Dto1D_shl_, sse2);
        ALL_LUMA_TU_S(ssd_s, pixel_ssd_s_, sse2);
        ALL_LUMA_TU_S(calcresidual, getResidual, sse2);
        ALL_LUMA_TU_S(transpose, transpose, sse2);

        ALL_LUMA_TU_S(intra_pred[PLANAR_IDX], intra_pred_planar, sse2);
        ALL_LUMA_TU_S(intra_pred[DC_IDX], intra_pred_dc, sse2);

        p.cu[BLOCK_4x4].intra_pred[2] = x265_intra_pred_ang4_2_sse2;
        p.cu[BLOCK_4x4].intra_pred[3] = x265_intra_pred_ang4_3_sse2;
        p.cu[BLOCK_4x4].intra_pred[4] = x265_intra_pred_ang4_4_sse2;
        p.cu[BLOCK_4x4].intra_pred[5] = x265_intra_pred_ang4_5_sse2;
        p.cu[BLOCK_4x4].intra_pred[6] = x265_intra_pred_ang4_6_sse2;
        p.cu[BLOCK_4x4].intra_pred[7] = x265_intra_pred_ang4_7_sse2;
        p.cu[BLOCK_4x4].intra_pred[8] = x265_intra_pred_ang4_8_sse2;
        p.cu[BLOCK_4x4].intra_pred[9] = x265_intra_pred_ang4_9_sse2;
        p.cu[BLOCK_4x4].intra_pred[10] = x265_intra_pred_ang4_10_sse2;
        p.cu[BLOCK_4x4].intra_pred[11] = x265_intra_pred_ang4_11_sse2;
        p.cu[BLOCK_4x4].intra_pred[12] = x265_intra_pred_ang4_12_sse2;
        p.cu[BLOCK_4x4].intra_pred[13] = x265_intra_pred_ang4_13_sse2;
        p.cu[BLOCK_4x4].intra_pred[14] = x265_intra_pred_ang4_14_sse2;
        p.cu[BLOCK_4x4].intra_pred[15] = x265_intra_pred_ang4_15_sse2;
        p.cu[BLOCK_4x4].intra_pred[16] = x265_intra_pred_ang4_16_sse2;
        p.cu[BLOCK_4x4].intra_pred[17] = x265_intra_pred_ang4_17_sse2;
        p.cu[BLOCK_4x4].intra_pred[18] = x265_intra_pred_ang4_18_sse2;
        p.cu[BLOCK_4x4].intra_pred[19] = x265_intra_pred_ang4_17_sse2;
        p.cu[BLOCK_4x4].intra_pred[20] = x265_intra_pred_ang4_16_sse2;
        p.cu[BLOCK_4x4].intra_pred[21] = x265_intra_pred_ang4_15_sse2;
        p.cu[BLOCK_4x4].intra_pred[22] = x265_intra_pred_ang4_14_sse2;
        p.cu[BLOCK_4x4].intra_pred[23] = x265_intra_pred_ang4_13_sse2;
        p.cu[BLOCK_4x4].intra_pred[24] = x265_intra_pred_ang4_12_sse2;
        p.cu[BLOCK_4x4].intra_pred[25] = x265_intra_pred_ang4_11_sse2;
        p.cu[BLOCK_4x4].intra_pred[26] = x265_intra_pred_ang4_26_sse2;
        p.cu[BLOCK_4x4].intra_pred[27] = x265_intra_pred_ang4_9_sse2;
        p.cu[BLOCK_4x4].intra_pred[28] = x265_intra_pred_ang4_8_sse2;
        p.cu[BLOCK_4x4].intra_pred[29] = x265_intra_pred_ang4_7_sse2;
        p.cu[BLOCK_4x4].intra_pred[30] = x265_intra_pred_ang4_6_sse2;
        p.cu[BLOCK_4x4].intra_pred[31] = x265_intra_pred_ang4_5_sse2;
        p.cu[BLOCK_4x4].intra_pred[32] = x265_intra_pred_ang4_4_sse2;
        p.cu[BLOCK_4x4].intra_pred[33] = x265_intra_pred_ang4_3_sse2;

        p.cu[BLOCK_4x4].sse_ss = x265_pixel_ssd_ss_4x4_mmx2;
        ALL_LUMA_CU(sse_ss, pixel_ssd_ss, sse2);

        p.chroma[X265_CSP_I422].cu[BLOCK_422_4x8].sse_pp = (pixelcmp_t)x265_pixel_ssd_ss_4x8_mmx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_8x16].sse_pp = (pixelcmp_t)x265_pixel_ssd_ss_8x16_sse2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].sse_pp = (pixelcmp_t)x265_pixel_ssd_ss_16x32_sse2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].sse_pp = (pixelcmp_t)x265_pixel_ssd_ss_32x64_sse2;

        p.cu[BLOCK_4x4].dct = x265_dct4_sse2;
        p.cu[BLOCK_8x8].dct = x265_dct8_sse2;
        p.cu[BLOCK_4x4].idct = x265_idct4_sse2;
        p.cu[BLOCK_8x8].idct = x265_idct8_sse2;

        p.idst4x4 = x265_idst4_sse2;

        LUMA_VSS_FILTERS(sse2);

        p.frameInitLowres = x265_frame_init_lowres_core_sse2;
    }
    if (cpuMask & X265_CPU_SSSE3)
    {
        p.scale1D_128to64 = x265_scale1D_128to64_ssse3;
        p.scale2D_64to32 = x265_scale2D_64to32_ssse3;

        // p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_ssse3; this one is broken
        ALL_LUMA_PU(satd, pixel_satd, ssse3);
        ASSIGN_SA8D(ssse3);
        INTRA_ANG_SSSE3(ssse3);

        p.dst4x4 = x265_dst4_ssse3;
        p.cu[BLOCK_8x8].idct = x265_idct8_ssse3;
        p.cu[BLOCK_4x4].count_nonzero = x265_count_nonzero_4x4_ssse3;
        p.cu[BLOCK_8x8].count_nonzero = x265_count_nonzero_8x8_ssse3;
        p.cu[BLOCK_16x16].count_nonzero = x265_count_nonzero_16x16_ssse3;
        p.cu[BLOCK_32x32].count_nonzero = x265_count_nonzero_32x32_ssse3;
        p.frameInitLowres = x265_frame_init_lowres_core_ssse3;

        p.pu[LUMA_4x4].convert_p2s = x265_filterPixelToShort_4x4_ssse3;
        p.pu[LUMA_4x8].convert_p2s = x265_filterPixelToShort_4x8_ssse3;
        p.pu[LUMA_4x16].convert_p2s = x265_filterPixelToShort_4x16_ssse3;
        p.pu[LUMA_8x4].convert_p2s = x265_filterPixelToShort_8x4_ssse3;
        p.pu[LUMA_8x8].convert_p2s = x265_filterPixelToShort_8x8_ssse3;
        p.pu[LUMA_8x16].convert_p2s = x265_filterPixelToShort_8x16_ssse3;
        p.pu[LUMA_8x32].convert_p2s = x265_filterPixelToShort_8x32_ssse3;
        p.pu[LUMA_16x4].convert_p2s = x265_filterPixelToShort_16x4_ssse3;
        p.pu[LUMA_16x8].convert_p2s = x265_filterPixelToShort_16x8_ssse3;
        p.pu[LUMA_16x12].convert_p2s = x265_filterPixelToShort_16x12_ssse3;
        p.pu[LUMA_16x16].convert_p2s = x265_filterPixelToShort_16x16_ssse3;
        p.pu[LUMA_16x32].convert_p2s = x265_filterPixelToShort_16x32_ssse3;
        p.pu[LUMA_16x64].convert_p2s = x265_filterPixelToShort_16x64_ssse3;
        p.pu[LUMA_32x8].convert_p2s = x265_filterPixelToShort_32x8_ssse3;
        p.pu[LUMA_32x16].convert_p2s = x265_filterPixelToShort_32x16_ssse3;
        p.pu[LUMA_32x24].convert_p2s = x265_filterPixelToShort_32x24_ssse3;
        p.pu[LUMA_32x32].convert_p2s = x265_filterPixelToShort_32x32_ssse3;
        p.pu[LUMA_32x64].convert_p2s = x265_filterPixelToShort_32x64_ssse3;
        p.pu[LUMA_64x16].convert_p2s = x265_filterPixelToShort_64x16_ssse3;
        p.pu[LUMA_64x32].convert_p2s = x265_filterPixelToShort_64x32_ssse3;
        p.pu[LUMA_64x48].convert_p2s = x265_filterPixelToShort_64x48_ssse3;
        p.pu[LUMA_64x64].convert_p2s = x265_filterPixelToShort_64x64_ssse3;
        p.pu[LUMA_24x32].convert_p2s = x265_filterPixelToShort_24x32_ssse3;
        p.pu[LUMA_12x16].convert_p2s = x265_filterPixelToShort_12x16_ssse3;
        p.pu[LUMA_48x64].convert_p2s = x265_filterPixelToShort_48x64_ssse3;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].p2s = x265_filterPixelToShort_4x4_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].p2s = x265_filterPixelToShort_4x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].p2s = x265_filterPixelToShort_4x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].p2s = x265_filterPixelToShort_8x4_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].p2s = x265_filterPixelToShort_8x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].p2s = x265_filterPixelToShort_8x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].p2s = x265_filterPixelToShort_8x32_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].p2s = x265_filterPixelToShort_16x4_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].p2s = x265_filterPixelToShort_16x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].p2s = x265_filterPixelToShort_16x12_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].p2s = x265_filterPixelToShort_16x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].p2s = x265_filterPixelToShort_16x32_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].p2s = x265_filterPixelToShort_32x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].p2s = x265_filterPixelToShort_32x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].p2s = x265_filterPixelToShort_32x24_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].p2s = x265_filterPixelToShort_32x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].p2s = x265_filterPixelToShort_4x4_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].p2s = x265_filterPixelToShort_4x8_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].p2s = x265_filterPixelToShort_4x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].p2s = x265_filterPixelToShort_4x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].p2s = x265_filterPixelToShort_8x4_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].p2s = x265_filterPixelToShort_8x8_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].p2s = x265_filterPixelToShort_8x12_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].p2s = x265_filterPixelToShort_8x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].p2s = x265_filterPixelToShort_8x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].p2s = x265_filterPixelToShort_8x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].p2s = x265_filterPixelToShort_12x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].p2s = x265_filterPixelToShort_16x8_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].p2s = x265_filterPixelToShort_16x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].p2s = x265_filterPixelToShort_16x24_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].p2s = x265_filterPixelToShort_16x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].p2s = x265_filterPixelToShort_16x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].p2s = x265_filterPixelToShort_24x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].p2s = x265_filterPixelToShort_32x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].p2s = x265_filterPixelToShort_32x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].p2s = x265_filterPixelToShort_32x48_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].p2s = x265_filterPixelToShort_32x64_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].p2s = x265_filterPixelToShort_4x2_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].p2s = x265_filterPixelToShort_8x2_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].p2s = x265_filterPixelToShort_8x6_ssse3;
        p.findPosFirstLast = x265_findPosFirstLast_ssse3;
    }
    if (cpuMask & X265_CPU_SSE4)
    {
        LUMA_ADDAVG(sse4);
        CHROMA_420_ADDAVG(sse4);
        CHROMA_422_ADDAVG(sse4);

        LUMA_FILTERS(sse4);
        CHROMA_420_HORIZ_FILTERS(sse4);
        CHROMA_420_VERT_FILTERS_SSE4(_sse4);
        CHROMA_422_HORIZ_FILTERS(_sse4);
        CHROMA_422_VERT_FILTERS_SSE4(_sse4);
        CHROMA_444_HORIZ_FILTERS(sse4);

        p.cu[BLOCK_8x8].dct = x265_dct8_sse4;
        p.quant = x265_quant_sse4;
        p.nquant = x265_nquant_sse4;
        p.dequant_normal = x265_dequant_normal_sse4;

        // p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_sse4; fails tests
        ALL_LUMA_PU(satd, pixel_satd, sse4);
        ASSIGN_SA8D(sse4);

        ALL_LUMA_TU_S(intra_pred[PLANAR_IDX], intra_pred_planar, sse4);
        ALL_LUMA_TU_S(intra_pred[DC_IDX], intra_pred_dc, sse4);
        INTRA_ANG_SSE4_COMMON(sse4);
        INTRA_ANG_SSE4_HIGH(sse4);

        p.planecopy_cp = x265_upShift_8_sse4;
        p.weight_pp = x265_weight_pp_sse4;
        p.weight_sp = x265_weight_sp_sse4;

        p.cu[BLOCK_4x4].psy_cost_pp = x265_psyCost_pp_4x4_sse4;
        p.cu[BLOCK_4x4].psy_cost_ss = x265_psyCost_ss_4x4_sse4;

        // TODO: check POPCNT flag!
        ALL_LUMA_TU_S(copy_cnt, copy_cnt_, sse4);
        ALL_LUMA_CU(psy_cost_pp, psyCost_pp, sse4);
        ALL_LUMA_CU(psy_cost_ss, psyCost_ss, sse4);

        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].p2s = x265_filterPixelToShort_2x4_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].p2s = x265_filterPixelToShort_2x8_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].p2s = x265_filterPixelToShort_6x8_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].p2s = x265_filterPixelToShort_2x8_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].p2s = x265_filterPixelToShort_2x16_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].p2s = x265_filterPixelToShort_6x16_sse4;
    }
    if (cpuMask & X265_CPU_AVX)
    {
        // p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_avx; fails tests
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].satd = x265_pixel_satd_16x24_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].satd = x265_pixel_satd_32x48_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].satd = x265_pixel_satd_24x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].satd = x265_pixel_satd_8x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].satd = x265_pixel_satd_8x12_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].satd = x265_pixel_satd_12x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].satd = x265_pixel_satd_4x32_avx;

        ALL_LUMA_PU(satd, pixel_satd, avx);
        ASSIGN_SA8D(avx);
        LUMA_VAR(avx);
        p.ssim_4x4x2_core = x265_pixel_ssim_4x4x2_core_avx;
        p.ssim_end_4 = x265_pixel_ssim_end4_avx;

        // copy_pp primitives
        // 16 x N
        p.pu[LUMA_64x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x64_avx;
        p.pu[LUMA_16x4].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x4_avx;
        p.pu[LUMA_16x8].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x8_avx;
        p.pu[LUMA_16x12].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x12_avx;
        p.pu[LUMA_16x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x16_avx;
        p.pu[LUMA_16x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x32_avx;
        p.pu[LUMA_16x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x64_avx;
        p.pu[LUMA_64x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x16_avx;
        p.pu[LUMA_64x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x32_avx;
        p.pu[LUMA_64x48].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x48_avx;
        p.pu[LUMA_64x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x64_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x4_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x12_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x24_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_16x32_avx;

        // 24 X N
        p.pu[LUMA_24x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_24x32_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_24x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_24x64_avx;

        // 32 x N
        p.pu[LUMA_32x8].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x8_avx;
        p.pu[LUMA_32x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x16_avx;
        p.pu[LUMA_32x24].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x24_avx;
        p.pu[LUMA_32x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x32_avx;
        p.pu[LUMA_32x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x64_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x24_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x48_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_32x64_avx;

        // 48 X 64
        p.pu[LUMA_48x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_48x64_avx;

        // copy_ss primitives
        // 16 X N
        p.cu[BLOCK_16x16].copy_ss = x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].copy_ss = x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].copy_ss = x265_blockcopy_ss_16x32_avx;

        // 32 X N
        p.cu[BLOCK_32x32].copy_ss = x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].copy_ss = x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].copy_ss = x265_blockcopy_ss_32x64_avx;

        // 64 X N
        p.cu[BLOCK_64x64].copy_ss = x265_blockcopy_ss_64x64_avx;

        // copy_ps primitives
        // 16 X N
        p.cu[BLOCK_16x16].copy_ps = (copy_ps_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].copy_ps = (copy_ps_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].copy_ps = (copy_ps_t)x265_blockcopy_ss_16x32_avx;

        // 32 X N
        p.cu[BLOCK_32x32].copy_ps = (copy_ps_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].copy_ps = (copy_ps_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].copy_ps = (copy_ps_t)x265_blockcopy_ss_32x64_avx;

        // 64 X N
        p.cu[BLOCK_64x64].copy_ps = (copy_ps_t)x265_blockcopy_ss_64x64_avx;

        // copy_sp primitives
        // 16 X N
        p.cu[BLOCK_16x16].copy_sp = (copy_sp_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].copy_sp = (copy_sp_t)x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].copy_sp = (copy_sp_t)x265_blockcopy_ss_16x32_avx;

        // 32 X N
        p.cu[BLOCK_32x32].copy_sp = (copy_sp_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].copy_sp = (copy_sp_t)x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].copy_sp = (copy_sp_t)x265_blockcopy_ss_32x64_avx;

        // 64 X N
        p.cu[BLOCK_64x64].copy_sp = (copy_sp_t)x265_blockcopy_ss_64x64_avx;

        p.frameInitLowres = x265_frame_init_lowres_core_avx;

        p.pu[LUMA_64x16].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x16_avx;
        p.pu[LUMA_64x32].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x32_avx;
        p.pu[LUMA_64x48].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x48_avx;
        p.pu[LUMA_64x64].copy_pp = (copy_pp_t)x265_blockcopy_ss_64x64_avx;
    }
    if (cpuMask & X265_CPU_XOP)
    {
        //p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_xop; this one is broken
        ALL_LUMA_PU(satd, pixel_satd, xop);
        ASSIGN_SA8D(xop);
        LUMA_VAR(xop);
        p.frameInitLowres = x265_frame_init_lowres_core_xop;
    }
    if (cpuMask & X265_CPU_AVX2)
    {
        p.pu[LUMA_48x64].satd = x265_pixel_satd_48x64_avx2;

        p.pu[LUMA_64x16].satd = x265_pixel_satd_64x16_avx2;
        p.pu[LUMA_64x32].satd = x265_pixel_satd_64x32_avx2;
        p.pu[LUMA_64x48].satd = x265_pixel_satd_64x48_avx2;
        p.pu[LUMA_64x64].satd = x265_pixel_satd_64x64_avx2;

        p.pu[LUMA_32x8].satd = x265_pixel_satd_32x8_avx2;
        p.pu[LUMA_32x16].satd = x265_pixel_satd_32x16_avx2;
        p.pu[LUMA_32x24].satd = x265_pixel_satd_32x24_avx2;
        p.pu[LUMA_32x32].satd = x265_pixel_satd_32x32_avx2;
        p.pu[LUMA_32x64].satd = x265_pixel_satd_32x64_avx2;

        p.pu[LUMA_16x4].satd = x265_pixel_satd_16x4_avx2;
        p.pu[LUMA_16x8].satd = x265_pixel_satd_16x8_avx2;
        p.pu[LUMA_16x12].satd = x265_pixel_satd_16x12_avx2;
        p.pu[LUMA_16x16].satd = x265_pixel_satd_16x16_avx2;
        p.pu[LUMA_16x32].satd = x265_pixel_satd_16x32_avx2;
        p.pu[LUMA_16x64].satd = x265_pixel_satd_16x64_avx2;

        p.cu[BLOCK_32x32].ssd_s = x265_pixel_ssd_s_32_avx2;
        p.cu[BLOCK_16x16].sse_ss = x265_pixel_ssd_ss_16x16_avx2;

        p.quant = x265_quant_avx2;
        p.nquant = x265_nquant_avx2;
        p.dequant_normal  = x265_dequant_normal_avx2;

        p.scale1D_128to64 = x265_scale1D_128to64_avx2;
        p.scale2D_64to32 = x265_scale2D_64to32_avx2;
        // p.weight_pp = x265_weight_pp_avx2; fails tests

        p.cu[BLOCK_16x16].calcresidual = x265_getResidual16_avx2;
        p.cu[BLOCK_32x32].calcresidual = x265_getResidual32_avx2;

        p.cu[BLOCK_16x16].blockfill_s = x265_blockfill_s_16x16_avx2;
        p.cu[BLOCK_32x32].blockfill_s = x265_blockfill_s_32x32_avx2;

        ALL_LUMA_TU(count_nonzero, count_nonzero, avx2);
        ALL_LUMA_TU_S(cpy1Dto2D_shl, cpy1Dto2D_shl_, avx2);
        ALL_LUMA_TU_S(cpy1Dto2D_shr, cpy1Dto2D_shr_, avx2);

        p.cu[BLOCK_8x8].copy_cnt = x265_copy_cnt_8_avx2;
        p.cu[BLOCK_16x16].copy_cnt = x265_copy_cnt_16_avx2;
        p.cu[BLOCK_32x32].copy_cnt = x265_copy_cnt_32_avx2;

        p.cu[BLOCK_8x8].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_8_avx2;
        p.cu[BLOCK_16x16].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_16_avx2;
        p.cu[BLOCK_32x32].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_32_avx2;

        p.cu[BLOCK_8x8].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_8_avx2;
        p.cu[BLOCK_16x16].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_16_avx2;
        p.cu[BLOCK_32x32].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_32_avx2;

        ALL_LUMA_TU_S(dct, dct, avx2);
        ALL_LUMA_TU_S(idct, idct, avx2);
        ALL_LUMA_CU_S(transpose, transpose, avx2);

        ALL_LUMA_PU(luma_vpp, interp_8tap_vert_pp, avx2);
        ALL_LUMA_PU(luma_vps, interp_8tap_vert_ps, avx2);
        ALL_LUMA_PU(luma_vsp, interp_8tap_vert_sp, avx2);
        ALL_LUMA_PU(luma_vss, interp_8tap_vert_ss, avx2);

        p.cu[BLOCK_16x16].add_ps = x265_pixel_add_ps_16x16_avx2;
        p.cu[BLOCK_32x32].add_ps = x265_pixel_add_ps_32x32_avx2;
        p.cu[BLOCK_64x64].add_ps = x265_pixel_add_ps_64x64_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].add_ps = x265_pixel_add_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].add_ps = x265_pixel_add_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].add_ps = x265_pixel_add_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].add_ps = x265_pixel_add_ps_32x64_avx2;

        p.cu[BLOCK_16x16].sub_ps = x265_pixel_sub_ps_16x16_avx2;
        p.cu[BLOCK_32x32].sub_ps = x265_pixel_sub_ps_32x32_avx2;
        p.cu[BLOCK_64x64].sub_ps = x265_pixel_sub_ps_64x64_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].sub_ps = x265_pixel_sub_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].sub_ps = x265_pixel_sub_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].sub_ps = x265_pixel_sub_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].sub_ps = x265_pixel_sub_ps_32x64_avx2;

        p.pu[LUMA_16x4].sad = x265_pixel_sad_16x4_avx2;
        p.pu[LUMA_16x8].sad = x265_pixel_sad_16x8_avx2;
        p.pu[LUMA_16x12].sad = x265_pixel_sad_16x12_avx2;
        p.pu[LUMA_16x16].sad = x265_pixel_sad_16x16_avx2;
        p.pu[LUMA_16x32].sad = x265_pixel_sad_16x32_avx2;

        p.pu[LUMA_16x4].convert_p2s = x265_filterPixelToShort_16x4_avx2;
        p.pu[LUMA_16x8].convert_p2s = x265_filterPixelToShort_16x8_avx2;
        p.pu[LUMA_16x12].convert_p2s = x265_filterPixelToShort_16x12_avx2;
        p.pu[LUMA_16x16].convert_p2s = x265_filterPixelToShort_16x16_avx2;
        p.pu[LUMA_16x32].convert_p2s = x265_filterPixelToShort_16x32_avx2;
        p.pu[LUMA_16x64].convert_p2s = x265_filterPixelToShort_16x64_avx2;
        p.pu[LUMA_32x8].convert_p2s = x265_filterPixelToShort_32x8_avx2;
        p.pu[LUMA_32x16].convert_p2s = x265_filterPixelToShort_32x16_avx2;
        p.pu[LUMA_32x24].convert_p2s = x265_filterPixelToShort_32x24_avx2;
        p.pu[LUMA_32x32].convert_p2s = x265_filterPixelToShort_32x32_avx2;
        p.pu[LUMA_32x64].convert_p2s = x265_filterPixelToShort_32x64_avx2;
        p.pu[LUMA_64x16].convert_p2s = x265_filterPixelToShort_64x16_avx2;
        p.pu[LUMA_64x32].convert_p2s = x265_filterPixelToShort_64x32_avx2;
        p.pu[LUMA_64x48].convert_p2s = x265_filterPixelToShort_64x48_avx2;
        p.pu[LUMA_64x64].convert_p2s = x265_filterPixelToShort_64x64_avx2;
        p.pu[LUMA_24x32].convert_p2s = x265_filterPixelToShort_24x32_avx2;
        p.pu[LUMA_48x64].convert_p2s = x265_filterPixelToShort_48x64_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].p2s = x265_filterPixelToShort_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].p2s = x265_filterPixelToShort_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].p2s = x265_filterPixelToShort_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].p2s = x265_filterPixelToShort_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].p2s = x265_filterPixelToShort_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].p2s = x265_filterPixelToShort_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].p2s = x265_filterPixelToShort_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].p2s = x265_filterPixelToShort_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].p2s = x265_filterPixelToShort_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].p2s = x265_filterPixelToShort_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].p2s = x265_filterPixelToShort_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].p2s = x265_filterPixelToShort_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].p2s = x265_filterPixelToShort_16x24_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].p2s = x265_filterPixelToShort_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].p2s = x265_filterPixelToShort_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].p2s = x265_filterPixelToShort_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].p2s = x265_filterPixelToShort_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].p2s = x265_filterPixelToShort_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].p2s = x265_filterPixelToShort_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].p2s = x265_filterPixelToShort_32x64_avx2;

        p.pu[LUMA_4x4].luma_hps = x265_interp_8tap_horiz_ps_4x4_avx2;
        p.pu[LUMA_4x8].luma_hps = x265_interp_8tap_horiz_ps_4x8_avx2;
        p.pu[LUMA_4x16].luma_hps = x265_interp_8tap_horiz_ps_4x16_avx2;

        if (cpuMask & X265_CPU_BMI2)
            p.scanPosLast = x265_scanPosLast_avx2_bmi2;
    }
}
#else // if HIGH_BIT_DEPTH

void setupAssemblyPrimitives(EncoderPrimitives &p, int cpuMask) // 8bpp
{
#if X86_64
    p.scanPosLast = x265_scanPosLast_x64;
#endif

    if (cpuMask & X265_CPU_SSE2)
    {
        /* We do not differentiate CPUs which support MMX and not SSE2. We only check
         * for SSE2 and then use both MMX and SSE2 functions */
        AVC_LUMA_PU(sad, mmx2);
        AVC_LUMA_PU(sad_x3, mmx2);
        AVC_LUMA_PU(sad_x4, mmx2);

        p.pu[LUMA_16x16].sad = x265_pixel_sad_16x16_sse2;
        p.pu[LUMA_16x16].sad_x3 = x265_pixel_sad_x3_16x16_sse2;
        p.pu[LUMA_16x16].sad_x4 = x265_pixel_sad_x4_16x16_sse2;
        p.pu[LUMA_16x8].sad  = x265_pixel_sad_16x8_sse2;
        p.pu[LUMA_16x8].sad_x3  = x265_pixel_sad_x3_16x8_sse2;
        p.pu[LUMA_16x8].sad_x4  = x265_pixel_sad_x4_16x8_sse2;
        HEVC_SAD(sse2);

        p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_mmx2;
        ALL_LUMA_PU(satd, pixel_satd, sse2);

        p.cu[BLOCK_4x4].sse_pp = x265_pixel_ssd_4x4_mmx;
        p.cu[BLOCK_8x8].sse_pp = x265_pixel_ssd_8x8_mmx;
        p.cu[BLOCK_16x16].sse_pp = x265_pixel_ssd_16x16_mmx;

        PIXEL_AVG_W4(mmx2);
        PIXEL_AVG(sse2);
        LUMA_VAR(sse2);

        ASSIGN_SA8D(sse2);
        p.chroma[X265_CSP_I422].cu[BLOCK_422_4x8].sse_pp = x265_pixel_ssd_4x8_mmx;
        ASSIGN_SSE_PP(sse2);
        ASSIGN_SSE_SS(sse2);

        LUMA_PU_BLOCKCOPY(pp, sse2);
        CHROMA_420_PU_BLOCKCOPY(pp, sse2);
        CHROMA_422_PU_BLOCKCOPY(pp, sse2);

        LUMA_CU_BLOCKCOPY(ss, sse2);
        LUMA_CU_BLOCKCOPY(sp, sse2);
        CHROMA_420_CU_BLOCKCOPY(ss, sse2);
        CHROMA_422_CU_BLOCKCOPY(ss, sse2);
        CHROMA_420_CU_BLOCKCOPY(sp, sse2);
        CHROMA_422_CU_BLOCKCOPY(sp, sse2);

        LUMA_VSS_FILTERS(sse2);
        CHROMA_420_VSS_FILTERS(_sse2);
        CHROMA_422_VSS_FILTERS(_sse2);
        CHROMA_444_VSS_FILTERS(sse2);
        CHROMA_420_VSP_FILTERS(_sse2);
        CHROMA_422_VSP_FILTERS(_sse2);
        CHROMA_444_VSP_FILTERS(_sse2);
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_vpp = x265_interp_4tap_vert_pp_2x4_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_vpp = x265_interp_4tap_vert_pp_2x8_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_vpp = x265_interp_4tap_vert_pp_4x2_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_vpp = x265_interp_4tap_vert_pp_2x16_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].filter_vpp = x265_interp_4tap_vert_pp_4x32_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_sse2;
#if X86_64
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_vpp = x265_interp_4tap_vert_pp_6x8_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_vpp = x265_interp_4tap_vert_pp_8x2_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_vpp = x265_interp_4tap_vert_pp_8x6_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_sse2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].filter_vpp = x265_interp_4tap_vert_pp_6x16_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_vpp = x265_interp_4tap_vert_pp_8x12_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_sse2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_vpp = x265_interp_4tap_vert_pp_8x64_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_sse2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_sse2;
#endif

        ALL_LUMA_PU(luma_hpp, interp_8tap_horiz_pp, sse2);
        p.pu[LUMA_4x4].luma_hpp = x265_interp_8tap_horiz_pp_4x4_sse2;
        ALL_LUMA_PU(luma_hps, interp_8tap_horiz_ps, sse2);
        p.pu[LUMA_4x4].luma_hps = x265_interp_8tap_horiz_ps_4x4_sse2;
        p.pu[LUMA_8x8].luma_hvpp = x265_interp_8tap_hv_pp_8x8_sse3;

        //p.frameInitLowres = x265_frame_init_lowres_core_mmx2;
        p.frameInitLowres = x265_frame_init_lowres_core_sse2;

        ALL_LUMA_TU(blockfill_s, blockfill_s, sse2);
        ALL_LUMA_TU_S(cpy2Dto1D_shl, cpy2Dto1D_shl_, sse2);
        ALL_LUMA_TU_S(cpy2Dto1D_shr, cpy2Dto1D_shr_, sse2);
        ALL_LUMA_TU_S(cpy1Dto2D_shl, cpy1Dto2D_shl_, sse2);
        ALL_LUMA_TU_S(cpy1Dto2D_shr, cpy1Dto2D_shr_, sse2);
        ALL_LUMA_TU_S(ssd_s, pixel_ssd_s_, sse2);

        ALL_LUMA_TU_S(intra_pred[PLANAR_IDX], intra_pred_planar, sse2);
        ALL_LUMA_TU_S(intra_pred[DC_IDX], intra_pred_dc, sse2);

        p.cu[BLOCK_4x4].intra_pred[2] = x265_intra_pred_ang4_2_sse2;
        p.cu[BLOCK_4x4].intra_pred[3] = x265_intra_pred_ang4_3_sse2;
        p.cu[BLOCK_4x4].intra_pred[4] = x265_intra_pred_ang4_4_sse2;
        p.cu[BLOCK_4x4].intra_pred[5] = x265_intra_pred_ang4_5_sse2;
        p.cu[BLOCK_4x4].intra_pred[6] = x265_intra_pred_ang4_6_sse2;
        p.cu[BLOCK_4x4].intra_pred[7] = x265_intra_pred_ang4_7_sse2;
        p.cu[BLOCK_4x4].intra_pred[8] = x265_intra_pred_ang4_8_sse2;
        p.cu[BLOCK_4x4].intra_pred[9] = x265_intra_pred_ang4_9_sse2;
        p.cu[BLOCK_4x4].intra_pred[10] = x265_intra_pred_ang4_10_sse2;
        p.cu[BLOCK_4x4].intra_pred[11] = x265_intra_pred_ang4_11_sse2;
        p.cu[BLOCK_4x4].intra_pred[12] = x265_intra_pred_ang4_12_sse2;
        p.cu[BLOCK_4x4].intra_pred[13] = x265_intra_pred_ang4_13_sse2;
        p.cu[BLOCK_4x4].intra_pred[14] = x265_intra_pred_ang4_14_sse2;
        p.cu[BLOCK_4x4].intra_pred[15] = x265_intra_pred_ang4_15_sse2;
        p.cu[BLOCK_4x4].intra_pred[16] = x265_intra_pred_ang4_16_sse2;
        p.cu[BLOCK_4x4].intra_pred[17] = x265_intra_pred_ang4_17_sse2;
        p.cu[BLOCK_4x4].intra_pred[18] = x265_intra_pred_ang4_18_sse2;
        p.cu[BLOCK_4x4].intra_pred[19] = x265_intra_pred_ang4_17_sse2;
        p.cu[BLOCK_4x4].intra_pred[20] = x265_intra_pred_ang4_16_sse2;
        p.cu[BLOCK_4x4].intra_pred[21] = x265_intra_pred_ang4_15_sse2;
        p.cu[BLOCK_4x4].intra_pred[22] = x265_intra_pred_ang4_14_sse2;
        p.cu[BLOCK_4x4].intra_pred[23] = x265_intra_pred_ang4_13_sse2;
        p.cu[BLOCK_4x4].intra_pred[24] = x265_intra_pred_ang4_12_sse2;
        p.cu[BLOCK_4x4].intra_pred[25] = x265_intra_pred_ang4_11_sse2;
        p.cu[BLOCK_4x4].intra_pred[26] = x265_intra_pred_ang4_26_sse2;
        p.cu[BLOCK_4x4].intra_pred[27] = x265_intra_pred_ang4_9_sse2;
        p.cu[BLOCK_4x4].intra_pred[28] = x265_intra_pred_ang4_8_sse2;
        p.cu[BLOCK_4x4].intra_pred[29] = x265_intra_pred_ang4_7_sse2;
        p.cu[BLOCK_4x4].intra_pred[30] = x265_intra_pred_ang4_6_sse2;
        p.cu[BLOCK_4x4].intra_pred[31] = x265_intra_pred_ang4_5_sse2;
        p.cu[BLOCK_4x4].intra_pred[32] = x265_intra_pred_ang4_4_sse2;
        p.cu[BLOCK_4x4].intra_pred[33] = x265_intra_pred_ang4_3_sse2;

        p.cu[BLOCK_4x4].intra_pred_allangs = x265_all_angs_pred_4x4_sse2;

        p.cu[BLOCK_4x4].calcresidual = x265_getResidual4_sse2;
        p.cu[BLOCK_8x8].calcresidual = x265_getResidual8_sse2;

        ALL_LUMA_TU_S(transpose, transpose, sse2);
        p.cu[BLOCK_64x64].transpose = x265_transpose64_sse2;

        p.ssim_4x4x2_core = x265_pixel_ssim_4x4x2_core_sse2;
        p.ssim_end_4 = x265_pixel_ssim_end4_sse2;

        p.cu[BLOCK_4x4].dct = x265_dct4_sse2;
        p.cu[BLOCK_8x8].dct = x265_dct8_sse2;
        p.cu[BLOCK_4x4].idct = x265_idct4_sse2;
#if X86_64
        p.cu[BLOCK_8x8].idct = x265_idct8_sse2;
#endif
        p.idst4x4 = x265_idst4_sse2;

        p.planecopy_sp = x265_downShift_16_sse2;
    }
    if (cpuMask & X265_CPU_SSE3)
    {
        ALL_CHROMA_420_PU(filter_hpp, interp_4tap_horiz_pp, sse3);
        ALL_CHROMA_422_PU(filter_hpp, interp_4tap_horiz_pp, sse3);
        ALL_CHROMA_444_PU(filter_hpp, interp_4tap_horiz_pp, sse3);
    }
    if (cpuMask & X265_CPU_SSSE3)
    {
        p.pu[LUMA_8x16].sad_x3 = x265_pixel_sad_x3_8x16_ssse3;
        p.pu[LUMA_8x32].sad_x3 = x265_pixel_sad_x3_8x32_ssse3;
        p.pu[LUMA_12x16].sad_x3 = x265_pixel_sad_x3_12x16_ssse3;
        HEVC_SAD_X3(ssse3);

        p.pu[LUMA_8x4].sad_x4  = x265_pixel_sad_x4_8x4_ssse3;
        p.pu[LUMA_8x8].sad_x4  = x265_pixel_sad_x4_8x8_ssse3;
        p.pu[LUMA_8x16].sad_x4 = x265_pixel_sad_x4_8x16_ssse3;
        p.pu[LUMA_8x32].sad_x4 = x265_pixel_sad_x4_8x32_ssse3;
        p.pu[LUMA_12x16].sad_x4 = x265_pixel_sad_x4_12x16_ssse3;
        HEVC_SAD_X4(ssse3);

        p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_ssse3;
        ALL_LUMA_PU(satd, pixel_satd, ssse3);

        ASSIGN_SA8D(ssse3);
        PIXEL_AVG(ssse3);
        PIXEL_AVG_W4(ssse3);
        INTRA_ANG_SSSE3(ssse3);

        ASSIGN_SSE_PP(ssse3);
        p.cu[BLOCK_4x4].sse_pp = x265_pixel_ssd_4x4_ssse3;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_4x8].sse_pp = x265_pixel_ssd_4x8_ssse3;

        p.dst4x4 = x265_dst4_ssse3;
        p.cu[BLOCK_8x8].idct = x265_idct8_ssse3;

        ALL_LUMA_TU(count_nonzero, count_nonzero, ssse3);

        // MUST be done after LUMA_FILTERS() to overwrite default version
        p.pu[LUMA_8x8].luma_hvpp = x265_interp_8tap_hv_pp_8x8_ssse3;

        p.frameInitLowres = x265_frame_init_lowres_core_ssse3;
        p.scale1D_128to64 = x265_scale1D_128to64_ssse3;
        p.scale2D_64to32 = x265_scale2D_64to32_ssse3;

        p.pu[LUMA_8x4].convert_p2s = x265_filterPixelToShort_8x4_ssse3;
        p.pu[LUMA_8x8].convert_p2s = x265_filterPixelToShort_8x8_ssse3;
        p.pu[LUMA_8x16].convert_p2s = x265_filterPixelToShort_8x16_ssse3;
        p.pu[LUMA_8x32].convert_p2s = x265_filterPixelToShort_8x32_ssse3;
        p.pu[LUMA_16x4].convert_p2s = x265_filterPixelToShort_16x4_ssse3;
        p.pu[LUMA_16x8].convert_p2s = x265_filterPixelToShort_16x8_ssse3;
        p.pu[LUMA_16x12].convert_p2s = x265_filterPixelToShort_16x12_ssse3;
        p.pu[LUMA_16x16].convert_p2s = x265_filterPixelToShort_16x16_ssse3;
        p.pu[LUMA_16x32].convert_p2s = x265_filterPixelToShort_16x32_ssse3;
        p.pu[LUMA_16x64].convert_p2s = x265_filterPixelToShort_16x64_ssse3;
        p.pu[LUMA_32x8].convert_p2s = x265_filterPixelToShort_32x8_ssse3;
        p.pu[LUMA_32x16].convert_p2s = x265_filterPixelToShort_32x16_ssse3;
        p.pu[LUMA_32x24].convert_p2s = x265_filterPixelToShort_32x24_ssse3;
        p.pu[LUMA_32x32].convert_p2s = x265_filterPixelToShort_32x32_ssse3;
        p.pu[LUMA_32x64].convert_p2s = x265_filterPixelToShort_32x64_ssse3;
        p.pu[LUMA_64x16].convert_p2s = x265_filterPixelToShort_64x16_ssse3;
        p.pu[LUMA_64x32].convert_p2s = x265_filterPixelToShort_64x32_ssse3;
        p.pu[LUMA_64x48].convert_p2s = x265_filterPixelToShort_64x48_ssse3;
        p.pu[LUMA_64x64].convert_p2s = x265_filterPixelToShort_64x64_ssse3;
        p.pu[LUMA_12x16].convert_p2s = x265_filterPixelToShort_12x16_ssse3;
        p.pu[LUMA_24x32].convert_p2s = x265_filterPixelToShort_24x32_ssse3;
        p.pu[LUMA_48x64].convert_p2s = x265_filterPixelToShort_48x64_ssse3;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].p2s = x265_filterPixelToShort_8x2_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].p2s = x265_filterPixelToShort_8x4_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].p2s = x265_filterPixelToShort_8x6_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].p2s = x265_filterPixelToShort_8x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].p2s = x265_filterPixelToShort_8x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].p2s = x265_filterPixelToShort_8x32_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].p2s = x265_filterPixelToShort_16x4_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].p2s = x265_filterPixelToShort_16x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].p2s = x265_filterPixelToShort_16x12_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].p2s = x265_filterPixelToShort_16x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].p2s = x265_filterPixelToShort_16x32_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].p2s = x265_filterPixelToShort_32x8_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].p2s = x265_filterPixelToShort_32x16_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].p2s = x265_filterPixelToShort_32x24_ssse3;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].p2s = x265_filterPixelToShort_32x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].p2s = x265_filterPixelToShort_8x4_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].p2s = x265_filterPixelToShort_8x8_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].p2s = x265_filterPixelToShort_8x12_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].p2s = x265_filterPixelToShort_8x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].p2s = x265_filterPixelToShort_8x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].p2s = x265_filterPixelToShort_8x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].p2s = x265_filterPixelToShort_12x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].p2s = x265_filterPixelToShort_16x8_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].p2s = x265_filterPixelToShort_16x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].p2s = x265_filterPixelToShort_16x24_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].p2s = x265_filterPixelToShort_16x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].p2s = x265_filterPixelToShort_16x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].p2s = x265_filterPixelToShort_24x64_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].p2s = x265_filterPixelToShort_32x16_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].p2s = x265_filterPixelToShort_32x32_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].p2s = x265_filterPixelToShort_32x48_ssse3;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].p2s = x265_filterPixelToShort_32x64_ssse3;
        p.findPosFirstLast = x265_findPosFirstLast_ssse3;
    }
    if (cpuMask & X265_CPU_SSE4)
    {
        p.sign = x265_calSign_sse4;
        p.saoCuOrgE0 = x265_saoCuOrgE0_sse4;
        p.saoCuOrgE1 = x265_saoCuOrgE1_sse4;
        p.saoCuOrgE1_2Rows = x265_saoCuOrgE1_2Rows_sse4;
        p.saoCuOrgE2[0] = x265_saoCuOrgE2_sse4;
        p.saoCuOrgE2[1] = x265_saoCuOrgE2_sse4;
        p.saoCuOrgE3[0] = x265_saoCuOrgE3_sse4;
        p.saoCuOrgE3[1] = x265_saoCuOrgE3_sse4;
        p.saoCuOrgB0 = x265_saoCuOrgB0_sse4;

        LUMA_ADDAVG(sse4);
        CHROMA_420_ADDAVG(sse4);
        CHROMA_422_ADDAVG(sse4);

        // TODO: check POPCNT flag!
        ALL_LUMA_TU_S(copy_cnt, copy_cnt_, sse4);

        p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_sse4;
        ALL_LUMA_PU(satd, pixel_satd, sse4);
        ASSIGN_SA8D(sse4);
        ASSIGN_SSE_SS(sse4);
        p.cu[BLOCK_64x64].sse_pp = x265_pixel_ssd_64x64_sse4;

        LUMA_PIXELSUB(sse4);
        CHROMA_420_PIXELSUB_PS(sse4);
        CHROMA_422_PIXELSUB_PS(sse4);

        LUMA_FILTERS(sse4);
        CHROMA_420_FILTERS(sse4);
        CHROMA_422_FILTERS(sse4);
        CHROMA_444_FILTERS(sse4);
        CHROMA_420_VSS_FILTERS_SSE4(_sse4);
        CHROMA_422_VSS_FILTERS_SSE4(_sse4);
        CHROMA_420_VSP_FILTERS_SSE4(_sse4);
        CHROMA_422_VSP_FILTERS_SSE4(_sse4);
        CHROMA_444_VSP_FILTERS_SSE4(_sse4);

        // MUST be done after LUMA_FILTERS() to overwrite default version
        p.pu[LUMA_8x8].luma_hvpp = x265_interp_8tap_hv_pp_8x8_ssse3;

        LUMA_CU_BLOCKCOPY(ps, sse4);
        CHROMA_420_CU_BLOCKCOPY(ps, sse4);
        CHROMA_422_CU_BLOCKCOPY(ps, sse4);

        p.cu[BLOCK_16x16].calcresidual = x265_getResidual16_sse4;
        p.cu[BLOCK_32x32].calcresidual = x265_getResidual32_sse4;
        p.cu[BLOCK_8x8].dct = x265_dct8_sse4;
        p.denoiseDct = x265_denoise_dct_sse4;
        p.quant = x265_quant_sse4;
        p.nquant = x265_nquant_sse4;
        p.dequant_normal = x265_dequant_normal_sse4;

        p.weight_pp = x265_weight_pp_sse4;
        p.weight_sp = x265_weight_sp_sse4;

        ALL_LUMA_TU_S(intra_pred[PLANAR_IDX], intra_pred_planar, sse4);
        ALL_LUMA_TU_S(intra_pred[DC_IDX], intra_pred_dc, sse4);
        ALL_LUMA_TU(intra_pred_allangs, all_angs_pred, sse4);

        INTRA_ANG_SSE4_COMMON(sse4);
        INTRA_ANG_SSE4(sse4);

        p.cu[BLOCK_4x4].psy_cost_pp = x265_psyCost_pp_4x4_sse4;
        p.cu[BLOCK_4x4].psy_cost_ss = x265_psyCost_ss_4x4_sse4;

        p.pu[LUMA_4x4].convert_p2s = x265_filterPixelToShort_4x4_sse4;
        p.pu[LUMA_4x8].convert_p2s = x265_filterPixelToShort_4x8_sse4;
        p.pu[LUMA_4x16].convert_p2s = x265_filterPixelToShort_4x16_sse4;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].p2s = x265_filterPixelToShort_2x4_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].p2s = x265_filterPixelToShort_2x8_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].p2s = x265_filterPixelToShort_4x2_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].p2s = x265_filterPixelToShort_4x4_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].p2s = x265_filterPixelToShort_4x8_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].p2s = x265_filterPixelToShort_4x16_sse4;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].p2s = x265_filterPixelToShort_6x8_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].p2s = x265_filterPixelToShort_2x8_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].p2s = x265_filterPixelToShort_2x16_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].p2s = x265_filterPixelToShort_4x4_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].p2s = x265_filterPixelToShort_4x8_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].p2s = x265_filterPixelToShort_4x16_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].p2s = x265_filterPixelToShort_4x32_sse4;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].p2s = x265_filterPixelToShort_6x16_sse4;

#if X86_64
        ALL_LUMA_CU(psy_cost_pp, psyCost_pp, sse4);
        ALL_LUMA_CU(psy_cost_ss, psyCost_ss, sse4);
#endif
    }
    if (cpuMask & X265_CPU_AVX)
    {
        p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].satd = x265_pixel_satd_16x24_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].satd = x265_pixel_satd_32x48_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].satd = x265_pixel_satd_24x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].satd = x265_pixel_satd_8x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].satd = x265_pixel_satd_8x12_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].satd = x265_pixel_satd_12x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].satd = x265_pixel_satd_4x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].satd = x265_pixel_satd_16x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].satd = x265_pixel_satd_32x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].satd = x265_pixel_satd_16x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].satd = x265_pixel_satd_32x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].satd = x265_pixel_satd_16x64_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].satd = x265_pixel_satd_16x8_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].satd = x265_pixel_satd_32x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].satd = x265_pixel_satd_8x4_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].satd = x265_pixel_satd_8x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].satd = x265_pixel_satd_8x8_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].satd = x265_pixel_satd_8x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].satd = x265_pixel_satd_4x8_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].satd = x265_pixel_satd_4x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].satd = x265_pixel_satd_4x4_avx;
        ALL_LUMA_PU(satd, pixel_satd, avx);
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].satd = x265_pixel_satd_4x4_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].satd = x265_pixel_satd_8x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].satd = x265_pixel_satd_16x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].satd = x265_pixel_satd_32x32_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].satd = x265_pixel_satd_8x4_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].satd = x265_pixel_satd_4x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].satd = x265_pixel_satd_16x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].satd = x265_pixel_satd_8x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].satd = x265_pixel_satd_32x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].satd = x265_pixel_satd_16x32_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].satd = x265_pixel_satd_16x12_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].satd = x265_pixel_satd_12x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].satd = x265_pixel_satd_16x4_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].satd = x265_pixel_satd_4x16_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].satd = x265_pixel_satd_32x24_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].satd = x265_pixel_satd_24x32_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].satd = x265_pixel_satd_32x8_avx;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].satd = x265_pixel_satd_8x32_avx;
        ASSIGN_SA8D(avx);
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].sa8d = x265_pixel_sa8d_32x32_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].sa8d = x265_pixel_sa8d_16x16_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_8x8].sa8d = x265_pixel_sa8d_8x8_avx;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_4x4].sa8d = x265_pixel_satd_4x4_avx;
        ASSIGN_SSE_PP(avx);
        p.chroma[X265_CSP_I420].cu[BLOCK_420_8x8].sse_pp = x265_pixel_ssd_8x8_avx;
        ASSIGN_SSE_SS(avx);
        LUMA_VAR(avx);

        p.pu[LUMA_12x16].sad_x3 = x265_pixel_sad_x3_12x16_avx;
        p.pu[LUMA_16x4].sad_x3  = x265_pixel_sad_x3_16x4_avx;
        HEVC_SAD_X3(avx);

        p.pu[LUMA_12x16].sad_x4 = x265_pixel_sad_x4_12x16_avx;
        p.pu[LUMA_16x4].sad_x4  = x265_pixel_sad_x4_16x4_avx;
        HEVC_SAD_X4(avx);

        p.ssim_4x4x2_core = x265_pixel_ssim_4x4x2_core_avx;
        p.ssim_end_4 = x265_pixel_ssim_end4_avx;

        p.cu[BLOCK_16x16].copy_ss = x265_blockcopy_ss_16x16_avx;
        p.cu[BLOCK_32x32].copy_ss = x265_blockcopy_ss_32x32_avx;
        p.cu[BLOCK_64x64].copy_ss = x265_blockcopy_ss_64x64_avx;
        p.chroma[X265_CSP_I420].cu[CHROMA_420_16x16].copy_ss = x265_blockcopy_ss_16x16_avx;
        p.chroma[X265_CSP_I420].cu[CHROMA_420_32x32].copy_ss = x265_blockcopy_ss_32x32_avx;
        p.chroma[X265_CSP_I422].cu[CHROMA_422_16x32].copy_ss = x265_blockcopy_ss_16x32_avx;
        p.chroma[X265_CSP_I422].cu[CHROMA_422_32x64].copy_ss = x265_blockcopy_ss_32x64_avx;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].copy_pp = x265_blockcopy_pp_32x8_avx;
        p.pu[LUMA_32x8].copy_pp = x265_blockcopy_pp_32x8_avx;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].copy_pp = x265_blockcopy_pp_32x16_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].copy_pp = x265_blockcopy_pp_32x16_avx;
        p.pu[LUMA_32x16].copy_pp = x265_blockcopy_pp_32x16_avx;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].copy_pp = x265_blockcopy_pp_32x24_avx;
        p.pu[LUMA_32x24].copy_pp = x265_blockcopy_pp_32x24_avx;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].copy_pp = x265_blockcopy_pp_32x32_avx;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].copy_pp = x265_blockcopy_pp_32x32_avx;
        p.pu[LUMA_32x32].copy_pp  = x265_blockcopy_pp_32x32_avx;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].copy_pp = x265_blockcopy_pp_32x48_avx;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].copy_pp = x265_blockcopy_pp_32x64_avx;
        p.pu[LUMA_32x64].copy_pp = x265_blockcopy_pp_32x64_avx;

        p.pu[LUMA_64x16].copy_pp = x265_blockcopy_pp_64x16_avx;
        p.pu[LUMA_64x32].copy_pp = x265_blockcopy_pp_64x32_avx;
        p.pu[LUMA_64x48].copy_pp = x265_blockcopy_pp_64x48_avx;
        p.pu[LUMA_64x64].copy_pp = x265_blockcopy_pp_64x64_avx;

        p.pu[LUMA_48x64].copy_pp = x265_blockcopy_pp_48x64_avx;

        p.frameInitLowres = x265_frame_init_lowres_core_avx;
    }
    if (cpuMask & X265_CPU_XOP)
    {
        //p.pu[LUMA_4x4].satd = p.cu[BLOCK_4x4].sa8d = x265_pixel_satd_4x4_xop; this one is broken
        ALL_LUMA_PU(satd, pixel_satd, xop);
        ASSIGN_SA8D(xop);
        LUMA_VAR(xop);
        p.cu[BLOCK_8x8].sse_pp = x265_pixel_ssd_8x8_xop;
        p.cu[BLOCK_16x16].sse_pp = x265_pixel_ssd_16x16_xop;
        p.frameInitLowres = x265_frame_init_lowres_core_xop;
    }
#if X86_64
    if (cpuMask & X265_CPU_AVX2)
    {
        p.planecopy_sp = x265_downShift_16_avx2;

        p.cu[BLOCK_32x32].intra_pred[DC_IDX] = x265_intra_pred_dc32_avx2;

        p.cu[BLOCK_16x16].intra_pred[PLANAR_IDX] = x265_intra_pred_planar16_avx2;
        p.cu[BLOCK_32x32].intra_pred[PLANAR_IDX] = x265_intra_pred_planar32_avx2;

        p.idst4x4 = x265_idst4_avx2;
        p.dst4x4 = x265_dst4_avx2;
        p.scale2D_64to32 = x265_scale2D_64to32_avx2;
        p.saoCuOrgE0 = x265_saoCuOrgE0_avx2;
        p.saoCuOrgE1 = x265_saoCuOrgE1_avx2;
        p.saoCuOrgE1_2Rows = x265_saoCuOrgE1_2Rows_avx2;
        p.saoCuOrgE2[0] = x265_saoCuOrgE2_avx2;
        p.saoCuOrgE2[1] = x265_saoCuOrgE2_32_avx2;
        p.saoCuOrgE3[0] = x265_saoCuOrgE3_avx2;
        p.saoCuOrgE3[1] = x265_saoCuOrgE3_32_avx2;
        p.saoCuOrgB0 = x265_saoCuOrgB0_avx2;
        p.sign = x265_calSign_avx2;

        p.cu[BLOCK_4x4].psy_cost_ss = x265_psyCost_ss_4x4_avx2;
        p.cu[BLOCK_8x8].psy_cost_ss = x265_psyCost_ss_8x8_avx2;
        p.cu[BLOCK_16x16].psy_cost_ss = x265_psyCost_ss_16x16_avx2;
        p.cu[BLOCK_32x32].psy_cost_ss = x265_psyCost_ss_32x32_avx2;
        p.cu[BLOCK_64x64].psy_cost_ss = x265_psyCost_ss_64x64_avx2;

        p.cu[BLOCK_4x4].psy_cost_pp = x265_psyCost_pp_4x4_avx2;
        p.cu[BLOCK_8x8].psy_cost_pp = x265_psyCost_pp_8x8_avx2;
        p.cu[BLOCK_16x16].psy_cost_pp = x265_psyCost_pp_16x16_avx2;
        p.cu[BLOCK_32x32].psy_cost_pp = x265_psyCost_pp_32x32_avx2;
        p.cu[BLOCK_64x64].psy_cost_pp = x265_psyCost_pp_64x64_avx2;

        p.pu[LUMA_8x4].addAvg = x265_addAvg_8x4_avx2;
        p.pu[LUMA_8x8].addAvg = x265_addAvg_8x8_avx2;
        p.pu[LUMA_8x16].addAvg = x265_addAvg_8x16_avx2;
        p.pu[LUMA_8x32].addAvg = x265_addAvg_8x32_avx2;

        p.pu[LUMA_12x16].addAvg = x265_addAvg_12x16_avx2;

        p.pu[LUMA_16x4].addAvg = x265_addAvg_16x4_avx2;
        p.pu[LUMA_16x8].addAvg = x265_addAvg_16x8_avx2;
        p.pu[LUMA_16x12].addAvg = x265_addAvg_16x12_avx2;
        p.pu[LUMA_16x16].addAvg = x265_addAvg_16x16_avx2;
        p.pu[LUMA_16x32].addAvg = x265_addAvg_16x32_avx2;
        p.pu[LUMA_16x64].addAvg = x265_addAvg_16x64_avx2;

        p.pu[LUMA_24x32].addAvg = x265_addAvg_24x32_avx2;

        p.pu[LUMA_32x8].addAvg = x265_addAvg_32x8_avx2;
        p.pu[LUMA_32x16].addAvg = x265_addAvg_32x16_avx2;
        p.pu[LUMA_32x24].addAvg = x265_addAvg_32x24_avx2;
        p.pu[LUMA_32x32].addAvg = x265_addAvg_32x32_avx2;
        p.pu[LUMA_32x64].addAvg = x265_addAvg_32x64_avx2;

        p.pu[LUMA_48x64].addAvg = x265_addAvg_48x64_avx2;

        p.pu[LUMA_64x16].addAvg = x265_addAvg_64x16_avx2;
        p.pu[LUMA_64x32].addAvg = x265_addAvg_64x32_avx2;
        p.pu[LUMA_64x48].addAvg = x265_addAvg_64x48_avx2;
        p.pu[LUMA_64x64].addAvg = x265_addAvg_64x64_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].addAvg = x265_addAvg_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].addAvg = x265_addAvg_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].addAvg = x265_addAvg_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].addAvg = x265_addAvg_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].addAvg = x265_addAvg_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].addAvg = x265_addAvg_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].addAvg = x265_addAvg_12x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].addAvg = x265_addAvg_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].addAvg = x265_addAvg_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].addAvg = x265_addAvg_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].addAvg = x265_addAvg_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].addAvg = x265_addAvg_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].addAvg = x265_addAvg_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].addAvg = x265_addAvg_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].addAvg = x265_addAvg_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].addAvg = x265_addAvg_32x32_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].addAvg = x265_addAvg_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].addAvg = x265_addAvg_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].addAvg = x265_addAvg_8x12_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].addAvg = x265_addAvg_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].addAvg = x265_addAvg_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].addAvg = x265_addAvg_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].addAvg = x265_addAvg_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].addAvg = x265_addAvg_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].addAvg = x265_addAvg_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].addAvg = x265_addAvg_16x24_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].addAvg = x265_addAvg_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].addAvg = x265_addAvg_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].addAvg = x265_addAvg_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].addAvg = x265_addAvg_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].addAvg = x265_addAvg_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].addAvg = x265_addAvg_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].addAvg = x265_addAvg_32x64_avx2;

        p.cu[BLOCK_16x16].add_ps = x265_pixel_add_ps_16x16_avx2;
        p.cu[BLOCK_32x32].add_ps = x265_pixel_add_ps_32x32_avx2;
        p.cu[BLOCK_64x64].add_ps = x265_pixel_add_ps_64x64_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].add_ps = x265_pixel_add_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].add_ps = x265_pixel_add_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].add_ps = x265_pixel_add_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].add_ps = x265_pixel_add_ps_32x64_avx2;

        p.cu[BLOCK_16x16].sub_ps = x265_pixel_sub_ps_16x16_avx2;
        p.cu[BLOCK_32x32].sub_ps = x265_pixel_sub_ps_32x32_avx2;
        p.cu[BLOCK_64x64].sub_ps = x265_pixel_sub_ps_64x64_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].sub_ps = x265_pixel_sub_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].sub_ps = x265_pixel_sub_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].sub_ps = x265_pixel_sub_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].sub_ps = x265_pixel_sub_ps_32x64_avx2;

        p.pu[LUMA_16x4].pixelavg_pp = x265_pixel_avg_16x4_avx2;
        p.pu[LUMA_16x8].pixelavg_pp = x265_pixel_avg_16x8_avx2;
        p.pu[LUMA_16x12].pixelavg_pp = x265_pixel_avg_16x12_avx2;
        p.pu[LUMA_16x16].pixelavg_pp = x265_pixel_avg_16x16_avx2;
        p.pu[LUMA_16x32].pixelavg_pp = x265_pixel_avg_16x32_avx2;
        p.pu[LUMA_16x64].pixelavg_pp = x265_pixel_avg_16x64_avx2;

        p.pu[LUMA_32x64].pixelavg_pp = x265_pixel_avg_32x64_avx2;
        p.pu[LUMA_32x32].pixelavg_pp = x265_pixel_avg_32x32_avx2;
        p.pu[LUMA_32x24].pixelavg_pp = x265_pixel_avg_32x24_avx2;
        p.pu[LUMA_32x16].pixelavg_pp = x265_pixel_avg_32x16_avx2;
        p.pu[LUMA_32x8].pixelavg_pp = x265_pixel_avg_32x8_avx2;

        p.pu[LUMA_64x64].pixelavg_pp = x265_pixel_avg_64x64_avx2;
        p.pu[LUMA_64x48].pixelavg_pp = x265_pixel_avg_64x48_avx2;
        p.pu[LUMA_64x32].pixelavg_pp = x265_pixel_avg_64x32_avx2;
        p.pu[LUMA_64x16].pixelavg_pp = x265_pixel_avg_64x16_avx2;

        p.pu[LUMA_16x16].satd = x265_pixel_satd_16x16_avx2;
        p.pu[LUMA_16x8].satd  = x265_pixel_satd_16x8_avx2;
        p.pu[LUMA_8x16].satd  = x265_pixel_satd_8x16_avx2;
        p.pu[LUMA_8x8].satd   = x265_pixel_satd_8x8_avx2;

        p.pu[LUMA_16x4].satd  = x265_pixel_satd_16x4_avx2;
        p.pu[LUMA_16x12].satd = x265_pixel_satd_16x12_avx2;
        p.pu[LUMA_16x32].satd = x265_pixel_satd_16x32_avx2;
        p.pu[LUMA_16x64].satd = x265_pixel_satd_16x64_avx2;

        p.pu[LUMA_32x8].satd   = x265_pixel_satd_32x8_avx2;
        p.pu[LUMA_32x16].satd   = x265_pixel_satd_32x16_avx2;
        p.pu[LUMA_32x24].satd   = x265_pixel_satd_32x24_avx2;
        p.pu[LUMA_32x32].satd   = x265_pixel_satd_32x32_avx2;
        p.pu[LUMA_32x64].satd   = x265_pixel_satd_32x64_avx2;
        p.pu[LUMA_48x64].satd   = x265_pixel_satd_48x64_avx2;
        p.pu[LUMA_64x16].satd   = x265_pixel_satd_64x16_avx2;
        p.pu[LUMA_64x32].satd   = x265_pixel_satd_64x32_avx2;
        p.pu[LUMA_64x48].satd   = x265_pixel_satd_64x48_avx2;
        p.pu[LUMA_64x64].satd   = x265_pixel_satd_64x64_avx2;

        p.pu[LUMA_32x8].sad = x265_pixel_sad_32x8_avx2;
        p.pu[LUMA_32x16].sad = x265_pixel_sad_32x16_avx2;
        p.pu[LUMA_32x24].sad = x265_pixel_sad_32x24_avx2;
        p.pu[LUMA_32x32].sad = x265_pixel_sad_32x32_avx2;
        p.pu[LUMA_32x64].sad = x265_pixel_sad_32x64_avx2;
        p.pu[LUMA_48x64].sad = x265_pixel_sad_48x64_avx2;
        p.pu[LUMA_64x16].sad = x265_pixel_sad_64x16_avx2;
        p.pu[LUMA_64x32].sad = x265_pixel_sad_64x32_avx2;
        p.pu[LUMA_64x48].sad = x265_pixel_sad_64x48_avx2;
        p.pu[LUMA_64x64].sad = x265_pixel_sad_64x64_avx2;

        p.pu[LUMA_8x4].sad_x3 = x265_pixel_sad_x3_8x4_avx2;
        p.pu[LUMA_8x8].sad_x3 = x265_pixel_sad_x3_8x8_avx2;
        p.pu[LUMA_8x16].sad_x3 = x265_pixel_sad_x3_8x16_avx2;

        p.pu[LUMA_8x8].sad_x4 = x265_pixel_sad_x4_8x8_avx2;
        p.pu[LUMA_16x8].sad_x4  = x265_pixel_sad_x4_16x8_avx2;
        p.pu[LUMA_16x12].sad_x4 = x265_pixel_sad_x4_16x12_avx2;
        p.pu[LUMA_16x16].sad_x4 = x265_pixel_sad_x4_16x16_avx2;
        p.pu[LUMA_16x32].sad_x4 = x265_pixel_sad_x4_16x32_avx2;

        p.cu[BLOCK_16x16].sse_pp = x265_pixel_ssd_16x16_avx2;
        p.cu[BLOCK_32x32].sse_pp = x265_pixel_ssd_32x32_avx2;
        p.cu[BLOCK_64x64].sse_pp = x265_pixel_ssd_64x64_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].sse_pp = x265_pixel_ssd_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].sse_pp = x265_pixel_ssd_32x32_avx2;

        p.cu[BLOCK_16x16].ssd_s = x265_pixel_ssd_s_16_avx2;
        p.cu[BLOCK_32x32].ssd_s = x265_pixel_ssd_s_32_avx2;

        p.cu[BLOCK_8x8].copy_cnt = x265_copy_cnt_8_avx2;
        p.cu[BLOCK_16x16].copy_cnt = x265_copy_cnt_16_avx2;
        p.cu[BLOCK_32x32].copy_cnt = x265_copy_cnt_32_avx2;

        p.cu[BLOCK_16x16].blockfill_s = x265_blockfill_s_16x16_avx2;
        p.cu[BLOCK_32x32].blockfill_s = x265_blockfill_s_32x32_avx2;

        ALL_LUMA_TU_S(cpy1Dto2D_shl, cpy1Dto2D_shl_, avx2);
        ALL_LUMA_TU_S(cpy1Dto2D_shr, cpy1Dto2D_shr_, avx2);

        p.cu[BLOCK_8x8].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_8_avx2;
        p.cu[BLOCK_16x16].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_16_avx2;
        p.cu[BLOCK_32x32].cpy2Dto1D_shl = x265_cpy2Dto1D_shl_32_avx2;

        p.cu[BLOCK_8x8].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_8_avx2;
        p.cu[BLOCK_16x16].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_16_avx2;
        p.cu[BLOCK_32x32].cpy2Dto1D_shr = x265_cpy2Dto1D_shr_32_avx2;

        ALL_LUMA_TU(count_nonzero, count_nonzero, avx2);
        p.denoiseDct = x265_denoise_dct_avx2;
        p.quant = x265_quant_avx2;
        p.nquant = x265_nquant_avx2;
        p.dequant_normal = x265_dequant_normal_avx2;

        p.cu[BLOCK_16x16].calcresidual = x265_getResidual16_avx2;
        p.cu[BLOCK_32x32].calcresidual = x265_getResidual32_avx2;

        p.scale1D_128to64 = x265_scale1D_128to64_avx2;
        p.weight_pp = x265_weight_pp_avx2;
        p.weight_sp = x265_weight_sp_avx2;

        // intra_pred functions
        p.cu[BLOCK_4x4].intra_pred[3] = x265_intra_pred_ang4_3_avx2;
        p.cu[BLOCK_4x4].intra_pred[4] = x265_intra_pred_ang4_4_avx2;
        p.cu[BLOCK_4x4].intra_pred[5] = x265_intra_pred_ang4_5_avx2;
        p.cu[BLOCK_4x4].intra_pred[6] = x265_intra_pred_ang4_6_avx2;
        p.cu[BLOCK_4x4].intra_pred[7] = x265_intra_pred_ang4_7_avx2;
        p.cu[BLOCK_4x4].intra_pred[8] = x265_intra_pred_ang4_8_avx2;
        p.cu[BLOCK_4x4].intra_pred[9] = x265_intra_pred_ang4_9_avx2;
        p.cu[BLOCK_4x4].intra_pred[11] = x265_intra_pred_ang4_11_avx2;
        p.cu[BLOCK_4x4].intra_pred[12] = x265_intra_pred_ang4_12_avx2;
        p.cu[BLOCK_4x4].intra_pred[13] = x265_intra_pred_ang4_13_avx2;
        p.cu[BLOCK_4x4].intra_pred[14] = x265_intra_pred_ang4_14_avx2;
        p.cu[BLOCK_4x4].intra_pred[15] = x265_intra_pred_ang4_15_avx2;
        p.cu[BLOCK_4x4].intra_pred[16] = x265_intra_pred_ang4_16_avx2;
        p.cu[BLOCK_4x4].intra_pred[17] = x265_intra_pred_ang4_17_avx2;
        p.cu[BLOCK_4x4].intra_pred[19] = x265_intra_pred_ang4_19_avx2;
        p.cu[BLOCK_4x4].intra_pred[20] = x265_intra_pred_ang4_20_avx2;
        p.cu[BLOCK_4x4].intra_pred[21] = x265_intra_pred_ang4_21_avx2;
        p.cu[BLOCK_4x4].intra_pred[22] = x265_intra_pred_ang4_22_avx2;
        p.cu[BLOCK_4x4].intra_pred[23] = x265_intra_pred_ang4_23_avx2;
        p.cu[BLOCK_4x4].intra_pred[24] = x265_intra_pred_ang4_24_avx2;
        p.cu[BLOCK_4x4].intra_pred[25] = x265_intra_pred_ang4_25_avx2;
        p.cu[BLOCK_4x4].intra_pred[27] = x265_intra_pred_ang4_27_avx2;
        p.cu[BLOCK_4x4].intra_pred[28] = x265_intra_pred_ang4_28_avx2;
        p.cu[BLOCK_4x4].intra_pred[29] = x265_intra_pred_ang4_29_avx2;
        p.cu[BLOCK_4x4].intra_pred[30] = x265_intra_pred_ang4_30_avx2;
        p.cu[BLOCK_4x4].intra_pred[31] = x265_intra_pred_ang4_31_avx2;
        p.cu[BLOCK_4x4].intra_pred[32] = x265_intra_pred_ang4_32_avx2;
        p.cu[BLOCK_4x4].intra_pred[33] = x265_intra_pred_ang4_33_avx2;
        p.cu[BLOCK_8x8].intra_pred[3] = x265_intra_pred_ang8_3_avx2;
        p.cu[BLOCK_8x8].intra_pred[33] = x265_intra_pred_ang8_33_avx2;
        p.cu[BLOCK_8x8].intra_pred[4] = x265_intra_pred_ang8_4_avx2;
        p.cu[BLOCK_8x8].intra_pred[32] = x265_intra_pred_ang8_32_avx2;
        p.cu[BLOCK_8x8].intra_pred[5] = x265_intra_pred_ang8_5_avx2;
        p.cu[BLOCK_8x8].intra_pred[31] = x265_intra_pred_ang8_31_avx2;
        p.cu[BLOCK_8x8].intra_pred[30] = x265_intra_pred_ang8_30_avx2;
        p.cu[BLOCK_8x8].intra_pred[6] = x265_intra_pred_ang8_6_avx2;
        p.cu[BLOCK_8x8].intra_pred[7] = x265_intra_pred_ang8_7_avx2;
        p.cu[BLOCK_8x8].intra_pred[29] = x265_intra_pred_ang8_29_avx2;
        p.cu[BLOCK_8x8].intra_pred[8] = x265_intra_pred_ang8_8_avx2;
        p.cu[BLOCK_8x8].intra_pred[28] = x265_intra_pred_ang8_28_avx2;
        p.cu[BLOCK_8x8].intra_pred[9] = x265_intra_pred_ang8_9_avx2;
        p.cu[BLOCK_8x8].intra_pred[27] = x265_intra_pred_ang8_27_avx2;
        p.cu[BLOCK_8x8].intra_pred[25] = x265_intra_pred_ang8_25_avx2;
        p.cu[BLOCK_8x8].intra_pred[12] = x265_intra_pred_ang8_12_avx2;
        p.cu[BLOCK_8x8].intra_pred[24] = x265_intra_pred_ang8_24_avx2;
        p.cu[BLOCK_8x8].intra_pred[11] = x265_intra_pred_ang8_11_avx2;
        p.cu[BLOCK_8x8].intra_pred[13] = x265_intra_pred_ang8_13_avx2;
        p.cu[BLOCK_8x8].intra_pred[20] = x265_intra_pred_ang8_20_avx2;
        p.cu[BLOCK_8x8].intra_pred[21] = x265_intra_pred_ang8_21_avx2;
        p.cu[BLOCK_8x8].intra_pred[22] = x265_intra_pred_ang8_22_avx2;
        p.cu[BLOCK_8x8].intra_pred[23] = x265_intra_pred_ang8_23_avx2;
        p.cu[BLOCK_8x8].intra_pred[14] = x265_intra_pred_ang8_14_avx2;
        p.cu[BLOCK_8x8].intra_pred[15] = x265_intra_pred_ang8_15_avx2;
        p.cu[BLOCK_8x8].intra_pred[16] = x265_intra_pred_ang8_16_avx2;
        p.cu[BLOCK_16x16].intra_pred[3] = x265_intra_pred_ang16_3_avx2;
        p.cu[BLOCK_16x16].intra_pred[4] = x265_intra_pred_ang16_4_avx2;
        p.cu[BLOCK_16x16].intra_pred[5] = x265_intra_pred_ang16_5_avx2;
        p.cu[BLOCK_16x16].intra_pred[6] = x265_intra_pred_ang16_6_avx2;
        p.cu[BLOCK_16x16].intra_pred[7] = x265_intra_pred_ang16_7_avx2;
        p.cu[BLOCK_16x16].intra_pred[8] = x265_intra_pred_ang16_8_avx2;
        p.cu[BLOCK_16x16].intra_pred[9] = x265_intra_pred_ang16_9_avx2;
        p.cu[BLOCK_16x16].intra_pred[12] = x265_intra_pred_ang16_12_avx2;
        p.cu[BLOCK_16x16].intra_pred[11] = x265_intra_pred_ang16_11_avx2;
        p.cu[BLOCK_16x16].intra_pred[13] = x265_intra_pred_ang16_13_avx2;
        p.cu[BLOCK_16x16].intra_pred[25] = x265_intra_pred_ang16_25_avx2;
        p.cu[BLOCK_16x16].intra_pred[28] = x265_intra_pred_ang16_28_avx2;
        p.cu[BLOCK_16x16].intra_pred[27] = x265_intra_pred_ang16_27_avx2;
        p.cu[BLOCK_16x16].intra_pred[29] = x265_intra_pred_ang16_29_avx2;
        p.cu[BLOCK_16x16].intra_pred[30] = x265_intra_pred_ang16_30_avx2;
        p.cu[BLOCK_16x16].intra_pred[31] = x265_intra_pred_ang16_31_avx2;
        p.cu[BLOCK_16x16].intra_pred[32] = x265_intra_pred_ang16_32_avx2;
        p.cu[BLOCK_16x16].intra_pred[33] = x265_intra_pred_ang16_33_avx2;
        p.cu[BLOCK_16x16].intra_pred[24] = x265_intra_pred_ang16_24_avx2;
        p.cu[BLOCK_16x16].intra_pred[23] = x265_intra_pred_ang16_23_avx2;
        p.cu[BLOCK_16x16].intra_pred[22] = x265_intra_pred_ang16_22_avx2;
        p.cu[BLOCK_32x32].intra_pred[34] = x265_intra_pred_ang32_34_avx2;
        p.cu[BLOCK_32x32].intra_pred[2] = x265_intra_pred_ang32_2_avx2;
        p.cu[BLOCK_32x32].intra_pred[26] = x265_intra_pred_ang32_26_avx2;
        p.cu[BLOCK_32x32].intra_pred[27] = x265_intra_pred_ang32_27_avx2;
        p.cu[BLOCK_32x32].intra_pred[28] = x265_intra_pred_ang32_28_avx2;
        p.cu[BLOCK_32x32].intra_pred[29] = x265_intra_pred_ang32_29_avx2;
        p.cu[BLOCK_32x32].intra_pred[30] = x265_intra_pred_ang32_30_avx2;
        p.cu[BLOCK_32x32].intra_pred[31] = x265_intra_pred_ang32_31_avx2;
        p.cu[BLOCK_32x32].intra_pred[32] = x265_intra_pred_ang32_32_avx2;
        p.cu[BLOCK_32x32].intra_pred[33] = x265_intra_pred_ang32_33_avx2;
        p.cu[BLOCK_32x32].intra_pred[25] = x265_intra_pred_ang32_25_avx2;
        p.cu[BLOCK_32x32].intra_pred[24] = x265_intra_pred_ang32_24_avx2;
        p.cu[BLOCK_32x32].intra_pred[23] = x265_intra_pred_ang32_23_avx2;
        p.cu[BLOCK_32x32].intra_pred[22] = x265_intra_pred_ang32_22_avx2;
        p.cu[BLOCK_32x32].intra_pred[21] = x265_intra_pred_ang32_21_avx2;
        p.cu[BLOCK_32x32].intra_pred[18] = x265_intra_pred_ang32_18_avx2;

        // all_angs primitives
        p.cu[BLOCK_4x4].intra_pred_allangs = x265_all_angs_pred_4x4_avx2;

        // copy_sp primitives
        p.cu[BLOCK_16x16].copy_sp = x265_blockcopy_sp_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_16x16].copy_sp = x265_blockcopy_sp_16x16_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_16x32].copy_sp = x265_blockcopy_sp_16x32_avx2;

        p.cu[BLOCK_32x32].copy_sp = x265_blockcopy_sp_32x32_avx2;
        p.chroma[X265_CSP_I420].cu[BLOCK_420_32x32].copy_sp = x265_blockcopy_sp_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[BLOCK_422_32x64].copy_sp = x265_blockcopy_sp_32x64_avx2;

        p.cu[BLOCK_64x64].copy_sp = x265_blockcopy_sp_64x64_avx2;

        // copy_ps primitives
        p.cu[BLOCK_16x16].copy_ps = x265_blockcopy_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].cu[CHROMA_420_16x16].copy_ps = x265_blockcopy_ps_16x16_avx2;
        p.chroma[X265_CSP_I422].cu[CHROMA_422_16x32].copy_ps = x265_blockcopy_ps_16x32_avx2;

        p.cu[BLOCK_32x32].copy_ps = x265_blockcopy_ps_32x32_avx2;
        p.chroma[X265_CSP_I420].cu[CHROMA_420_32x32].copy_ps = x265_blockcopy_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].cu[CHROMA_422_32x64].copy_ps = x265_blockcopy_ps_32x64_avx2;

        p.cu[BLOCK_64x64].copy_ps = x265_blockcopy_ps_64x64_avx2;

        ALL_LUMA_TU_S(dct, dct, avx2);
        ALL_LUMA_TU_S(idct, idct, avx2);
        ALL_LUMA_CU_S(transpose, transpose, avx2);

        ALL_LUMA_PU(luma_vpp, interp_8tap_vert_pp, avx2);
        ALL_LUMA_PU(luma_vps, interp_8tap_vert_ps, avx2);
        ALL_LUMA_PU(luma_vsp, interp_8tap_vert_sp, avx2);
        ALL_LUMA_PU(luma_vss, interp_8tap_vert_ss, avx2);

        // missing 4x8, 4x16, 24x32, 12x16 for the fill set of luma PU
        p.pu[LUMA_4x4].luma_hpp = x265_interp_8tap_horiz_pp_4x4_avx2;
        p.pu[LUMA_4x8].luma_hpp = x265_interp_8tap_horiz_pp_4x8_avx2;
        p.pu[LUMA_4x16].luma_hpp = x265_interp_8tap_horiz_pp_4x16_avx2;
        p.pu[LUMA_8x4].luma_hpp = x265_interp_8tap_horiz_pp_8x4_avx2;
        p.pu[LUMA_8x8].luma_hpp = x265_interp_8tap_horiz_pp_8x8_avx2;
        p.pu[LUMA_8x16].luma_hpp = x265_interp_8tap_horiz_pp_8x16_avx2;
        p.pu[LUMA_8x32].luma_hpp = x265_interp_8tap_horiz_pp_8x32_avx2;
        p.pu[LUMA_16x4].luma_hpp = x265_interp_8tap_horiz_pp_16x4_avx2;
        p.pu[LUMA_16x8].luma_hpp = x265_interp_8tap_horiz_pp_16x8_avx2;
        p.pu[LUMA_16x12].luma_hpp = x265_interp_8tap_horiz_pp_16x12_avx2;
        p.pu[LUMA_16x16].luma_hpp = x265_interp_8tap_horiz_pp_16x16_avx2;
        p.pu[LUMA_16x32].luma_hpp = x265_interp_8tap_horiz_pp_16x32_avx2;
        p.pu[LUMA_16x64].luma_hpp = x265_interp_8tap_horiz_pp_16x64_avx2;
        p.pu[LUMA_32x8].luma_hpp  = x265_interp_8tap_horiz_pp_32x8_avx2;
        p.pu[LUMA_32x16].luma_hpp = x265_interp_8tap_horiz_pp_32x16_avx2;
        p.pu[LUMA_32x24].luma_hpp = x265_interp_8tap_horiz_pp_32x24_avx2;
        p.pu[LUMA_32x32].luma_hpp = x265_interp_8tap_horiz_pp_32x32_avx2;
        p.pu[LUMA_32x64].luma_hpp = x265_interp_8tap_horiz_pp_32x64_avx2;
        p.pu[LUMA_64x64].luma_hpp = x265_interp_8tap_horiz_pp_64x64_avx2;
        p.pu[LUMA_64x48].luma_hpp = x265_interp_8tap_horiz_pp_64x48_avx2;
        p.pu[LUMA_64x32].luma_hpp = x265_interp_8tap_horiz_pp_64x32_avx2;
        p.pu[LUMA_64x16].luma_hpp = x265_interp_8tap_horiz_pp_64x16_avx2;
        p.pu[LUMA_48x64].luma_hpp = x265_interp_8tap_horiz_pp_48x64_avx2;
        p.pu[LUMA_24x32].luma_hpp = x265_interp_8tap_horiz_pp_24x32_avx2;
        p.pu[LUMA_12x16].luma_hpp = x265_interp_8tap_horiz_pp_12x16_avx2;

        p.pu[LUMA_4x4].luma_hps = x265_interp_8tap_horiz_ps_4x4_avx2;
        p.pu[LUMA_4x8].luma_hps = x265_interp_8tap_horiz_ps_4x8_avx2;
        p.pu[LUMA_4x16].luma_hps = x265_interp_8tap_horiz_ps_4x16_avx2;
        p.pu[LUMA_8x4].luma_hps = x265_interp_8tap_horiz_ps_8x4_avx2;
        p.pu[LUMA_8x8].luma_hps = x265_interp_8tap_horiz_ps_8x8_avx2;
        p.pu[LUMA_8x16].luma_hps = x265_interp_8tap_horiz_ps_8x16_avx2;
        p.pu[LUMA_8x32].luma_hps = x265_interp_8tap_horiz_ps_8x32_avx2;
        p.pu[LUMA_16x8].luma_hps = x265_interp_8tap_horiz_ps_16x8_avx2;
        p.pu[LUMA_16x16].luma_hps = x265_interp_8tap_horiz_ps_16x16_avx2;
        p.pu[LUMA_16x12].luma_hps = x265_interp_8tap_horiz_ps_16x12_avx2;
        p.pu[LUMA_16x4].luma_hps = x265_interp_8tap_horiz_ps_16x4_avx2;
        p.pu[LUMA_16x32].luma_hps = x265_interp_8tap_horiz_ps_16x32_avx2;
        p.pu[LUMA_16x64].luma_hps = x265_interp_8tap_horiz_ps_16x64_avx2;

        p.pu[LUMA_32x32].luma_hps = x265_interp_8tap_horiz_ps_32x32_avx2;
        p.pu[LUMA_32x16].luma_hps = x265_interp_8tap_horiz_ps_32x16_avx2;
        p.pu[LUMA_32x24].luma_hps = x265_interp_8tap_horiz_ps_32x24_avx2;
        p.pu[LUMA_32x8].luma_hps = x265_interp_8tap_horiz_ps_32x8_avx2;
        p.pu[LUMA_32x64].luma_hps = x265_interp_8tap_horiz_ps_32x64_avx2;
        p.pu[LUMA_48x64].luma_hps = x265_interp_8tap_horiz_ps_48x64_avx2;
        p.pu[LUMA_64x64].luma_hps = x265_interp_8tap_horiz_ps_64x64_avx2;
        p.pu[LUMA_64x48].luma_hps = x265_interp_8tap_horiz_ps_64x48_avx2;
        p.pu[LUMA_64x32].luma_hps = x265_interp_8tap_horiz_ps_64x32_avx2;
        p.pu[LUMA_64x16].luma_hps = x265_interp_8tap_horiz_ps_64x16_avx2;
        p.pu[LUMA_12x16].luma_hps = x265_interp_8tap_horiz_ps_12x16_avx2;
        p.pu[LUMA_24x32].luma_hps = x265_interp_8tap_horiz_ps_24x32_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_hpp = x265_interp_4tap_horiz_pp_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_hpp = x265_interp_4tap_horiz_pp_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_hpp = x265_interp_4tap_horiz_pp_32x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_hpp = x265_interp_4tap_horiz_pp_16x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_hpp = x265_interp_4tap_horiz_pp_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_hpp = x265_interp_4tap_horiz_pp_2x8_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_hpp = x265_interp_4tap_horiz_pp_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_hpp = x265_interp_4tap_horiz_pp_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_hpp = x265_interp_4tap_horiz_pp_4x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_hpp = x265_interp_4tap_horiz_pp_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_hpp = x265_interp_4tap_horiz_pp_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_hpp = x265_interp_4tap_horiz_pp_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_hpp = x265_interp_4tap_horiz_pp_16x32_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_hpp = x265_interp_4tap_horiz_pp_6x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].filter_hpp = x265_interp_4tap_horiz_pp_6x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_hpp = x265_interp_4tap_horiz_pp_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_hpp = x265_interp_4tap_horiz_pp_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_hpp = x265_interp_4tap_horiz_pp_32x8_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_hpp = x265_interp_4tap_horiz_pp_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_hpp = x265_interp_4tap_horiz_pp_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_hpp = x265_interp_4tap_horiz_pp_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_hpp = x265_interp_4tap_horiz_pp_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_hpp = x265_interp_4tap_horiz_pp_8x32_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].filter_hpp = x265_interp_4tap_horiz_pp_12x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_hps = x265_interp_4tap_horiz_ps_32x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_hps = x265_interp_4tap_horiz_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_hps = x265_interp_4tap_horiz_ps_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_hps = x265_interp_4tap_horiz_ps_8x8_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_hps = x265_interp_4tap_horiz_ps_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_hps = x265_interp_4tap_horiz_ps_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_hps = x265_interp_4tap_horiz_ps_4x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_hps = x265_interp_4tap_horiz_ps_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_hps = x265_interp_4tap_horiz_ps_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_hps = x265_interp_4tap_horiz_ps_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_hps = x265_interp_4tap_horiz_ps_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_hps = x265_interp_4tap_horiz_ps_8x16_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_hps = x265_interp_4tap_horiz_ps_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_hps = x265_interp_4tap_horiz_ps_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_hps = x265_interp_4tap_horiz_ps_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_hps = x265_interp_4tap_horiz_ps_16x4_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_hps = x265_interp_4tap_horiz_ps_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_hps = x265_interp_4tap_horiz_ps_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_hps = x265_interp_4tap_horiz_ps_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_hps = x265_interp_4tap_horiz_ps_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_hps = x265_interp_4tap_horiz_ps_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_hps = x265_interp_4tap_horiz_ps_2x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_hps = x265_interp_4tap_horiz_ps_6x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_hpp = x265_interp_4tap_horiz_pp_24x32_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_vpp = x265_interp_4tap_vert_pp_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_vpp = x265_interp_4tap_vert_pp_2x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_vpp = x265_interp_4tap_vert_pp_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_vpp = x265_interp_4tap_vert_pp_6x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_vpp = x265_interp_4tap_vert_pp_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_vpp = x265_interp_4tap_vert_pp_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].filter_vpp = x265_interp_4tap_vert_pp_12x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_vpp = x265_interp_4tap_vert_pp_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_vpp = x265_interp_4tap_vert_pp_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_vpp = x265_interp_4tap_vert_pp_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_vpp = x265_interp_4tap_vert_pp_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_vpp = x265_interp_4tap_vert_pp_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_vpp = x265_interp_4tap_vert_pp_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_vpp = x265_interp_4tap_vert_pp_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_vpp = x265_interp_4tap_vert_pp_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_vpp = x265_interp_4tap_vert_pp_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_vpp = x265_interp_4tap_vert_pp_32x32_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_vps = x265_interp_4tap_vert_ps_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_vps = x265_interp_4tap_vert_ps_2x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_vps = x265_interp_4tap_vert_ps_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_vps = x265_interp_4tap_vert_ps_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_vps = x265_interp_4tap_vert_ps_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_vps = x265_interp_4tap_vert_ps_6x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_vps = x265_interp_4tap_vert_ps_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_vps = x265_interp_4tap_vert_ps_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_vps = x265_interp_4tap_vert_ps_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_vps = x265_interp_4tap_vert_ps_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_vps = x265_interp_4tap_vert_ps_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_vps = x265_interp_4tap_vert_ps_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].filter_vps = x265_interp_4tap_vert_ps_12x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_vps = x265_interp_4tap_vert_ps_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_vps = x265_interp_4tap_vert_ps_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_vps = x265_interp_4tap_vert_ps_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_vps = x265_interp_4tap_vert_ps_4x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_vps = x265_interp_4tap_vert_ps_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_vps = x265_interp_4tap_vert_ps_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_vps = x265_interp_4tap_vert_ps_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_vps = x265_interp_4tap_vert_ps_32x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_vps = x265_interp_4tap_vert_ps_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_vps = x265_interp_4tap_vert_ps_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_vps = x265_interp_4tap_vert_ps_32x8_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_vsp = x265_interp_4tap_vert_sp_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_vsp = x265_interp_4tap_vert_sp_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_vsp = x265_interp_4tap_vert_sp_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_vsp = x265_interp_4tap_vert_sp_32x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_vsp = x265_interp_4tap_vert_sp_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_vsp = x265_interp_4tap_vert_sp_2x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_vsp = x265_interp_4tap_vert_sp_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_vsp = x265_interp_4tap_vert_sp_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_vsp = x265_interp_4tap_vert_sp_4x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_vsp = x265_interp_4tap_vert_sp_6x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_vsp = x265_interp_4tap_vert_sp_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_vsp = x265_interp_4tap_vert_sp_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_vsp = x265_interp_4tap_vert_sp_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_vsp = x265_interp_4tap_vert_sp_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_vsp = x265_interp_4tap_vert_sp_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].filter_vsp = x265_interp_4tap_vert_sp_12x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_vsp = x265_interp_4tap_vert_sp_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_vsp = x265_interp_4tap_vert_sp_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_vsp = x265_interp_4tap_vert_sp_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_vsp = x265_interp_4tap_vert_sp_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_vsp = x265_interp_4tap_vert_sp_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_vsp = x265_interp_4tap_vert_sp_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_vsp = x265_interp_4tap_vert_sp_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_vsp = x265_interp_4tap_vert_sp_32x24_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x4].filter_vss = x265_interp_4tap_vert_ss_4x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x8].filter_vss = x265_interp_4tap_vert_ss_8x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x16].filter_vss = x265_interp_4tap_vert_ss_16x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].filter_vss = x265_interp_4tap_vert_ss_32x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x4].filter_vss = x265_interp_4tap_vert_ss_2x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_2x8].filter_vss = x265_interp_4tap_vert_ss_2x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x2].filter_vss = x265_interp_4tap_vert_ss_4x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x8].filter_vss = x265_interp_4tap_vert_ss_4x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_4x16].filter_vss = x265_interp_4tap_vert_ss_4x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_6x8].filter_vss = x265_interp_4tap_vert_ss_6x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x2].filter_vss = x265_interp_4tap_vert_ss_8x2_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x4].filter_vss = x265_interp_4tap_vert_ss_8x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x6].filter_vss = x265_interp_4tap_vert_ss_8x6_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x16].filter_vss = x265_interp_4tap_vert_ss_8x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_8x32].filter_vss = x265_interp_4tap_vert_ss_8x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_12x16].filter_vss = x265_interp_4tap_vert_ss_12x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x4].filter_vss = x265_interp_4tap_vert_ss_16x4_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x8].filter_vss = x265_interp_4tap_vert_ss_16x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x12].filter_vss = x265_interp_4tap_vert_ss_16x12_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_16x32].filter_vss = x265_interp_4tap_vert_ss_16x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].filter_vss = x265_interp_4tap_vert_ss_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].filter_vss = x265_interp_4tap_vert_ss_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].filter_vss = x265_interp_4tap_vert_ss_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].filter_vss = x265_interp_4tap_vert_ss_32x24_avx2;

        //i422 for chroma_vss
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_vss = x265_interp_4tap_vert_ss_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_vss = x265_interp_4tap_vert_ss_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_vss = x265_interp_4tap_vert_ss_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vss = x265_interp_4tap_vert_ss_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_vss = x265_interp_4tap_vert_ss_2x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_vss = x265_interp_4tap_vert_ss_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_vss = x265_interp_4tap_vert_ss_4x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_vss = x265_interp_4tap_vert_ss_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_vss = x265_interp_4tap_vert_ss_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_vss = x265_interp_4tap_vert_ss_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vss = x265_interp_4tap_vert_ss_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_vss = x265_interp_4tap_vert_ss_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_vss = x265_interp_4tap_vert_ss_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_vss = x265_interp_4tap_vert_ss_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].filter_vss = x265_interp_4tap_vert_ss_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_vss = x265_interp_4tap_vert_ss_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_vss = x265_interp_4tap_vert_ss_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_vss = x265_interp_4tap_vert_ss_8x12_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].filter_vss = x265_interp_4tap_vert_ss_6x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_vss = x265_interp_4tap_vert_ss_2x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_vss = x265_interp_4tap_vert_ss_16x24_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].filter_vss = x265_interp_4tap_vert_ss_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].filter_vss = x265_interp_4tap_vert_ss_4x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x4].filter_vss = x265_interp_4tap_vert_ss_2x4_avx2;

        //i444 for chroma_vss
        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_vss = x265_interp_4tap_vert_ss_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_vss = x265_interp_4tap_vert_ss_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_vss = x265_interp_4tap_vert_ss_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_vss = x265_interp_4tap_vert_ss_32x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x64].filter_vss = x265_interp_4tap_vert_ss_64x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_vss = x265_interp_4tap_vert_ss_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_vss = x265_interp_4tap_vert_ss_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_vss = x265_interp_4tap_vert_ss_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_vss = x265_interp_4tap_vert_ss_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_vss = x265_interp_4tap_vert_ss_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_vss = x265_interp_4tap_vert_ss_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_vss = x265_interp_4tap_vert_ss_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_12x16].filter_vss = x265_interp_4tap_vert_ss_12x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_vss = x265_interp_4tap_vert_ss_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_vss = x265_interp_4tap_vert_ss_4x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_vss = x265_interp_4tap_vert_ss_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_vss = x265_interp_4tap_vert_ss_24x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_vss = x265_interp_4tap_vert_ss_32x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_vss = x265_interp_4tap_vert_ss_8x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x32].filter_vss = x265_interp_4tap_vert_ss_64x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_vss = x265_interp_4tap_vert_ss_32x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x48].filter_vss = x265_interp_4tap_vert_ss_64x48_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_48x64].filter_vss = x265_interp_4tap_vert_ss_48x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x16].filter_vss = x265_interp_4tap_vert_ss_64x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_vss = x265_interp_4tap_vert_ss_16x64_avx2;

        p.pu[LUMA_16x16].luma_hvpp = x265_interp_8tap_hv_pp_16x16_avx2;

        p.pu[LUMA_32x8].convert_p2s = x265_filterPixelToShort_32x8_avx2;
        p.pu[LUMA_32x16].convert_p2s = x265_filterPixelToShort_32x16_avx2;
        p.pu[LUMA_32x24].convert_p2s = x265_filterPixelToShort_32x24_avx2;
        p.pu[LUMA_32x32].convert_p2s = x265_filterPixelToShort_32x32_avx2;
        p.pu[LUMA_32x64].convert_p2s = x265_filterPixelToShort_32x64_avx2;
        p.pu[LUMA_64x16].convert_p2s = x265_filterPixelToShort_64x16_avx2;
        p.pu[LUMA_64x32].convert_p2s = x265_filterPixelToShort_64x32_avx2;
        p.pu[LUMA_64x48].convert_p2s = x265_filterPixelToShort_64x48_avx2;
        p.pu[LUMA_64x64].convert_p2s = x265_filterPixelToShort_64x64_avx2;
        p.pu[LUMA_48x64].convert_p2s = x265_filterPixelToShort_48x64_avx2;
        p.pu[LUMA_24x32].convert_p2s = x265_filterPixelToShort_24x32_avx2;

        p.chroma[X265_CSP_I420].pu[CHROMA_420_24x32].p2s = x265_filterPixelToShort_24x32_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x8].p2s = x265_filterPixelToShort_32x8_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x16].p2s = x265_filterPixelToShort_32x16_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x24].p2s = x265_filterPixelToShort_32x24_avx2;
        p.chroma[X265_CSP_I420].pu[CHROMA_420_32x32].p2s = x265_filterPixelToShort_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].p2s = x265_filterPixelToShort_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].p2s = x265_filterPixelToShort_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].p2s = x265_filterPixelToShort_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].p2s = x265_filterPixelToShort_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].p2s = x265_filterPixelToShort_32x64_avx2;

        //i422 for chroma_hpp
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].filter_hpp = x265_interp_4tap_horiz_pp_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].filter_hpp = x265_interp_4tap_horiz_pp_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_hpp = x265_interp_4tap_horiz_pp_2x16_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_hpp = x265_interp_4tap_horiz_pp_2x16_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_hpp = x265_interp_4tap_horiz_pp_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_hpp = x265_interp_4tap_horiz_pp_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_hpp = x265_interp_4tap_horiz_pp_4x16_avx2;
        
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_hpp = x265_interp_4tap_horiz_pp_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_hpp = x265_interp_4tap_horiz_pp_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_hpp = x265_interp_4tap_horiz_pp_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_hpp = x265_interp_4tap_horiz_pp_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_hpp = x265_interp_4tap_horiz_pp_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_hpp = x265_interp_4tap_horiz_pp_8x12_avx2;
        
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].filter_hpp = x265_interp_4tap_horiz_pp_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_hpp = x265_interp_4tap_horiz_pp_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_hpp = x265_interp_4tap_horiz_pp_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_hpp = x265_interp_4tap_horiz_pp_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_hpp = x265_interp_4tap_horiz_pp_16x24_avx2;
        
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_hpp = x265_interp_4tap_horiz_pp_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_hpp = x265_interp_4tap_horiz_pp_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_hpp = x265_interp_4tap_horiz_pp_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_hpp = x265_interp_4tap_horiz_pp_32x48_avx2;
        
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_hpp = x265_interp_4tap_horiz_pp_2x8_avx2;

        //i444 filters hpp

        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_hpp = x265_interp_4tap_horiz_pp_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_hpp = x265_interp_4tap_horiz_pp_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_hpp = x265_interp_4tap_horiz_pp_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_hpp = x265_interp_4tap_horiz_pp_32x32_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_hpp = x265_interp_4tap_horiz_pp_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_hpp = x265_interp_4tap_horiz_pp_4x16_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_hpp = x265_interp_4tap_horiz_pp_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_hpp = x265_interp_4tap_horiz_pp_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_hpp = x265_interp_4tap_horiz_pp_8x32_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_hpp = x265_interp_4tap_horiz_pp_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_hpp = x265_interp_4tap_horiz_pp_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_hpp = x265_interp_4tap_horiz_pp_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_hpp = x265_interp_4tap_horiz_pp_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_hpp = x265_interp_4tap_horiz_pp_16x64_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_12x16].filter_hpp = x265_interp_4tap_horiz_pp_12x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_hpp = x265_interp_4tap_horiz_pp_24x32_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_hpp = x265_interp_4tap_horiz_pp_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_hpp = x265_interp_4tap_horiz_pp_32x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_hpp = x265_interp_4tap_horiz_pp_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_hpp = x265_interp_4tap_horiz_pp_32x8_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_64x64].filter_hpp = x265_interp_4tap_horiz_pp_64x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x32].filter_hpp = x265_interp_4tap_horiz_pp_64x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x48].filter_hpp = x265_interp_4tap_horiz_pp_64x48_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x16].filter_hpp = x265_interp_4tap_horiz_pp_64x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_48x64].filter_hpp = x265_interp_4tap_horiz_pp_48x64_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_hps = x265_interp_4tap_horiz_ps_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_hps = x265_interp_4tap_horiz_ps_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_hps = x265_interp_4tap_horiz_ps_4x16_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_hps = x265_interp_4tap_horiz_ps_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_hps = x265_interp_4tap_horiz_ps_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_hps = x265_interp_4tap_horiz_ps_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_hps = x265_interp_4tap_horiz_ps_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_hps = x265_interp_4tap_horiz_ps_8x64_avx2; //adding macro call
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_hps = x265_interp_4tap_horiz_ps_8x12_avx2; //adding macro call

        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].filter_hps = x265_interp_4tap_horiz_ps_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_hps = x265_interp_4tap_horiz_ps_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_hps = x265_interp_4tap_horiz_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_hps = x265_interp_4tap_horiz_ps_16x64_avx2;//adding macro call
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_hps = x265_interp_4tap_horiz_ps_16x24_avx2;//adding macro call

        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_hps = x265_interp_4tap_horiz_ps_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_hps = x265_interp_4tap_horiz_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_hps = x265_interp_4tap_horiz_ps_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_hps = x265_interp_4tap_horiz_ps_32x48_avx2;

        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_hps = x265_interp_4tap_horiz_ps_2x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].filter_hps = x265_interp_4tap_horiz_ps_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_hps = x265_interp_4tap_horiz_ps_2x16_avx2;

        //i444 chroma_hps
        p.chroma[X265_CSP_I444].pu[LUMA_64x32].filter_hps = x265_interp_4tap_horiz_ps_64x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x48].filter_hps = x265_interp_4tap_horiz_ps_64x48_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x16].filter_hps = x265_interp_4tap_horiz_ps_64x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x64].filter_hps = x265_interp_4tap_horiz_ps_64x64_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_hps = x265_interp_4tap_horiz_ps_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_hps = x265_interp_4tap_horiz_ps_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_hps = x265_interp_4tap_horiz_ps_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_hps = x265_interp_4tap_horiz_ps_32x32_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_hps = x265_interp_4tap_horiz_ps_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_hps = x265_interp_4tap_horiz_ps_4x16_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_hps = x265_interp_4tap_horiz_ps_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_hps = x265_interp_4tap_horiz_ps_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_hps = x265_interp_4tap_horiz_ps_8x32_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_hps = x265_interp_4tap_horiz_ps_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_hps = x265_interp_4tap_horiz_ps_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_hps = x265_interp_4tap_horiz_ps_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_hps = x265_interp_4tap_horiz_ps_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_hps = x265_interp_4tap_horiz_ps_16x64_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_hps = x265_interp_4tap_horiz_ps_24x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_48x64].filter_hps = x265_interp_4tap_horiz_ps_48x64_avx2;

        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_hps = x265_interp_4tap_horiz_ps_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_hps = x265_interp_4tap_horiz_ps_32x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_hps = x265_interp_4tap_horiz_ps_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_hps = x265_interp_4tap_horiz_ps_32x8_avx2;

        //i422 for chroma_vsp
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_vsp = x265_interp_4tap_vert_sp_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_vsp = x265_interp_4tap_vert_sp_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_vsp = x265_interp_4tap_vert_sp_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vsp = x265_interp_4tap_vert_sp_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_vsp = x265_interp_4tap_vert_sp_2x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_vsp = x265_interp_4tap_vert_sp_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_vsp = x265_interp_4tap_vert_sp_4x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_vsp = x265_interp_4tap_vert_sp_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_vsp = x265_interp_4tap_vert_sp_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_vsp = x265_interp_4tap_vert_sp_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vsp = x265_interp_4tap_vert_sp_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].filter_vsp = x265_interp_4tap_vert_sp_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_vsp = x265_interp_4tap_vert_sp_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_vsp = x265_interp_4tap_vert_sp_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_vsp = x265_interp_4tap_vert_sp_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_24x64].filter_vsp = x265_interp_4tap_vert_sp_24x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_vsp = x265_interp_4tap_vert_sp_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_vsp = x265_interp_4tap_vert_sp_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_vsp = x265_interp_4tap_vert_sp_8x12_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_6x16].filter_vsp = x265_interp_4tap_vert_sp_6x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x16].filter_vsp = x265_interp_4tap_vert_sp_2x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_vsp = x265_interp_4tap_vert_sp_16x24_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].filter_vsp = x265_interp_4tap_vert_sp_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x32].filter_vsp = x265_interp_4tap_vert_sp_4x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x4].filter_vsp = x265_interp_4tap_vert_sp_2x4_avx2;

        //i444 for chroma_vsp
        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_vsp = x265_interp_4tap_vert_sp_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_vsp = x265_interp_4tap_vert_sp_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_vsp = x265_interp_4tap_vert_sp_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_vsp = x265_interp_4tap_vert_sp_32x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x64].filter_vsp = x265_interp_4tap_vert_sp_64x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_vsp = x265_interp_4tap_vert_sp_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_vsp = x265_interp_4tap_vert_sp_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_vsp = x265_interp_4tap_vert_sp_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_vsp = x265_interp_4tap_vert_sp_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_vsp = x265_interp_4tap_vert_sp_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_vsp = x265_interp_4tap_vert_sp_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_vsp = x265_interp_4tap_vert_sp_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_12x16].filter_vsp = x265_interp_4tap_vert_sp_12x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_vsp = x265_interp_4tap_vert_sp_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_vsp = x265_interp_4tap_vert_sp_4x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_vsp = x265_interp_4tap_vert_sp_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_vsp = x265_interp_4tap_vert_sp_24x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_vsp = x265_interp_4tap_vert_sp_32x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_vsp = x265_interp_4tap_vert_sp_8x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x32].filter_vsp = x265_interp_4tap_vert_sp_64x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_vsp = x265_interp_4tap_vert_sp_32x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x48].filter_vsp = x265_interp_4tap_vert_sp_64x48_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_48x64].filter_vsp = x265_interp_4tap_vert_sp_48x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_64x16].filter_vsp = x265_interp_4tap_vert_sp_64x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_vsp = x265_interp_4tap_vert_sp_16x64_avx2;

        //i422 for chroma_vps
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_vps = x265_interp_4tap_vert_ps_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_vps = x265_interp_4tap_vert_ps_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_vps = x265_interp_4tap_vert_ps_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vps = x265_interp_4tap_vert_ps_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_vps = x265_interp_4tap_vert_ps_2x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_vps = x265_interp_4tap_vert_ps_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_vps = x265_interp_4tap_vert_ps_4x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_vps = x265_interp_4tap_vert_ps_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_vps = x265_interp_4tap_vert_ps_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_vps = x265_interp_4tap_vert_ps_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vps = x265_interp_4tap_vert_ps_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].filter_vps = x265_interp_4tap_vert_ps_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_vps = x265_interp_4tap_vert_ps_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_vps = x265_interp_4tap_vert_ps_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_vps = x265_interp_4tap_vert_ps_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_vps = x265_interp_4tap_vert_ps_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_vps = x265_interp_4tap_vert_ps_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].filter_vps = x265_interp_4tap_vert_ps_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_vps = x265_interp_4tap_vert_ps_8x12_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x4].filter_vps = x265_interp_4tap_vert_ps_2x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_vps = x265_interp_4tap_vert_ps_16x24_avx2;

        //i444 for chroma_vps
        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_vps = x265_interp_4tap_vert_ps_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_vps = x265_interp_4tap_vert_ps_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_vps = x265_interp_4tap_vert_ps_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_vps = x265_interp_4tap_vert_ps_32x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_vps = x265_interp_4tap_vert_ps_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_vps = x265_interp_4tap_vert_ps_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_vps = x265_interp_4tap_vert_ps_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_vps = x265_interp_4tap_vert_ps_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_vps = x265_interp_4tap_vert_ps_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_vps = x265_interp_4tap_vert_ps_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_vps = x265_interp_4tap_vert_ps_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_12x16].filter_vps = x265_interp_4tap_vert_ps_12x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_vps = x265_interp_4tap_vert_ps_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_vps = x265_interp_4tap_vert_ps_4x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_vps = x265_interp_4tap_vert_ps_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_vps = x265_interp_4tap_vert_ps_24x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_vps = x265_interp_4tap_vert_ps_32x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_vps = x265_interp_4tap_vert_ps_8x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_vps = x265_interp_4tap_vert_ps_16x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_vps = x265_interp_4tap_vert_ps_32x64_avx2;

        //i422 for chroma_vpp
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x32].filter_vpp = x265_interp_4tap_vert_pp_16x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x8].filter_vpp = x265_interp_4tap_vert_pp_2x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x16].filter_vpp = x265_interp_4tap_vert_pp_16x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x32].filter_vpp = x265_interp_4tap_vert_pp_32x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x8].filter_vpp = x265_interp_4tap_vert_pp_16x8_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x16].filter_vpp = x265_interp_4tap_vert_pp_32x16_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x64].filter_vpp = x265_interp_4tap_vert_pp_16x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x64].filter_vpp = x265_interp_4tap_vert_pp_8x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x64].filter_vpp = x265_interp_4tap_vert_pp_32x64_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_32x48].filter_vpp = x265_interp_4tap_vert_pp_32x48_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_12x32].filter_vpp = x265_interp_4tap_vert_pp_12x32_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_8x12].filter_vpp = x265_interp_4tap_vert_pp_8x12_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_2x4].filter_vpp = x265_interp_4tap_vert_pp_2x4_avx2;
        p.chroma[X265_CSP_I422].pu[CHROMA_422_16x24].filter_vpp = x265_interp_4tap_vert_pp_16x24_avx2;

        //i444 for chroma_vpp
        p.chroma[X265_CSP_I444].pu[LUMA_4x4].filter_vpp = x265_interp_4tap_vert_pp_4x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x8].filter_vpp = x265_interp_4tap_vert_pp_8x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x16].filter_vpp = x265_interp_4tap_vert_pp_16x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x32].filter_vpp = x265_interp_4tap_vert_pp_32x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x4].filter_vpp = x265_interp_4tap_vert_pp_8x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x8].filter_vpp = x265_interp_4tap_vert_pp_4x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x8].filter_vpp = x265_interp_4tap_vert_pp_16x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x16].filter_vpp = x265_interp_4tap_vert_pp_8x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x16].filter_vpp = x265_interp_4tap_vert_pp_32x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x32].filter_vpp = x265_interp_4tap_vert_pp_16x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x12].filter_vpp = x265_interp_4tap_vert_pp_16x12_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_12x16].filter_vpp = x265_interp_4tap_vert_pp_12x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x4].filter_vpp = x265_interp_4tap_vert_pp_16x4_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_4x16].filter_vpp = x265_interp_4tap_vert_pp_4x16_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x24].filter_vpp = x265_interp_4tap_vert_pp_32x24_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_24x32].filter_vpp = x265_interp_4tap_vert_pp_24x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x8].filter_vpp = x265_interp_4tap_vert_pp_32x8_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_8x32].filter_vpp = x265_interp_4tap_vert_pp_8x32_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_16x64].filter_vpp = x265_interp_4tap_vert_pp_16x64_avx2;
        p.chroma[X265_CSP_I444].pu[LUMA_32x64].filter_vpp = x265_interp_4tap_vert_pp_32x64_avx2;

        if (cpuMask & X265_CPU_BMI2)
            p.scanPosLast = x265_scanPosLast_avx2_bmi2;
    }
#endif
}
#endif // if HIGH_BIT_DEPTH

} // namespace x265

extern "C" {
#ifdef __INTEL_COMPILER

/* Agner's patch to Intel's CPU dispatcher from pages 131-132 of
 * http://agner.org/optimize/optimizing_cpp.pdf (2011-01-30)
 * adapted to x265's cpu schema. */

// Global variable indicating cpu
int __intel_cpu_indicator = 0;
// CPU dispatcher function
void x265_intel_cpu_indicator_init(void)
{
    uint32_t cpu = x265::cpu_detect();

    if (cpu & X265_CPU_AVX)
        __intel_cpu_indicator = 0x20000;
    else if (cpu & X265_CPU_SSE42)
        __intel_cpu_indicator = 0x8000;
    else if (cpu & X265_CPU_SSE4)
        __intel_cpu_indicator = 0x2000;
    else if (cpu & X265_CPU_SSSE3)
        __intel_cpu_indicator = 0x1000;
    else if (cpu & X265_CPU_SSE3)
        __intel_cpu_indicator = 0x800;
    else if (cpu & X265_CPU_SSE2 && !(cpu & X265_CPU_SSE2_IS_SLOW))
        __intel_cpu_indicator = 0x200;
    else if (cpu & X265_CPU_SSE)
        __intel_cpu_indicator = 0x80;
    else if (cpu & X265_CPU_MMX2)
        __intel_cpu_indicator = 8;
    else
        __intel_cpu_indicator = 1;
}

/* __intel_cpu_indicator_init appears to have a non-standard calling convention that
 * assumes certain registers aren't preserved, so we'll route it through a function
 * that backs up all the registers. */
void __intel_cpu_indicator_init(void)
{
    x265_safe_intel_cpu_indicator_init();
}

#else // ifdef __INTEL_COMPILER
void x265_intel_cpu_indicator_init(void) {}

#endif // ifdef __INTEL_COMPILER
}
