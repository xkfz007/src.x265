/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
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

#ifndef X265_MC_H
#define X265_MC_H

#define LOWRES(cpu) \
    void x265_frame_init_lowres_core_ ## cpu(const pixel* src0, pixel* dst0, pixel* dsth, pixel* dstv, pixel* dstc, \
                                             intptr_t src_stride, intptr_t dst_stride, int width, int height);
LOWRES(mmx2)
LOWRES(sse2)
LOWRES(ssse3)
LOWRES(avx)
LOWRES(xop)

#define DECL_SUF(func, args) \
    void func ## _mmx2 args; \
    void func ## _sse2 args; \
    void func ## _ssse3 args;
DECL_SUF(x265_pixel_avg_64x64, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_64x48, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_64x32, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_64x16, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_48x64, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_32x64, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_32x32, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_32x24, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_32x16, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_32x8,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_24x32, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x64, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x32, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x16, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x12, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x8,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_16x4,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_12x16, (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_8x32,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_8x16,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_8x8,   (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_8x4,   (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_4x16,  (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_4x8,   (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))
DECL_SUF(x265_pixel_avg_4x4,   (pixel*, intptr_t, const pixel*, intptr_t, const pixel*, intptr_t, int))

#undef LOWRES
#undef DECL_SUF

#endif // ifndef X265_MC_H
