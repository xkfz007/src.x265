/*
 * HEVC video decoder
 *
 * Copyright (C) 2012 - 2013 Guillaume Martres
 * Copyright (C) 2013 - 2014 Pierre-Edouard Lepere
 *
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef AVCODEC_X86_HEVCDSP_H
#define AVCODEC_X86_HEVCDSP_H

struct SAOParams;
struct AVFrame;
struct UpsamplInf;
struct HEVCWindow;

//#define OPTI_ASM

#define idct_dc_proto(size, bitd, opt) \
		void ff_hevc_idct##size##_dc_add_##bitd##_##opt(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride)

#define PEL_LINK2(dst, idx1, idx2, idx3, name, D, opt) \
dst[idx1][idx2][idx3] = ff_hevc_put_hevc_ ## name ## _ ## D ## _ ## opt; \
dst ## _bi[idx1][idx2][idx3] = ff_hevc_put_hevc_bi_ ## name ## _ ## D ## _ ## opt; \
dst ## _uni[idx1][idx2][idx3] = ff_hevc_put_hevc_uni_ ## name ## _ ## D ## _ ## opt; \
dst ## _uni_w[idx1][idx2][idx3] = ff_hevc_put_hevc_uni_w_ ## name ## _ ## D ## _ ## opt; \
dst ## _bi_w[idx1][idx2][idx3] = ff_hevc_put_hevc_bi_w_ ## name ## _ ## D ## _ ## opt

#ifdef OPTI_ASM
#define PEL_LINK(dst, idx1, idx2, idx3, name, D, opt) \
PEL_LINK2(dst, idx1, idx2, idx3, name, D, opt)
#else
#define PEL_LINK(dst, idx1, idx2, idx3, name, D, opt) \
PEL_LINK2(dst, idx1, idx2, idx3, name, D, sse)
#endif

#define PEL_PROTOTYPE2(name, D, opt) \
void ff_hevc_put_hevc_ ## name ## _ ## D ## _##opt(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,int width); \
void ff_hevc_put_hevc_bi_ ## name ## _ ## D ## _##opt(uint8_t *_dst, ptrdiff_t _dststride, uint8_t *_src, ptrdiff_t _srcstride, int16_t *src2, ptrdiff_t src2stride, int height, intptr_t mx, intptr_t my, int width); \
void ff_hevc_put_hevc_uni_ ## name ## _ ## D ## _##opt(uint8_t *_dst, ptrdiff_t _dststride, uint8_t *_src, ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my, int width); \
void ff_hevc_put_hevc_uni_w_ ## name ## _ ## D ## _##opt(uint8_t *_dst, ptrdiff_t _dststride, uint8_t *_src, ptrdiff_t _srcstride, int height, int denom, int wx, int ox, intptr_t mx, intptr_t my, int width); \
void ff_hevc_put_hevc_bi_w_ ## name ## _ ## D ## _##opt(uint8_t *_dst, ptrdiff_t _dststride, uint8_t *_src, ptrdiff_t _srcstride, int16_t *src2, ptrdiff_t src2stride, int height, int denom, int wx0, int wx1, int ox0, int ox1, intptr_t mx, intptr_t my, int width)


#define WEIGHTING_PROTOTYPE2(width, bitd, opt) \
void ff_hevc_put_hevc_uni_w##width##_##bitd##_##opt(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride, int height, int denom,  int _wx, int _ox); \
void ff_hevc_put_hevc_bi_w##width##_##bitd##_##opt(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride, int16_t *_src2, ptrdiff_t _src2stride, int height, int denom,  int _wx0,  int _wx1, int _ox0, int _ox1)

#ifdef OPTI_ASM
#define WEIGHTING_PROTOTYPE(width, bitd, opt) \
		WEIGHTING_PROTOTYPE2(width, bitd, opt)
#else
#define WEIGHTING_PROTOTYPE(width, bitd, opt) \
		WEIGHTING_PROTOTYPE2(width, bitd, sse)
#endif

#ifdef OPTI_ASM
#define PEL_PROTOTYPE(name, D, opt) \
PEL_PROTOTYPE2(name, D, opt)
#else
#define PEL_PROTOTYPE(name, D, opt) \
PEL_PROTOTYPE2(name, D, sse)
#endif

///////////////////////////////////////////////////////////////////////////////
//IDCT functions
///////////////////////////////////////////////////////////////////////////////
void ff_hevc_transform_skip_8_sse(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);

void ff_hevc_transform_4x4_luma_8_sse4(int16_t *coeffs);
void ff_hevc_transform_4x4_luma_10_sse4(int16_t *coeffs);
void ff_hevc_transform_4x4_luma_12_sse4(int16_t *coeffs);

#define IDCT_FUNC(s, b) void ff_hevc_transform_ ## s ## x ## s ##_## b ##_sse4\
            (int16_t *coeffs, int col_limit);

IDCT_FUNC(4, 8)
IDCT_FUNC(4, 10)
IDCT_FUNC(4, 12)
IDCT_FUNC(8, 8)
IDCT_FUNC(8, 10)
IDCT_FUNC(8, 12)
IDCT_FUNC(16, 8)
IDCT_FUNC(16, 10)
IDCT_FUNC(16, 12)
IDCT_FUNC(32, 8)
IDCT_FUNC(32, 10)
IDCT_FUNC(32, 12)

void ff_hevc_transform_4x4_add_8_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_8x8_add_8_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_16x16_add_8_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_32x32_add_8_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);

void ff_hevc_transform_4x4_add_10_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_8x8_add_10_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_16x16_add_10_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_32x32_add_10_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);

void ff_hevc_transform_4x4_add_12_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_8x8_add_12_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_16x16_add_12_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);
void ff_hevc_transform_32x32_add_12_sse4(uint8_t *dst, int16_t *coeffs, ptrdiff_t stride);

///////////////////////////////////////////////////////////////////////////////
// MC functions
///////////////////////////////////////////////////////////////////////////////
#define EPEL_PROTOTYPES(fname, bitd, opt) \
        PEL_PROTOTYPE(fname##4,  bitd, opt); \
        PEL_PROTOTYPE(fname##6,  bitd, opt); \
        PEL_PROTOTYPE(fname##8,  bitd, opt); \
        PEL_PROTOTYPE(fname##12, bitd, opt); \
        PEL_PROTOTYPE(fname##16, bitd, opt); \
        PEL_PROTOTYPE(fname##24, bitd, opt); \
        PEL_PROTOTYPE(fname##32, bitd, opt); \
        PEL_PROTOTYPE(fname##48, bitd, opt); \
        PEL_PROTOTYPE(fname##64, bitd, opt)

#define QPEL_PROTOTYPES(fname, bitd, opt) \
        PEL_PROTOTYPE(fname##4,  bitd, opt); \
        PEL_PROTOTYPE(fname##8,  bitd, opt); \
        PEL_PROTOTYPE(fname##12, bitd, opt); \
        PEL_PROTOTYPE(fname##16, bitd, opt); \
        PEL_PROTOTYPE(fname##24, bitd, opt); \
        PEL_PROTOTYPE(fname##32, bitd, opt); \
        PEL_PROTOTYPE(fname##48, bitd, opt); \
        PEL_PROTOTYPE(fname##64, bitd, opt)

#define WEIGHTING_PROTOTYPES(bitd, opt) \
        WEIGHTING_PROTOTYPE(2, bitd, opt); \
        WEIGHTING_PROTOTYPE(4, bitd, opt); \
        WEIGHTING_PROTOTYPE(6, bitd, opt); \
        WEIGHTING_PROTOTYPE(8, bitd, opt); \
        WEIGHTING_PROTOTYPE(12, bitd, opt); \
        WEIGHTING_PROTOTYPE(16, bitd, opt); \
        WEIGHTING_PROTOTYPE(24, bitd, opt); \
        WEIGHTING_PROTOTYPE(32, bitd, opt); \
        WEIGHTING_PROTOTYPE(48, bitd, opt); \
        WEIGHTING_PROTOTYPE(64, bitd, opt)

///////////////////////////////////////////////////////////////////////////////
// QPEL_PIXELS EPEL_PIXELS
///////////////////////////////////////////////////////////////////////////////
EPEL_PROTOTYPES(pel_pixels ,  8, sse4);
EPEL_PROTOTYPES(pel_pixels , 10, sse4);
EPEL_PROTOTYPES(pel_pixels , 12, sse4);
///////////////////////////////////////////////////////////////////////////////
// EPEL
///////////////////////////////////////////////////////////////////////////////
EPEL_PROTOTYPES(epel_h ,  8, sse4);
EPEL_PROTOTYPES(epel_h , 10, sse4);
EPEL_PROTOTYPES(epel_h , 12, sse4);

EPEL_PROTOTYPES(epel_v ,  8, sse4);
EPEL_PROTOTYPES(epel_v , 10, sse4);
EPEL_PROTOTYPES(epel_v , 12, sse4);

EPEL_PROTOTYPES(epel_hv ,  8, sse4);
EPEL_PROTOTYPES(epel_hv , 10, sse4);
EPEL_PROTOTYPES(epel_hv , 12, sse4);

///////////////////////////////////////////////////////////////////////////////
// QPEL
///////////////////////////////////////////////////////////////////////////////
QPEL_PROTOTYPES(qpel_h ,  8, sse4);
QPEL_PROTOTYPES(qpel_h , 10, sse4);
QPEL_PROTOTYPES(qpel_h , 12, sse4);

QPEL_PROTOTYPES(qpel_v,  8, sse4);
QPEL_PROTOTYPES(qpel_v, 10, sse4);
QPEL_PROTOTYPES(qpel_v, 12, sse4);

QPEL_PROTOTYPES(qpel_hv,  8, sse4);
QPEL_PROTOTYPES(qpel_hv, 10, sse4);
QPEL_PROTOTYPES(qpel_hv, 12, sse4);


WEIGHTING_PROTOTYPES(8, sse4);
WEIGHTING_PROTOTYPES(10, sse4);
WEIGHTING_PROTOTYPES(12, sse4);

///////////////////////////////////////////////////////////////////////////////
// IDCT
///////////////////////////////////////////////////////////////////////////////
idct_dc_proto(4, 8,mmxext);
idct_dc_proto(8, 8,mmxext);
idct_dc_proto(16,8,  sse2);
idct_dc_proto(32,8,  sse2);

idct_dc_proto(32,8,  avx2);


idct_dc_proto(4, 10,mmxext);
idct_dc_proto(8, 10,  sse2);
idct_dc_proto(16,10,  sse2);
idct_dc_proto(32,10,  sse2);
idct_dc_proto(8, 10,   avx);
idct_dc_proto(16,10,   avx);
idct_dc_proto(32,10,   avx);

idct_dc_proto(16,10,  avx2);
idct_dc_proto(32,10,  avx2);

///////////////////////////////////////////////////////////////////////////////
// SAO functions
///////////////////////////////////////////////////////////////////////////////
//#ifndef OPTI_ASM
void ff_hevc_sao_edge_filter_0_8_sse(uint8_t *_dst, uint8_t *_src,
                                     ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                     int *borders, int _width, int _height, int c_idx,
                                     uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);
void ff_hevc_sao_edge_filter_1_8_sse(uint8_t *_dst, uint8_t *_src,
                                     ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                     int *borders, int _width, int _height, int c_idx,
                                     uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);
//#endif
void ff_hevc_sao_edge_filter_0_10_sse(uint8_t *_dst, uint8_t *_src,
                                      ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                      int *borders, int _width, int _height, int c_idx,
                                      uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);
void ff_hevc_sao_edge_filter_1_10_sse(uint8_t *_dst, uint8_t *_src,
                                      ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                      int *borders, int _width, int _height, int c_idx,
                                      uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);
void ff_hevc_sao_edge_filter_0_12_sse(uint8_t *_dst, uint8_t *_src,
                                      ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                      int *borders, int _width, int _height, int c_idx,
                                      uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);
void ff_hevc_sao_edge_filter_1_12_sse(uint8_t *_dst, uint8_t *_src,
                                      ptrdiff_t _stride_dst, ptrdiff_t _stride_src, struct SAOParams *sao,
                                      int *borders, int _width, int _height, int c_idx,
                                      uint8_t *vert_edge, uint8_t *horiz_edge, uint8_t *diag_edge);

//#ifndef OPTI_ASM
void ff_hevc_sao_band_filter_0_8_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                     struct SAOParams *sao, int *borders, int width, int height, int c_idx);
//#endif
void ff_hevc_sao_band_filter_0_10_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                      struct SAOParams *sao, int *borders, int width, int height, int c_idx);
void ff_hevc_sao_band_filter_0_12_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                      struct SAOParams *sao, int *borders, int width, int height, int c_idx);

//#ifdef SVC_EXTENSION

    void ff_upsample_filter_block_luma_h_all_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
            int x_EL, int x_BL, int block_w, int block_h, int widthEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_h_all_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
                int x_EL, int x_BL, int block_w, int block_h, int widthEL,
                const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_luma_v_all_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_v_all_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);

    void ff_upsample_filter_block_luma_h_x2_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
            int x_EL, int x_BL, int block_w, int block_h, int widthEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_h_x2_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
                int x_EL, int x_BL, int block_w, int block_h, int widthEL,
                const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_luma_v_x2_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_v_x2_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);

    void ff_upsample_filter_block_luma_h_x1_5_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
            int x_EL, int x_BL, int block_w, int block_h, int widthEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_h_x1_5_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
                int x_EL, int x_BL, int block_w, int block_h, int widthEL,
                const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_luma_v_x1_5_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_cr_v_x1_5_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);

    void ff_upsample_filter_block_luma_h_8_8_sse( int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
            int x_EL, int x_BL, int block_w, int block_h, int widthEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
    void ff_upsample_filter_block_luma_v_8_8_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
            int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
   void ff_upsample_filter_block_cr_h_8_8_sse( int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,
            int x_EL, int x_BL, int block_w, int block_h, int widthEL,
            const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
   void ff_upsample_filter_block_cr_v_8_8_sse(uint8_t *dst, ptrdiff_t dststride, int16_t *_src, ptrdiff_t _srcstride,
           int y_BL, int x_EL, int y_EL, int block_w, int block_h, int widthEL, int heightEL,
           const struct HEVCWindow *Enhscal, struct UpsamplInf *up_info);
//#endif

#endif // AVCODEC_X86_HEVCDSP_H
