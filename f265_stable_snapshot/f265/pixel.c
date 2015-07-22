// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Pixel-processing functions.

#include "f265/enc.h"

// Normalize distortion (sad/satd) to 8 bit depth.
#ifdef F265_LBD
#define DIST_BD_NORM(sum, shift) (sum)
#else
#define DIST_BD_NORM(sum, shift) ((sum) >> (shift))
#endif

// Store the block specified at the destination specified. Both the destination
// and the source can be unaligned.
void venc_copy_block(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int width, int height)
{
    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = src[x];
}

// Same for 16-bit pixels.
void venc_copy_block_s16(int16_t *dst, int dst_stride, int16_t *src, int src_stride, int width, int height)
{
    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = src[x];
}

// Transpose and store the block specified at the destination specified.
void venc_transpose_block(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int width, int height)
{
    for (int y = 0; y < height; y++, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x*dst_stride + y] = src[x];
}

// Pad a plane horizontally. Note that the assembly implementation accesses
// pixels as uint32_t in the plane.
static void noinline venc_pad_horizontal(f265_pix *dst, f265_pix *src, int32_t stride, int32_t nb_rows, int32_t len)
{
    // The early return is needed both for performance and the assembly
    // implementation.
    if (unlikely(!nb_rows|!len)) return;

    // Assembly hook needed for aligned values.
    for (int i = 0; i < nb_rows; i++, dst += stride, src += stride)
        for (int j = 0; j < len; j++) dst[j] = *src;
}

// Pad a plane vertically.
static void noinline venc_pad_vertical(f265_pix *dst, f265_pix *src, int32_t stride, int32_t nb_rows, int32_t len)
{
    for (int i = 0; i < nb_rows; i++, dst += stride) memcpy(dst, src, len*F265_PSIZE);
}

// Pad all sides of a plane (left, right, top, bottom) with its boundary pixels.
// 'bX' are the boundary offsets, 'pX' are the padding sizes (from the
// boundary).
void venc_pad_plane(f265_pix *plane, int32_t stride,
                    int32_t bl, int32_t br, int32_t bt, int32_t bb,
                    int32_t pl, int32_t pr, int32_t pt, int32_t pb)
{
    // Boundary top row aligned on the left padding.
    f265_pix *top_left = plane + bt*stride + bl - pl;

    // Boundary bottom row aligned on the left padding.
    f265_pix *bot_left = plane + bb*stride + bl - pl;

    // Boundary top row aligned on the right boundary.
    f265_pix *top_right = plane + bt*stride + br;

    // Number of rows within the plane boundaries.
    int32_t nb_plane_rows = bb - bt + 1;

    // Size of padded line.
    int32_t padded_line_size = br - bl + 1 + pr + pl;

    // Pad left, right, top, bottom.
    venc_pad_horizontal(top_left, top_left + pl, stride, nb_plane_rows, pl);
    venc_pad_horizontal(top_right + 1, top_right, stride, nb_plane_rows, pr);
    venc_pad_vertical(top_left - pt*stride, top_left, stride, pt, padded_line_size);
    venc_pad_vertical(bot_left + stride, bot_left, stride, pb, padded_line_size);
}

// Pad the reconstructed YUV planes.
void venc_pad_rec_yuv_planes(f265_enc_thread *t)
{
    f265_frame *f = t->src_frame;
    int32_t br = t->enc->gd.pix_dim[0]-1, bb = t->enc->gd.pix_dim[1]-1;
    int32_t stride = t->me.ref_stride;
    f265_pix *planes[3] = { f->rec_planes[0], f->rec_planes[4], f->rec_planes[5] };
    for (int i = 0; i < 3; i++)
    {
        int32_t pad = F265_LUMA_PLANE_PADDING>>!!i;
        venc_pad_plane(planes[i], stride, 0, br>>!!i, 0, bb>>!!i, pad, pad, pad, pad);
    }
}

// Compute the luma halfpel values of the block specified. The block height is
// arbitrary. The block width must be aligned on a 16-byte boundary.
//
// The 8-taps filter processes 4 pixels on each side of a halfpel. The position
// of a halfpel in a plane corresponds to the position of the fullpel to the
// left/top of the halfpel.
//
// The function correctly computes all the halfpels on the top and bottom sides
// of the block. However, the function incorrectly computes the three left-most
// columns and the four right-most columns, due to optimizations.
//
// WARNING: this function does not yet work for HBD!
void venc_compute_hpel(f265_pix *planes[4], int stride, int width, int height, uint8_t *spill)
{
    #define FILTER_8TAPS(src, stride)\
        *(src-3*stride)*-1 +\
        *(src-2*stride)*4 +\
        *(src-1*stride)*-11 +\
        *(src+0*stride)*40 +\
        *(src+1*stride)*40 +\
        *(src+2*stride)*-11 +\
        *(src+3*stride)*4 +\
        *(src+4*stride)*-1

    #ifdef F265_LBD
    f265_pix *f_plane = planes[0], *h_plane = planes[1], *v_plane = planes[2], *d_plane = planes[3];

    // Store the vertical halfpels of the current row to compute the diagonal
    // halfpels. We are not following the specification here. We interpolate the
    // vertical halfpels horizontally to get the diagonal halfpels. In LBD this
    // yields the correct result because there is no shifting involved. In HBD
    // we have to follow the specification.
    int16_t *row_buf = (int16_t*)spill;

    // Process each row.
    for (int y = 0; y < height; y++, f_plane += stride, h_plane += stride, v_plane += stride, d_plane += stride)
    {
        // Compute the horizontal and vertical halfpels (correctly).
        for (int x = 0; x < width; x++)
        {
            int hpel = FILTER_8TAPS(f_plane + x, 1);
            int vpel = FILTER_8TAPS(f_plane + x, stride);
            h_plane[x] = F265_CLAMP((hpel+32)>>6, 0, 255);
            v_plane[x] = F265_CLAMP((vpel+32)>>6, 0, 255);
            row_buf[x] = vpel;
        }

        // Compute the diagonal halfpels (incorrect values on boundaries).
        for (int x = 3; x < width - 4; x++)
        {
            int dpel = FILTER_8TAPS(row_buf + x, 1);
            d_plane[x] = F265_CLAMP((dpel+2048)>>12, 0, 255);
        }
    }
    #endif

    #undef FILTER_8TAPS
}

// Compute the halfpel planes of the reconstructed frame.
void venc_compute_hpel_planes(f265_enc_thread *t)
{
    f265_frame *f = t->src_frame;
    int pad = F265_LUMA_PLANE_PADDING;
    int width = t->enc->gd.pix_dim[0], height = t->enc->gd.pix_dim[1];
    int aligned_width = F265_ALIGN_VAL(width, 16);
    int stride = t->me.ref_stride;

    // Plane processing bounds:
    // - 4 halfpels outside the top side (all correct).
    // - 3 halfpels outside the bottom side (all correct).
    // - 16 halfpels outside the left side (3 incorrect).
    // - 16 halfpels outside the right side (4 incorrect).
    //
    // The last pixel we correctly interpolate on a side is the pixel we must
    // replicate on the side, alignment concerns notwithstanding. This pixel is
    // equal to the fullpel on the frame boundary.
    f265_pix *planes[4];
    for (int i = 0; i < 4; i++) planes[i] = f->rec_planes[i] - 4*stride - 16;
    venc_compute_hpel(planes, stride, aligned_width + 2*16, 4 + height + 3, t->store);

    // Plane padding bounds:
    // - top:    duplicate row at -4.
    // - bottom: duplicate row at height + 2.
    // - left:   duplicate column at -8 (for alignment).
    // - right:  duplicate column at width + 7 (for alignment).
    //
    // We do more work than necessary if the frame width is not aligned on 16,
    // but the padding amount is the same on both sides that way. Leave it be.
    for (int i = 1; i < 4; i++)
        venc_pad_plane(f->rec_planes[i], stride, -8, width + 7, -4, height + 2, pad - 8, pad - 8, pad - 4, pad - 3);
}

// Interpolation functions implementation.
// The interpolation functions assumes that the reference pictures are
// padded at the borders (by duplicating the pixels). No checks are performed.
// Luma needs 3 pixels at the left/top and 4 at right/bottom borders.
// Chroma needs 1 pixel at left/top and 2 at right/bottom borders.

// Perform filtering at integer position. No real filtering, just scale up and copy.
void venc_filter_int(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                     int32_t width, int32_t height, int32_t shift)
{
    int x, y;
    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            dst[x] = src[x] << shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// The following functions are deprecated.

// Perform filtering at fractional positions horizontally.
//
// Arguments:
// src, dst                 :           source and destination buffers
// src_stride, dst_stide    :           source and destination strides
// width, height            :           block width and height
// coef                     :           filter coefficients
void venc_filter_hor_luma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                          int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int32_t x, y, tmp;

    src -= 3;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+1]*coef[1] + src[x+2]*coef[2] + src[x+3]*coef[3] +
                  src[x+4]*coef[4] + src[x+5]*coef[5] + src[x+6]*coef[6] + src[x+7]*coef[7];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Perform filtering at fractional positions vertically.
void venc_filter_ver_luma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                          int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int x, y;
    int32_t tmp;

    src -= 3*src_stride;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+src_stride]*coef[1] + src[x+2*src_stride]*coef[2] + src[x+3*src_stride]*coef[3] +
                  src[x+4*src_stride]*coef[4] + src[x+5*src_stride]*coef[5] + src[x+6*src_stride]*coef[6] + src[x+7*src_stride]*coef[7];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Same as above except src type is int16_t.
void venc_filter_ver_luma_s16(int16_t *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                              int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int x, y;
    int32_t tmp;

    src -= 3*src_stride;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+src_stride]*coef[1] + src[x+2*src_stride]*coef[2] + src[x+3*src_stride]*coef[3] +
                  src[x+4*src_stride]*coef[4] + src[x+5*src_stride]*coef[5] + src[x+6*src_stride]*coef[6] + src[x+7*src_stride]*coef[7];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Apply a 4-tap horizontal filtering on chroma component
void venc_filter_hor_chroma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                            int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int x, y;
    int32_t tmp;
    src -= 1;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+1]*coef[1] + src[x+2]*coef[2] + src[x+3]*coef[3];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Apply a 4-tap vertical filtering on chroma component
void venc_filter_ver_chroma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                            int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int x, y;
    int32_t tmp;

    src -= src_stride;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+src_stride]*coef[1] + src[x+2*src_stride]*coef[2] + src[x+3*src_stride]*coef[3];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Same as above except src type is int16_t.
void venc_filter_ver_chroma_s16(int16_t *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                                int32_t width, int32_t height, int8_t *coef, int32_t shift)
{
    int x, y;
    int32_t tmp;

    src -= src_stride;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = src[x]*coef[0] + src[x+src_stride]*coef[1] + src[x+2*src_stride]*coef[2] + src[x+3*src_stride]*coef[3];
            dst[x] = tmp >> shift;
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Perform luma interpolation.
void venc_interpol_luma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                        int32_t width, int32_t height, int32_t xfrac, int32_t yfrac, int32_t bitdepth)
{
    if (xfrac == 0)
    {
        if (yfrac == 0)
        {
            // Integer position. Just copy and scale.
            venc_filter_int(src, src_stride, dst, dst_stride, width, height, 14-bitdepth);
        }
        else
        {
            // Integer position on vertical direction. Perform horizontal filtering.
            venc_filter_ver_luma(src, src_stride, dst, dst_stride, width, height, (int8_t*)f265_if_luma[yfrac],
                                 bitdepth-8);
        }
    }
    else
    {
        if (yfrac == 0)
        {
            // Integer position on horizontal direction. Perform vertical filtering.
            venc_filter_hor_luma(src, src_stride, dst, dst_stride, width, height, (int8_t*)f265_if_luma[xfrac],
                                 bitdepth-8);
        }
        else
        {

            // Fractional position on both directions.
            // We first perform horizontal then vertical filtering.
            // In the horizontal filtering, we generate 3 more lines on top border and
            // 4 on the bottom. Those pixels are needed by the vertical filtering.
            int16_t tmp[(64+7)*64];
            venc_filter_hor_luma(src-3*src_stride, src_stride, tmp, 64, width, height+7,
                                 (int8_t*)f265_if_luma[xfrac], bitdepth-8);
            venc_filter_ver_luma_s16(tmp + 3*64, 64, dst, dst_stride, width, height,
                                     (int8_t*)f265_if_luma[yfrac], 6);
        }
    }
}

// Perform chroma interpolation.
void venc_interpol_chroma(f265_pix *src, int32_t src_stride, int16_t *dst, int32_t dst_stride,
                          int32_t width, int32_t height, int32_t xfrac, int32_t yfrac, int32_t bitdepth)
{
    if (xfrac == 0)
    {
        if (yfrac == 0)
        {
            // Integer position. Just copy and scale.
            venc_filter_int(src, src_stride, dst, dst_stride, width, height, 14-bitdepth);
        }
        else
        {
            // Integer position on vertical direction. Perform horizontal filtering.
            venc_filter_ver_chroma(src, src_stride, dst, dst_stride, width, height, (int8_t*)f265_if_chroma[yfrac],
                                   bitdepth-8);
        }
    }
    else
    {
        if (yfrac == 0)
        {
            // Integer position on horizontal direction. Perform vertical filtering.
            venc_filter_hor_chroma(src, src_stride, dst, dst_stride, width, height, (int8_t*)f265_if_chroma[xfrac],
                                   bitdepth-8);
        }
        else
        {
            // Fractional position on both directions.
            // We first perform horizontal then vertical filtering.
            // In the horizontal filtering, we generate 1 more line on top border and
            // 2 on the bottom. Those pixels are needed by the vertical filtering.
            int16_t tmp[(32+3)*32];
            venc_filter_hor_chroma(src-src_stride, src_stride, tmp, 32, width, height+3,
                                   (int8_t*)f265_if_chroma[xfrac], bitdepth-8);
            venc_filter_ver_chroma_s16(tmp + 32, 32, dst, dst_stride, width, height,
                                       (int8_t*)f265_if_chroma[yfrac], 6);
        }
    }
}

// Scale down pixels and saturate (default weighted prediction).
void venc_scale_pix_uni(int16_t *src, int32_t src_stride, f265_pix *dst, int32_t dst_stride,
                        int32_t width, int32_t height, int32_t bitdepth)
{
    int x, y;
    int32_t tmp;
    int8_t shift = 14 - bitdepth;
    int16_t offset = shift ? (1 << (shift-1)) : 0;
    f265_pix maxval = (1 << bitdepth) - 1;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = (src[x] + offset) >> shift;
            dst[x] = F265_CLAMP(tmp, 0, maxval);
        }

        src += src_stride;
        dst += dst_stride;
    }
}

// Scale down pixels and saturate (default weighted prediction).
void venc_scale_pix_bi(int16_t *src0, int16_t *src1, int32_t src_stride, f265_pix *dst, int32_t dst_stride,
                       int32_t width, int32_t height, int32_t bitdepth)
{
    int x, y;
    int32_t tmp;
    int8_t shift = 15 - bitdepth;
    int16_t offset = 1 << (shift-1);

    f265_pix maxval = (1 << bitdepth) - 1;

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = (src0[x] + src1[x] + offset) >> shift;
            dst[x] = F265_CLAMP(tmp, 0, maxval);
        }

        src0 += src_stride;
        src1 += src_stride;
        dst += dst_stride;
    }
}

// Apply explicit weighted prediction for uni-direction.
// FIXME: argument order violates convention (destination first).
void venc_weight_pix_uni(int16_t *src, int32_t src_stride, f265_pix *dst, int32_t dst_stride,
                         int32_t width, int32_t height, f265_weight *weights, int32_t bitdepth)
{

    int x, y;
    int32_t tmp;
    int16_t round;
    int16_t w = weights->num;
    int16_t off = weights->off;
    int16_t log2_wd = weights->log2_denom + 14 - bitdepth;
    f265_pix maxval = (1<<bitdepth)-1;

    assert(log2_wd >= 0);

    if (log2_wd)
    {
        round = 1<<(log2_wd-1);
        for (y=0; y<height; y++)
        {
            for (x=0; x<width; x++)
            {
                tmp = ((src[x] * w + round) >> log2_wd) + off;
                dst[x] = F265_CLAMP(tmp, 0, maxval);
            }
            src += src_stride;
            dst += dst_stride;
        }
    }
    else
    {
        for (y=0; y<height; y++)
        {
            for (x=0; x<width; x++)
            {
                tmp = src[x] * w + off;
                dst[x] = F265_CLAMP(tmp, 0, maxval);
            }
            src += src_stride;
            dst += dst_stride;
        }
    }
}

// Apply explicit weighted prediction for bi-direction.
void venc_weight_pix_bi(int16_t *src0, int16_t *src1, int32_t src_stride, f265_pix *dst, int32_t dst_stride,
                        int32_t width, int32_t height, f265_weight *weights0, f265_weight *weights1, int32_t bitdepth)
{
    int x, y;
    int32_t tmp;
    int16_t w0 = weights0->num;
    int16_t w1 = weights1->num;
    int16_t log2_wd = weights0->log2_denom + 14 - bitdepth;
    int32_t off = (weights0->off + weights1->off + 1) << log2_wd;
    int16_t shift  = log2_wd + 1;
    f265_pix maxval = (1 << bitdepth) - 1;

    assert(log2_wd >= 0);

    for (y=0; y<height; y++)
    {
        for (x=0; x<width; x++)
        {
            tmp = (src0[x] * w0 + src1[x] * w1 + off) >> shift;
            dst[x] = F265_CLAMP(tmp, 0, maxval);
        }

        src0 += src_stride;
        src1 += src_stride;
        dst += dst_stride;
    }
}

// Apply explicit weighted prediction if enabled, otherwise apply default weighted prediction.
void venc_ref_weighted_pred(int16_t *src0, int16_t *src1, int32_t src_stride, f265_pix *dst, int32_t dst_stride,
                            int32_t bitdepth, f265_frac_ref_block *rb)
{
    f265_weight *w0, *w1;

    // Get the weight pointers for the current image component.
    w0 = rb->ctx[0]->weights + rb->comp;

    if (rb->bipred_flag)
    {
        // Get the weight pointers for the current image component.
        w1 = rb->ctx[1]->weights + rb->comp;

        if (w0->used_flag||w1->used_flag)
        {
            // We do weighted prediction if one or both of the references are
            // weighted. Non-weighted references have default weights for that
            // purpose.
            venc_weight_pix_bi(src0, src1, src_stride, dst, dst_stride, rb->size[0], rb->size[1], w0, w1, bitdepth);
        }
        else
        {
            venc_scale_pix_bi(src0, src1, src_stride, dst, dst_stride, rb->size[0], rb->size[1], bitdepth);
        }
    }
    else
    {
        if (w0->used_flag)
        {
            venc_weight_pix_uni(src0, src_stride, dst, dst_stride, rb->size[0], rb->size[1], w0, bitdepth);
        }
        else
        {
            venc_scale_pix_uni(src0, src_stride, dst, dst_stride, rb->size[0], rb->size[1], bitdepth);
        }
    }
}

// Perform the fractional interpolation of the luma reference block.
void venc_get_frac_ref_luma(f265_enc *enc, f265_pix *dst, int32_t dst_stride, f265_frac_ref_block *rb)
{
    int8_t bitdepth;
    int8_t xfrac, yfrac;
    int32_t ref_stride, xmv, ymv;
    f265_pix *plane, *ref;
    int16_t pred[2][64*64];

    // Image component must be set to luma.
    assert (rb->comp == 0);

    // Get bit-depth.
    bitdepth = enc->gd.bit_depth[0];
    ref_stride = enc->gd.stride;

    for (int reflist=0; reflist<rb->bipred_flag+1; reflist++)
    {
        // Get integer part of the MVs.
        xmv = rb->mv[reflist].x >> 2;
        ymv = rb->mv[reflist].y >> 2;
        // Get fractional part of the MVs.
        xfrac = rb->mv[reflist].x & 3;
        yfrac = rb->mv[reflist].y & 3;

        // Get reference block position using block position and MVs.
        plane = rb->ctx[reflist]->planes[0];
        ref = plane + rb->pos[0] + xmv + (rb->pos[1] + ymv) * ref_stride;
        venc_interpol_luma(ref, ref_stride, pred[reflist], 64, rb->size[0], rb->size[1], xfrac, yfrac, bitdepth);
    }
    venc_ref_weighted_pred(pred[0], pred[1], 64, dst, dst_stride, bitdepth, rb);
}

// Perform the fractional interpolation of the chroma reference block.
// Block position and size must be already scaled to chroma dimensions.
// MV integer and fractional parts computation is made within this function according to chroma idc.
// As of version JCTVC-L1003_v34 of HEVC spec (Jan 2013), 4:2:0 is assumed for
// chroma interpolation (see clause 8.5.3.3.3.1). This function must be checked when other
// chroma formats are fully specified.
void venc_get_frac_ref_chroma(f265_enc *enc, f265_pix *dst, int32_t dst_stride, f265_frac_ref_block *rb)
{
    int8_t bitdepth, plane_idx;
    int8_t xfrac, yfrac;
    int32_t ref_stride, xmv, ymv;
    int8_t *ref_csf;
    f265_pix *plane, *ref;

    // Temporary prediction buffer.
    int16_t pred[2][64*64];

    // Get bit-depth.
    bitdepth = enc->gd.bit_depth[1];
    ref_stride = enc->gd.stride;
    // Get chroma scale factor.
    // For chroma_format_idc=1, i.e. (4:2:0) csf[0]=csf[1]=1, thus, the right shift to
    // get integer part of MVs is 3 and the mask for fractional part of MVs is 7.
    ref_csf = enc->gd.csf;

    // Get the image plane (U or V).
    plane_idx = 3 + rb->comp;

    for (int reflist=0; reflist<rb->bipred_flag+1; reflist++)
    {
        // Get integer part of the MVs.
        xmv = rb->mv[reflist].x >> (ref_csf[0]+2);
        ymv = rb->mv[reflist].y >> (ref_csf[1]+2);
        // Get fractional part of the MVs.
        xfrac = rb->mv[reflist].x & ((1<<(ref_csf[0]+2))-1);
        yfrac = rb->mv[reflist].y & ((1<<(ref_csf[1]+2))-1);
        // Get reference block position using block position and MVs.
        plane = rb->ctx[reflist]->planes[plane_idx];
        ref = plane + rb->pos[0] + xmv + (rb->pos[1] + ymv) * ref_stride;
        venc_interpol_chroma(ref, ref_stride, pred[reflist], 64, rb->size[0], rb->size[1],
                             xfrac, yfrac, bitdepth);
    }

    venc_ref_weighted_pred(pred[0], pred[1], 64, dst, dst_stride, bitdepth, rb);
}

// Perform the fractional interpolation of the reference block specified.
void venc_get_frac_ref(f265_enc *enc, f265_pix *dst, int32_t dst_stride, f265_frac_ref_block *rb)
{
    if (!rb->comp) venc_get_frac_ref_luma(enc, dst, dst_stride, rb);
    else venc_get_frac_ref_chroma(enc, dst, dst_stride, rb);
}

// End of interpolation functions.


///////////////////////////////////////////////////////////////////////////////
// Luma/chroma quarterpel interpolation.

// Each interpolation function has two versions. The "pix" version clips the
// final interpolation result to the pixel size. The "s16" version stores the
// intermediate interpolation result as 16-bit.

// Luma 8-taps filter.
#define FILTER_8TAPS(src, stride, factors)\
    (*((src)-3*(stride))*(factors)[0] +\
     *((src)-2*(stride))*(factors)[1] +\
     *((src)-1*(stride))*(factors)[2] +\
     *((src)+0*(stride))*(factors)[3] +\
     *((src)+1*(stride))*(factors)[4] +\
     *((src)+2*(stride))*(factors)[5] +\
     *((src)+3*(stride))*(factors)[6] +\
     *((src)+4*(stride))*(factors)[7])

// Chroma 4-taps filter.
#define FILTER_4TAPS(src, stride, factors)\
    (*((src)-1*(stride))*(factors)[0] +\
     *((src)+0*(stride))*(factors)[1] +\
     *((src)+1*(stride))*(factors)[2] +\
     *((src)+2*(stride))*(factors)[3])

// The following functions were auto-generated by snippets/interpol.py.
void venc_interpol_luma_qpel_pix_h_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_luma[frac&3];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_8TAPS(src+x, 1, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_luma_qpel_s16_h_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = bd - 8;
    const int8_t *factors = f265_if_luma[frac&3];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_8TAPS(src+x, 1, factors)>>shift;
}

void venc_interpol_luma_qpel_pix_v_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_luma[frac>>2];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_8TAPS(src+x, src_stride, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_luma_qpel_s16_v_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = bd - 8;
    const int8_t *factors = f265_if_luma[frac>>2];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_8TAPS(src+x, src_stride, factors)>>shift;
}

void venc_interpol_luma_qpel_pix_d_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 20 - bd;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_luma[frac>>2];

    // Interpolate 3 rows above the block and 4 rows below.
    int16_t *h_buf = (int16_t*)spill;
    venc_interpol_luma_qpel_s16_h_c(h_buf, width, src - 3*src_stride, src_stride, frac,
                                    ((bd-8)<<16)|(width<<8)|(height+7), NULL);
    h_buf += 3*width;

    for (int y = 0; y < height; y++, dst += dst_stride, h_buf += width)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_8TAPS(h_buf+x, width, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_luma_qpel_s16_d_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                     int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    const int8_t *factors = f265_if_luma[frac>>2];

    // Interpolate 3 rows above the block and 4 rows below.
    int16_t *h_buf = (int16_t*)spill;
    venc_interpol_luma_qpel_s16_h_c(h_buf, width, src - 3*src_stride, src_stride, frac,
                                    ((bd-8)<<16)|(width<<8)|(height+7), NULL);
    h_buf += 3*width;

    for (int y = 0; y < height; y++, dst += dst_stride, h_buf += width)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_8TAPS(h_buf+x, width, factors)>>shift;
}

void venc_interpol_chroma_qpel_pix_h_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_chroma[frac&7];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_4TAPS(src+x, 1, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_chroma_qpel_s16_h_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = bd - 8;
    const int8_t *factors = f265_if_chroma[frac&7];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_4TAPS(src+x, 1, factors)>>shift;
}

void venc_interpol_chroma_qpel_pix_v_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_chroma[frac>>3];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_4TAPS(src+x, src_stride, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_chroma_qpel_s16_v_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = bd - 8;
    const int8_t *factors = f265_if_chroma[frac>>3];

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_4TAPS(src+x, src_stride, factors)>>shift;
}

void venc_interpol_chroma_qpel_pix_d_c(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 20 - bd;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;
    const int8_t *factors = f265_if_chroma[frac>>3];

    // Interpolate 1 rows above the block and 2 rows below.
    int16_t *h_buf = (int16_t*)spill;
    venc_interpol_chroma_qpel_s16_h_c(h_buf, width, src - 1*src_stride, src_stride, frac,
                                      ((bd-8)<<16)|(width<<8)|(height+3), NULL);
    h_buf += 1*width;

    for (int y = 0; y < height; y++, dst += dst_stride, h_buf += width)
        for (int x = 0; x < width; x++)
        {
            int pix = (FILTER_4TAPS(h_buf+x, width, factors) + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}

void venc_interpol_chroma_qpel_s16_d_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
                                       int packed_dims, uint8_t *spill)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 6;
    const int8_t *factors = f265_if_chroma[frac>>3];

    // Interpolate 1 rows above the block and 2 rows below.
    int16_t *h_buf = (int16_t*)spill;
    venc_interpol_chroma_qpel_s16_h_c(h_buf, width, src - 1*src_stride, src_stride, frac,
                                      ((bd-8)<<16)|(width<<8)|(height+3), NULL);
    h_buf += 1*width;

    for (int y = 0; y < height; y++, dst += dst_stride, h_buf += width)
        for (int x = 0; x < width; x++)
            dst[x] = FILTER_4TAPS(h_buf+x, width, factors)>>shift;
}
// End of auto-generated code.

#undef FILTER_8TAPS
#undef FILTER_4TAPS

// Scale the pixels without interpolating.
void venc_scale_qpel_c(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int packed_dims)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 14 - bd;

    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)
        for (int x = 0; x < width; x++)
            dst[x] = src[x]<<shift;
}

// Compute the average of two pixel blocks. The destination is unstriden since
// this is the most common use case. Only valid AMP luma block sizes are
// supported in the assembly dispatch version.
void venc_avg_pix_c(f265_pix *dst, f265_pix *src0, int src0_stride, f265_pix *src1, int src1_stride, int packed_dims)
{
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    for (int y = 0; y < height; y++, dst += width, src0 += src0_stride, src1 += src1_stride)
        for (int x = 0; x < width; x++)
            dst[x] = (src0[x] + src1[x] + 1) >> 1;
}

// Unweighted bi-prediction.
void venc_avg_pix_s16_c(f265_pix *dst, int dst_stride, int16_t *src0, int16_t *src1, int src_stride, int packed_dims)
{
    #ifdef F265_LBD
    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #else
    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    #endif
    int shift = 15 - bd;
    int bias = 1<<(shift-1), max_val = (1<<bd)-1;

    for (int y = 0; y < height; y++, dst += dst_stride, src0 += src_stride, src1 += src_stride)
        for (int x = 0; x < width; x++)
        {
            int pix = (src0[x] + src1[x] + bias)>>shift;
            dst[x] = F265_CLAMP(pix, 0, max_val);
        }
}


///////////////////////////////////////////////////////////////////////////////
// C implementation of the distortion functions.

// Compute the sum of absolute differences (SAD) of two blocks.
// The subshift argument controls the SAD sub-sampling.
// 0 : no sub-sampling, 1 : one line out of two is used.
int32_t venc_sad(f265_pix *src0, int32_t stride0, f265_pix *src1, int32_t stride1,
                 int32_t width, int32_t height, int32_t subshift, int32_t bitdepth)
{
    int32_t sad = 0;
    int32_t step = 1 << subshift;

    stride0 <<= subshift;
    stride1 <<= subshift;

    for (int32_t y = 0; y < height; y += step)
    {
        for (int32_t x = 0; x < width; x++)
        {
            sad += F265_ABS(src0[x] - src1[x]);
        }
        src0 += stride0;
        src1 += stride1;
    }

    // Scale up to match non-sampled SAD.
    sad <<= subshift;

    // Normalize distortion to 8-bit depth if needed.
    sad = DIST_BD_NORM(sad, bitdepth-8);

    return sad;
}

// Fast SAD version (rename once the symbol above no longer conflicts).
int venc_fsad_c(f265_pix *src0, int stride0, f265_pix *src1, int stride1, int packed_dims)
{
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    int sad = 0;
    for (int y = 0; y < height; y++, src0 += stride0, src1 += stride1)
        for (int x = 0; x < width; x++)
            sad += F265_ABS(src0[x] - src1[x]);
    return sad;
}

// 3-or-4 block SAD distortion function. Requires less loads and loops than the
// single block function. The costs and references are stored in arrays. Costs
// in the array are replaced. All the references have the same stride, but the
// source has an independent stride. The width and height are packed together
// ((width<<8)|height) to avoid using the stack on amd64.
//
// Note: costs[3] is written by the assembly implementation for SAD3. Be
// careful.
void venc_sad3_c(int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims)
{
    for (int i = 0; i < 3; i++)
        costs[i] = venc_fsad_c(src, src_stride, refs[i], ref_stride, packed_dims);
}

void venc_sad4_c(int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims)
{
    for (int i = 0; i < 4; i++)
        costs[i] = venc_fsad_c(src, src_stride, refs[i], ref_stride, packed_dims);
}

// Same as the SAD function, but compute the sum of square differences (SSD).
int32_t venc_ssd(f265_pix *src0, int32_t stride0, f265_pix *src1, int32_t stride1,
                 int32_t width, int32_t height, int32_t subshift, int32_t bitdepth)
{
    int32_t ssd = 0;
    int32_t step = 1 << subshift;

    stride0 <<= subshift;
    stride1 <<= subshift;

    for (int32_t y = 0; y < height; y += step)
    {
        for (int32_t x = 0; x < width; x++)
        {
            int tmp = src0[x] - src1[x];
            ssd += tmp*tmp;
        }
        src0 += stride0;
        src1 += stride1;
    }

    // Scale up to match non-sampled SAD.
    ssd <<= subshift;

    // Normalize distortion to 8-bit depth if needed.
    ssd = DIST_BD_NORM(ssd, bitdepth-8);

    return ssd;
}

// Same as the SAD function, but compute the sum of square differences (SSD).
int32_t venc_ssd16(int16_t *src0, int32_t stride0, int16_t *src1, int32_t stride1,
                   int32_t width, int32_t height, int32_t subshift, int32_t bitdepth)
{
    int32_t ssd = 0;
    int32_t step = 1 << subshift;

    stride0 <<= subshift;
    stride1 <<= subshift;

    for (int32_t y = 0; y < height; y += step)
    {
        for (int32_t x = 0; x < width; x++)
        {
            int tmp = src0[x] - src1[x];
            ssd += tmp*tmp;
        }
        src0 += stride0;
        src1 += stride1;
    }

    // Scale up to match non-sampled SAD.
    ssd <<= subshift;

    // Normalize distortion to 8-bit depth if needed.
    ssd = DIST_BD_NORM(ssd, bitdepth-8);

    return ssd;
}

// Compute the sum of absolute differences after a Hadamard transform.
// The transform is a straightforward application of matrix multiplication :
// T = H * D * Transpose(H)
int32_t venc_satd_NxN(f265_pix *src0, int32_t stride0, f265_pix *src1, int32_t stride1, int32_t dim, int32_t bitdepth)
{
    int32_t d[F265_H_DIM][F265_H_DIM], m[F265_H_DIM][F265_H_DIM];
    int32_t shift = dim >> 2;
    int32_t sum=0, tmp;
    int x, y, i;

    // Compute the difference.
    for (y=0; y<dim; y++)
    {
        for (x=0; x<dim; x++)
        {
            d[y][x] = src0[x] - src1[x];
        }
        src0 += stride0;
        src1 += stride1;
    }

    // First matrix multiplication.
    for (y=0; y<dim; y++)
    {
        for (x=0; x<dim; x++)
        {
            for (i=0, tmp=0; i<dim; i++)
            {
                tmp += f265_h8x8[y][i] * d[i][x];
            }
            m[y][x] = tmp;
        }
    }

    // Second matrix multiplication and accumulation.
    for (y=0; y<dim; y++)
    {
        for (x=0; x<dim; x++)
        {
            for (i=0, tmp=0; i < dim; i++)
            {
                tmp += m[y][i] * f265_h8x8[i][x];
            }
            sum += F265_ABS(tmp);
        }
    }

     return (sum + shift) >> shift;
}

// Compute the satd of a block. Only transform sizes 2x2, 4x4, 8x8 are needed considering HM implementation.
// For the blocks that have sizes not equal to 2x2, 4x4 or 8x8, we determine the largest possible transform size
// to use then we loop over the sub-blocks of that size.
int32_t venc_satd(f265_pix *src0, int32_t stride0, f265_pix *src1, int32_t stride1,
                  int32_t width, int32_t height, int32_t bitdepth)
{
    // Block sizes must be multiple of 2 at least.
    assert (( width % 2 == 0) && (height % 2 == 0));

    // Get transform size.
    int8_t tr_size = 2;
    if (( width % 4 == 0) && (height % 4 == 0)) tr_size = 4;
    if (( width % 8 == 0) && (height % 8 == 0)) tr_size = 8;

    // Get source offsets.
    int32_t off_src0 = stride0 * tr_size;
    int32_t off_src1 = stride1 * tr_size;

    uint32_t sum = 0;

    for (int32_t y=0; y<height; y+=tr_size)
    {
        for (int32_t x=0; x<width; x+=tr_size)
        {
            sum += venc_satd_NxN(src0+x, stride0, src1+x, stride1, tr_size, bitdepth);
        }
        src0 += off_src0;
        src1 += off_src1;
    }

    // Normalize distortion to 8 bit depth if needed.
    sum = DIST_BD_NORM(sum, bitdepth-8);

    return sum;
}

// End of sad/satd functions.


