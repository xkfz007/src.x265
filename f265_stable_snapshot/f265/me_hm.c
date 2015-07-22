// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// HM-compatible motion estimation functions.

#include "f265/enc.h"

#ifdef VAN_USE_HM_ME

// Estimate the number of bits needed to encode a value.
int32_t venc_hm_bits_count(int32_t val)
{
    int32_t length = 1;
    int32_t tmp = (val <= 0) ? (-val << 1) + 1: (val << 1);

    while (tmp != 1)
    {
        tmp >>= 1;
        length += 2;
    }

    return length;
}

// Return cost.
int32_t venc_hm_bits_cost(int32_t bits, int32_t lambda)
{
    return bits * lambda >> 16;
}

// Estimate MV cost.
// Arguments:
// x,y:      motion vector components.
// mvp:      motion vector predictor.
// lambda:   mv lambda.
// mv_scale: mv scale factor. 0: quarter pel, 1: half pel, 2: full pel.
int32_t venc_hm_mv_cost(int32_t x, int32_t y, f265_mv *mvp, int32_t lambda, int32_t mv_scale)
{
    int32_t bits;
    bits = venc_hm_bits_count((x << mv_scale) - mvp->x) + venc_hm_bits_count((y << mv_scale) - mvp->y);
    return venc_hm_bits_cost(bits, lambda);
}

// Clip the motion vector.
// Arguments:
// mv:   MV to clip (input/output).
// off:  extra offset for interpolation filter, possible values (4 or 7).
//       This is needed for bit-exactness with HM, because HM uses
//       16 pixels when padding, versus 8 in f265.
void venc_hm_mv_clip(f265_enc_thread *t, f265_cb *cb, f265_mv *mv, int32_t off)
{
    // Get encoder pointer.
    f265_enc *enc = t->enc;

    int32_t cb_ox = cb->cb_off[0] + t->ctb_off[0];
    int32_t cb_oy = cb->cb_off[1] + t->ctb_off[1];

    // Set picture data (follow HM behavior).
    int32_t ctb_size = t->enc->gd.ctb_size;
    int32_t pic_width = enc->gd.clip_dim[0];
    int32_t pic_height = enc->gd.clip_dim[1];

    // Compute MV min and max.
    int32_t hor_min = (-ctb_size - off - cb_ox) << 2;
    int32_t hor_max = (pic_width + off - cb_ox) << 2;
    int32_t ver_min = (-ctb_size - off - cb_oy) << 2;
    int32_t ver_max = (pic_height + off - cb_oy) << 2;

    // Clip and store MV.
    mv->x = F265_MIN(hor_max, F265_MAX(hor_min, mv->x));
    mv->y = F265_MIN(ver_max, F265_MAX(ver_min, mv->y));
}

// Set motion estimation search range.
void venc_hm_me_range(f265_enc_thread *t, f265_cb *cb, f265_mv *mv, int32_t range)
{
    f265_me_ctx *mex = &t->me;
    f265_mv v;

    // Scale range to quarter pixel resolution.
    range <<= 2;

    // Set left-top position.
    v.x = mv->x - range;
    v.y = mv->y - range;
    venc_hm_mv_clip(t, cb, &v, 7);
    mex->me_bounds[0] = v.x >> 2;
    mex->me_bounds[1] = v.y >> 2;

    // Set right-bottom position.
    v.x = mv->x + range;
    v.y = mv->y + range;
    venc_hm_mv_clip(t, cb, &v, 7);
    mex->me_bounds[2] = v.x >> 2;
    mex->me_bounds[3] = v.y >> 2;
}

// Perform diamond pattern search.
void venc_hm_tz_diamond(f265_enc_thread *t, f265_frac_ref_block *rb, int32_t start_x, int32_t start_y, int32_t dist)
{
    f265_me_ctx *mex = &t->me;
    mex->hm_best_round++;

    // Test points 2, 4, 5 and 7 for distance equals to 1.
    if (dist == 1)
    {
        if (start_y - 1 >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, start_x, start_y - 1, 2, 1);
        if (start_x - 1 >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - 1, start_y, 4, 1);
        if (start_x + 1 <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + 1, start_y, 5, 1);
        if (start_y + 1 <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, start_x, start_y + 1, 7, 1);
        return;
    }

    // Test points 1 to 8 for distances in the range [2, 8].
    if (dist < 9)
    {
        int8_t d = dist >> 1;
        if (start_y - dist >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, start_x, start_y - dist, 2, dist);
        if (start_y - d >= mex->me_bounds[1])
        {
            if (start_x - d >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - d, start_y - d, 1, d);
            if (start_x + d <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + d, start_y - d, 3, d);
        }
        if (start_x - dist >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - dist, start_y, 4, dist);
        if (start_x + dist <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + dist, start_y, 5, dist);
        if (start_y + d <= mex->me_bounds[3])
        {
            if (start_x - d >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - d, start_y + d, 6, d);
            if (start_x + d <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + d, start_y + d, 8, d);
        }
        if (start_y + dist <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, start_x, start_y + dist, 7, dist);

        return;
    }

    // Test points 1 to 8 plus 8 extra points for distances greater than 8.

    // Test top, left, right and bottom points.
    if (start_y - dist >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, start_x, start_y - dist, 0, dist);
    if (start_x - dist >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - dist, start_y, 0, dist);
    if (start_x + dist <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + dist, start_y, 0, dist);
    if (start_y + dist <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, start_x, start_y + dist, 0, dist);

    // Test points that are on the diamond sides.
    int8_t d = dist >> 2, dd;
    for (int i = d; i < dist; i += d)
    {
        dd = dist - i;
        if (start_y - dd >= mex->me_bounds[1])
        {
            if (start_x - i >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - i, start_y - dd, 0, dist);
            if (start_x + i <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + i, start_y - dd, 0, dist);
        }
        if (start_y + dd <= mex->me_bounds[3])
        {
            if (start_x - i >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, start_x - i, start_y + dd, 0, dist);
            if (start_x + i <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, start_x + i, start_y + dd, 0, dist);
        }
    }
}

// Test two additional outer points.
void venc_hm_tz_2points(f265_enc_thread *t, f265_frac_ref_block *rb)
{
    // Two points search table.
    int8_t tp_tab[9][4] =
    {
         {  0,  0,  0,  0},
         { -1,  0,  0, -1},
         { -1, -1,  1, -1},
         {  0, -1,  1,  0},
         { -1,  1, -1, -1},
         {  1, -1,  1,  1},
         { -1,  0,  0,  1},
         { -1,  1,  1,  1},
         {  1,  0,  0,  1},
    };

    f265_me_ctx *mex = &t->me;
    int16_t start_x = mex->ref[0].hm_best_mv.x;
    int16_t start_y = mex->ref[0].hm_best_mv.y;
    int8_t dir = mex->hm_best_dir;
    int16_t x0 = start_x + tp_tab[dir][0];
    int16_t y0 = start_y + tp_tab[dir][1];
    int16_t x1 = start_x + tp_tab[dir][2];
    int16_t y1 = start_y + tp_tab[dir][3];

    switch (dir)
    {
    case 1:
        if (x0 >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
        if (y1 >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        break;
    case 2:
        if (y0 >= mex->me_bounds[1])
        {
            if (x0 >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
            if (x1 <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        }
        break;
    case 3:
        if (y0 >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
        if (x1 <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        break;
    case 4:
        if (x0 >= mex->me_bounds[0])
        {
            if (y0 <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
            if (y1 >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        }
        break;
    case 5:
        if (x0 <= mex->me_bounds[2])
        {
            if (y0 >= mex->me_bounds[1]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
            if (y1 <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        }
        break;
    case 6:
        if (x0 >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
        if (y1 <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        break;
    case 7:
        if (y0 <= mex->me_bounds[3])
        {
            if (x0 >= mex->me_bounds[0]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
            if (x1 <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        }
        break;
    case 8:
        if (x0 <= mex->me_bounds[2]) venc_hm_tz_eval_mv(t, rb, x0, y0, 0, 2);
        if (y1 <= mex->me_bounds[3]) venc_hm_tz_eval_mv(t, rb, x1, y1, 0, 2);
        break;
    }
}

// Copy 8u pixels to 16s pixels (needed by the bi-directional prediction search).
void venc_hm_copy_8u_to_16s(f265_pix *dst, int32_t dst_stride, uint8_t *src, int32_t src_stride,
                            int32_t width, int32_t height)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            dst[x] = (int16_t)src[x];
        }
        src += src_stride;
        dst += dst_stride;
    }
}

// Test both MVPs and select the best.
int32_t venc_hm_test_mvp(f265_enc_thread *t, f265_cb *cb, f265_mv mvp_list[], f265_frac_ref_block *rb)
{
    // Set the first MVP as the best.
    int32_t bd = t->enc->gd.bit_depth[0];
    f265_me_ctx *mex = &t->me;
    mex->ref[0].pmv = mvp_list[0];

    // Quit if the MVPs are equal.
    if (mvp_list[0].p == mvp_list[1].p) return 0;

    // Interpolated pixel buffer.
    f265_pix us_ref[64 * 64];

    // Set source plane.
    f265_pix *src = mex->src_planes[0] + rb->pos[1] * mex->hm_src_stride + rb->pos[0];

    // Get the up-sampled reference for the first MVP.
    venc_hm_upsample(t, cb, rb, &mvp_list[0], us_ref, 64);

    // Get the cost.
    int32_t cost = venc_sad(src, mex->hm_src_stride, us_ref, 64, rb->size[0], rb->size[1], 0, bd);

    // Get the up-sampled reference for the second MVP.
    venc_hm_upsample(t, cb, rb, &mvp_list[1], us_ref, 64);

    // Select the second MVP if better.
    if (venc_sad(src, mex->hm_src_stride, us_ref, 64, rb->size[0], rb->size[1], 0, bd) < cost)
    {
        mex->ref[0].pmv = mvp_list[1];
        return 1;
    }

    return 0;
}

// Perform full search motion estimation at full pixel resolution.
void venc_hm_me_full_search(f265_enc_thread *t, f265_frac_ref_block *rb)
{
    f265_me_ctx *mex = &t->me;

    // Set search range.
    int16_t x_min = mex->me_bounds[0];
    int16_t y_min = mex->me_bounds[1];
    int16_t x_max = mex->me_bounds[2];
    int16_t y_max = mex->me_bounds[3];
    int32_t stride_ref = mex->ref_stride;
    int32_t src_stride = mex->hm_src_stride;

    // Get the reference pixels.
    f265_pix *ref = mex->ref[0].ref_planes[0] + rb->pos[1] * mex->ref_stride + rb->pos[0];

    // Get the source pixels.
    f265_pix *src = mex->src_planes[0] + rb->pos[1] * mex->hm_src_stride + rb->pos[0];

    int32_t width = rb->size[0];
    int32_t height = rb->size[1];
    int32_t shift = mex->hm_subshift;
    int32_t bd = t->enc->gd.bit_depth[0];
    f265_mv *mvp = &mex->ref[0].pmv;
    int32_t lambda = mex->hm_lambda_me;
    int16_t x_best = 0, y_best = 0;
    int32_t cost, cost_best = F265_MAX_SAD;
    int32_t x, y;

    for (y = y_min, ref += y_min * stride_ref; y <= y_max; y++, ref += stride_ref)
    {
        for (x = x_min; x <= x_max; x++)
        {
            cost = venc_sad(src, src_stride, ref + x, stride_ref, width, height, shift, bd);
            #ifdef VAN_TRACE_FS_MV_COST
            printf("mv(%d,%d) SAD=%d ", x, y, cost);
            #endif
            cost += venc_hm_mv_cost(x, y, mvp, lambda, 2);
            #ifdef VAN_TRACE_FS_MV_COST
            printf("cost=%d\n", cost);
            #endif
            if (cost < cost_best)
            {
                cost_best = cost;
                x_best = x;
                y_best = y;
            }
        }
    }
    #ifdef VAN_TRACE_FS_MV_COST
    printf("best mv = (%d,%d)\n", x_best, y_best);
    #endif

    // Save best MV.
    mex->ref[0].hm_best_mv.x = x_best;
    mex->ref[0].hm_best_mv.y = y_best;
}

// Perform fractional pixel search refinement.
// Arguments:
// us_ref: array of up-sampled reference buffers.
// frac: fractional part of the MV.
// frac_step: 2 for half and 1 for quarter pixel.
int32_t venc_hm_frac_refine(f265_enc_thread *t, f265_frac_ref_block *rb, f265_pix us_ref[][64 * 64],
                            int32_t stride_ref, f265_mv *frac, int32_t frac_step)
{
    // A word on HM implementation of fractional pixel motion estimation:

    // Surrounding half pixel positions relative to central pixel (0,0) (at quarter resolution).
    //    (-2,-2)       ( 0,-2)        ( 2,-2)
    //
    //    (-2, 0)       ( 0, 0)        ( 2, 0)
    //
    //    (-2, 2)       ( 0, 2)        ( 2, 2)

    // Half pixel interpolated planes.
    //    d    v    d
    //
    //    h    i    h
    //
    //    d    v    d

    // In the fractional pixel resolution refinement, 9 pixels (the 8 surrounding pixels plus the
    // central pixel) are tested.
    // When testing for half pixel, 4 planes are needed. Some pixels share the same plane:
    // e.g. (-2,0) and (2,0) are on the same plane (plane h). Thus, the planes may need to be
    // extended by one line or one column depending on the case: e.g. plane h needs to be extended
    // by one column on the left.
    // The interpolated buffer pointers point at the origin of the plane : e.g. for plane h, the
    // buffer point to the pixel at the relative position (-2,0).
    // The plane index is given by the absolute value of the MV frac :
    // index_x = abs (frac_x) and index_y = abs(frac_y)
    // The pointer offset is computed according to the position of the pixel relative to the interpolated
    // buffer origin. For example, pixel (-2,0) is on the plane (2,0) (or plane h) and the offset is 0.
    // Pixel (2,0) is on the plane (2,0) (or plane h) and the offset is 1. Pixel (0, 2) is on the
    // plane (0,2) (or plane v) and the offset is the buffer stride.
    // For quarter pixel search, all the surrounding pixels lay on different planes, thus 8 planes need to
    // be generated. The ninth plane where the central pixel lays is already available. Note that HM re-test
    // the central pixel (best half pixel position) but it needs not to.

    // Refine pattern at half pel resolution.
    int16_t const f265_refine_half[9*2] =
    {
         0,  0,
         0, -1,
         0,  1,
        -1,  0,
         1,  0,
        -1, -1,
         1, -1,
        -1,  1,
         1,  1
    };

    // Refine pattern at quarter pel resolution.
    int16_t const f265_refine_quart[9*2] =
    {
         0,  0,
         0, -1,
         0,  1,
        -1, -1,
         1, -1,
        -1,  0,
         1,  0,
        -1,  1,
         1,  1
    };

    int16_t frac_x, frac_y;
    int32_t cost, best_cost = F265_MAX_SAD;
    int8_t idx_best = 0;
    f265_me_ctx *mex =  &t->me;
    int32_t bd = t->enc->gd.bit_depth[0];

    // Get refinement table.
    // HM uses two different pattern tables for half and quarter.
    // The pattern is the same (all surrounding pixel are tested) but the order differs.
    // We need to use these tables for bit exactness.
    int16_t const *ref_tab = frac_step == 2 ? f265_refine_half : f265_refine_quart;

    // Loop over the surrounding pixels.
    for (int16_t i = 0; i < 9; i++)
    {
        // Get the fractional MV.
        // When testing half pixel, me->frac MV will contain (0,0).
        // When testing quarter pixel, me->frac MV will contain best half pel MV.
        frac_x = frac->x + ref_tab[i << 1] * frac_step;
        frac_y = frac->y + ref_tab[(i << 1) + 1] * frac_step;

        int16_t buf_idx = (ref_tab[(i << 1) + 1] + 1) * 3 + ref_tab[i << 1] + 1;
        f265_pix *ref = us_ref[buf_idx];

        // Get the source block.
        f265_pix *src = mex->src_planes[0] + rb->pos[1] * mex->hm_src_stride + rb->pos[0];

        // Compute the distortion using satd/sad.
        if (mex->dist_func_id)
            cost = venc_satd(src, mex->hm_src_stride, ref, stride_ref, rb->size[0], rb->size[1], bd);
        else
            cost = venc_sad(src, mex->hm_src_stride, ref, stride_ref, rb->size[0], rb->size[1], 0, bd);

        // Add MV cost.
        cost += venc_hm_mv_cost(mex->ref[0].hm_best_mv.x + frac_x, mex->ref[0].hm_best_mv.y + frac_y,
                                &mex->ref[0].pmv, mex->hm_lambda_me, 0);

        // Store the best MV.
        if (cost < best_cost)
        {
            best_cost = cost;
            idx_best = i;
        }
    }

    // Save best MV.
    frac->x += ref_tab[(idx_best << 1)] * frac_step;
    frac->y += ref_tab[(idx_best << 1) + 1] * frac_step;

    // Save MV best cost.
    mex->best_cost = best_cost;

    // Return cost.
    return best_cost;
}

// Interpolate a block of pixels.
void venc_hm_upsample(f265_enc_thread *t, f265_cb *cb, f265_frac_ref_block *rb, f265_mv *mv,
                      f265_pix *us_ref, int32_t stride)
{
    f265_me_ctx *mex = &t->me;

    // Set to 0 even for bi-predicted partitions.
    rb->bipred_flag = 0;

    // Set color component to luma (default).
    rb->comp = 0;

    // Set the reference context parameters.
    rb->ctx[0]->planes[0] = mex->ref[0].ref_planes[0];

    // FIXME. We need to set the WP.
    rb->ctx[0]->weights = mex->ref[0].weights;
    rb->ctx[1]->weights = mex->ref[0].weights;

    // Set the MV and clip.
    rb->mv[0].p = mv->p;
    venc_hm_mv_clip(t, cb, &rb->mv[0], 4);

    // Get up-sampled reference.
    venc_get_frac_ref_luma(t->enc, us_ref, stride, rb);
}

// Perform the up-sampling of the reference block for fractional search for the 9 possible positions.
void venc_hm_upsample_all(f265_enc_thread *t, f265_cb *cb, f265_frac_ref_block *rb,
                          f265_pix us_ref[][64 * 64], f265_mv *frac, int32_t frac_step)
{
    // Optimize this, interpolation is redundant for some MVs. See documentation.
    for (int y = -1; y < 2; y++)
        for (int x = -1; x < 2; x++)
        {
            rb->mv[0].x = t->me.ref[0].hm_best_mv.x + frac->x + x * frac_step;
            rb->mv[0].y = t->me.ref[0].hm_best_mv.y + frac->y + y * frac_step;
            venc_hm_mv_clip(t, cb, &rb->mv[0], 4);
            venc_get_frac_ref_luma(t->enc, us_ref[(y+1)*3+(x+1)], 64, rb);
        }
}

// Perform fractional search.
// Arguments:
// frac: the output fractional part of the MV.
int32_t venc_hm_frac_search(f265_enc_thread *t, f265_cb *cb, f265_frac_ref_block *rb, f265_mv *frac)
{
    int32_t cost;
    f265_me_ctx *mex = &t->me;

    // Up-sampled planes.
    f265_pix us_ref[9][64 * 64];

    // Set to 0 even for bi-predicted partitions.
    rb->bipred_flag = 0;

    // Set color component to luma (default).
    rb->comp = 0;

    // Set up the reference context.
    rb->ctx[0]->planes[0] = mex->ref[0].ref_planes[0];
    rb->ctx[0]->weights = mex->ref[0].weights;
    rb->ctx[1]->weights = mex->ref[0].weights;

    // Reset fractional MV.
    frac->p = 0;

    // Search the fractional MV at half pixel resolution.

    // Generate up-sampled planes.
    venc_hm_upsample_all(t, cb, rb, us_ref, frac, 2);

    // Half pixel search.
    venc_hm_frac_refine(t, rb, us_ref, 64, frac, 2);

    // Search the fractional MV at quarter pixel resolution.

    // Generate up-sampled planes.
    venc_hm_upsample_all(t, cb, rb, us_ref, frac, 1);

    // Quarter pixel search.
    cost = venc_hm_frac_refine(t, rb, us_ref, 64, frac, 1);

    return cost;
}

// Perform TZ search algorithm.
void venc_hm_tz_search(f265_enc_thread *t, f265_frac_ref_block *rb, f265_cb *cb)
{
    f265_me_ctx *mex = &t->me;
    f265_mv *best_mv = &mex->ref[0].hm_best_mv;

    // Set search range. TZ implies uni-directional ME.
    int16_t range = mex->hm_range[0];

    // Clip best MV.
    venc_hm_mv_clip(t, cb, best_mv, 7);

    // Test MVP.
    venc_hm_tz_eval_mv(t, rb, best_mv->x, best_mv->y, 0, 0);

    // Test zero MV.
    venc_hm_tz_eval_mv(t, rb, 0, 0, 0, 0);

    // Set start point.
    int16_t start_x = best_mv->x;
    int16_t start_y = best_mv->y;

    int dist;

    // Perform diamond search.
    for (dist = 1; dist <= range; dist *= 2)
    {
        venc_hm_tz_diamond(t, rb, start_x, start_y, dist);

        // Break if max round is reached.
        if (mex->hm_best_round >= 3) break;
    }

    // Test two additional points.
    if (mex->hm_best_dist == 1)
    {
        mex->hm_best_dist = 0;
        venc_hm_tz_2points(t, rb);
    }

    // Test in raster mode if the best distance is large.
    int const raster = 5;
    if (mex->hm_best_dist > raster)
    {
        mex->hm_best_dist = raster;
        for (int y = mex->me_bounds[1]; y <= mex->me_bounds[3]; y += raster)
            for (int x = mex->me_bounds[0]; x <= mex->me_bounds[2]; x += raster)
            {
                venc_hm_tz_eval_mv(t, rb, x, y, 0, dist);
            }
    }

    // Refine.
    while (mex->hm_best_dist > 0)
    {
        start_x = best_mv->x;
        start_y = best_mv->y;
        mex->hm_best_dist = 0;
        mex->hm_best_dir =  0;

        for (dist = 1; dist < range + 1; dist *= 2)
        {
            venc_hm_tz_diamond(t, rb, start_x, start_y, dist);
        }

        // Two additional points.
        if (mex->hm_best_dist == 1)
        {
            mex->hm_best_dist = 0;
            if (mex->hm_best_dir != 0)
            {
                venc_hm_tz_2points(t, rb);
            }
        }
    }
}

// Evaluate one MV and update the best MV if a better MV is found.
void venc_hm_tz_eval_mv(f265_enc_thread *t, f265_frac_ref_block *rb, int32_t x, int32_t y, int32_t dir, int32_t dist)
{
    f265_me_ctx *mex = &t->me;

    int32_t stride_ref = mex->ref_stride;
    int32_t src_stride = mex->hm_src_stride;
    f265_pix *ref = mex->ref[0].ref_planes[0] + (rb->pos[1] + y) * mex->ref_stride + rb->pos[0] + x;

    // Get source pixels.
    f265_pix *src = mex->src_planes[0] + rb->pos[1] * mex->hm_src_stride + rb->pos[0];

    int32_t width = rb->size[0];
    int32_t height = rb->size[1];
    int32_t shift = mex->hm_subshift;
    int32_t bd = t->enc->gd.bit_depth[0];
    int32_t lambda = mex->hm_lambda_me;

    // Compute MV cost.
    int32_t mv_cost = venc_hm_mv_cost(x, y, &mex->ref[0].pmv, lambda, 2);

    // Compute total cost.
    int32_t cost = venc_sad(src, src_stride, ref, stride_ref, width, height, shift, bd) + mv_cost;

    // Update if a new minimum is found.
    if (cost < mex->best_cost)
    {
        mex->ref[0].hm_best_mv.x = x;
        mex->ref[0].hm_best_mv.y = y;
        mex->best_cost = cost;
        mex->hm_best_dir = dir;
        mex->hm_best_dist = dist;
        mex->hm_best_round = 0;
    }
}

// Pre-process source frame in bi-prediction search by removing the contribution of the dual reference frame,
// i.e. new_src(x,y) = 2 * src(x,y) - ref_dual(x,y).
void venc_hm_preproc_src(f265_enc_thread *t, f265_frac_ref_block *rb, f265_pix new_src[], int32_t stride_new_src)
{
    f265_me_ctx *mex = &t->me;

    f265_pix *refd = mex->hm_refd + rb->pos[1] * mex->hm_src_stride + rb->pos[0];
    f265_pix *src = mex->src_planes[0] + rb->pos[1] * mex->hm_src_stride + rb->pos[0];

    int16_t width = rb->size[0];
    int16_t height = rb->size[1];

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            new_src[x] = (src[x] << 1) - refd[x];
        }
        new_src += stride_new_src;
        src += mex->hm_src_stride;
        refd += mex->ref_stride;
    }
}

// Perform motion estimation.
// Arguments:
// mv_init: central point. For uni-directional, the central point will be equal to best predictor.
//          For bi-directional, it will be equal to best MV of the uni-directional search.
// mv_bits: motion vector bits.
// bipred:  bi-directional prediction. Set to 1 when performing the search in
//          the second reference. Set to 0 otherwise.
// bits and cost: cumulative bits and cost (to include reference list and index).
void venc_hm_me_search(f265_enc_thread *t, f265_cb *cb, f265_frac_ref_block *rb, int32_t bipred,
                       f265_mv *mv_init, int32_t *mv_bits, int32_t *bits, int32_t *cost)
{
    f265_me_ctx *mex = &t->me;
    f265_mv *best_mv = &mex->ref[0].hm_best_mv;

    // Fractional MV.
    f265_mv frac;

    // Use a local variable to avoid clipping mv_init.
    f265_mv mv = *mv_init;

    // Clip init MV.
    venc_hm_mv_clip(t, cb, &mv, 7);

    // Set search range.
    venc_hm_me_range(t, cb, &mv, mex->hm_range[bipred]);

    // Set initial MV for full pixel search.
    best_mv->x = mv.x >> 2;
    best_mv->y = mv.y >> 2;

    // Set cost weight to uni-prediction case.
    double weight = 1.0;

    if (bipred)
    {
        // Remove the contribution of the dual reference.
        mex->src_planes[0] = mex->hm_new_src;

        // Set weight to bi-prediction case.
        weight = 0.5;
    }

    // Init best cost.
    mex->best_cost = F265_MAX_SAD;

    //  Perform full pixel search.
    if (mex->hm_algo)
    {
        // TZ algorithm.
        venc_hm_tz_search(t, rb, cb);
    }
    else
    {
        // Full search.
        venc_hm_me_full_search(t, rb);
    }

    // Convert MV to quarter pixel resolution. Need to be done before the
    // fractional search because of the cost computation.
    best_mv->x <<= 2;
    best_mv->y <<= 2;

    // Perform fractional search.
    *cost = venc_hm_frac_search(t, cb, rb, &frac);

    // Set the best fractional MV.
    best_mv->x += frac.x;
    best_mv->y += frac.y;

    // Update MV bits and MV cost.
    f265_mv *mvp = &mex->ref[0].pmv;
    *mv_bits = venc_hm_bits_count(best_mv->x - mvp->x) + venc_hm_bits_count(best_mv->y - mvp->y);
    int32_t mv_cost = venc_hm_bits_cost(*mv_bits, mex->hm_lambda_me);

    // Update ME bits.
    *bits += *mv_bits;

    // Update ME cost.
    *cost = (int32_t)(floor(weight * ((double)*cost - (double)mv_cost))
            + (double)venc_hm_bits_cost(*bits, mex->hm_lambda_me));
}

// Compute inter prediction mode bits.
void venc_hm_pred_mode_bits(int32_t part_mode, int32_t is_p, int32_t part_idx,
                            int32_t prev_mode, int32_t pred_mode_bits[])
{
    switch (part_mode)
    {
        case F265_PART_UN:
            pred_mode_bits[0] = (!is_p) ? 3 : 1;
            pred_mode_bits[1] = 3;
            pred_mode_bits[2] = 5;
            break;

        case F265_PART_H1:
        case F265_PART_H2:
        case F265_PART_H3:
            if (is_p)
            {
                pred_mode_bits[0] = 3;
                pred_mode_bits[1] = 0;
                pred_mode_bits[2] = 0;
            }
            else
            {
                int32_t the_table[2][3][3] = {{{0,0,3}, {0,0,0}, {0,0,0}}, {{5,7,7}, {7,5,7}, {6,6,6}}};
                memcpy(pred_mode_bits, the_table[part_idx][prev_mode], 3*sizeof(int32_t));
            }
            break;

        case F265_PART_V1:
        case F265_PART_V2:
        case F265_PART_V3:
            if (is_p)
            {
                pred_mode_bits[0] = 3;
                pred_mode_bits[1] = 0;
                pred_mode_bits[2] = 0;
            }
            else
            {
                int32_t the_table[2][3][3] = {{{0,2,3}, {0,0,0}, {0,0,0} }, {{5,7,7}, {5,5,7}, {6,6,6}}};
                memcpy(pred_mode_bits, the_table[part_idx][prev_mode], 3 * sizeof(int32_t));
            }
            break;

        case F265_PART_HV:
            printf("F265 error: inter NxN partitioning not supported in f265 !\n");
            MC_ASSERT(0);
            break;

        default:
            printf("F265 error: invalid partition mode: %d !\n", part_mode);
            MC_ASSERT(0);
            break;
    }
}

// Select best MV predictor index (for MVD computation).
// Arguments:
// bits: bit needed to date to encode MV, MVP index, reference list and reference index.
// cost: best cost to date.
void venc_hm_select_mvp(f265_enc_thread *t, f265_mv mvp_list[], int32_t *mvp_idx,
                        int32_t *mv_bits, int32_t *bits, int32_t *cost)
{
    f265_me_ctx *mex = &t->me;

    // Note that the MVP list must contain exactly 2 MV, even if they are duplicate.

    // Quit if both MVP are equal.
    if (mvp_list[0].p == mvp_list[1].p) return;

    // Set the MVP index.
    int idx = 1 - *mvp_idx;

    // Get the MV cost in bits.
    // Note that in theory, 1 bit for the MVP index should be added. This bit is already
    // included in the function argument "bits", and is always equal to "1" regardless of
    // the selected MVP.
    int32_t tmp = venc_hm_bits_count(mex->ref[0].hm_best_mv.x - mvp_list[idx].x)
                + venc_hm_bits_count(mex->ref[0].hm_best_mv.y - mvp_list[idx].y);

    // If the current MVP is still the best return.
    if (tmp >= *mv_bits) return;

    // MVP index changed, so update.
    int32_t old_bits = *bits;

    // Update bits and cost.
    *bits = *bits - *mv_bits + tmp;
    *cost = *cost - venc_hm_bits_cost(old_bits, mex->hm_lambda_me) + venc_hm_bits_cost(*bits, mex->hm_lambda_me);

    // Save new best MVP and its index.
    mex->ref[0].pmv = mvp_list[idx];
    *mvp_idx =  idx;
}

// Set up the motion estimation.
// Most of the processing needs not be done for every reference/CB/partition.
// We put it here to keep the HM compatibility code localized.
//
// WARNING!
//
// This code changes the semantics of t->me.src_planes and co. In this code, the
// planes are not centered on the source block location.
//
// This code is not safe to use with frame multithreading (tile multithreading
// should work).
void venc_hm_setup_me(f265_enc_thread *t, f265_cb *cb, int32_t part_height,
                      int32_t bipred_flag, int32_t list, int32_t ref_idx, f265_weight *wp)
{
    f265_main_data *md = &t->enc->md;
    f265_frame *src_frame = t->src_frame;

    // Hadamard flag (to enable SATD for sub-pixel search).
    int had_flag = 1;

    // Fast ME (to enable TZ algorithm).
    #ifndef VAN_HM_RDO_USE_FS
    int fast_me_flag = 1;
    #else
    int fast_me_flag = 0;
    #endif

    // Fast encoding decision (to enable sub-sampled SAD).
    #if 0
    int fast_enc_decision = 1;
    #else
    int fast_enc_decision = 0;
    #endif

    // Get the ME context.
    f265_me_ctx *mex =  &t->me;

    // Chroma ME search.
    mex->chroma_flag = 0;

    // Distortion function ID (for sub-pixel search). 0:SAD, 1:HAD.
    mex->dist_func_id = had_flag;

    // Set SAD sub-sampling step.
    mex->hm_subshift = fast_enc_decision && part_height > 8 ? 1 : 0;

    // ME range.
    mex->hm_range[0] = 64;
    mex->hm_range[1] = 4;

    // Weighted prediction.
    wp->used_flag = 0;
    mex->ref[0].weights = mex->ref[1].weights = wp;

    // Search algorithm (set to fast uni and full for bi).
    mex->hm_algo = fast_me_flag && !bipred_flag;

    // MV Lambda.
    mex->hm_lambda_me = (uint32_t)floor(65536.0 * sqrt(t->hm_lambda[0]));

    // QP factor for B-slices is assumed to be equal to 1.0 when computing hm_lambda in
    // venc_hm_calc_lambda function. In case of mismatch, check the HM configuration file
    // for the QP factor and also for the GOP size.
    #ifdef VAN_LOAD_HM_ME
    MC_ASSERT(cb->hm_lambda_me == mex->hm_lambda_me);
    #endif

    // Set source frame.
    if (t->src_frame->frame_type == F265_FRAME_P)
    {
        // For uni-directional search, BD conversion is not needed.

        f265_ref_ctx *refx = &t->ref_ctx[list][ref_idx];
        f265_frame *ref_frame = refx->frame;
        mex->src_planes[0] = t->src_frame->src_planes[0];
        mex->hm_src_stride = t->enc->gd.stride;
        mex->ref[0].ref_planes[0] = ref_frame->rec_planes[0];
        mex->ref_stride = t->enc->gd.stride;
        return;
    }

    // Convert source and reconstructed reference frames from
    // 8u bit depth to 16s bit depth.
    int stride = t->enc->gd.stride;
    int *pix_dim = t->enc->gd.pix_dim;

    // Reference frame offset (from buffer origin to picture origin).
    int32_t ref_off = F265_LUMA_PLANE_PADDING * stride + F265_LUMA_PLANE_PADDING;

    // Set the plane pointers.
    mex->src_planes[0] = md->hm_src_plane_16s;
    mex->hm_new_src = md->hm_new_hm_src_plane_16s;
    mex->hm_refd = md->hm_ref_dual_plane_16s;
    mex->ref[0].ref_planes[0] = md->hm_ref_plane_16s[list][ref_idx] + ref_off;

    // Set stride.
    mex->hm_src_stride = mex->ref_stride = stride;

    // If a new frame, convert the whole DPB. Otherwise, quit.
    if (src_frame->abs_poc == md->hm_last_poc) return;

    // Update POC.
    md->hm_last_poc = src_frame->abs_poc;

    // Convert source frame.
    venc_hm_copy_8u_to_16s(md->hm_src_plane_16s, stride, (uint8_t *)t->src_frame->src_planes[0],
                           stride, stride, pix_dim[1]);

    // Convert reference frames (the whole DPB).
    int num_list = t->src_frame->frame_type;
    for (int l = 0; l < num_list; l++)
    {
        int num_ref = t->nb_ref_idx[l];
        for (int r = 0; r < num_ref; r++)
        {
            f265_ref_ctx *refx = &t->ref_ctx[l][r];
            f265_frame *ref_frame = refx->frame;
            venc_hm_copy_8u_to_16s(md->hm_ref_plane_16s[l][r], stride, (uint8_t *)ref_frame->rec_planes[0] - ref_off,
                                   stride, stride, pix_dim[1] + 2 * F265_LUMA_PLANE_PADDING);
        }
    }
}

// Map the index of the current reference frame from l1 to l0.
int venc_hm_map_l1_to_l0(f265_enc_thread *t, int32_t ref_idx)
{
    f265_frame *ref_l0;
    f265_frame *ref_l1;

    ref_l1 = t->ref_ctx[1][ref_idx].frame;

    for (int i = 0; i < t->nb_ref_idx[1]; i++)
    {
        ref_l0 = t->ref_ctx[0][i].frame;
        if (ref_l0->abs_poc == ref_l1->abs_poc) return i;
    }
    return -1;
}

// Temp function to test HM compatible inter search.
// Encode a CTB. Only inter is enabled. Dispatch lbd/hbd functions
// according to frame type.
#if 0
void venc_hm_analyze_inter_ctb(f265_enc_thread *t)
{
    void f265_lbd_hm_analyze_inter_cb(f265_enc_thread *, f265_cb *);
    void f265_hbd_hm_analyze_inter_cb(f265_enc_thread *, f265_cb *);

    if (t->src_frame->frame_type == F265_FRAME_P)
        f265_lbd_hm_analyze_inter_cb(t, t->cb);
    else f265_hbd_hm_analyze_inter_cb(t, t->cb);
}
#endif

#if 0
// Temp function to test HM compatible inter search.
// Encode a CB. Only inter is enabled.
void venc_hm_analyze_inter_cb(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return;

    // Handle split CB.
    if (cb->flags&F265_CB_SPLIT)
    {
        for (int i = 0; i < 4; i++) venc_hm_analyze_inter_cb(t, t->cb + cb->child_idx + i);
        return;
    }

    // Inter.
    int32_t part_mode = cb->inter_part;
    if (!(cb->flags&F265_CB_INTRA) && !(cb->flags&F265_CB_SKIP))
        venc_hm_search_inter_cb(t, cb, part_mode);
}
#endif

// Search for the best inter prediction mode, reference list, reference indices,
// motion vectors and motion vector predictor indices.
// Arguments:
// part_idx: partition index to be tested.
// prev_part_mode: best mode selected of previous partition in same CB.
int32_t venc_hm_search_inter_cb(f265_enc_thread *t, f265_cb *cb, int32_t part_idx, int32_t *prev_pred_mode)
{
    // Do not perform ME if it's a merge mode.
    int merge_flag = (cb->inter_bf[part_idx] & 7) != 5;
    if (merge_flag) return 0;

    // Partition mode.
    int32_t part_mode = cb->inter_part;

    // Get the ME context.
    f265_me_ctx *mex =  &t->me;

    // Weighted prediction.
    f265_weight wp;

    // Bi-directional prediction flag.
    int bi_flag;

    // Variables for cost comparison.
    int mv_bits[2], bits, cost;
    int32_t mvp_idx;
    int best_cost[3] = {F265_MAX_SAD, F265_MAX_SAD, F265_MAX_SAD};
    int best_bits[3];

    // Prediction mode bits.
    int32_t pred_mode_bits[3];

    // Best reference index for each list.
    int best_ref_idx[2];

    // Best motion vector for each list.
    f265_mv best_mv[2];

    // Best motion vector predictor index for each list.
    int best_mvp_idx[2] = {-1,-1};

    // Buffers to save best MV and MVP index for each reference frame.
    f265_mv mv_buf[2][16];
    int32_t mvp_idx_buf[2][16];

    // Buffers to save cost, bits and MV to map l0 search result to l1.
    int cost_l0[16], bits_l0[16];
    f265_mv mv_l0[16];

    // Map reference index from l1 to l0.
    int l0_idx;

    // Enable/disable l1 search result.
    int cost_validl1 = F265_MAX_SAD;
    int bits_validl1 = F265_MAX_SAD;
    int ref_idx_validl1 = 0;
    int mvp_idx_validl1 = 0;
    f265_mv mv_validl1; mv_validl1.p = 0;

    // Get the neighbours.
    f265_inter_neighbour_mv neighbours[6];
    venc_get_neighbour_mvs(t, cb, part_idx, neighbours);

    // Get partition size and offset.
    int part_off[2];
    int cb_size = 1 << cb->lg_bs;
    int part_size[2];
    venc_get_part_loc(part_size, part_off, part_mode, part_idx, cb_size);

    f265_frac_ref_block rb;
    f265_ref_ctx ctx[2];
    rb.ctx[0] = &ctx[0];
    rb.ctx[1] = &ctx[1];

    // Set partition size and position.
    rb.size[0] = part_size[0];
    rb.size[1] = part_size[1];
    rb.pos[0] = cb->cb_off[0] + t->ctb_off[0] + part_off[0];
    rb.pos[1] = cb->cb_off[1] + t->ctb_off[1] + part_off[1];
    rb.bipred_flag = bi_flag = 0;

    // Get prediction mode (reference list selection) bits.
    int32_t is_p = t->src_frame->frame_type == F265_FRAME_P;
    venc_hm_pred_mode_bits(part_mode, is_p, part_idx, *prev_pred_mode, pred_mode_bits);

    // Get the number of lists.
    int list, num_list = t->src_frame->frame_type;

    // Loop over available lists.
    for (list = 0; list < num_list; list++)
    {
        int num_ref = t->nb_ref_idx[list];
        for (int8_t ref_idx = 0; ref_idx < num_ref; ref_idx++)
        {
            // Get number of bits used to encode the partition mode.
            bits = pred_mode_bits[list];

            // Add reference bits.
            if (num_ref > 1)
            {
                // Add bits for reference index.
                bits += ref_idx + 1;

                // Remove 1 bit if reference index points to the last reference frame.
                if (ref_idx == num_ref - 1) bits--;
            }

            // Add 1 bit for MVP index.
            bits += 1;

            // Setup ME.
            venc_hm_setup_me(t, cb, part_size[1], bi_flag, list, ref_idx, &wp);

            // Get the MV predictors.
            f265_mv mvp_list[2];
            venc_get_pmv(t, part_idx, neighbours, ref_idx, list, mvp_list);

            // Test MV predictors and select the best.
            mvp_idx = venc_hm_test_mvp(t, cb, mvp_list, &rb);

            // Do not perform ME if the current reference was tested in l0.
            // Reuse the MVs from l0 and update cost and bits count.

            // A word about l0 to l1 cost inference.
            // If a reference frame, when searching in l1, was already tested
            // in l0, the motion estimation is not performed. In this case, the
            // MV is inherited (and the cost and bits count adjusted accordingly)
            // from the search in l0 even if the best MV could be different (because
            // of the MVP being different).
            // The inherited best MV (and its cost) is only used to determine the best
            // reference when performing bi-directional search. It will be ignored when
            // performing uni-directional prediction, meaning that the reference
            // from l1 is taken only if a full ME search was performed.

            l0_idx = -1;
            if (list) l0_idx = venc_hm_map_l1_to_l0(t, ref_idx);
            if (l0_idx >= 0)
            {
                // Set the best MV.
                mex->ref[0].hm_best_mv = mv_l0[l0_idx];

                // Set total cost.
                cost = cost_l0[l0_idx];

                // Remove cost due to bits consumption.
                cost -= venc_hm_bits_cost(bits_l0[l0_idx], mex->hm_lambda_me);

                // Get MV bits
                mv_bits[list] = venc_hm_bits_count(mv_l0[l0_idx].x - mex->ref[0].pmv.x)
                              + venc_hm_bits_count(mv_l0[l0_idx].y - mex->ref[0].pmv.y);

                // Add MV bits to bit count.
                bits += mv_bits[list];

                // Add cost due to bits consumption.
                cost += venc_hm_bits_cost(bits, mex->hm_lambda_me);
            }
            else
            {
                // Call ME search.
                venc_hm_me_search(t, cb, &rb, bi_flag, &mex->ref[0].pmv, &mv_bits[list], &bits, &cost);
            }

            // Select the best MVP index for MVD coding. Not performed if the MVPs are equal.
            venc_hm_select_mvp(t, mvp_list, &mvp_idx, &mv_bits[list], &bits, &cost);

            // Save best MV of the current reference frame.
            mv_buf[list][ref_idx] = mex->ref[0].hm_best_mv;
            mvp_idx_buf[list][ref_idx] = mvp_idx;

            // Save result from l0 to use in l1.
            if (!list)
            {
                cost_l0[ref_idx] = cost;
                bits_l0[ref_idx] = bits;
                mv_l0[ref_idx] = mex->ref[0].hm_best_mv;
            }

            // Check best and update.
            if (cost < best_cost[list])
            {
                best_cost[list] = cost;
                best_bits[list] = bits;
                best_mv[list] = mex->ref[0].hm_best_mv;
                best_ref_idx[list] = ref_idx;
                best_mvp_idx[list] = mvp_idx;
            }

            if (list && cost < cost_validl1 && l0_idx < 0)
            {
                // For l1, update only if the reference frame was actually
                // searched (not inferred from l0 search). See comment
                // above on l0/l1 cost mapping.
                mv_validl1 = mex->ref[0].hm_best_mv;
                cost_validl1 = cost;
                bits_validl1 = bits;
                ref_idx_validl1 = ref_idx;
                mvp_idx_validl1 =  mvp_idx;
            }
        } // Loop over ref_idx
    } // Loop over list_idx

    // Perform bi-directional inter prediction if required, otherwise finalize
    // uni-directional search.
    if (t->src_frame->frame_type != F265_FRAME_B) goto check_best;

    // Partitions 8x4 and 4x8 not allowed in bi-directional inter.
    if (cb_size == 8 && (part_size[0] < 8 || part_size[1] < 8)) goto check_best;

    // The following code assumes that only one bi-prediction iteration is
    // used (fast encoder setting).

    // Set bi-directional inter to true.
    rb.bipred_flag = bi_flag = 1;

    // Working variables for bi-directional inter search.
    f265_mv mv_bi[2];
    int8_t ref_idx_bi[2];
    int32_t me_bits_bi[2];
    int mvp_idx_bi[2];

    for (int l = 0; l < 2; l++)
    {
        // Make a copy of MV and ref_idx from l0/l1 uni-directional search.
        mv_bi[l] = best_mv[l];
        mvp_idx_bi[l] = best_mvp_idx[l];
        ref_idx_bi[l] = best_ref_idx[l];

        // Compute ME bits.
        me_bits_bi[l] = best_bits[l] - pred_mode_bits[l];
    }

    // Set the reference to search according to best cost.
    // The reference with best cost will have its contribution
    // removed from the source, and the search will performed on
    // the dual reference (from the other list).
    list = 0;

    #if 0
    // This assumes fast encoder decision is on.
    // It's not the case in the tests I'm running.
    // FC 2014/01/10
    if (best_cost[0] <= best_cost[1]) list = 1;
    #endif

    // Set the dual reference list.
    int list_d = 1 - list;

    // Setup ME.
    venc_hm_setup_me(t, cb, part_size[1], bi_flag, list_d, ref_idx_bi[list_d], &wp);

    // Perform motion compensation using best MV and best reference index from uni-directional search.
    f265_pix *tmp  = mex->hm_refd + rb.pos[1]*mex->ref_stride + rb.pos[0];
    venc_hm_upsample(t, cb, &rb, &mv_bi[list_d], tmp, mex->ref_stride);

    // Remove the contribution of the dual reference.
    tmp = mex->hm_new_src + rb.pos[1]*mex->hm_src_stride + rb.pos[0];
    venc_hm_preproc_src(t, &rb, tmp, mex->hm_src_stride);

    // Get the MV predictors.
    f265_mv mvp_list[2];
    venc_get_pmv(t, part_idx, neighbours, ref_idx_bi[list_d], list_d, mvp_list);

    // Set MV bits for the dual reference.
    mv_bits[list_d] = venc_hm_bits_count(mv_bi[list_d].x - mvp_list[mvp_idx_bi[list_d]].x)
                    + venc_hm_bits_count(mv_bi[list_d].y - mvp_list[mvp_idx_bi[list_d]].y);

    // Loop over the reference frames.
    int num_ref = t->nb_ref_idx[list];
    for (int8_t ref_idx = 0; ref_idx < num_ref; ref_idx++)
    {
        bits = pred_mode_bits[2] + me_bits_bi[list_d];
        if (num_ref > 1)
        {
            bits += ref_idx + 1;
            if (ref_idx == num_ref - 1) bits -= 1;
        }

        // Add 1 bit for MVP index.
        bits += 1;

        // Setup ME.
        venc_hm_setup_me(t, cb, part_size[1], bi_flag, list, ref_idx, &wp);

        // Get the MV predictors.
        venc_get_pmv(t, part_idx, neighbours, ref_idx, list, mvp_list);

        // Set the predictor for the current reference.
        mex->ref[0].pmv = mvp_list[mvp_idx_buf[list][ref_idx]];

        // Perform ME.
        venc_hm_me_search(t, cb, &rb, bi_flag, &mv_buf[list][ref_idx], &mv_bits[list], &bits, &cost);

        // Select the best MVP index for MVD coding.
        venc_hm_select_mvp(t, mvp_list, &mvp_idx_buf[list][ref_idx], &mv_bits[list], &bits, &cost);

        // Check cost, update if better.
        if (cost < best_cost[2])
        {
            best_cost[2] = cost;
            best_bits[2] = bits;
            me_bits_bi[list] = bits - pred_mode_bits[2] - me_bits_bi[list_d];
            mv_bi[list] = mex->ref[0].hm_best_mv;
            ref_idx_bi[list] = ref_idx;
            mvp_idx_bi[list] = mvp_idx_buf[list][ref_idx];
        }
    }  // Loop over ref_idx.

    // Finalize best selection.
    check_best:

    // This is used to invalidate the cost from l1 if the reference frame
    // was not actually searched. See comment above on l0/l1 mapping.
    best_mv[1] = mv_validl1;
    best_bits[1] = bits_validl1;
    best_cost[1] = cost_validl1;
    best_ref_idx[1] = ref_idx_validl1;
    best_mvp_idx[1] = mvp_idx_validl1;

    // Final check best prediction mode (l0, l1 or bi).
    if (best_cost[2] <= best_cost[0] && best_cost[2] <= best_cost[1])
    {
        // Transfer the best values to the CB.
        cb->inter_bf[part_idx] = 5;
        for (int i = 0; i < 2; i++, list=!list)
        {
            cb->ref_idx[part_idx][i] = ref_idx_bi[i];
            if (cb->ref_idx[part_idx][i] != -1) cb->mv[part_idx][i] = mv_bi[i];
            else cb->mv[part_idx][i].p = 0;
            cb->inter_bf[part_idx] |= mvp_idx_bi[i]<<(3+i);
        }

        *prev_pred_mode = 2;
        return best_bits[2];
    }
    else if (best_cost[0] <= best_cost[1])
    {
        // Transfer the best values to the CB.
        cb->mv[part_idx][0] = best_mv[0];
        cb->mv[part_idx][1].p = 0;
        cb->ref_idx[part_idx][0] = best_ref_idx[0];
        cb->ref_idx[part_idx][1] = -1;
        cb->inter_bf[part_idx] = (best_mvp_idx[0]<<3)|5;

        *prev_pred_mode = 0;
        return best_bits[0];
    }
    else
    {
        // Transfer the best values to the CB.
        cb->mv[part_idx][0].p = 0;
        cb->mv[part_idx][1] = best_mv[1];
        cb->ref_idx[part_idx][0] = -1;
        cb->ref_idx[part_idx][1] = best_ref_idx[1];
        cb->inter_bf[part_idx] = (best_mvp_idx[1]<<4)|5;

        *prev_pred_mode = 1;
        return best_bits[1];
    }
}

#endif // VAN_USE_HM_ME

