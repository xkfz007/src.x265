// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Intra prediction.

#include "f265/enc.h"

// Get the number of available neighbour pixels on the top/left edge of the
// specified intra block.
static finline void venc_get_intra_nb_avail(f265_enc_thread *t, int avail[2], int comp, int ct_ox, int ct_oy)
{
    f265_gen_data *gd = &t->enc->gd;
    int csf = !!comp;

    // Convert the chroma case to the luma case.
    ct_ox<<=csf;
    ct_oy<<=csf;

    // 4x4 block stride in a CTB.
    int b4_stride = gd->ctb_size>>2;

    // Raster scan address of the current 4x4 block.
    int b4_xy = (ct_oy*b4_stride + ct_ox)>>2;

    // Query the intra availability map.
    uint32_t packed_avail = gd->intra_avail_map[b4_xy];
    int availx = 1<<(packed_avail&0xf);
    int availy = 1<<(packed_avail>>4);

    // Apply the intra cutoff.
    availx = F265_MIN(availx, t->intra_avail_cutoff[0] - ct_ox);
    availy = F265_MIN(availy, t->intra_avail_cutoff[1] - ct_oy);

    // Handle the availability of CTB A and B.
    int ctb_a_avail_flag = !!(t->ctb_neighbours & 8);
    int ctb_b_avail_flag = !!(t->ctb_neighbours & 2);
    if (!(ct_oy|ctb_b_avail_flag)) availx = 0;
    if (!(ct_ox|ctb_a_avail_flag)) availy = 0;

    avail[0] = availx>>csf;
    avail[1] = availy>>csf;
}

// Get the intra filtering flags.
// - filter_edge_flag: true if the edges are filtered for DC/vertical/horizontal.
// - filter_neighbour_flag: true if the neighbours are filtered.
// - neighbour_bilinear_flag: true if the neighbours are filtered bilinearly.
void venc_get_intra_filter_flags(int *filter_edge_flag, int *filter_neighbour_flag, int *neighbour_bilinear_flag,
                                 int comp, int lg_bs, int mode, int smooth_intra_flag)
{
    int bs = 1<<lg_bs;
    *filter_edge_flag = !comp && bs != 32;
    *filter_neighbour_flag = !comp && f265_intra_mode_dist[mode] > f265_intra_dist_thresholds[lg_bs-2];
    *neighbour_bilinear_flag = bs == 32 && smooth_intra_flag;
}

// Get the intra encoding flags.
// - dst_flag: true if the 4x4 DST is used.
// - order: coefficient scan order.
void venc_get_intra_encode_flags(int *dst_flag, int *order, int comp, int lg_bs, int mode)
{
    *dst_flag = !comp && lg_bs == 2;
    *order = 0;
    if (lg_bs == 2 || !comp && lg_bs == 3)
        *order = mode >= 22 && mode <= 30 ? 1 : mode >= 6 && mode <= 14 ? 2 : 0;
}

// Extract the unfiltered neighbour pixels of the specified intra block for one
// image component. 'src' points to the top-left neighbour pixel of the block,
// 'nx' and 'ny' are the number of pixels to predict in each direction. 'bd' is
// the bit depth.
static finline void venc_extract_intra_neighbours(f265_pix dst[129], f265_pix *src, int src_stride,
                                                  int availx, int availy, int nx, int ny, int bd)
{
    // The following logic relies on the slice layout restrictions.

    // Copy top-left tentatively.
    dst[0] = src[0];

    // Left is fully available, copy.
    if (likely(availy >= ny))
    {
        for (int i = 0; i < ny; i++) dst[65+i] = src[(1+i)*src_stride];
    }

    // Left is partially available, copy and broadcast.
    else if (likely(availy > 0))
    {
        for (int i = 0; i < availy; i++) dst[65+i] = src[(1+i)*src_stride];
        f265_pix p = dst[64+availy];
        for (int i = availy; i < ny; i++) dst[65+i] = p;
    }

    // Left and top-left are not available but top is. Broadcast the first
    // pixel directly above the block.
    else if (likely(availx > 0))
    {
        f265_pix p = src[1];
        dst[0] = p;
        for (int i = 0; i < ny; i++) dst[65+i] = p;
    }

    // Nothing is available, perform DC prediction.
    else
    {
        f265_pix p = 1<<(bd-1);
        for (int i = 0; i < nx+1; i++) dst[i] = p;
        for (int i = 0; i < ny; i++) dst[65+i] = p;
        return;
    }

    // Top is fully available, copy.
    if (likely(availx >= nx))
    {
        for (int i = 0; i < nx; i++) dst[1+i] = src[1+i];
    }

    // Top is partially available, copy and broadcast.
    else if (likely(availx > 0))
    {
        for (int i = 0; i < availx; i++) dst[1+i] = src[1+i];
        f265_pix p = dst[availx];
        for (int i = availx; i < nx; i++) dst[1+i] = p;
    }

    // Top-left, top, top-right are not available. Broadcast the first pixel
    // directly left of the block.
    else
    {
        f265_pix p = dst[65];
        for (int i = 0; i < nx+1; i++) dst[i] = p;
    }
}

// Predict the unfiltered neighbour pixels of the specified intra block at
// pixel offset (ct_ox, ct_oy) in the CTB with block size 'bs'. 'rec_flag' is
// true if the reconstructed pixels are used for the prediction, false if the
// source pixels are used as approximation. This function assumes that
// constrained intra prediction is not used.
//
// Layout of the destination array, by offset: 0 (top-left), 1 (top and
// top-right), 65 (left and bottom left).
void venc_predict_intra_neighbours(f265_enc_thread *t, f265_pix dst[129], int rec_flag,
                                   int comp, int bs, int ct_ox, int ct_oy)
{
    int avail[2];
    venc_get_intra_nb_avail(t, avail, comp, ct_ox, ct_oy);
    int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox-1, ct_oy-1);
    f265_pix *src = rec_flag ? t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off :
                               t->src_frame->src_planes[comp] + plane_off;
    venc_extract_intra_neighbours(dst, src, t->me.ref_stride, avail[0], avail[1], bs<<1, bs<<1,
                                  t->enc->gd.bit_depth[!!comp]);
}

// Filter the neighbour pixels of the block specified.
void venc_filter_intra_neighbours(f265_pix *dst, f265_pix *src, int bs, int bd, int bilinear_flag)
{
    int bs2 = bs<<1;
    int top_left = src[0], top_last = src[bs2], left_last = src[64+bs2];

    // Check for bilinear filtering.
    if (bilinear_flag)
    {
        int top_middle = src[32], left_middle = src[64+32];
        int threshold = 1<<(bd-5);
        bilinear_flag = F265_ABS(top_left + top_last - (top_middle<<1)) < threshold &&
                        F265_ABS(top_left + left_last - (left_middle<<1)) < threshold;
        if (bilinear_flag)
        {
            dst[0] = top_left;
            dst[64] = top_last;
            dst[64+64] = left_last;
            for (int i = 0; i < 63; i++)
            {
                dst[1+i] =  ((63-i)*top_left + (i+1)*top_last  + 32)>>6;
                dst[65+i] = ((63-i)*top_left + (i+1)*left_last + 32)>>6;
            }
            return;
        }
    }

    // Regular filtering.
    dst[0] = ((top_left<<1) + src[1] + src[65] + 2)>>2;
    dst[bs2] = top_last;
    dst[64+bs2] = left_last;
    for (int i = 1; i < bs2; i++) dst[i] = ((src[i]<<1) + src[i-1] + src[i+1] + 2)>>2;
    dst[65] = ((src[65]<<1) + top_left + src[66] + 2)>>2;
    for (int i = 66; i < 64+bs2; i++) dst[i] = ((src[i]<<1) + src[i-1] + src[i+1] + 2)>>2;
}

// Intra planar prediction.
void venc_predict_intra_planar(f265_pix *dst, f265_pix *nbuf, int lg_bs)
{
    int bs = 1<<lg_bs;
    int top_right = nbuf[1+bs];
    int bottom_left = nbuf[65+bs];
    for (int y = 0; y < bs; y++)
        for (int x = 0; x < bs; x++)
            dst[y*bs+x] = ((bs-1-x)*nbuf[65+y] + (bs-1-y)*nbuf[1+x] + (x+1)*top_right + (y+1)*bottom_left + bs)
                          >>(lg_bs+1);
}

// Intra DC prediction.
void venc_predict_intra_dc(f265_pix *dst, f265_pix *nbuf, int lg_bs, int filter_edge_flag)
{
    int bs = 1<<lg_bs;
    int dc_val = bs;
    for (int i = 0; i < bs; i++) dc_val += nbuf[1+i] + nbuf[65+i];
    dc_val = dc_val>>(lg_bs+1);
    for (int i = 0; i < bs*bs; i++) dst[i] = dc_val;

    if (filter_edge_flag)
    {
        dst[0] = ((dc_val<<1) + nbuf[1] + nbuf[65] + 2)>>2;
        for (int i = 1; i < bs; i++)
        {
            dst[i] = (nbuf[i+1] + 3*dc_val + 2)>>2;
            dst[i*bs] = (nbuf[65+i] + 3*dc_val + 2)>>2;
        }
    }
}

// Intra angular prediction.
void venc_predict_intra_angular(f265_pix *dst, f265_pix *nbuf, int lg_bs, int bd, int filter_edge_flag, int mode)
{
    int bs = 1<<lg_bs;

    // Flip the neighbours in the horizontal case.
    int hor_flag = mode < 18;
    f265_pix ntmp[129];
    if (hor_flag)
    {
        ntmp[0] = nbuf[0];
        for (int i = 0; i < (bs<<1); i++) { ntmp[1+i] = nbuf[65+i]; ntmp[65+i] = nbuf[1+i]; }
        nbuf = ntmp;
    }

    // Get the prediction angle.
    int angle_off = hor_flag ? 10-mode : mode-26;
    int angle = f265_intra_angle_table[8+angle_off];

    // Vertical prediction.
    if (!angle)
    {
        for (int y = 0; y < bs; y++)
            for (int x = 0; x < bs; x++)
                dst[y*bs+x] = nbuf[1+x];

        if (filter_edge_flag)
        {
            int top_left = nbuf[0], top = nbuf[1];
            for (int y = 0; y < bs; y++)
                dst[y*bs] = F265_CLAMP(top + ((nbuf[65+y] - top_left)>>1), 0, (1<<bd)-1);
        }
    }

    // Angular prediction.
    else
    {
        // Get the reference pixels. The reference base is the first pixel to the
        // top (nbuf[1]).
        f265_pix ref_buf[64], *ref;

        // Use the projected left neighbours and the top neighbours.
        if (angle < 0)
        {
            // Number of neighbours projected. Mind the negative shift properties
            // (-1>>5 == -32>>5 == -1). Note that the specification projects one
            // more pixel than necessary (off-by-one error).
            int nb_projected = -((bs*angle)>>5) - 1;
            ref = ref_buf + nb_projected + 1;

            // Project the neighbours.
            int inv_angle = f265_intra_inv_angle_table[-angle_off-1];
            int inv_angle_sum = 128;
            for (int i = 0; i < nb_projected; i++)
            {
                inv_angle_sum += inv_angle;
                ref[-2-i] = nbuf[64+(inv_angle_sum>>8)];
            }

            // Copy the top-left and top pixels.
            for (int i = 0; i < bs+1; i++) ref[-1+i] = nbuf[i];
        }

        // Use the top and top-right neighbours.
        else ref = nbuf+1;

        // Pass every row.
        int angle_sum = 0;
        for (int y = 0; y < bs; y++)
        {
            angle_sum += angle;
            int off = angle_sum>>5;
            int frac = angle_sum&31;

            // Interpolate.
            if (frac) for (int x = 0; x < bs; x++) dst[y*bs+x] = ((32-frac)*ref[off+x] + frac*ref[off+x+1] + 16)>>5;

            // Copy.
            else for (int x = 0; x < bs; x++) dst[y*bs+x] = ref[off+x];
        }
    }

    // Flip for horizontal.
    if (hor_flag)
        for (int y = 0; y < bs; y++)
            for (int x = y+1; x < bs; x++)
            {
                F265_SWAP(int, dst[y*bs+x], dst[x*bs+y]);
            }
}

// Predict the pixels of the intra mode specified.
void venc_predict_intra_mode(f265_pix *dst, f265_pix *neighbours, int lg_bs, int bd, int mode, int filter_edge_flag)
{
    if (mode == 0) venc_predict_intra_planar(dst, neighbours, lg_bs);
    else if (mode == 1) venc_predict_intra_dc(dst, neighbours, lg_bs, filter_edge_flag);
    else venc_predict_intra_angular(dst, neighbours, lg_bs, bd, filter_edge_flag, mode);
}

// Predict the intra block with the mode and the CTB offset specified.
void venc_predict_intra(f265_enc_thread *t, f265_pix *dst, int comp, int lg_bs, int mode, int ct_ox, int ct_oy)
{
    int chroma_flag = !!comp;
    int bs = 1<<lg_bs;
    int bd = t->enc->gd.bit_depth[chroma_flag];
    int smooth_intra_flag = F265_GET_FLAG(t->enc->gd.eflags, F265_PF_SMOOTH_INTRA);
    int filter_edge_flag, filter_neighbour_flag, neighbour_bilinear_flag;

    // Unfiltered/filtered neighbours.
    F265_ALIGN64 f265_pix nbuf[2][129];

    // Predict the unfiltered neighbours. Assuming 4:2:0.
    venc_predict_intra_neighbours(t, nbuf[0], 1, comp, bs, ct_ox, ct_oy);

    // Filter the neighbours.
    venc_get_intra_filter_flags(&filter_edge_flag, &filter_neighbour_flag, &neighbour_bilinear_flag,
                                comp, lg_bs, mode, smooth_intra_flag);
    if (filter_neighbour_flag) venc_filter_intra_neighbours(nbuf[1], nbuf[0], bs, bd, neighbour_bilinear_flag);
    f265_pix *neighbours = nbuf[filter_neighbour_flag];

    // Do the prediction.
    venc_predict_intra_mode(dst, neighbours, lg_bs, bd, mode, filter_edge_flag);
}

// Load the neighbour's intra mode at the (X,Y) 4x4 block offset from the CTB
// origin.
static int venc_get_intra_neighbour_at(f265_enc_thread *t, int x, int y)
{
    // The block is outside the CTB, predict DC. FIXME: set the CTB neighbours
    // correctly to avoid testing for this.
    if (y < 0) return 1;

    // Get the coding block mode.
    int cb_idx;
    int part_idx;
    venc_lookup_pmap(t, x, y, &cb_idx, &part_idx);
    f265_cb *cb = t->cb + cb_idx;

    // Return DC if not intra.
    if ((cb->flags & (F265_CB_INTRA | F265_CB_PRESENT)) != (F265_CB_INTRA | F265_CB_PRESENT))
        return 1;

    int luma_mode = cb->intra_luma_mode[part_idx];

    // Return DC if PCM.
    if (luma_mode == 35) return 1;

    return luma_mode;
}

// Generate the list of three luma intra MPMs.
void venc_get_intra_pred_mode(f265_enc_thread *t, f265_cb *cb, int partition_idx, int *mpm_list)
{
    // Get the current CB position.
    int part_size[2], part_off[2];
    venc_get_part_loc4(part_size, part_off, cb->intra_luma_mode[1] == -1 ? 0 : 3, partition_idx, 1<<cb->lg_bs);
    int off[2];
    for (int i = 0; i < 2; i++)
        off[i] = ((cb->cb_off[i]) >> 2) + part_off[i];

    // Get the neighbour modes.
    mpm_list[0] = venc_get_intra_neighbour_at(t, off[0] - 1, off[1]);
    mpm_list[1] = venc_get_intra_neighbour_at(t, off[0],     off[1] - 1);

    // Duplicate mode.
    if (mpm_list[0] == mpm_list[1])
    {
        if (mpm_list[0] < 2)
        {
            mpm_list[0] = 0;
            mpm_list[1] = 1;
            mpm_list[2] = 26;
        }
        else
        {
            mpm_list[1] = 2 + ((mpm_list[0] + 29) % 32);
            mpm_list[2] = 2 + ((mpm_list[0] - 1) % 32);
        }
    }

    // Non-duplicate.
    else
    {
        int plane_present = (mpm_list[0] == 0) | (mpm_list[1] == 0);
        int dc_present = (mpm_list[0] == 1) | (mpm_list[1] == 1);
        int mask = (dc_present << 1) | plane_present;
        uint8_t default_table[4] = { 0, 1, 0, 26 };
        mpm_list[2] = default_table[mask];
    }
}

