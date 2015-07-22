// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Inter prediction.

#include "f265/enc.h"

// Perform the motion compensation for the coding block and image component
// specified.
void venc_mc_cb(f265_enc_thread *t, f265_cb *cb, f265_pix *dst, int dst_stride, int comp)
{
    int plane_stride = t->plane_stride;
    int chroma_flag = !!comp;
    int part_mode = cb->inter_part;
    int nb_parts = f265_nb_parts[part_mode];
    int pred_size = 1<<(cb->lg_bs-chroma_flag);
    int cb_plane_off = ((t->ctb_off[1] + cb->cb_off[1])>>chroma_flag)*plane_stride +
                       ((t->ctb_off[0] + cb->cb_off[0])>>chroma_flag);

    // Pass each partition.
    for (int part_idx = 0; part_idx < nb_parts; part_idx++)
    {
        // Set the block location.
        int part_size[2], part_off[2];
        venc_get_part_loc(part_size, part_off, part_mode, part_idx, pred_size);
        int plane_off = cb_plane_off + part_off[1]*plane_stride + part_off[0];
        int awi = f265_width_to_awi[part_size[0]];
        int packed_dims = awi<<24|part_size[0]<<8|part_size[1];
        f265_pix *dst_part = dst + part_off[1]*dst_stride + part_off[0];

        // Pass each list.
        f265_ref_ctx *rcs[2];
        f265_mv mvs[2];
        int list_count = 0;
        for (int list = 0; list < 2; list++)
        {
            int ref_idx = cb->ref_idx[part_idx][list];
            if (ref_idx == -1) continue;
            rcs[list_count] = t->ref_ctx[list] + ref_idx;
            mvs[list_count] = f265_clip_mv(cb->mv[part_idx][list], t->mc_bounds64);
            list_count++;
        }

        // Bi-prediction.
        if (list_count == 2)
        {
            if (chroma_flag) venc_mc_chroma_b(t, dst_part, dst_stride, rcs, mvs, packed_dims, plane_off, comp);
            else venc_mc_luma_b(t, dst_part, dst_stride, rcs, mvs, packed_dims, plane_off);
        }

        // Uni-prediction.
        else
        {
            if (chroma_flag) venc_mc_chroma_p(t, dst_part, dst_stride, rcs[0], mvs[0], packed_dims, plane_off, comp);
            else venc_mc_luma_p(t, dst_part, dst_stride, rcs[0], mvs[0], packed_dims, plane_off);
        }
    }
}

// Eventually auto-generate the following code.

// Do the motion compensation for the current uni-predicted block.
void venc_mc_luma_p(f265_enc_thread *t, f265_pix *dst, int dst_stride, f265_ref_ctx *rc, f265_mv mv,
                    int packed_dims, int plane_off)
{
    int plane_stride = t->plane_stride;
    int mx = mv.x, my = mv.y;
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff, awi = packed_dims>>24;
    int ref_off = plane_off + (my>>2)*plane_stride + (mx>>2);

    // No weighted prediction.
    if (!rc->weights->used_flag)
    {
        int hpel_idx = ((my&3)<<2) + (mx&3);

        // Interpolate.
        if (likely(hpel_idx&5))
        {
            int frac_x = mx&3, frac_y = my&3;
            int frac = (frac_y<<2)|frac_x, func_idx = 3*awi + (!!frac_y<<1) + !!frac_x - 1;
            venc_interpol_luma_qpel_pix[func_idx](dst, dst_stride, rc->planes[0] + ref_off, plane_stride, frac,
                                                  packed_dims, t->store);
        }

        // Copy.
        else
            venc_copy_block(dst, dst_stride, rc->planes[f265_hpel_src0[hpel_idx]] + ref_off, plane_stride,
                            width, height);
    }

    // Weighted prediction.
    else
    {
        // FIXME.
    }
}

void venc_mc_chroma_p(f265_enc_thread *t, f265_pix *dst, int dst_stride, f265_ref_ctx *rc, f265_mv mv,
                      int packed_dims, int plane_off, int comp)
{
    int plane_stride = t->plane_stride;
    int mx = mv.x, my = mv.y;
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    int ref_off = plane_off + (my>>3)*plane_stride + (mx>>3);
    int frac_x = mx&7, frac_y = my&7;
    int frac = (frac_y<<3)|frac_x, frac_idx = (!!frac_y<<1) + !!frac_x - 1;

    // Temporary. FIXME.
    void (*frac_funcs[3])(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac, int packed_dims,
                          uint8_t *spill) =
        { venc_interpol_chroma_qpel_pix_h_c, venc_interpol_chroma_qpel_pix_v_c, venc_interpol_chroma_qpel_pix_d_c };

    // No weighted prediction.
    if (!rc->weights[comp].used_flag)
    {
        // Interpolate.
        if (likely(frac))
            frac_funcs[frac_idx](dst, dst_stride, rc->planes[3+comp] + ref_off, plane_stride, frac, packed_dims,
                                 t->store);

        // Copy.
        else
            venc_copy_block(dst, dst_stride, rc->planes[3+comp] + ref_off, plane_stride, width, height);
    }

    // Weighted prediction.
    else
    {
        // FIXME.
    }
}

// Do the motion compensation for the current bi-predicted block.
void venc_mc_luma_b(f265_enc_thread *t, f265_pix *dst, int dst_stride, f265_ref_ctx *rc[2], f265_mv mv[2],
                    int packed_dims, int plane_off)
{
    int plane_stride = t->plane_stride;
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    int aligned_block_size = F265_ALIGN_VAL(width*height, 64);
    int alloc_size = 2*2*aligned_block_size;
    int16_t *ref_buf = (int16_t*)t->store;
    t->store += alloc_size;

    // Temporary. FIXME.
    void (*frac_funcs[3])(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac, int packed_dims,
                          uint8_t *spill) =
        { venc_interpol_luma_qpel_s16_h_c, venc_interpol_luma_qpel_s16_v_c, venc_interpol_luma_qpel_s16_d_c };

    // Interpolate each reference in a temporary buffer.
    for (int list = 0; list < 2; list++)
    {
        int mx = mv[list].x, my = mv[list].y;
        int ref_off = plane_off + (my>>2)*plane_stride + (mx>>2);
        int frac_x = mx&3, frac_y = my&3;
        int frac = (frac_y<<2)|frac_x, frac_idx = (!!frac_y<<1) + !!frac_x - 1;
        int16_t *ref_dst = ref_buf + aligned_block_size*list;

        // Interpolate.
        if (likely(frac))
            frac_funcs[frac_idx](ref_dst, width, rc[list]->planes[0] + ref_off, plane_stride, frac,
                                 packed_dims, t->store);

        // Scale.
        else
            venc_scale_qpel_c(ref_dst, width, rc[list]->planes[0] + ref_off, plane_stride, packed_dims);
    }

    // Unweighted bi-prediction.
    if (!rc[0]->weights->used_flag && !rc[1]->weights->used_flag)
        venc_avg_pix_s16_c(dst, dst_stride, ref_buf, ref_buf + aligned_block_size, width, packed_dims);

    // Weighted bi-prediction.
    else
    {
        // FIXME.
    }

    t->store -= alloc_size;
}

void venc_mc_chroma_b(f265_enc_thread *t, f265_pix *dst, int dst_stride, f265_ref_ctx *rc[2], f265_mv mv[2],
                      int packed_dims, int plane_off, int comp)
{
    int plane_stride = t->plane_stride;
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff;
    int aligned_block_size = F265_ALIGN_VAL(width*height, 64);
    int alloc_size = 2*2*aligned_block_size;
    int16_t *ref_buf = (int16_t*)t->store;
    t->store += alloc_size;

    // Temporary. FIXME.
    void (*frac_funcs[3])(int16_t *dst, int dst_stride, f265_pix *src, int src_stride, int frac, int packed_dims,
                          uint8_t *spill) =
        { venc_interpol_chroma_qpel_s16_h_c, venc_interpol_chroma_qpel_s16_v_c, venc_interpol_chroma_qpel_s16_d_c };

    // Interpolate each reference in a temporary buffer.
    for (int list = 0; list < 2; list++)
    {
        int mx = mv[list].x, my = mv[list].y;
        int ref_off = plane_off + (my>>3)*plane_stride + (mx>>3);
        int frac_x = mx&7, frac_y = my&7;
        int frac = (frac_y<<3)|frac_x, frac_idx = (!!frac_y<<1) + !!frac_x - 1;
        int16_t *ref_dst = ref_buf + aligned_block_size*list;

        // Interpolate.
        if (likely(frac))
            frac_funcs[frac_idx](ref_dst, width, rc[list]->planes[3+comp] + ref_off, plane_stride, frac,
                                 packed_dims, t->store);

        // Scale.
        else
            venc_scale_qpel_c(ref_dst, width, rc[list]->planes[3+comp] + ref_off, plane_stride, packed_dims);
    }

    // Unweighted bi-prediction.
    if (!rc[0]->weights[comp].used_flag && !rc[1]->weights[comp].used_flag)
        venc_avg_pix_s16_c(dst, dst_stride, ref_buf, ref_buf + aligned_block_size, width, packed_dims);

    // Weighted bi-prediction.
    else
    {
        // FIXME.
    }

    t->store -= alloc_size;
}


///////////////////////

// Neighbour motion vector order.
typedef enum f265_inter_neighbour_direction
{
    F265_ND_TOP_LEFT = 0,
    F265_ND_TOP,
    F265_ND_TOP_RIGHT,
    F265_ND_LEFT,
    F265_ND_BOTTOM_LEFT,
    F265_ND_COLLOCATED,
} f265_inter_neighbour_direction;


static int get_tmv_list(f265_enc_thread *t);
static int get_tmv_ref(f265_enc_thread *t);
static int get_tmv_force_list_flag(f265_enc_thread *t);
static f265_ref_ctx* get_tmv_ref_ctx(f265_enc_thread *t);
static void venc_save_temporal_ref_poc(f265_enc_thread *t);
static void venc_find_temp_ref_idx(f265_enc_thread *t, f265_frame *f);
static void venc_calculate_temporal_scaling(f265_enc_thread *t);
static void venc_prepare_temporal_ref(f265_enc_thread *t, f265_frame *f);
static void venc_get_inter_neighbour_mv_at(f265_enc_thread *t, int cur_enc_idx, int x, int y,
                                           f265_inter_neighbour_mv *neighbour);
static int venc_get_collocated_mv(f265_enc_thread *t, int x, int y, f265_inter_neighbour_mv *neighbour);
static f265_mv venc_scale_mv(f265_mv mv, int tb, int tx);
static int venc_load_direct_pmv(f265_enc_thread *t, f265_mv *pmv, f265_inter_neighbour_mv *neighbours,
                               uint32_t reference_list, f265_frame *ref_frame);
static int venc_load_indirect_pmv(f265_enc_thread *t, f265_mv *pmv, f265_inter_neighbour_mv *neighbours,
                                  uint32_t reference_list, uint32_t ref_idx);
static inline int venc_is_same_inter_neighbour(f265_inter_neighbour_mv *cand1, f265_inter_neighbour_mv *cand2);


// Return the collocated reference list.
static int get_tmv_list(f265_enc_thread *t)
{
    return F265_GET_FLAG(t->temp_ref_info, 1 << 4);
}

// Return the collocated reference index.
static int get_tmv_ref(f265_enc_thread *t)
{
    return t->temp_ref_info & 0xf;
}

// Return the force candidate list flag.
static int get_tmv_force_list_flag(f265_enc_thread *t)
{
    return F265_GET_FLAG(t->temp_ref_info, 1 << 5);
}

// Return the reference context associated to the reference index.
static f265_ref_ctx* get_tmv_ref_ctx(f265_enc_thread *t)
{
    return t->ref_ctx[get_tmv_list(t)] + get_tmv_ref(t);
}

// Save the reference POC in the frame. Required for temporal MV scaling.
static void venc_save_temporal_ref_poc(f265_enc_thread *t)
{
    for (int l = 0; l < 2; l++)
    {
        for (int i = 0; i < 16; i++)
        {
            if (i < t->src_frame->nb_ref_idx[l])
                t->src_frame->ref_poc[l][i] = t->ref_ctx[l][i].frame->abs_poc;
            else
                t->src_frame->ref_poc[l][i] = -1;
        }
    }
}

// This function imitates the HM way of deciding the temporal candidate list.
// Note that HM always sets collocated_ref_idx to 0. Only the lists are searched.
// This code assumes that long term references are unsupported.
static void venc_find_temp_ref_idx(f265_enc_thread *t, f265_frame *f)
{
    // Choose list 0 and reference index 0 by default.
    f->temp_ref_info = 0;

    if (f->frame_type != F265_FRAME_B) return;

    int gop_size = t->enc->gd.hm_gop_size;
    if (gop_size)
    {
        f265_hm_gop_entry *gop_entry = f->gop_entry;

        // Search for the index of closest reference in L0 (left frame),
        // and search for the index of the closest frame in L1 (right frame).
        // The index of the reference is the index set the configuration file and
        // computed as POC_off(curr_frame) - POC_off(ref_frame).
        // The POC offset set in the configuration file is the wrap-around POC
        // of a frame using gop_size as period.

        int left=1, right=-1;
        for (int i = 0; i< gop_entry->nb_refs; i++)
        {
            int ref = gop_entry->refs[i];
            if (ref > 0 && (ref < right || right == -1))
            {
                right = ref;
            }
            else if (ref < 0 && (ref > left || left == 1))
            {
                left = ref;
            }
        }

        // Compute the POC offsets of the left and right frames using the
        // POC offset of the current frame.
        if (right > -1)
        {
            right = right + gop_entry->poc_off - 1;
        }
        if (left < 1)
        {
            left = left + gop_entry->poc_off - 1;
            while (left < 0)
            {
                left += gop_size;
            }
        }

        // Search in the gop structure for the QP offsets of the left and right frames.
        int left_qp = 0, right_qp = 0;
        for (int i = 0; i < gop_size; i++)
        {
            gop_entry = t->enc->gd.hm_gop + i;
            if (gop_entry->poc_off == (left % gop_size) + 1)
            {
                left_qp = gop_entry->qp_offset;
            }
            if (gop_entry->poc_off == (right % gop_size) + 1)
            {
                right_qp = gop_entry->qp_offset;
            }
        }

        // Set the reference list (the collocated_from_l0_flag syntax element).
        f->temp_ref_info = 1 << 4;
        if (right > -1 && right_qp < left_qp) f->temp_ref_info = 0;
    }

    // If the dbp contains a picture that comes after the current picture, the spec requires
    // that the TMV is loaded from the opposite list if the reference list is available.
    if (f->dpb_pos > 0) f->temp_ref_info |= 1 << 5;

    t->temp_ref_info = f->temp_ref_info;
}

// Calculate the temporal scaling factors.
static void venc_calculate_temporal_scaling(f265_enc_thread *t)
{
    f265_frame *temp_frame = get_tmv_ref_ctx(t)->frame;
    int64_t col_poc = temp_frame->abs_poc;

    for (int list = 0; list < 2; list++)
    {
        for (int i = 0; i < temp_frame->nb_ref_idx[list]; i++)
        {
            int delta_poc = F265_CLAMP(col_poc - temp_frame->ref_poc[list][i], -128, 127);
            assert(delta_poc != 0);
            t->temp_dsf[i][list] = (0x4000 + (F265_ABS(delta_poc)>>1)) / delta_poc;
        }
    }
}

// Setup all info required for temporal prediction in the frame.
static void venc_prepare_temporal_ref(f265_enc_thread *t, f265_frame *f)
{
    // Clear the temporal reference.
    memset(t->temp_ref_idx, -1, 5*5*2);

    // Only do this if we support temporal MV.
    if (!(t->enc->gd.eflags & F265_PF_TMV)) return;

    // Find the temporal reference list and index.
    venc_find_temp_ref_idx(t, f);

    // Save the ref POC for later use.
    venc_save_temporal_ref_poc(t);

    // Calculate temporal scaling delta POC.
    venc_calculate_temporal_scaling(t);
}

// Calculate the scaling offset for spatial and temporal PMVs. Should be called
// every time the reference frame or temporal reference changes.
void venc_frame_set_up_inter_pred(f265_enc_thread *t)
{
    if (t->src_frame->frame_type == F265_FRAME_I) return;

    int curr_poc = t->src_frame->abs_poc;

    // Calculate spatial scaling delta POC.
    for (int list = 0; list < t->nb_lists; list++)
    {
        for (int i = 0; i < t->nb_ref_idx[list]; i++)
        {
            int delta_poc = F265_CLAMP(curr_poc - t->ref_ctx[list][i].frame->abs_poc, -128, 127);
            t->spatial_dsf[i][list][0] = delta_poc;
            t->spatial_dsf[i][list][1] = (0x4000 + (F265_ABS(delta_poc)>>1)) / delta_poc;
        }
    }

    venc_prepare_temporal_ref(t, t->src_frame);
}

// Load the temporal PMVs. Should be called at every CTB.
void venc_ctb_set_up_inter_pred(f265_enc_thread *t)
{
    // Skip this function if we don't use TMV or this is intra.
    if (!(t->enc->gd.eflags & F265_PF_TMV) || !t->nb_lists) return;

    int b4_stride = t->enc->gd.pix_dim[0] >> 2;
    int b4_x = (t->ctb_off[0] >> 2);
    int b4_y = (t->ctb_off[1] >> 2);
    int b4_off = b4_y * b4_stride + b4_x;

    // Find the number of 16x16 blocks used for temporal prediction in the CTB
    // in each direction within the frame borders.
    int max[2];
    for (int i = 0; i < 2; i++)
    {
        // Since we may use a TMV from a neighbour on the right of the current
        // CTB, allocate an additional block horizontally.
        max[i] = ((1 << t->cb[0].lg_bs) >> 4) + !i;

        // Clip to the frame border.
        int delta_past_border = (t->ctb_off[i] + (1 << t->cb[0].lg_bs)) - t->enc->gd.pix_dim[i];
        if (unlikely(delta_past_border >= 0))
            max[i] -= (delta_past_border >> 4) + !i;
    }

    f265_frame *temp_frame = get_tmv_ref_ctx(t)->frame;
    f265_mv (*mv)[2] = temp_frame->mv + b4_off;
    int8_t (*ref)[2] = temp_frame->ref_idx + b4_off;

    // Mark all temporal references as unavailable.
    memset(t->temp_ref_idx, -1, sizeof(t->temp_ref_idx));

    // Load available TMV and reference.
    for (int row = 0; row < max[1]; row++)
    {
        for (int col = 0; col < max[0]; col++)
        {
            t->temp_mv[(row * 5) + col][0] = mv[col << 2][0];
            t->temp_mv[(row * 5) + col][1] = mv[col << 2][1];
            t->temp_ref_idx[(row * 5) + col][0] = ref[col << 2][0];
            t->temp_ref_idx[(row * 5) + col][1] = ref[col << 2][1];
        }

        // The scaling factors cancel out here.
        mv += t->enc->gd.pix_dim[0];
        ref += t->enc->gd.pix_dim[0];
    }
}

// Cache a spatial neighbour at the (X,Y) 4x4 block offset from the CTB origin.
static void venc_get_inter_neighbour_mv_at(f265_enc_thread *t, int cur_enc_idx, int x, int y,
                                           f265_inter_neighbour_mv *neighbour)
{
    int cb_idx;
    // Get the CB and partition index.
    int part_idx;
    venc_lookup_pmap(t, x, y, &cb_idx, &part_idx);

    f265_cb *cb = t->cb + cb_idx;
    if ((cb->flags & (F265_CB_INTRA | F265_CB_PRESENT)) != F265_CB_PRESENT || cb->enc_idx > cur_enc_idx)
    {
        // Mark the PMV as unavailable.
        for (int i = 0; i < 2; i++)
            neighbour->ref_idx[i] = -1;
        return;
    }

    for (int i = 0; i < 2; i++)
    {
        // Get the MV.
        neighbour->mv[i] = cb->mv[part_idx][i];

        // Get the reference idx.
        neighbour->ref_idx[i] = cb->ref_idx[part_idx][i];
    }
}

// Cache the collocated MV. Same arguments as above. Return true if the motion
// vector is available.
static int venc_get_collocated_mv(f265_enc_thread *t, int x, int y, f265_inter_neighbour_mv *neighbour)
{
    // Import the collocated MV.
    f265_mv *mv = t->temp_mv[x + (y * 5)] + 0;
    int8_t *ref = t->temp_ref_idx[x + (y * 5)] + 0;

    // Since special rules apply when loading the PMV, we need to load the list one at a time.
    for (int list = 0; list < 2; list++)
    {
        // By default, we want to load from the list we are setting.
        int loaded_list = list;

        // If a reference uses a POC greater than the current POC, force both lists to use the same MV.
        if (get_tmv_force_list_flag(t))
            loaded_list = !get_tmv_list(t);

        // Load the required MV.
        neighbour->ref_idx[list] = ref[loaded_list];
        neighbour->mv[list] = mv[loaded_list];

        // If the MV is invalid, load the MV from the other list.
        if (neighbour->ref_idx[list] == -1)
        {
            loaded_list = !loaded_list;
            neighbour->ref_idx[list] = ref[loaded_list];
            neighbour->mv[list] = mv[loaded_list];
        }

        // Set the reference list being used.
        neighbour->ref_idx[list] |= loaded_list<<4;
    }

    // If both reference index are negatives, the MV is invalid.
    return neighbour->unified_ref != -1;
}

// Cache all neighbours MVs. Should be done once per partition.
void venc_get_neighbour_mvs(f265_enc_thread *t, f265_cb *cb, uint32_t partition_idx,
                            f265_inter_neighbour_mv *neighbours)
{
    // Get the partition sizes and offsets in 4x4 blocks.
    int part_size[2], part_off[2], mv_off[2];
    venc_get_part_loc4(part_size, part_off, cb->inter_part, partition_idx, 1 << cb->lg_bs);

    for (int i = 0; i < 2; i++)
        mv_off[i] = part_off[i] + (cb->cb_off[i] >> 2);

    // Neighbour offset table. The index is the neighbour index. The value is a
    // bitfield. There are four bits per neighbour, 2 bits for X and 2 bits for
    // Y. The first bit is the block offset (0 or -1). The second bit is true if
    // the partition size must be added.
    // 1001 1101 0110 0111 0101 => 0x9d675;
    int neighbour_offset_bitfield = 0x9d675;

    // Set the temporal neighbour to unavailable by default.
    neighbours[F265_ND_COLLOCATED].unified_ref = -1;

    // Load the spatial MVs.
    int cur_enc_idx = cb->enc_idx;
    for (int i = 0; i < 5; i++, neighbour_offset_bitfield >>=4)
    {
        int offset_by_one_x = F265_GET_FLAG(neighbour_offset_bitfield, 1<<0);
        int part_size_x_mask = F265_GET_FLAG(neighbour_offset_bitfield, 1<<1);
        int offset_by_one_y = F265_GET_FLAG(neighbour_offset_bitfield, 1<<2);
        int part_size_y_mask = F265_GET_FLAG(neighbour_offset_bitfield, 1<<3);

        venc_get_inter_neighbour_mv_at(t,
                                       cur_enc_idx,
                                       mv_off[0] + part_size[0] * part_size_x_mask - offset_by_one_x,
                                       mv_off[1] + part_size[1] * part_size_y_mask - offset_by_one_y,
                                       neighbours + i);
    }

    // Get the temporal MV, if available.
    if (t->enc->gd.eflags & F265_PF_TMV)
    {
        // There is 2 possible neighbours of interest: bottom right and middle.
        // We use the bottom right MV if it's available and has a valid
        // candidate. If it is unavailable, we use the middle MV.

        // Check if the bottom right position is inside the frame.
        // FIXME: refactor this code.
        int bottom_right_inside_flag = 1;
        for (int i = 0; i < 2; i++)
        {
            if (t->ctb_off[i] + ((mv_off[i] + part_size[i])<<2) >= t->enc->gd.pix_dim[i])
            {
                bottom_right_inside_flag = 0;
                break;
            }
        }

        // Get the offset of the bottom right neighbour in 16x16 blocks.
        int tmv_off[2];
        for (int i = 0; i < 2; i++)
            tmv_off[i] = (mv_off[i] + part_size[i]) >> 2;

        // Fallback to middle if bottom-right is not available.
        if (!bottom_right_inside_flag ||
            !venc_get_collocated_mv(t, tmv_off[0], tmv_off[1], neighbours + F265_ND_COLLOCATED))
        {
            // Get the 16x16 aligned offset centered on the current partition.
            for (int i = 0; i < 2; i++)
                tmv_off[i] = (mv_off[i] + (part_size[i] >> 1)) >> 2;

            venc_get_collocated_mv(t, tmv_off[0], tmv_off[1], neighbours + F265_ND_COLLOCATED);
        }
    }
}

// Scale a MV by the POC distance.
static f265_mv venc_scale_mv(f265_mv mv, int tb, int tx)
{
    int scale = F265_CLAMP((tb * tx + 32) >> 6, -4096, 4095);
    for (int i = 0; i < 2; i++)
    {
        // Scale while handling the negative rounding.
        int s = (scale * mv.v[i] + 127 + (scale * mv.v[i] < 0)) >> 8;
        mv.v[i] = F265_CLAMP(s, -32768, 32767);
    }

    return mv;
}

// Load the direct spatial MV if available. Return true if a valid MV was found.
static int venc_load_direct_pmv(f265_enc_thread *t, f265_mv *pmv, f265_inter_neighbour_mv *neighbours,
                                uint32_t reference_list, f265_frame *ref_frame)
{
    for (int i = 0; i < 2; i++, reference_list = !reference_list)
    {
        if (neighbours->ref_idx[reference_list] != -1 &&
            t->ref_ctx[reference_list][neighbours->ref_idx[reference_list]].frame == ref_frame)
        {
            *pmv = neighbours->mv[reference_list];
            return 1;
        }
    }

    return 0;
}

// Load the indirect spatial MV if available. An indirect MV has to be scaled.
// Return true if a valid MV was found.
static int venc_load_indirect_pmv(f265_enc_thread *t, f265_mv *pmv, f265_inter_neighbour_mv *neighbours,
                                  uint32_t reference_list, uint32_t ref_idx)
{
    int ref_mv_ref_list = reference_list;
    for (int i = 0; i < 2; i++, ref_mv_ref_list = !ref_mv_ref_list)
    {
        if (neighbours->ref_idx[ref_mv_ref_list] != -1)
        {
            *pmv = venc_scale_mv(neighbours->mv[ref_mv_ref_list], t->spatial_dsf[ref_idx][reference_list][0],
                                 t->spatial_dsf[neighbours->ref_idx[ref_mv_ref_list]][ref_mv_ref_list][1]);
            return 1;
        }
    }

    return 0;
}

// Create the PMV list. Requires a f265_inter_neighbour_mv array populated with
// venc_get_neighbour_mvs(). Should be called for each reference index.
void venc_get_pmv(f265_enc_thread *t, int part_idx, f265_inter_neighbour_mv *neighbours,
                  uint32_t ref_idx, uint32_t reference_list, f265_mv *pmv)
{
    // FIXME: it's the reverse declaration order, can be optimized.
    int8_t left_neigh_offset[2] = { F265_ND_BOTTOM_LEFT, F265_ND_LEFT };
    int8_t top_neigh_offset[3] = { F265_ND_TOP_RIGHT, F265_ND_TOP, F265_ND_TOP_LEFT };

    // Set the first candidate to 0 in case nothing matches.
    pmv[0].p = 0;

    // Get the referenced frame.
    f265_frame *ref_frame = t->ref_ctx[reference_list][ref_idx].frame;

    // Number of motion vectors found.
    uint32_t index;

    // Flag indicating if the left neighbour was found.
    int left_found;

    // Left candidate.
    do
    {
        index = 1; // Suppose a candidate will be found.

        for (int neigh = 0; neigh < 2; neigh++)
            if (venc_load_direct_pmv(t, pmv, neighbours + left_neigh_offset[neigh], reference_list, ref_frame))
                goto LEFT_NEIGHBOUR_FOUND;

        for (int neigh = 0; neigh < 2; neigh++)
            if (venc_load_indirect_pmv(t, pmv, neighbours + left_neigh_offset[neigh], reference_list, ref_idx))
                goto LEFT_NEIGHBOUR_FOUND;

        index = 0; // No valid candidate was found.
    } while (0);
    LEFT_NEIGHBOUR_FOUND:
    left_found = index;

    // Top candidate.
    do
    {
        f265_mv *curr_pmv = pmv + index;
        index++;

        for (int neigh = 0; neigh < 3; neigh++)
            if (venc_load_direct_pmv(t, curr_pmv, neighbours + top_neigh_offset[neigh], reference_list, ref_frame))
                goto DIRECT_TOP_NEIGHBOUR_FOUND;

        index--;
    } while (0);
    DIRECT_TOP_NEIGHBOUR_FOUND:

    // If the left neighbour is absent, search for another PMV in the top neighbour.
    if (!left_found)
    {
        f265_mv *curr_pmv = pmv + index;
        index++;

        for (int neigh = 0; neigh < 3; neigh++)
            if (venc_load_indirect_pmv(t, curr_pmv, neighbours + top_neigh_offset[neigh], reference_list, ref_idx))
                goto INDIRECT_TOP_NEIGHBOUR_FOUND;

        index--;
    }
    INDIRECT_TOP_NEIGHBOUR_FOUND:

    // If two motion vectors were found, test for duplicate.
    if (index == 2) index -= pmv[0].p == pmv[1].p;

    // Get the collocated candidate. At this step, either the first candidate
    // was found or its value is 0.
    if (index < 2)
    {
        // Set the second candidate to 0 in case it doesn't get loaded.
        pmv[1].p = 0;
        int temp_ref_idx = neighbours[F265_ND_COLLOCATED].ref_idx[reference_list];
        if (temp_ref_idx != -1)
            pmv[index] =
                venc_scale_mv(neighbours[F265_ND_COLLOCATED].mv[reference_list],
                              t->spatial_dsf[ref_idx][reference_list][0],
                              t->temp_dsf[temp_ref_idx&0xf][temp_ref_idx >> 4]);
    }
}

// Test if two neighbours are equal. FIXME. Replace this by memcmp (make sure
// all candidates are initialized with zero MVs when a list is unused to avoid
// Valgrind errors and mismatches).
static int venc_is_same_inter_neighbour(f265_inter_neighbour_mv *cand1, f265_inter_neighbour_mv *cand2)
{
    int same = 1;
    for (int list_idx = 0; list_idx < 2; list_idx++)
    {
        same &= cand1->ref_idx[list_idx] == cand2->ref_idx[list_idx];
        same &= (cand1->ref_idx[list_idx] == -1) |
                (cand1->mv[list_idx].p == cand2->mv[list_idx].p);
    }
    return same;
}

// Generate the required number of merge candidates. This function does not
// handle the 4x8/8x4 case.
static void venc_gen_merge_candidate(f265_enc_thread *t, f265_cb *cb, uint32_t partition_idx,
                                     f265_inter_neighbour_mv *neighbours, f265_inter_neighbour_mv *merge_candidate)
{
    f265_inter_neighbour_mv unsplit_neighbours[6];
    int partition_type = cb->inter_part;

    // Test the single merge candidate list.
    if ((t->enc->gd.parallel_merge_level != 2) & (cb->lg_bs == 3) & !!partition_type)
    {
        // The merge candidate must be generated as if it is an unsplit CB.
        partition_idx = 0;
        partition_type = F265_PART_UN;

        // Load the unsplit neighbours.
        int real_inter_part = cb->inter_part;
        cb->inter_part = F265_PART_UN;
        venc_get_neighbour_mvs(t, cb, partition_idx, unsplit_neighbours);
        cb->inter_part = real_inter_part;
        neighbours = unsplit_neighbours;
    }

    // Setup the test order.
    uint8_t spatial_order[5] = { F265_ND_LEFT, F265_ND_TOP, F265_ND_TOP_RIGHT, F265_ND_BOTTOM_LEFT, F265_ND_TOP_LEFT };

    // B frame flag.
    int b_flag = t->src_frame->frame_type == F265_FRAME_B;

    // Number of merge candidates required.
    int req_cand_num = t->enc->gd.merge_cand;

    // Get the partition size and offset.
    int part_size[2]; int part_off[2];
    venc_get_part_loc(part_size, part_off, partition_type, partition_idx, 1<<cb->lg_bs);

    // Get the position of the partition top-left pixel.
    int pos[2];
    for (int i = 0; i < 2; ++i)
        pos[i] = t->ctb_off[i] + cb->cb_off[i] + part_off[i];

    // Test if each neighbour is aligned on the parallel_merge_level boundary.
    // If it is not, then the neighbour is in the same parallel block as the
    // current block and cannot be used for prediction.
    int parallel_merge_mask = (1 << t->enc->gd.parallel_merge_level) - 1;
    int is_left_aligned = (pos[0] & parallel_merge_mask) == 0;
    int is_top_aligned = (pos[1] & parallel_merge_mask) == 0;
    int is_right_aligned = ((pos[0] + part_size[0]) & parallel_merge_mask) == 0;
    int is_bottom_aligned = ((pos[1] + part_size[1]) & parallel_merge_mask) == 0;

    // We do not support the inter HV mode. Hence, there are two cases. Either
    // we're predicting for the first partition or the second partition. The
    // first partition uses the neighbour blocks unconditionally. The second
    // partition uses the neighbour blocks unless they correspond to the first
    // partition. We assume that we're predicting for the second prediction and
    // mask with the following value, which accounts for the partition index.
    int ignore_part_test = !partition_idx;

    // Check if the left/top partition is in the CB.
    // Test for H1, H2 and H3.
    int left_block_not_in_cb_flag = ((0x3b >> partition_type)&1) | ignore_part_test;

    // Test for V1, V2 and V3.
    int top_block_not_in_cb_flag = ((0xcd >> partition_type)&1) | ignore_part_test;

    // Bitfield indicating the neighbour availability requirements.
    // 0- Check left_block_not_in_cb_flag.
    // 1- Check top_block_not_in_cb_flag.
    // 2- Check is_left_aligned.
    // 3- Check is_top_aligned.
    // 4- Check is_right_aligned.
    // 5- Check is_bottom_aligned.
    // 6- Check for equality with the left candidate.
    // 7- Check for equality with the top candidate.
    // 8-9 Index in ref_avail[] below.
    uint16_t neighbour_avail_check_mask[5] = { 0x005, 0x14a, 0x298, 0x264, 0x2cc };

    // Merge all flags together.
    int flag_mask = (left_block_not_in_cb_flag  << 0) |
                    (top_block_not_in_cb_flag   << 1) |
                    (is_left_aligned            << 2) |
                    (is_top_aligned             << 3) |
                    (is_right_aligned           << 4) |
                    (is_bottom_aligned          << 5);

    // Indicates if the left/top neighbours are available for comparison.
    // Index: 0- left, 1- top, 2- trash.
    int ref_avail[3] = { 0, 0 };

    // Count of candidates found.
    int cand_idx = 0;

    // Pass all spatial neighbours.
    for (int i = 0; i < 5; i++)
    {
        int ni = spatial_order[i];
        uint32_t mask = neighbour_avail_check_mask[i];

        // Check if the MV is valid.
        int available = neighbours[ni].unified_ref != -1;

        // Test if the neighbour is in the same CB.
        available &= ((flag_mask | (~mask)) & 0x3) == 0x3;

        // Test if the neighbour is in another parallel merge block. The neighbour is valid if
        // it's in a different block in any relevant direction.
        available &= !!((flag_mask & mask) & 0x3c);

        // Save the availability of the current neighbour at the appropriate
        // location.
        ref_avail[mask>>8] = available;

        // Test if the current MV is a disallowed duplicate.
        if (!!(mask & 64) & ref_avail[0])
            available &= !venc_is_same_inter_neighbour(neighbours + F265_ND_LEFT, neighbours + ni);

        if (!!(mask & 128) & ref_avail[1])
            available &= !venc_is_same_inter_neighbour(neighbours + F265_ND_TOP, neighbours + ni);

        // Keep the candidate.
        if (available)
        {
            memcpy(merge_candidate + cand_idx, neighbours + ni, sizeof(f265_inter_neighbour_mv));

            // Return if we reached the required number of candidates.
            if (++cand_idx == req_cand_num) return;

            // Only 4 candidates may come from the spatial neighbours.
            if (cand_idx == 4) break;
        }
    }

    // Load the temporal MV candidate.
    if (neighbours[F265_ND_COLLOCATED].unified_ref != -1)
    {
        f265_inter_neighbour_mv *mrg_cand = merge_candidate + cand_idx;

        // Set the second list as invalid. Will be overwritten if it's a B
        // frame.
        mrg_cand->ref_idx[1] = -1;

        for (int list = 0; list < (1+b_flag); list++)
        {
            int ref_list = neighbours[F265_ND_COLLOCATED].ref_idx[list] >> 4;
            int ref_idx = neighbours[F265_ND_COLLOCATED].ref_idx[list] & 0xf;

            // Scale the MV. The temporal candidate uses reference index 0.
            mrg_cand->mv[list] = venc_scale_mv(neighbours[F265_ND_COLLOCATED].mv[list],
                                               t->spatial_dsf[0][list][0],
                                               t->temp_dsf[ref_idx][ref_list]);
            mrg_cand->ref_idx[list] = 0;
        }

        if (++cand_idx == req_cand_num) return;
    }

    // Create new candidates by mixing the previous ones.
    if (b_flag)
    {
        uint32_t priority_list0 = 0xedc984;  // { 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3 };
        uint32_t priority_list1 = 0xb73621;  // { 1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2 };

        int max_idx = cand_idx * (cand_idx - 1);
        for (int idx = 0; idx < max_idx; idx++, priority_list0 >>= 2, priority_list1 >>= 2)
        {
            f265_inter_neighbour_mv *l0_cand = merge_candidate + (priority_list0 & 3);
            f265_inter_neighbour_mv *l1_cand = merge_candidate + (priority_list1 & 3);

            // Load the reference indices of the candidates.
            int l0_ref_idx = l0_cand->ref_idx[0];
            int l1_ref_idx = l1_cand->ref_idx[1];

            // If a reference is invalid, skip the current mix.
            if ((l0_ref_idx == -1) | (l1_ref_idx == -1))
                continue;

            // Make sure both references do not point to the same place in the same frame.
            if ((t->ref_ctx[0][l0_ref_idx].frame != t->ref_ctx[1][l1_ref_idx].frame) |
                (l0_cand->mv[0].p != l1_cand->mv[1].p))
            {
                merge_candidate[cand_idx].mv[0] = l0_cand->mv[0];
                merge_candidate[cand_idx].mv[1] = l1_cand->mv[1];
                merge_candidate[cand_idx].ref_idx[0] = l0_ref_idx;
                merge_candidate[cand_idx].ref_idx[1] = l1_ref_idx;
                if (++cand_idx == req_cand_num) return;
            }
        }
    }

    // Generate null motion vector with different reference indices.
    int max_idx = b_flag ? F265_MIN(t->src_frame->nb_ref_idx[0], t->src_frame->nb_ref_idx[1]) :
                           t->src_frame->nb_ref_idx[0];
    for (int i = 0; i < max_idx; i++)
    {
        merge_candidate[cand_idx].mv[0].p = 0;
        merge_candidate[cand_idx].mv[1].p = 0;
        merge_candidate[cand_idx].ref_idx[0] = i;
        merge_candidate[cand_idx].ref_idx[1] = b_flag ? i : -1;
        if (++cand_idx == req_cand_num) return;
    }

    // Generate null motion vector with reference index 0.
    while (cand_idx < req_cand_num)
    {
        merge_candidate[cand_idx].mv[0].p = 0;
        merge_candidate[cand_idx].mv[1].p = 0;
        merge_candidate[cand_idx].ref_idx[0] = 0;
        merge_candidate[cand_idx].ref_idx[1] = b_flag ? 0 : -1;
        cand_idx++;
    }
}

// Create the merge candidate list. Requires a f265_inter_neighbours_mv array
// populated with venc_get_neighbour_mvs().
void venc_get_merge_candidate(f265_enc_thread *t, f265_cb *cb, uint32_t partition_idx,
                              f265_inter_neighbour_mv *neighbours, f265_inter_neighbour_mv *merge_candidate)
{
    venc_gen_merge_candidate(t, cb, partition_idx, neighbours, merge_candidate);

    // Prevent bi-prediction for 8x4/4x8 blocks.
    if (cb->lg_bs == 3 && cb->inter_part != F265_PART_UN)
    {
        int nb_cand = t->enc->gd.merge_cand;
        for (int i = 0; i < nb_cand; i++)
        {
            f265_inter_neighbour_mv *cand = merge_candidate + i;
            if (cand->ref_idx[0] != -1 && cand->ref_idx[1] != -1)
            {
                cand->ref_idx[1] = -1;
                cand->mv[1].p = 0;
            }
        }
    }
}

