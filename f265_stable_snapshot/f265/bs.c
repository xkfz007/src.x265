// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Bitstream encoding functions, excluding the CTB data.
#include "f265/enc.h"

#ifdef VAN_TRACE_SYNTAX
// True if the syntax element values are printed.
int venc_trace_syntax_flag;

// Print a CABAC bin of a syntax element as needed.
static void venc_trace_bin(int type, int ctx, int bin)
{
  if (!venc_trace_syntax_flag) return;
  if (type == 0 || type == 2) printf("    %d/%02x\n", bin, ctx);
  else printf("    %d\n", bin);
}
#endif

// Write the current byte in the bitstream without modifying bits_left.
static void venc_vlc_write_byte(f265_vlc_bs *vbs, int byte)
{
    // Direct write.
    if (byte > 3)
    {
        *vbs->cur++ = byte;
        vbs->trail_zeroes = 0;
        return;
    }

    // Add an emulation prevention code.
    if (vbs->trail_zeroes == 2)
    {
        *vbs->cur++ = 3;
        vbs->trail_zeroes = 0;
    }

    // Write the current byte.
    *vbs->cur++ = byte;
    vbs->trail_zeroes += !byte;
}

// Add a single bit (flag) to the bitstream.
void venc_vlc_put_bit(f265_vlc_bs *vbs, uint32_t code)
{
    assert(code <= 1);

    // Append the bit to the current byte without completing it.
    if (vbs->bits_left > 1)
    {
        vbs->bits_left--;
        vbs->cur[0] |= code<<vbs->bits_left;
        return;
    }

    // Complete and write the current byte.
    venc_vlc_write_byte(vbs, vbs->cur[0]|code);
    vbs->cur[0] = 0;
    vbs->bits_left = 8;
}

// Add up to 64 bits to the bitstream.
void venc_vlc_put_bits(f265_vlc_bs *vbs, uint64_t code, int length)
{
    // Append the bits to the current byte without completing it.
    if (length < vbs->bits_left)
    {
        vbs->bits_left -= length;
        vbs->cur[0] |= code<<vbs->bits_left;
        return;
    }

    // Complete and write the current byte.
    length -= vbs->bits_left;
    venc_vlc_write_byte(vbs, vbs->cur[0]|(code>>length));

    // Write the next complete bytes.
    while (length >= 8)
    {
        length -= 8;
        venc_vlc_write_byte(vbs, (code>>length)&0xff);
    }

    // Begin the last partial byte.
    vbs->bits_left = 8-length;
    *vbs->cur = (code&((1<<length)-1))<<vbs->bits_left;
}

// Add a VLC code to the bitstream.
void venc_vlc_put_vlc(f265_vlc_bs *vbs, uint32_t code)
{
    // If the code is small enough, use the precomputed table.
    if (likely(code < 64))
    {
        venc_vlc_put_bits(vbs, f265_vlc_table[code][0], f265_vlc_table[code][1]);
    }

    // Compute the value.
    else
    {
        int32_t i = 0;
        uint32_t nn = code + 1;

        // Find the number of bits needed to encode the value.
        while (nn >>= 1) i++;
        i += i + 1;

        venc_vlc_put_bits(vbs, code+1, i);
    }
}

// Write the stop bit and align the bitstream on the next byte boundary.
void venc_vlc_flush(f265_vlc_bs *vbs)
{
    venc_vlc_put_bits(vbs, 1<<(vbs->bits_left-1), vbs->bits_left);
}

// Convert a signed VLC code to an unsigned VLC code.
#define F265_SIGNED_VLC_CODE(code) (2*F265_ABS(code) - ((code) > 0))

// VLC writing support. Uncomment to debug.
#if 0
#define VENC_ANNOUNCE_VLC(name) printf("\n"name":\n");
#define VENC_PUT_FLAG(BS, CODE, NAME) venc_vlc_put_bit_debug(BS, CODE, NAME)
#define VENC_PUT_BITS(BS, CODE, LENGTH, NAME) venc_vlc_put_bits_debug(BS, CODE, LENGTH, NAME)
#define VENC_PUT_UE_V(BS, CODE, NAME) venc_vlc_put_vlc_debug(BS, CODE, 0, NAME)
#define VENC_PUT_SE_V(BS, CODE, NAME) venc_vlc_put_vlc_debug(BS, CODE, 1, NAME)
static void venc_vlc_put_bit_debug(f265_vlc_bs *vbs, uint32_t code, char *name)
{
    printf("%-50s flag  %d\n", name, code);
    venc_vlc_put_bit(vbs, code);
}
static void venc_vlc_put_bits_debug(f265_vlc_bs *vbs, uint64_t code, int length, char *name)
{
    printf("%-50s u(%02d) ", name, length);
    for (int i = 0; i < length; i++) printf("%d", (int)(code>>(length-1-i)&1));
    printf("\n");
    venc_vlc_put_bits(vbs, code, length);
}
static void venc_vlc_put_vlc_debug(f265_vlc_bs *vbs, int32_t code, int se_flag, char *name)
{
    printf("%-50s %s(v) %d\n", name,  se_flag ? "se" : "ue", code);
    if (se_flag) code = F265_SIGNED_VLC_CODE(code);
    venc_vlc_put_vlc(vbs, code);
}
#else
#define VENC_ANNOUNCE_VLC(name)
#define VENC_PUT_FLAG(BS, CODE, NAME) venc_vlc_put_bit(BS, CODE)
#define VENC_PUT_BITS(BS, CODE, LENGTH, NAME) venc_vlc_put_bits(BS, CODE, LENGTH)
#define VENC_PUT_UE_V(BS, CODE, NAME) venc_vlc_put_vlc(BS, CODE)
#define VENC_PUT_SE_V(BS, CODE, NAME) venc_vlc_put_vlc(BS, F265_SIGNED_VLC_CODE(CODE))
#endif

// Return the identifiers for VPS/SPS. Identify the stream origin.
static int venc_get_vps_id(f265_enc *enc) { return enc->gd.hm_gop_compat_flag ? 0 : 11; }
static int venc_get_sps_id(f265_enc *enc) { return enc->gd.hm_gop_compat_flag ? 0 : 2; }

// Write the video parameter set in the bitstream.
void venc_write_vps(f265_vlc_bs *vbs, f265_enc *enc)
{
    f265_gen_data *gd = &enc->gd;

    VENC_ANNOUNCE_VLC("VPS");
    VENC_PUT_BITS(vbs, venc_get_vps_id(enc), 4, "vps_video_parameter_set_id");
    VENC_PUT_BITS(vbs, 3, 2, "vps_reserved_three_2bits");
    VENC_PUT_BITS(vbs, 0, 6, "vps_reserved_zero_6bits");
    VENC_PUT_BITS(vbs, 0, 3, "vps_max_sub_layers_minus1");
    VENC_PUT_FLAG(vbs, 1, "vps_temporal_id_nesting_flag");
    VENC_PUT_BITS(vbs, 0xffff, 16, "vps_reserved_ffff_16bits");
    VENC_PUT_BITS(vbs, 0, 2, "XXX_profile_space"); // Must be 0.
    VENC_PUT_BITS(vbs, 1, 1, "XXX_tier_flag");     // Main tier = 0, High tier = 1.
    VENC_PUT_BITS(vbs, 1, 5, "XXX_profile_idc");   // Main = 1, Main10 = 2, MainStill = 3.
    VENC_PUT_BITS(vbs, 0x60000000, 32, "XXX_profile_compatibility_flag");
    VENC_PUT_FLAG(vbs, 0, "general_progressive_source_flag");
    VENC_PUT_FLAG(vbs, 0, "general_interlaced_source_flag");
    VENC_PUT_FLAG(vbs, 0, "general_non_packed_constraint_flag");
    VENC_PUT_FLAG(vbs, 0, "general_frame_only_constraint_flag");
    VENC_PUT_BITS(vbs, 0, 44, "XXX_reserved_zero_44bits");
    VENC_PUT_BITS(vbs, 153, 8, "general_level_idc"); // 30* Level. 5.1*30=153.
    VENC_PUT_FLAG(vbs, 1, "vps_sub_layer_ordering_info_present_flag");
    VENC_PUT_UE_V(vbs, gd->nb_refs, "vps_max_dec_pic_buffering_minus1");
    VENC_PUT_UE_V(vbs, gd->nb_reordered_frames, "vps_num_reorder_pics");
    VENC_PUT_UE_V(vbs, 0, "vps_max_latency_increase_plus1");
    VENC_PUT_BITS(vbs, 0, 6, "vps_max_nuh_reserved_zero_layer_id");
    VENC_PUT_UE_V(vbs, 0, "vps_max_op_sets_minus1");
    VENC_PUT_FLAG(vbs, 0, "vps_timing_info_present_flag");
    VENC_PUT_FLAG(vbs, 0, "vps_extension_flag");
    venc_vlc_flush(vbs);
}

// Write the sequence parameter set in the bitstream.
void venc_write_sps(f265_vlc_bs *vbs, f265_enc *enc)
{
    f265_gen_data *gd = &enc->gd;

    VENC_ANNOUNCE_VLC("SPS");
    VENC_PUT_BITS(vbs, venc_get_vps_id(enc), 4, "sps_video_parameter_set_id");
    VENC_PUT_BITS(vbs, 0, 3, "sps_max_sub_layers_minus1");
    VENC_PUT_FLAG(vbs, 1, "sps_temporal_id_nesting_flag");
    VENC_PUT_BITS(vbs, 0, 2, "XXX_profile_space");
    VENC_PUT_BITS(vbs, 1, 1, "XXX_tier_flag");
    VENC_PUT_BITS(vbs, 1, 5, "XXX_profile_idc");
    VENC_PUT_BITS(vbs, 0x60000000, 32, "XXX_profile_compatibility_flag");
    VENC_PUT_FLAG(vbs, 0, "general_progressive_source_flag");
    VENC_PUT_FLAG(vbs, 0, "general_interlaced_source_flag");
    VENC_PUT_FLAG(vbs, 0, "general_non_packed_constraint_flag");
    VENC_PUT_FLAG(vbs, 0, "general_frame_only_constraint_flag");
    VENC_PUT_BITS(vbs, 0, 44, "XXX_reserved_zero_44bits");
    VENC_PUT_BITS(vbs, 153, 8, "general_level_idc");
    VENC_PUT_UE_V(vbs, venc_get_sps_id(enc), "sps_seq_parameter_set_id");
    VENC_PUT_UE_V(vbs, gd->chroma_format, "chroma_format_idc");
    VENC_PUT_UE_V(vbs, gd->pix_dim[0], "pic_width_in_luma_samples");
    VENC_PUT_UE_V(vbs, gd->pix_dim[1], "pic_height_in_luma_samples");

    // Handle clipping. Assumes chroma format == 1.
    int clip_flag = gd->hm_gop_compat_flag || gd->pix_dim[0] != gd->clip_dim[0] || gd->pix_dim[1] != gd->clip_dim[1];
    VENC_PUT_FLAG(vbs, clip_flag, "conformance_window_flag");
    if (clip_flag)
    {
        VENC_PUT_UE_V(vbs, 0, "conf_win_left_offset");
        VENC_PUT_UE_V(vbs, (gd->pix_dim[0] - gd->clip_dim[0])>>1, "conf_win_right_offset");
        VENC_PUT_UE_V(vbs, 0, "conf_win_top_offset");
        VENC_PUT_UE_V(vbs, (gd->pix_dim[1] - gd->clip_dim[1])>>1, "conf_win_bottom_offset");
    }

    VENC_PUT_UE_V(vbs, gd->bit_depth[0] - 8, "bit_depth_luma_minus8");
    VENC_PUT_UE_V(vbs, gd->bit_depth[1] - 8, "bit_depth_chroma_minus8");
    VENC_PUT_UE_V(vbs, gd->poc_bits - 4, "log2_max_pic_order_cnt_lsb_minus4");
    VENC_PUT_FLAG(vbs, 1, "sps_sub_layer_ordering_info_present_flag");
    VENC_PUT_UE_V(vbs, gd->nb_refs, "sps_max_dec_pic_buffering_minus1");
    VENC_PUT_UE_V(vbs, gd->nb_reordered_frames, "sps_num_reorder_pics");
    VENC_PUT_UE_V(vbs, 0, "sps_max_latency_increase_plus1");
    VENC_PUT_UE_V(vbs, gd->cb_range[0]-3, "log2_min_coding_block_size_minus3");
    VENC_PUT_UE_V(vbs, gd->cb_range[1] - gd->cb_range[0], "log2_diff_max_min_coding_block_size");
    VENC_PUT_UE_V(vbs, gd->tb_range[0]-2, "log2_min_transform_block_size_minus2");
    VENC_PUT_UE_V(vbs, gd->tb_range[1] - gd->tb_range[0], "log2_diff_max_min_transform_block_size");
    VENC_PUT_UE_V(vbs, gd->tb_depth[1], "max_transform_hierarchy_depth_inter");
    VENC_PUT_UE_V(vbs, gd->tb_depth[0], "max_transform_hierarchy_depth_intra");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_SCALING), "scaling_list_enabled_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_AMP), "amp_enabled_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_SAO), "sample_adaptive_offset_enabled_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_PCM), "pcm_enabled_flag");
    VENC_PUT_UE_V(vbs, 0, "num_short_term_ref_pic_sets");
    VENC_PUT_FLAG(vbs, 0, "long_term_ref_pics_present_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_TMV), "sps_temporal_mvp_enable_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_SMOOTH_INTRA), "sps_strong_intra_smoothing_enable_flag");
    VENC_PUT_FLAG(vbs, 0, "vui_parameters_present_flag");
    VENC_PUT_FLAG(vbs, 0, "sps_extension_flag");

    venc_vlc_flush(vbs);
}

// Write the picture parameter set in the bitstream.
void venc_write_pps(f265_vlc_bs *vbs, f265_enc *enc)
{
    f265_gen_data *gd = &enc->gd;
    int deblock_flag = F265_GET_FLAG(gd->eflags, F265_PF_DEBLOCK);

    VENC_ANNOUNCE_VLC("PPS");
    VENC_PUT_UE_V(vbs, 0, "pps_pic_parameter_set_id");
    VENC_PUT_UE_V(vbs, venc_get_sps_id(enc), "pps_seq_parameter_set_id");
    VENC_PUT_FLAG(vbs, 0, "dependent_slice_segments_enabled_flag");
    VENC_PUT_FLAG(vbs, 0, "output_flag_present_flag");
    VENC_PUT_BITS(vbs, 0, 3, "num_extra_slice_header_bits");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_SIGN_HIDING), "sign_data_hiding_flag");
    VENC_PUT_FLAG(vbs, 0, "cabac_init_present_flag");
    VENC_PUT_UE_V(vbs, gd->default_nb_ref_idx[0]-1, "num_ref_idx_l0_default_active_minus1");
    VENC_PUT_UE_V(vbs, gd->default_nb_ref_idx[1]-1, "num_ref_idx_l1_default_active_minus1");

    // Suboptimal for HM compatibility.
    VENC_PUT_SE_V(vbs, 0, "init_qp_minus26");

    VENC_PUT_FLAG(vbs, 0, "constrained_intra_pred_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_TRANSFORM_SKIP), "transform_skip_enabled_flag");
    VENC_PUT_FLAG(vbs, 0, "cu_qp_delta_enabled_flag");
    VENC_PUT_SE_V(vbs, 0, "pps_cb_qp_offset");
    VENC_PUT_SE_V(vbs, 0, "pps_cr_qp_offset");
    VENC_PUT_FLAG(vbs, 0, "pps_slice_chroma_qp_offsets_present_flag");
    VENC_PUT_FLAG(vbs, 0, "weighted_pred_flag");
    VENC_PUT_FLAG(vbs, 0, "weighted_bipred_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_TRANSQUANT_BYPASS), "transquant_bypass_enable_flag");
    VENC_PUT_FLAG(vbs, 0, "tiles_enabled_flag");
    VENC_PUT_FLAG(vbs, F265_GET_FLAG(gd->eflags, F265_PF_WPP), "entropy_coding_sync_enabled_flag");
    VENC_PUT_FLAG(vbs, 1, "loop_filter_across_slices_enabled_flag");

    VENC_PUT_FLAG(vbs, !deblock_flag, "deblocking_filter_control_present_flag");
    if (!deblock_flag)
    {
        VENC_PUT_FLAG(vbs, 0, "deblocking_filter_override_enabled_flag");
        VENC_PUT_FLAG(vbs, 1, "pps_disable_deblocking_filter_flag");
    }

    VENC_PUT_FLAG(vbs, 0, "pps_scaling_list_data_present_flag");
    VENC_PUT_FLAG(vbs, 0, "lists_modification_present_flag");
    VENC_PUT_UE_V(vbs, gd->parallel_merge_level-2, "log2_parallel_merge_level_minus2");
    VENC_PUT_FLAG(vbs, 0, "slice_segment_header_extension_present_flag");
    VENC_PUT_FLAG(vbs, 0, "pps_extension_flag");
    venc_vlc_flush(vbs);
}

// Write the slice segment header in the bitstream.
void venc_write_slice_header(f265_vlc_bs *vbs, f265_enc *enc, f265_frame *f, int *chunk_sizes, int nb_chunks,
                             uint32_t seg_flags, int ctb_xy, int nal_type)
{
    f265_gen_data *gd = &enc->gd;
    int deblock_flag = F265_GET_FLAG(gd->eflags, F265_PF_DEBLOCK);
    int tmv_flag = F265_GET_FLAG(enc->gd.eflags, F265_PF_TMV);

    VENC_ANNOUNCE_VLC("Slice header");
    VENC_PUT_FLAG(vbs, 1, "first_slice_segment_in_pic_flag");
    if (nal_type >= 16 && nal_type <= 23) VENC_PUT_FLAG(vbs, 0, "no_output_of_prior_pics_flag");
    VENC_PUT_UE_V(vbs, 0, "slice_pic_parameter_set_id");
    VENC_PUT_UE_V(vbs, 2-f->frame_type, "slice_type");

    if (!(f->gen_flags&F265_FF_IDR))
    {
        VENC_PUT_BITS(vbs, f->h265_poc&((1<<gd->poc_bits)-1), gd->poc_bits, "slice_pic_order_cnt_lsb");
        VENC_PUT_FLAG(vbs, 0, "short_term_ref_pic_set_sps_flag");

        // Intra-predicted RPS. All the frames are used for reference, unless
        // it's a CRA frame.
        int use_ref_flag = !(f->gen_flags&F265_FF_CRA);
        VENC_PUT_UE_V(vbs, f->dpb_neg, "num_negative_pics");
        VENC_PUT_UE_V(vbs, f->dpb_pos, "num_positive_pics");
        int64_t last_poc = f->h265_poc;
        for (int i = 0; i < f->dpb_neg; i++)
        {
            int64_t next_poc = f->dpb[f->dpb_neg-1-i]->h265_poc;
            VENC_PUT_UE_V(vbs, last_poc-next_poc-1, "delta_poc_s0_minus1");
            VENC_PUT_FLAG(vbs, use_ref_flag, "used_by_curr_pic_s0_flag");
            last_poc = next_poc;
        }
        last_poc = f->h265_poc;
        for (int i = 0; i < f->dpb_pos; i++)
        {
            int64_t next_poc = f->dpb[f->dpb_neg+i]->h265_poc;
            VENC_PUT_UE_V(vbs, next_poc-last_poc-1, "delta_poc_s1_minus1");
            VENC_PUT_FLAG(vbs, use_ref_flag, "used_by_curr_pic_s1_flag");
            last_poc = next_poc;
        }

        // Present even for non-IDR I frames.
        if (tmv_flag) VENC_PUT_FLAG(vbs, 1, "slice_temporal_mvp_enabled_flag");
    }

    if (F265_GET_FLAG(gd->eflags, F265_PF_SAO))
    {
        VENC_PUT_FLAG(vbs, 1, "slice_sao_luma_flag");
        VENC_PUT_FLAG(vbs, 1, "slice_sao_chroma_flag");
    }

    if (f->frame_type != F265_FRAME_I)
    {
        int override_flag = gd->default_nb_ref_idx[0] != f->nb_ref_idx[0] ||
                            f->frame_type == F265_FRAME_B && gd->default_nb_ref_idx[1] != f->nb_ref_idx[1];
        VENC_PUT_FLAG(vbs, override_flag, "num_ref_idx_active_override_flag");
        if (override_flag)
        {
            VENC_PUT_UE_V(vbs, f->nb_ref_idx[0]-1, "num_ref_idx_l0_active_minus1");
            if (f->frame_type == F265_FRAME_B)
                VENC_PUT_UE_V(vbs, f->nb_ref_idx[1]-1, "num_ref_idx_l1_active_minus1");
        }

        if (f->frame_type == F265_FRAME_B)
            VENC_PUT_FLAG(vbs, 0, "mvd_l1_zero_flag");

        if (tmv_flag)
        {
            int col_list = F265_GET_FLAG(f->temp_ref_info, 1<<4);

            if (f->frame_type == F265_FRAME_B)
                VENC_PUT_FLAG(vbs, !col_list, "collocated_from_l0_flag");

            if (f->nb_ref_idx[col_list] > 1)
                VENC_PUT_UE_V(vbs, (f->temp_ref_info&0xf), "collocated_ref_idx");
        }

        VENC_PUT_UE_V(vbs, 5-gd->merge_cand, "five_minus_max_num_merge_cand");
    }

    // HM compatibility.
    VENC_PUT_SE_V(vbs, f->qp - 26, "slice_qp_delta");

    if (deblock_flag) VENC_PUT_FLAG(vbs, 1, "slice_loop_filter_across_slices_enabled_flag");

    if (gd->eflags&F265_PF_WPP)
    {
        VENC_PUT_UE_V(vbs, nb_chunks-1, "num_entry_point_offsets");

        if (nb_chunks > 1)
        {
            int max_chunk_size = 0;
            for (int i = 0; i < nb_chunks - 1; i++)
                max_chunk_size = F265_MAX(max_chunk_size, chunk_sizes[i]);
            int nb_bits = 32 - __builtin_clz(max_chunk_size-1);

            VENC_PUT_UE_V(vbs, nb_bits-1, "offset_len_minus1");

            for (int i = 0; i < nb_chunks - 1; i++)
                VENC_PUT_BITS(vbs, chunk_sizes[i]-1, nb_bits, "entry_point_offset_minus1");
        }
    }

    venc_vlc_flush(vbs);
}

// Return the NAL type of the frame.
static int venc_get_frame_nal_type(f265_frame *f)
{
    // True if the current frame can be used for reference.
    int ref_flag = f->gop_entry ? f->gop_entry->ref_flag : F265_GET_FLAG(f->gen_flags, F265_FF_REF);

    // IDR => IDR_W_RADL (HM compatibility, there are no leading frames).
    if (f->gen_flags&F265_FF_IDR) return 19;

    // CRA => CRA_NUT.
    if (f->gen_flags&F265_FF_CRA) return 21;

    // The current frame is displayed before a CRA frame => RASL.
    if (f->dpb_size && (f->dpb[f->dpb_size-1]->gen_flags&F265_FF_CRA) && f->abs_poc < f->dpb[f->dpb_size-1]->abs_poc)
        return 8+ref_flag;

    // Regular frame => TRAIL.
    return ref_flag;
}

// Begin packaging the NAL data.
void venc_begin_nal(f265_nal_bs *nbs, int nal_type, int zero_byte_flag)
{
    uint8_t *dst = nbs->nal_start;

    // Annex B encapsulation:
    // - leading_zero_8bits (0x00), optional.
    // - zero_byte (0x00), mandatory for VPS/SPS/PPS or the first NAL of the
    //   the access unit.
    // - start_code_prefix (0x000001).
    // - NAL header (2 bytes).
    // - NAL data (variable).
    // - trailing_zero_8bits (0x00), optional.

    // Annex B zero byte and start code prefix.
    if (nbs->format == F265_FORMAT_BYTE)
    {
        if (zero_byte_flag) *dst++ = 0;
        *dst++ = 0;
        *dst++ = 0;
        *dst++ = 1;
    }

    // Leave space for the NAL size.
    else dst += nbs->nal_size_len;

    // Write the NAL header.
    // +-----+------+----------+---------+
    // | fzb | type | layer id | temp id |
    // +-----+------+----------+---------+
    //      fzb: forbidden_zero_bit (MSB of first byte)
    //     type: nal_unit_type (Bits 1 to 6 of first byte)
    // layer id: nuh_layer_id (LSB of first byte, and first 5 bits of second byte)
    //  temp id: nuh_temporal_id_plus1 (last 3 bits of second byte)
    //
    // The spec currently says "nuh_layer_id shall be equal to 0", and
    // "The value of nuh_temporal_id_plus1 shall not be equal to 0.".
    //
    // Left shifting the type by 1 makes sure fzb is 0, type occupies the
    // *middle* bits of the first byte, and that the layer id's first bit
    // is 0 (the 6 bits have to be 0). Signaling 1 in the second byte
    // makes sure the layer id is 0, and that nuh_temporal_id_plus1 is
    // not 0.
    *dst++ = nal_type << 1;
    *dst++ = 1;

    // Set the data start position.
    nbs->data_start = dst;
}

// Finish packaging the NAL data.
void venc_finish_nal(f265_nal_bs *nbs, int data_size)
{
    uint8_t *nal_end = nbs->data_start + data_size;

    // Write the NAL size in the bytestream.
    if (nbs->format != F265_FORMAT_BYTE)
    {
        // The NAL length excludes the length field itself.
        int nal_size_len = nbs->nal_size_len;
        int len = nal_end - nbs->nal_start - nal_size_len;
        uint8_t *dst = nbs->nal_start + nal_size_len - 1;
        while (nal_size_len--)
        {
            *dst-- = len;
            len >>= 8;
        }
    }

    // Pass to the next NAL.
    nbs->nal_start = nal_end;
}

// Begin writing the VLC bitstream in a buffer.
void venc_begin_vlc(f265_vlc_bs *vbs, uint8_t *base)
{
    vbs->base = vbs->cur = base;
    vbs->bits_left = 8;
    vbs->trail_zeroes = 0;
    base[0] = 0;
}

// Write a NAL containing only VLC data.
void venc_write_pure_vlc_nal(f265_nal_bs *nbs, f265_enc *enc, int nal_type,
                             void (*data_writer)(f265_vlc_bs*, f265_enc*))
{
    f265_vlc_bs vbs;
    venc_begin_nal(nbs, nal_type, 1);
    venc_begin_vlc(&vbs, nbs->data_start);
    data_writer(&vbs, enc);
    venc_finish_nal(nbs, vbs.cur - vbs.base);
}

// Write a video parameter set NAL.
static void venc_write_vps_nal(f265_nal_bs *nbs, f265_enc *enc)
{
    venc_write_pure_vlc_nal(nbs, enc, 32, venc_write_vps);
}

// Write a sequence parameter set NAL.
static void venc_write_sps_nal(f265_nal_bs *nbs, f265_enc *enc)
{
    venc_write_pure_vlc_nal(nbs, enc, 33, venc_write_sps);
}

// Write a picture parameter set NAL.
static void venc_write_pps_nal(f265_nal_bs *nbs, f265_enc *enc)
{
    venc_write_pure_vlc_nal(nbs, enc, 34, venc_write_pps);
}

// Write a slice segment NAL. The segment pointers are updated.
void venc_write_segment_nal(f265_nal_bs *nbs, f265_enc *enc, f265_frame *f, uint32_t **seg, uint8_t **chunks,
                            int nal_type, int zero_byte_flag)
{
    uint32_t seg_flags = (*seg)[0];
    int ctb_xy = (*seg)[1];
    int nb_chunks = (*seg)[2];
    int *chunk_sizes = (int*)((*seg) + 3);

    // Begin the NAL.
    venc_begin_nal(nbs, nal_type, zero_byte_flag);

    // Write the segment header.
    f265_vlc_bs vbs;
    venc_begin_vlc(&vbs, nbs->data_start);
    venc_write_slice_header(&vbs, enc, f, chunk_sizes, nb_chunks, seg_flags, ctb_xy, nal_type);
    int seg_hdr_size = (vbs.cur - vbs.base);

    // Write the segment chunks.
    int seg_data_size = 0;
    for (int i = 0; i < nb_chunks; i++) seg_data_size += chunk_sizes[i];
    memcpy(nbs->data_start + seg_hdr_size, *chunks, seg_data_size);

    // Finish the NAL.
    venc_finish_nal(nbs, seg_hdr_size + seg_data_size);

    // Update the pointers.
    *seg += 3 + nb_chunks;
    *chunks += seg_data_size;
}

// Write the NAL units of the current frame.
void venc_encode_frame_nals(f265_enc *enc, f265_frame *f, uint8_t *chunk_buf)
{
    f265_enc_req *req = enc->md.req;
    uint32_t frame_flags = f->gen_flags;
    f265_frame_map *fmap = f->fmap;

    // Process each output frame.
    for (int output_idx = 0; output_idx < req->nb_outputs; output_idx++)
    {
        f265_output_frame *of = req->outputs + output_idx;

        // Begin the NAL bytestream.
        f265_nal_bs nbs;
        nbs.nal_start = of->bs;
        nbs.format = of->format;
        nbs.nal_size_len = of->nal_size_len;

        // Manage the zero_byte status.
        int zero_byte_flag = 1;

        // FIXME: HM compatibility.
        if (f->abs_poc) frame_flags = 0;

        f->last_idr_abs_poc = 0;
        if (enc->md.enc_frames[0])
            f->last_idr_abs_poc = enc->md.enc_frames[0]->last_idr_abs_poc;

        // Write the VPS/SPS/PPS NALs.
        if (frame_flags & F265_FF_IDR)
        {
            if (of->vps_flag) venc_write_vps_nal(&nbs, enc);
            if (of->sps_flag) venc_write_sps_nal(&nbs, enc);
            if (of->pps_flag) venc_write_pps_nal(&nbs, enc);
            zero_byte_flag = !(of->vps_flag|of->sps_flag|of->pps_flag);
            f->last_idr_abs_poc = f->abs_poc;
        }
        f->h265_poc = f->abs_poc - f->last_idr_abs_poc;

        // Write each section.
        int nal_type = venc_get_frame_nal_type(f);
        for (int section_idx = 0; section_idx < fmap->nb_sections; section_idx++)
        {
            int ctb_idx = fmap->section_map[section_idx][0];
            int nb_ctbs = fmap->section_map[section_idx][1];
            uint32_t *nb_segs, *seg;
            uint8_t *chunks;
            venc_get_chunk_buf_offsets(enc, chunk_buf, &nb_segs, &seg, &chunks, ctb_idx, nb_ctbs);

            // Write each segment NAL.
            for (int seg_idx = 0; seg_idx < *nb_segs; seg_idx++, zero_byte_flag = 0)
                venc_write_segment_nal(&nbs, enc, f, &seg, &chunks, nal_type, zero_byte_flag);
        }

        // Finish the NAL bytestream.
        of->bs_size = nbs.nal_start - of->bs;

        // Set the frame properties.
        of->timestamp = f->timestamp;
        of->duration = f->duration;
        of->frame_type = f->frame_type;
    }
}

// FIXME: when we actually support delayed SAO, handle cbs->mode in these
// functions.

// Initialize the CABAC engine. Make sure buf is aligned on 32-bits.
void venc_init_cabac_engine(f265_cabac_bs *cbs, uint8_t *buf)
{
    // Leave four bytes free to unify the snapshot (backup the last four bytes).
    cbs->raw.base = cbs->raw.cur = buf + 4;

    // The first bit is discarded, so we need to encode 32+1 bits before
    // outputting a 32-bit chunk.
    cbs->raw.range = 510;
    cbs->raw.low = 0;
    cbs->raw.bits_left = 32+1;
    cbs->raw.chunks_outstanding = 0;
}

// Initialize the CABAC contexts for a slice type and QP.
void venc_init_cabac_contexts(f265_cabac_bs *cbs, int slice_type, int qp)
{
    qp = F265_CLAMP(qp, 0, 51);
    const uint8_t *init_table = f265_cabac_ctx_init_table[slice_type];
    for (int i = 0; i < F265_NB_CABAC_CTX; i++)
    {
        int init_val = init_table[i];
        int slope = (init_val>>4)*5 - 45;
        int offset = ((init_val&15)<<3) - 16;
        int tmp = F265_CLAMP(((slope*qp)>>4) + offset, 1, 126);
        int mps = tmp > 63;
        int state = mps ? tmp - 64 : 63 - tmp;
        int ctx = (state<<1)|mps;
        cbs->contexts[i] = ctx;
    }
}

// Load the CABAC engine state from a snapshot.
void venc_load_cabac_engine(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    cbs->raw.cur = snapshot->cur;
    cbs->raw.low = snapshot->low;
    cbs->raw.range = snapshot->range;
    cbs->raw.bits_left = snapshot->bits_left;
    cbs->raw.chunks_outstanding = snapshot->chunks_outstanding;
    *(uint32_t*)(cbs->raw.cur-4) = snapshot->prev_chunk;
}

// Load the CABAC contexts from a snapshot.
void venc_load_cabac_contexts(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    memcpy(cbs->contexts, snapshot->contexts, sizeof(cbs->contexts));
}

// Load the full CABAC state from a snapshot.
void venc_load_cabac(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    venc_load_cabac_engine(cbs, snapshot);
    venc_load_cabac_contexts(cbs, snapshot);
}

// Save the CABAC engine state in a snapshot. We are guaranteed to be still
// aligned when we get a snapshot.
void venc_save_cabac_engine(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    snapshot->cur = cbs->raw.cur;
    snapshot->low = cbs->raw.low;
    snapshot->range = cbs->raw.range;
    snapshot->bits_left = cbs->raw.bits_left;
    snapshot->chunks_outstanding = cbs->raw.chunks_outstanding;
    snapshot->prev_chunk = *(uint32_t*)(cbs->raw.cur-4);
}

// Save the CABAC contexts in a snapshot.
void venc_save_cabac_contexts(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    memcpy(snapshot->contexts, cbs->contexts, sizeof(cbs->contexts));
}

// Save the full CABAC state in a snapshot.
void venc_save_cabac(f265_cabac_bs *cbs, f265_cabac_snapshot *snapshot)
{
    venc_save_cabac_engine(cbs, snapshot);
    venc_save_cabac_contexts(cbs, snapshot);
}

// The CABAC bitstream writing is optimized as follow.
//
// Conceptually, the output bitstream corresponds to the value of the low
// register, assuming that the register grows as large as needed. The value in
// the low register increases when the LPS symbol is selected. The value never
// decrease because it represents a minimum bound. Each step of the
// renormalization operation appends a zero at the right of the low register to
// add precision.
//
// The implementation simulates this logic with a 64-bit register. The bitstream
// is written left-to-right in 32-bit chunks. Once written, those bits are
// removed from the low register, which prevents the register from growing
// indefinitely. Picture:
//                  32-bit chunk 0                 32-bit chunk 1
// Low: [1001100101010101010001101110101][1001000101010101010001101110101][...]
//
// When the value of the low register increases, all the bitstream can
// potentially change by the propagation of the carry bit in the addition.
// However, by the nature of interval division, the propagation of the carry can
// occur at most once for a given bit, excluding the low 10 bits of the low
// register. Thus, we write a chunk when we have accumulated at least 32+10 bits
// in the low register and we apply the propagation of the carry to the
// previously-written chunks.
//
// Three cases can happen when we write a 32-bit chunk, based on its content
// (the chunk bits are within the brackets, the carry bit is to the left):
// 1) 0[11111111111111111111111111111111].
// 2) 0[bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb] (at least one 'b' is zero).
// 3) 1[bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb] (at least one 'b' is zero).
//
// In case 1, the chunk is not written immediately because a future increase in
// the low register could change the value of the chunk and the previous chunks.
// Instead, the outstanding chunk counter is incremented by 1 and the chunk bits
// are removed from the low register. The carry bit is necessarily 0 when the
// chunk bits are all 1, so we need not consider it.
//
// In cases 2 and 3, the chunk contains at least one zero bit. A future increase
// in the low register may change the chunk value, but the propagation will stop
// in that chunk. In case 2, there is no carry bit. Thus, the outstanding chunks
// are written (0xffffffff) in the bitstream, then the current chunk is written.
// In case 3, there is a carry bit which propagates in the outstanding chunks up
// to the previously-written chunk. The carry bit is added to the
// previously-written chunk, the outstanding chunks are written as 0x00000000
// and the current chunk is written. In both cases, the chunk bits and the carry
// bit are removed from the low register, and the outstanding chunk counter is
// cleared.
void venc_write_cabac_chunk(f265_cabac_bs *cbs)
{
    // Not enough accumulated bits, bail out.
    if (cbs->raw.bits_left > 0) return;

    // Extract the carry bit and the 32 chunk bits from the low register.
    uint64_t chunk = cbs->raw.low >> (-cbs->raw.bits_left+10);
    int carry = chunk>>32;
    chunk &= 0xffffffff;
    cbs->raw.low &= (0x400<<(-cbs->raw.bits_left)) - 1;
    cbs->raw.bits_left += 32;

    // Case 1: outstanding chunk.
    if (chunk == 0xffffffff) cbs->raw.chunks_outstanding++;

    // Case 2 or 3: actual bitstream writing. The bitstream is guaranteed to be
    // aligned on 32-bit at this point.
    else
    {
        // Add the carry bit to the previously-written chunk, big endian.
        uint32_t *cur = (uint32_t*)cbs->raw.cur;
        cur[-1] = F265_BYTESWAP32(F265_BYTESWAP32(cur[-1]) + carry);

        // Write the outstanding chunks.
        uint32_t outstanding_val = carry - 1;
        for (; cbs->raw.chunks_outstanding; cbs->raw.chunks_outstanding--) *cur++ = outstanding_val;

        // Write the current chunk.
        *cur++ = F265_BYTESWAP32((uint32_t)chunk);
        cbs->raw.cur = (uint8_t*)cur;
    }
}

// Renormalize the BAC after encoding a context bin.
void venc_renormalize_cabac(f265_cabac_bs *cbs)
{
    // Add precision bits until the range is no longer under 256.
    int shift = 8 - (63 - __builtin_clzll(cbs->raw.range));
    cbs->raw.low <<= shift;
    cbs->raw.range <<= shift;
    cbs->raw.bits_left -= shift;

    // Write the ready chunks.
    venc_write_cabac_chunk(cbs);
}

// Encode a bin with a probability context.
void venc_encode_context_bin_raw(f265_cabac_bs *cbs, int ctx_idx, int bin)
{
    // Get the context and LPS value.
    int context = cbs->contexts[ctx_idx];
    int range_lps = f265_cabac_range_table[context>>1][(cbs->raw.range>>6)&3];
    #ifdef VAN_TRACE_SYNTAX
    venc_trace_bin(0, context, bin);
    #endif

    // Update the context.
    cbs->contexts[ctx_idx] = f265_cabac_transit_table[context][bin];

    // Select the MPS tentatively.
    cbs->raw.range -= range_lps;

    // Select the LPS as needed.
    if (bin != (context & 1))
    {
        cbs->raw.low += cbs->raw.range;
        cbs->raw.range = range_lps;
    }

    // Renormalize.
    venc_renormalize_cabac(cbs);
}

// Encode a bin without a probability context (assume 50% probability).
void venc_encode_bypass_bin_raw(f265_cabac_bs *cbs, int bin)
{
    #ifdef VAN_TRACE_SYNTAX
    venc_trace_bin(1, 0, bin);
    #endif
    cbs->raw.low <<= 1;
    if (bin) cbs->raw.low += cbs->raw.range;
    cbs->raw.bits_left--;
    venc_write_cabac_chunk(cbs);
}

// Encode several bins in bypass.
void venc_encode_bypass_bins_raw(f265_cabac_bs *cbs, int bins, int nb_bins)
{
    for (int i = nb_bins - 1; i >= 0; i--) venc_encode_bypass_bin_raw(cbs, (bins>>i)&1);
}

// Encoding a terminating bin with value 0.
void venc_encode_term0_bin_raw(f265_cabac_bs *cbs)
{
    #ifdef VAN_TRACE_SYNTAX
    venc_trace_bin(2, 0, 0);
    #endif
    cbs->raw.range -= 2;
    venc_renormalize_cabac(cbs);
}

// Encode a terminating bin with value 1 and flush the bitstream.
void venc_encode_term1_bin(f265_cabac_bs *cbs)
{
    // We assume that there is no ready chunk.
    assert(cbs->raw.bits_left > 0);

    // Get the final bits of the low register. We encode the LPS symbol and set
    // the stop bit to 1.
    cbs->raw.low += cbs->raw.range - 2;
    cbs->raw.low |= 1;

    // We must put the final bits of the low register as-is into the bitstream.
    // There are (32-bits_left + 10) bits in the low register, plus a carry bit.
    // There might be outstanding chunks.
    //
    // We shift the low 10 bits left to promote them to chunk bits. To flush the
    // outstanding chunks, we add a dummy zero byte to the low register so that
    // the last chunk contains some zero bits.
    //
    // There is one or two chunks to output, depending on the number of bits
    // present in the low register. We output the first chunk if there is one.
    // Then, we shift the remaining bits to the left of the last chunk and
    // output the last chunk.
    //
    // We remove the extra bytes at the end (1 dummy byte plus up to four empty
    // bytes in the last chunk).
    cbs->raw.low <<= 18;
    cbs->raw.bits_left -= 18;
    venc_write_cabac_chunk(cbs);
    int nb_extra = 1 + 4 - ((32 - cbs->raw.bits_left + 7)>>3);
    cbs->raw.low <<= cbs->raw.bits_left;
    cbs->raw.bits_left = 0;
    venc_write_cabac_chunk(cbs);
    cbs->raw.cur -= nb_extra;
}

// Get the offsets in the chunk buffer for the section specified.
void venc_get_chunk_buf_offsets(f265_enc *enc, uint8_t *buf, uint32_t **nb_segs, uint32_t **seg, uint8_t **chunks,
                                int ctb_idx, int nb_ctbs)
{
    // Assuming the worst case of 1 segment per CTB.
    *nb_segs = (uint32_t*)(buf + ctb_idx*enc->gd.ctb_chunk_size);
    *seg = *nb_segs + 1;
    *chunks = (uint8_t*)(*seg + 4*nb_ctbs);
}

// Initialize the current segment.
void venc_init_segment(f265_enc_thread *t, int dependent_flag)
{
    t->seg[0] = dependent_flag;
    t->seg[1] = t->ctb_xy;
    t->seg[2] = 0;
    t->init_segment_flag = 0;
}

// Initialize the current segment chunk. This doesn't modify the chunk buffer.
void venc_init_seg_chunk(f265_seg_chunk *chunk, uint8_t *base)
{
    chunk->base = chunk->cur = base;
    chunk->trail_zeroes = 0;
    chunk->start_flag = 1;
}

// Write end_of_slice_segment_flag/end_of_sub_stream_one_bit.
void venc_write_stream_end(f265_cabac_bs *cbs, int segment_end_flag, int chunk_end_flag)
{
    if (segment_end_flag) venc_encode_term1_bin(cbs);
    else
    {
        venc_trace_syntax_element("end_of_slice_segment_flag");
        venc_encode_term0_bin_raw(cbs);
        if (chunk_end_flag) venc_encode_term1_bin(cbs);
    }
}

// Escape and append the CABAC bitstream to a segment chunk. If flushing, all
// the data is written. Otherwise, the last four bytes are skipped since they
// may still change. The chunk object is updated as the data is written. The
// CABAC state and the prior chunk bytes are unaffected.
void venc_write_seg_chunk(f265_cabac_bs *cbs, f265_seg_chunk *chunk, int chunk_end_flag)
{
    // Track the number of trailing zeroes without branches.
    const uint8_t mask_table_in[3] = { 0, 1, 3 };
    const uint8_t mask_table_out[4] = { 0, 1, 0, 2 };

    // Get the byte range to copy, if any.
    uint8_t *src = cbs->raw.base;
    uint8_t *src_end = cbs->raw.cur;
    if (likely(!chunk_end_flag))
    {
        if (unlikely(src_end - src <= 4)) return;
        src_end -= 4;
    }

    // Import the chunk state.
    uint8_t *dst = chunk->cur;
    int32_t prev_zero_mask = mask_table_in[chunk->trail_zeroes];

    // Escape the RBSP.
    while (src < src_end)
    {
        // Check for start code emulation.
        if (unlikely(((prev_zero_mask&3) == 3) && !(*src&0xfc)))
        {
            *dst++ = 3;
            prev_zero_mask <<= 1;
        }

        // Store bit 1 for each zero byte.
        prev_zero_mask <<= 1;
        prev_zero_mask |= !*src;

        // Write the byte.
        *dst++ = *src++;
    }

    // Update the chunk.
    chunk->cur = dst;
    chunk->trail_zeroes = mask_table_out[prev_zero_mask&3];
}

// Write the CTB and the end-of-stream flags in the segment chunk tentatively.
void venc_write_ctb_chunk(f265_enc_thread *t, f265_seg_chunk *chunk, int segment_end_flag, int chunk_end_flag)
{
    f265_cabac_bs *cbs = &t->cbs;

    // Write the end-of-stream flags of the previous CTB in the chunk.
    if (!chunk->start_flag) venc_write_stream_end(cbs, 0, 0);

    // Write the current CTB.
    venc_write_ctb(t);

    // Write the end-of-stream flags of the last CTB of the chunk.
    if (chunk_end_flag) venc_write_stream_end(cbs, segment_end_flag, 1);

    // Append the current CTB to the chunk.
    venc_write_seg_chunk(cbs, chunk, chunk_end_flag);
}

// Commit the data written in a segment chunk.
void venc_commit_seg_chunk(f265_enc_thread *t, f265_seg_chunk *chunk, int segment_end_flag, int chunk_end_flag)
{
    f265_cabac_bs *cbs = &t->cbs;

    // End the current chunk and pass to the next, if any.
    if (unlikely(chunk_end_flag))
    {
        t->seg[3+t->seg[2]++] = chunk->cur - chunk->base;

        venc_init_seg_chunk(chunk, chunk->cur);
        venc_init_cabac_engine(cbs, t->cabac_buf);

        // End the current segment and pass to the next (without initializing).
        if (segment_end_flag)
        {
            (*t->nb_segs)++;
            t->seg += 3 + t->seg[2];
            t->init_segment_flag = 1;
        }
    }

    // Not the end of the chunk.
    else
    {
        chunk->start_flag = 0;

        // Discard the bytes that were committed to the chunk (all but the last
        // four).
        if (likely(cbs->raw.cur != cbs->raw.base))
        {
            *(uint32_t*)cbs->raw.base = *(uint32_t*)(cbs->raw.cur-4);
            cbs->raw.cur = cbs->raw.base + 4;
        }
    }
}

