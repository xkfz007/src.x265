// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Block entropy coding. This file is compiled twice. If F265_RDO is defined, we
// measure the loss of entropy, otherwise we encode the block. Functions must
// either be static or have a different name to avoid symbol clashes.

#include "f265/enc.h"

// Dispatch calls.
#ifdef F265_RDO
#define venc_encode_context_bin(cbs, ctx_idx, bin) venc_encode_context_bin_rdo(cbs, ctx_idx, bin)
#define venc_encode_bypass_bin(cbs, bin) venc_encode_bypass_bin_rdo(cbs, bin)
#define venc_encode_bypass_bins(cbs, bins, nb_bins) venc_encode_bypass_bins_rdo(cbs, bins, nb_bins)
#define venc_encode_term0_bin(cbs) venc_encode_term0_bin_rdo(cbs)
#else
#define venc_encode_context_bin(cbs, ctx_idx, bin) (cbs)->encode_context_bin(cbs, ctx_idx, bin)
#define venc_encode_bypass_bin(cbs, bin) (cbs)->encode_bypass_bin(cbs, bin)
#define venc_encode_bypass_bins(cbs, bins, nb_bins) (cbs)->encode_bypass_bins(cbs, bins, nb_bins)
#define venc_encode_term0_bin(cbs) (cbs)->encode_term0_bin(cbs)
#endif


///////////////////////////////////////////////////////////////////////////////
// RDO and non-RDO encoding functions.

// Compute an unary binarization string.
static finline void venc_unary_bin(int val, int max, int *bins, int *nb_bins)
{
    if (val == max)
    {
        *bins = 0xffff>>(16-max);
        *nb_bins = max;
    }
    else
    {
        *bins = (0xffff>>(15-val))-1;
        *nb_bins = val+1;
    }
}

// Compute a k-order exponential golumb binarization string.
static finline void venc_egk_bin(int val, int k, int *bins, int *nb_bins)
{
    // Offset corresponding to the order parameter.
    int o = val + (1<<k);

    // Number of bits in the remainder.
    int n = 31 - __builtin_clz(o);

    // Value of the remainder.
    int r = o - (1<<n);

    // Number of bits in the prefix.
    int len = 1 + n - k;

    // The first part is the number of leading ones followed by the 0. The
    // second part is the remainder.
    *bins = (((1<<len)-2)<<n)|(r&((1<<n)-1));
    *nb_bins = len+n;
}

// Retrieve the CB indices of the top and left CBs.
static finline void venc_get_left_and_top_cb_idx(f265_enc_thread *t, f265_cb *cb, int *left_cb_idx, int *top_cb_idx)
{
    // Optimize the pmap lookups eventually.
    int b4_off[2] = { cb->cb_off[0]>>2, cb->cb_off[1]>>2 };
    int part_idx;
    venc_lookup_pmap(t, b4_off[0]-1, b4_off[1], left_cb_idx, &part_idx);
    venc_lookup_pmap(t, b4_off[0], b4_off[1]-1, top_cb_idx, &part_idx);
}

// Return the context index offset of the split_cu_flag.
static finline int venc_get_split_cu_ctx_off(f265_enc_thread *t, f265_cb *cb)
{
    int left_cb_idx, top_cb_idx;
    venc_get_left_and_top_cb_idx(t, cb, &left_cb_idx, &top_cb_idx);
    return (t->cb[left_cb_idx].lg_bs<cb->lg_bs) + (t->cb[top_cb_idx].lg_bs<cb->lg_bs);
}

// Return the context index offset of the cu_skip_flag.
static finline int venc_get_cu_skip_ctx_off(f265_enc_thread *t, f265_cb *cb)
{
    int left_cb_idx, top_cb_idx;
    venc_get_left_and_top_cb_idx(t, cb, &left_cb_idx, &top_cb_idx);
    return F265_GET_FLAG(t->cb[left_cb_idx].flags, F265_CB_SKIP) + F265_GET_FLAG(t->cb[top_cb_idx].flags, F265_CB_SKIP);
}

static finline void venc_write_split_cu_flag(f265_enc_thread *t, f265_cb *cb, int flag)
{
    // This condition could be optimized with a flag eventually.
    if (cb->lg_bs == t->enc->gd.cb_range[0] || cb->flags&F265_CB_FORBIDDEN) return;
    venc_trace_syntax_element("split_cu_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_CU + venc_get_split_cu_ctx_off(t, cb), flag);
}

static finline void venc_write_cu_skip_flag(f265_enc_thread *t, f265_cb *cb, int flag)
{
    if (t->src_frame->frame_type == F265_FRAME_I) return;
    venc_trace_syntax_element("cu_skip_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_CU_SKIP + venc_get_cu_skip_ctx_off(t, cb), flag);
}

static finline void venc_write_pred_mode_flag(f265_enc_thread *t, f265_cb *cb, int flag)
{
    if (t->src_frame->frame_type == F265_FRAME_I) return;
    venc_trace_syntax_element("pred_mode_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_PRED_MODE, flag);
}

static finline void venc_write_intra_part_mode(f265_enc_thread *t, f265_cb *cb)
{
    if (cb->lg_bs != t->enc->gd.cb_range[0]) return;
    // Partitioning mode: 1, 0 (UN, HV).
    venc_trace_syntax_element("part_mode");
    venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE, cb->intra_luma_mode[1] == -1);
}

// Write the value of prev_intra_luma_pred_flag. Return the value of the next
// syntax element. If the value is less than 3, then it is mpm_idx, otherwise
// it is rem_intra_luma_pred_mode plus 3.
static finline int venc_write_intra_luma_pred(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    int luma_mode = cb->intra_luma_mode[part_idx];

    // Consider caching the MPM list.
    int mpm_list[3];
    venc_get_intra_pred_mode(t, cb, part_idx, mpm_list);

    // Check if the luma mode matches a MPM. If not, count the number of
    // MPMs that are lower than the luma mode to compute the remainder.
    int mpm_idx = 0, rem_sub_val = 0;
    for (; mpm_idx < 3 && luma_mode != mpm_list[mpm_idx]; mpm_idx++)
        if (luma_mode > mpm_list[mpm_idx])
            rem_sub_val++;

    venc_trace_syntax_element("prev_intra_luma_pred_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_INTRA_LUMA_PRED, mpm_idx<3);
    return mpm_idx < 3 ? mpm_idx : luma_mode - rem_sub_val + 3;
}

static finline void venc_write_intra_mpm_idx(f265_enc_thread *t, int mpm_idx)
{
    // FIXME: convert that table to bitfield.
    int val = f265_mpm_bin_table[mpm_idx];
    venc_trace_syntax_element("mpm_idx");
    venc_encode_bypass_bins(&t->cbs, val&3, val>>2);
}

static finline void venc_write_intra_rem_mode(f265_enc_thread *t, int rem_mode)
{
    venc_trace_syntax_element("rem_intra_luma_pred_mode");
    venc_encode_bypass_bins(&t->cbs, rem_mode, 5);
}

static finline void venc_write_intra_chroma_mode(f265_enc_thread *t, f265_cb *cb)
{
    venc_trace_syntax_element("intra_chroma_pred_mode");
    int luma_mode = cb->intra_luma_mode[0], chroma_mode = cb->intra_chroma_mode;

    // Same mode as luma (intra_chroma_pred_mode == 4 -> binarization 0).
    if (luma_mode == chroma_mode)
        venc_encode_context_bin(&t->cbs, F265_CO_INTRA_CHROMA_PRED, 0);

    // One of the common chroma modes.
    else
    {
        // Common chroma mode wanted. If the chroma mode is 34, this is the luma
        // mode, otherwise this is the chroma mode. FIXME: replace by a larger
        // table to avoid the loop?
        int mode_wanted = (chroma_mode == 34) ? luma_mode : chroma_mode;
        int mode_idx = 0;
        while (f265_chroma_mode_table[mode_idx] != mode_wanted) mode_idx++;
        venc_encode_context_bin(&t->cbs, F265_CO_INTRA_CHROMA_PRED, 1);
        venc_encode_bypass_bins(&t->cbs, mode_idx, 2);
    }
}

static finline void venc_write_inter_part_mode(f265_enc_thread *t, f265_cb *cb)
{
    int part_mode = cb->inter_part;
    int lg_bs = cb->lg_bs;
    int min_bs_flag = lg_bs == t->enc->gd.cb_range[0];
    venc_trace_syntax_element("part_mode");

    // UN (1).
    if (part_mode == F265_PART_UN)
        venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE, 1);

    // !AMP && !min_bs_flag || 8x8: 01, 00 (H2, V2).
    // min_bs_flag:                 01, 001, 000 (H2, V2, HV).
    else if (min_bs_flag | !(t->enc->gd.eflags&F265_PF_AMP))
    {
        venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE, 0);
        venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE+1, part_mode == F265_PART_H2);
        // Note: we're swapping context indices 2 and 3 specified by the spec to
        // simplify the AMP case. We're assuming HV is not possible here
        // (hardcoded "1").
        if (min_bs_flag & lg_bs != 3 & part_mode != F265_PART_H2)
            venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE+3, 1);
    }

    // AMP: 011, 001, 0100, 0101, 0000, 0001 (H2, V2, H1, H3, V1, V3).
    else
    {
        // Binarization string for the contexts. There are 4 bits per
        // partitioning mode. The first three bits are the context bin values in
        // encoding order. The fourth bit is unused. We derive the bypass bin
        // value directly from the partitioning mode.
        uint32_t ctx_str = 0x00220460;

        for (int i = 0; i < 3; i++)
            venc_encode_context_bin(&t->cbs, F265_CO_PART_MODE+i, (ctx_str>>((part_mode<<2)+i))&1);
        if (part_mode > 3) venc_encode_bypass_bin(&t->cbs, part_mode&1);
    }
}

static finline void venc_write_merge_idx(f265_enc_thread *t, f265_cb *cb, int merge_idx)
{
    int max = t->enc->gd.merge_cand - 1;
    if (max == 0) return;
    venc_trace_syntax_element("merge_idx");
    venc_encode_context_bin(&t->cbs, F265_CO_MERGE_IDX, merge_idx>0);
    if (max == 1 || !merge_idx) return;
    int bins, nb_bins;
    venc_unary_bin(merge_idx-1, max-1, &bins, &nb_bins);
    venc_encode_bypass_bins(&t->cbs, bins, nb_bins);
}

static finline void venc_write_inter_pred_idc(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    if (t->src_frame->frame_type != F265_FRAME_B) return;
    venc_trace_syntax_element("inter_pred_idc");

    // Not a 8x4 or 4x8 partition.
    if ((cb->lg_bs != 3) | (cb->inter_part == F265_PART_UN))
    {
        int cb_depth = t->enc->gd.cb_range[1] - cb->lg_bs;
        int bi_flag = (cb->ref_idx[part_idx][0] != -1) & (cb->ref_idx[part_idx][1] != -1);
        venc_encode_context_bin(&t->cbs, F265_CO_INTER_PRED+cb_depth, bi_flag);
        if (bi_flag) return;
    }

    int l1_flag = cb->ref_idx[part_idx][1] != -1;
    venc_encode_context_bin(&t->cbs, F265_CO_INTER_PRED+4, l1_flag);
}

static finline void venc_write_ref_idx(f265_enc_thread *t, f265_cb *cb, int part_idx, int list)
{
    int max = t->nb_ref_idx[list]-1, bins, nb_bins;
    int ref_idx = cb->ref_idx[part_idx][list];
    venc_unary_bin(ref_idx, max, &bins, &nb_bins);
    int nb_ctx_bins = F265_MIN(nb_bins, 2);

    // This loop unifies the case max == 0.
    if (nb_ctx_bins) venc_trace_syntax_element("ref_idx");
    for (int i = 0; i < nb_ctx_bins; i++)
        venc_encode_context_bin(&t->cbs, F265_CO_REF_IDX+i, (bins>>(nb_bins-i-1))&1);

    if (nb_bins > 2)
    {
        int nb_bypass_bins = nb_bins-2;
        venc_encode_bypass_bins(&t->cbs, bins&((1<<nb_bypass_bins)-1), nb_bypass_bins);
    }
}

static finline void venc_write_mvd(f265_enc_thread *t, f265_cb *cb, int part_idx, int list)
{
    // Integrate condition with mvd_l1_zero_flag here, assuming we ever want
    // that.

    // Consider caching the MVP.
    f265_inter_neighbour_mv neighbours[6];
    venc_get_neighbour_mvs(t, cb, part_idx, neighbours);
    f265_mv mvp_list[2];
    venc_get_pmv(t, part_idx, neighbours, cb->ref_idx[part_idx][list], list, mvp_list);
    int mvp_idx = (cb->inter_bf[part_idx]>>(3+list))&1;
    f265_mv mvp = mvp_list[mvp_idx];
    f265_mv mv = cb->mv[part_idx][list];
    int16_t mvd[2] = { mv.x - mvp.x, mv.y - mvp.y };
    int16_t abs_mvd[2] = { F265_ABS(mvd[0]), F265_ABS(mvd[1]) };

    for (int i = 0; i < 2; i++)
    {
        venc_trace_syntax_element("abs_mvd_greater0_flag");
        venc_encode_context_bin(&t->cbs, F265_CO_MVD_GREATER0, abs_mvd[i]>0);
    }

    for (int i = 0; i < 2; i++)
        if (abs_mvd[i])
        {
            venc_trace_syntax_element("abs_mvd_greater1_flag");
            venc_encode_context_bin(&t->cbs, F265_CO_MVD_GREATER1, abs_mvd[i]>1);
        }

    for (int i = 0; i < 2; i++)
    {
        if (!abs_mvd[i]) continue;
        if (abs_mvd[i] > 1)
        {
            int bins, nb_bins;
            venc_egk_bin(abs_mvd[i]-2, 1, &bins, &nb_bins);
            venc_trace_syntax_element("abs_mvd_minus2");
            venc_encode_bypass_bins(&t->cbs, bins, nb_bins);
        }
        venc_trace_syntax_element("mvd_sign_flag");
        venc_encode_bypass_bin(&t->cbs, mvd[i]<0);
    }
}

static finline void venc_write_mvp(f265_enc_thread *t, f265_cb *cb, int part_idx, int list)
{
    int flag = (cb->inter_bf[part_idx]>>(3+list))&1;
    venc_trace_syntax_element("mvp_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_MVP, flag);
}

static finline void venc_write_inter_part(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    int merge_idx = cb->inter_bf[part_idx]&7;
    venc_trace_syntax_element("merge_flag");
    venc_encode_context_bin(&t->cbs, F265_CO_MERGE_FLAG, merge_idx<5);
    if (merge_idx < 5)
    {
        venc_write_merge_idx(t, cb, merge_idx);
        return;
    }

    venc_write_inter_pred_idc(t, cb, part_idx);
    for (int list = 0; list < 2; list++)
    {
        if (cb->ref_idx[part_idx][list] == -1) continue;
        venc_write_ref_idx(t, cb, part_idx, list);
        venc_write_mvd(t, cb, part_idx, list);
        venc_write_mvp(t, cb, part_idx, list);
    }
}

// Encode the root CBF value. Return true if if the transform tree is
// non-empty.
static finline int venc_write_rqt_root_cbf(f265_enc_thread *t, f265_cb *cb)
{
    if (cb->inter_part == F265_PART_UN && (cb->inter_bf[0]&7) != 5) return 1;
    int flag = !!(t->tt.tn[0]&7);
    venc_trace_syntax_element("rqt_root_cbf");
    venc_encode_context_bin(&t->cbs, F265_CO_ROOT_CBF, flag);
    return flag;
}

// Write the coefficients in a non-empty transform block.
// FIXME: check if forced inlining of chroma_flag is beneficial.
static finline void venc_write_tb(f265_cabac_bs *cbs, f265_tt_enc *tt, int chroma_flag)
{
    f265_tb_enc *tb = tt->tb;
    f265_sb_enc *sb = tt->sb;
    int16_t *levels = tt->levels;
    int lg_bs = tb->lg_bs;
    int lg_sb = lg_bs-2;
    int nb_sb = 1<<(lg_sb<<1);
    int order = tb->order;
    uint64_t sb_nz_flags = tb->nz_flags[0];
    uint64_t sb_neighbour_flags[2] = { tb->nz_flags[1], tb->nz_flags[2] };
    const uint8_t *sb_scan = f265_scan_map_data + f265_scan_map_idx[lg_sb][order];
    const uint8_t *coeff_scan = f265_scan_map_data + f265_scan_map_idx[2][order];

    // Set the non-zero coefficient context information.
    const uint8_t *coeff_nz_ctx_idx = f265_coeff_nz_ctx_table;
    int coeff_nz_base_ctx = 0;
    if (lg_bs == 2) coeff_nz_ctx_idx += 4*16;
    else if (lg_bs == 3) coeff_nz_base_ctx += (!chroma_flag && order) ? 15 : 9;
    else coeff_nz_base_ctx += chroma_flag ? 12 : 21;

    // Intialize the last greater-than-1 context counter.
    int gt1_ctx_counter = 1;

    // Cache the sign hiding flags.
    uint64_t sign_hiding_flags = tb->sign_hiding_flags;

    // Index of the last non-zero subblock in encoding order.
    int last_sb_idx = __builtin_ctzll(sb_nz_flags);

    // Index of the last coefficient in the current subblock whose non-zero flag
    // value cannot be inferred in encoding order.
    int last_coeff_nz_idx = __builtin_ctzll(sb->nz_flags) + 1;

    // Encode the last coefficient position.
    {
        // Get the position of the last coefficient in raster scan.
        int sb_pos = sb_scan[last_sb_idx];
        int coeff_pos = coeff_scan[last_coeff_nz_idx - 1];
        int pos[2];
        pos[0] = ((sb_pos&((1<<lg_sb)-1))<<2) + (coeff_pos&3);
        pos[1] = ((sb_pos>>lg_sb)<<2) + (coeff_pos>>2);

        // Swap the order for vertical.
        if (unlikely(tb->order == 2)) F265_SWAP(int, pos[0], pos[1]);

        // Encode each position.
        int ctx_idx = F265_CO_LAST_SIG_COEFF + (chroma_flag ? 15 : 3*(lg_bs-2) + ((lg_bs-1)>>2));
        int ctx_shift = chroma_flag ? lg_bs - 2 : (lg_bs+1)>>2;
        int packed_suffix_bits = 0, packed_suffix_len = 0;
        for (int i = 0; i < 2; i++, ctx_idx += 18)
        {
            int tmp = f265_last_coeff_table[pos[i]];
            int prefix_ones = tmp&15;
            int suffix_len = tmp>>4;

            // Prefix leading ones.
            venc_trace_syntax_element(i ? "last_sig_coeff_y_prefix" : "last_sig_coeff_x_prefix");
            for (int j = 0; j < prefix_ones; j++)
                venc_encode_context_bin(cbs, ctx_idx + (j>>ctx_shift), 1);

            // Prefix trailing zero if this is not the maximum TR value for the
            // transform block size (4x4: 3, 8x8: 5, 16x16: 7, 32x32: 9).
            if (likely(prefix_ones != (lg_bs<<1) - 1))
                venc_encode_context_bin(cbs, ctx_idx + (prefix_ones>>ctx_shift), 0);

            // Pack the suffix bits.
            packed_suffix_bits <<= suffix_len;
            packed_suffix_bits |= (pos[i]&((1<<suffix_len)-1));
            packed_suffix_len += suffix_len;
        }

        // Kludge to trace correctly.
        #ifdef VAN_TRACE_SYNTAX
        for (int i = 0; i < 2; i++)
        {
            int suffix_len = f265_last_coeff_table[pos[i]]>>4;
            if (suffix_len) venc_trace_syntax_element(i ? "last_sig_coeff_y_suffix" : "last_sig_coeff_x_suffix");
            venc_encode_bypass_bins(cbs, pos[i]&((1<<suffix_len)-1), suffix_len);
        }
        #else
        // Suffix bits.
        venc_encode_bypass_bins(cbs, packed_suffix_bits, packed_suffix_len);
        #endif
    }

    // Pretend that the first subblock is non-empty to force the processing of
    // its non-zero flags.
    sb_nz_flags |= 1llu<<(nb_sb-1);

    // Loop over the subblocks in encoding order.
    for (int sb_idx = last_sb_idx, middle_sb_flag = 0;
         sb_idx < nb_sb;
         sb_idx++, middle_sb_flag = 1, last_coeff_nz_idx = 0)
    {
        // Get the subblock position in raster scan.
        int sb_pos = sb_scan[sb_idx];

        // True if the current subblock is not the first or the last.
        middle_sb_flag &= !!sb_pos;

        // Compute the values of the subblock neighbour flags for the context
        // derivation of the non-zero flags.
        int neighbour_flags[2] = { (sb_neighbour_flags[0]>>sb_pos)&1, (sb_neighbour_flags[1]>>sb_pos)&1 };
        int sb_neighbour_idx = neighbour_flags[0]|neighbour_flags[1];
        int coeff_neighbour_idx = (neighbour_flags[0]|(neighbour_flags[1]<<1))<<4;

        // Encode the subblock non-zero flag unless it is inferred to be 1.
        int sb_enc_flag = (sb_nz_flags>>sb_idx)&1;
        if (middle_sb_flag)
        {
            venc_trace_syntax_element("coded_sub_block_flag");
            venc_encode_context_bin(cbs, F265_CO_CODED_SUB_BLOCK + 2*chroma_flag + sb_neighbour_idx, sb_enc_flag);
        }

        // Skip empty subblocks.
        if (!sb_enc_flag) continue;

        // Expand the subblock data and pass to the next.
        int nz_flags = sb->nz_flags;
        int signs = sb->signs;
        int remain_flags = sb->remain_flags;
        int gt1_flags = sb->gt1_flags;
        int gt2_flag = sb->packed_data>>5;
        int nb_nz = sb->packed_data&31;
        int nb_gt1 = F265_MIN(nb_nz, 8);
        sb++;

        // Encode the coefficient non-zero flags. The non-zero flag of the first
        // coefficient is present unless the subblock non-zero flag is present
        // and the other coefficients are zero.
        int skip_first_coeff_nz_flag = middle_sb_flag & (nz_flags == (1<<15));
        for (int i = last_coeff_nz_idx; i < 16-skip_first_coeff_nz_flag; i++)
        {
            int pos = coeff_scan[i];
            int ctx = coeff_nz_base_ctx + coeff_nz_ctx_idx[coeff_neighbour_idx + pos];
            if (sb_pos) ctx += chroma_flag ? 0 : 3;
            else if (!pos) ctx = 0;
            int flag = (nz_flags>>i)&1;
            venc_trace_syntax_element("sig_coeff_flag");
            venc_encode_context_bin(cbs, F265_CO_SIG_COEFF + 27*chroma_flag + ctx, flag);
        }

        // Update the context set for the greater-than-1 flags.
        int gt1_ctx_set = (((!chroma_flag && sb_pos)<<1) + !gt1_ctx_counter)<<2;
        gt1_ctx_counter = 1;

        // Encode the greater-than-1 flags. FIXME: also count the number of
        // coefficients inferred using the greater-than-1 flags. Use conditional
        // compilation to use popcnt below instead when available.
        int nb_gt1_inferred = 0;
        for (int i = 0; i < nb_gt1; i++)
        {
            int gt1_flag = (gt1_flags>>i)&1;
            int ctx = gt1_ctx_set + gt1_ctx_counter;
            gt1_ctx_counter = f265_gt1_ctx_counter_table[(gt1_flag<<2)|gt1_ctx_counter];
            nb_gt1_inferred += !gt1_flag;
            venc_trace_syntax_element("coeff_abs_level_greater1_flag");
            venc_encode_context_bin(cbs, F265_CO_COEFF_GREATER1 + 16*chroma_flag + ctx, gt1_flag);
        }

        // Encode the greater-than-2 flag.
        if (gt1_flags)
        {
            int ctx = gt1_ctx_set>>2;
            venc_trace_syntax_element("coeff_abs_level_greater2_flag");
            venc_encode_context_bin(cbs, F265_CO_COEFF_GREATER2 + 4*chroma_flag + ctx, gt2_flag);
        }

        // Encode the signs.
        int sign_hidden_flag = (sign_hiding_flags>>sb_pos)&1;
        if (nb_nz - sign_hidden_flag) venc_trace_syntax_element("coeff_sign_flag");
        venc_encode_bypass_bins(cbs, signs >> sign_hidden_flag, nb_nz - sign_hidden_flag);

        // Encode the coefficient remaining levels.
        if (remain_flags)
        {
            // The current base level decreases as we consume the last
            // greater-than-1 and greater-than-2 flags. Since we skip the
            // coefficients already inferred, the values of those flags are
            // necessarily 1 when they are present.
            nb_gt1 -= nb_gt1_inferred + !gt2_flag;

            // Current rice parameter.
            int rice = 0;

            do
            {
                // Get the coefficient position, update the bitmap.
                int idx = __builtin_ctz(remain_flags);
                remain_flags -= (1<<idx);

                // Encode the remaining level.
                int base_level = 1 + (nb_gt1>0) + gt2_flag;
                int level = levels[idx];
                int delta = level - base_level;
                nb_gt1--; gt2_flag = 0;

                // TR max corresponding to the rice parameter.
                int tr_max = 4<<rice;

                venc_trace_syntax_element("coeff_abs_level_remaining");

                // The prefix does not begin by "1111".
                if (likely(delta < tr_max))
                {
                    // Number of bits in the prefix.
                    int len = (delta>>rice) + 1;

                    // The first part is the number of leading ones followed by
                    // the 0. The second part is the remainder.
                    // Example: "110" + "1".
                    int res = (((1<<len)-2)<<rice)|(delta&((1<<rice)-1));
                    venc_encode_bypass_bins(cbs, res, len+rice);
                }

                // The prefix begins by "1111" (exponential golumb with rice
                // parameter).
                else
                {
                    // Value to encode.
                    int v = delta - tr_max;

                    // Offset corresponding to the rice parameter.
                    int o = v + (2<<rice);

                    // Number of bits in the remainder.
                    int n = 31 - __builtin_clz(o);

                    // Value of the remainder.
                    int r = o - (1<<n);

                    // Number of bits in the prefix.
                    int len = 4 + n - rice;

                    // Same as above.
                    int res = (((1<<len)-2)<<n)|(r&((1<<n)-1));

                    venc_encode_bypass_bins(cbs, res, len+n);
                }

                // Update the rice parameter if the level busts the threshold.
                rice += level > (3<<rice) && rice != 4;

            } while (remain_flags);

            // Pass to the next subblock levels.
            levels += 16;
        }
    }

    tt->tb++;
    tt->sb = sb;
    tt->levels = levels;
}

// Encode the current TT. Used in RDOQ to update the contexts.
F265_UNUSED static void venc_encode_rdoq_sub_part(f265_enc_thread *t, int yuv_flag, int comp, int lg_bs, int depth)
{
    // Chroma CBFs.
    if (comp)
    {
        venc_trace_syntax_element(comp == 1 ? "cbf_cb" : "cbf_cr");
        venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+depth, yuv_flag);
    }

    // Luma CBFs.
    else
    {
        venc_trace_syntax_element("cbf_luma");
        venc_encode_context_bin(&t->cbs, F265_CO_CBF_LUMA+!depth, yuv_flag);
    }

    if (yuv_flag)
        venc_write_tb(&t->cbs, &t->tt, !!comp);
}

// Advance the transform node pointer until all the transform tree has been
// explored. This is used to skip empty transform trees.
F265_UNUSED static void venc_skip_tt(uint8_t **tn)
{
    if (*(*tn)++ & 8)
        for (int i = 0; i < 4; i++)
            venc_skip_tt(tn);
}

// Write the transform tree recursively. The CBF mask bits are 0 if the
// corresponding YUV components are inferred.
F265_UNUSED static void venc_write_tt(f265_enc_thread *t, f265_cb *cb, int cbf_mask, int max_depth, int lg_bs,
                                      int block_idx)
{
    int depth = cb->lg_bs-lg_bs;
    int tn_val = *t->tt.tn++;
    int split_tt_flag = !!(tn_val&8);

    // Split transform flag. FIXME: optimize this by computing in caller and
    // passing the value.
    int intra_hv_flag = F265_GET_FLAG(cb->flags, F265_CB_INTRA) && cb->intra_luma_mode[1] != -1;
    if (lg_bs > t->enc->gd.tb_range[0] &&
        lg_bs <= t->enc->gd.tb_range[1] &&
        depth < max_depth+intra_hv_flag &&
        (!intra_hv_flag || depth))
    {
        venc_trace_syntax_element("split_transform_flag");
        venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, split_tt_flag);
    }

    // Chroma CBFs.
    if (lg_bs > 2)
        for (int c = 1; c < 3; c++)
            if (cbf_mask&(1<<c))
            {
                venc_trace_syntax_element(c == 1 ? "cbf_cb" : "cbf_cr");
                venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+depth, (tn_val>>c)&1);
            }

    // Split the transform tree.
    if (split_tt_flag)
    {
        // The luma CBF is uninferred in the sub tree.
        int sub_cbf_mask = (tn_val&6)|1;
        for (int i = 0; i < 4; i++) venc_write_tt(t, cb, sub_cbf_mask, max_depth, lg_bs-1, i);
        return;
    }

    // Luma CBF.
    if (cbf_mask&1)
    {
        venc_trace_syntax_element("cbf_luma");
        venc_encode_context_bin(&t->cbs, F265_CO_CBF_LUMA+!depth, tn_val&1);
    }

    // Handle delta QP here.

    // Encode the luma transform block.
    if (tn_val&1) venc_write_tb(&t->cbs, &t->tt, 0);

    // Encode the chroma transform blocks. If the current block is 4x4, then
    // the chroma CBFs are those of the parent 8x8 block. Otherwise, the chroma
    // CBFs are those of the current block.
    int chroma_cbf = tn_val;
    if (lg_bs == 2)
    {
        if (block_idx != 3) return;
        chroma_cbf = cbf_mask;
    }
    for (int c = 1; c < 3; c++) if (chroma_cbf&(1<<c)) venc_write_tb(&t->cbs, &t->tt, 1);
}

///////////////////////////////////////////////////////////////////////////////
// Non-RDO functions.

#ifndef F265_RDO

// Encode the CB recursively.
F265_UNUSED static void venc_write_cb(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return;

    #ifdef VAN_TRACE_SYNTAX
    printf("\n");
    venc_show_loc(t, cb->cb_off[0], cb->cb_off[1], 1<<cb->lg_bs);
    #endif

    // Reset the QP written flag if the current CB is the first CB of the
    // current QG (add a flag for this in cb->flags).

    // Split CB.
    int split_cu_flag = F265_GET_FLAG(cb->flags, F265_CB_SPLIT);
    venc_write_split_cu_flag(t, cb, split_cu_flag);
    if (split_cu_flag)
    {
        for (int i = 0; i < 4; i++) venc_write_cb(t, t->cb + cb->child_idx + i);
        return;
    }

    int intra_flag = F265_GET_FLAG(cb->flags, F265_CB_INTRA);

    // Late skip test. This is required since we do not detect all skip cases in
    // the analysis.
    if (!intra_flag && cb->inter_part == F265_PART_UN && (cb->inter_bf[0]&7) != 5 && (*t->tt.tn&7) == 0)
        F265_SET_FLAG(cb->flags, F265_CB_SKIP, 1);

    // Skip CB.
    int cu_skip_flag = F265_GET_FLAG(cb->flags, F265_CB_SKIP);
    venc_write_cu_skip_flag(t, cb, cu_skip_flag);
    if (cu_skip_flag)
    {
        venc_write_merge_idx(t, cb, cb->inter_bf[0]&7);
        venc_skip_tt(&t->tt.tn);
        return;
    }

    // Handle PCM here.

    // Prediction mode.
    venc_write_pred_mode_flag(t, cb, intra_flag);

    // Encode the intra partitioning modes.
    if (intra_flag)
    {
        venc_write_intra_part_mode(t, cb);
        int nb_parts = cb->intra_luma_mode[1] == -1 ? 1 : 4;
        int mode_vals[4];
        for (int i = 0; i < nb_parts; i++) mode_vals[i] = venc_write_intra_luma_pred(t, cb, i);
        for (int i = 0; i < nb_parts; i++)
        {
            int val = mode_vals[i];
            if (val < 3) venc_write_intra_mpm_idx(t, val);
            else venc_write_intra_rem_mode(t, val-3);
        }
        venc_write_intra_chroma_mode(t, cb);
    }

    // Encode the inter partitioning modes and the root CBF.
    else
    {
        venc_write_inter_part_mode(t, cb);
        int nb_parts = f265_nb_parts[cb->inter_part];
        for (int i = 0; i < nb_parts; i++) venc_write_inter_part(t, cb, i);
        if (!venc_write_rqt_root_cbf(t, cb))
        {
            venc_skip_tt(&t->tt.tn);
            return;
        }
    }

    // Encode the transform tree. At the root, the chroma CBFs are uninferred,
    // and the luma CBF is uninferred if the CB is intra or a chroma component
    // is non-zero.
    venc_write_tt(t, cb, 6|(intra_flag|!!(t->tt.tn[0]&6)), t->enc->gd.tb_depth[!intra_flag], cb->lg_bs, 0);
}

// Write the current CTB in the bitstream.
void venc_write_ctb(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_SYNTAX
    venc_trace_syntax_flag = 1;
    #endif

    // Reset the transform tree pointers.
    venc_init_transform_tree(t);

    // Encode the SAO parameters.

    // Encode the CBs recursively.
    venc_write_cb(t, t->cb);
}
#endif

