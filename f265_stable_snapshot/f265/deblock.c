// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Deblocking filter.

#include "f265/enc.h"

// Documentation.
//
// This deblocking filter implementation processes horizontal rows of 8x8
// blocks. This model is compatible with frame multithreading, where the CTBs
// need to be deblocked one line at a time. The implementation also benefits
// from locality of reference.
//
// In the following picture, the '***' block is the current 8x8 block. All the
// labeled edges (of 2*4 pixels) except 'K' and 'L' have been deblocked. The
// filtering order requires that the vertical edge to the right of the current
// block (K) be deblocked before the top edge (L).
//
//      +-C-+-E-+-G-+
//      A   B   D   F
//      |   |   |   |
//      +-J-+-L-+---+
//      H   I***K   |
//      |   |***|   |
//      +---+---+---+
//
// There are four logical steps in the deblocking implementation. All those
// steps are skipped if the frame is not deblocked.
//
// The first step occurs during the CTB encoding. The deblocking edges of each
// 8x8 block are identified based on the layout of the prediction and transform
// trees of the coding blocks. Following the specification convention, a 8x8
// block "owns" its left and top edges. The luma non-zero TBs are also marked.
//
// The remaining steps occur over the current horizontal row of 8x8 blocks.
//
// The second step computes the QP, the unfiltered flags and the boundary
// strength of each 4-pixel edge in parallel. The values obtained are stored in
// arrays in processing order for each direction. The boundary strength is set
// to 0 for the edges that lie on or outside the frame boundaries.
//
// The third steps performs the table lookups for the Beta and tC parameters in
// parallel. The parameters are stored in arrays like in the second step.
//
// The fourth step performs the pixel filtering. The vertical edges are
// processed first, followed by the horizontal edges. Several 8x8 blocks may be
// processed in one batch. The out-of-bounds edges have a null Beta value so
// that the pixels outside the frame boundaries can be referenced safely during
// the batch processing. The frame padding ensures that the out-of-bounds
// accesses do not cause problems. The chroma edges are processed alongside the
// luma edges.
//
// By convention, we identify the horizontal/vertical directions as 0/1.
//
// NOTE:
//
// The above describes the planned final implementation. The current
// implementation is naive and slow.

#ifdef VAN_DUMP_DEBLOCK_FILTER

// Dump file.
#define VAN_DUMP_DEBLOCK_FILTER_PATH "/tmp/hevc_df.dump"
FILE *venc_df_dump_file;

// Size of the frame in 4x4 blocks.
int venc_df_b4_dim[2];

// Deblocking filter map for the frame. The first index is the direction. The
// second is the 4x4 block index in raster scan. The third index is as follow:
// deblock_flag|bs|qpl|qpc|beta|lTc|cTc|dE|dEp|dEq.
int8_t (*venc_f265_df_map[2])[10];
int8_t (*venc_hm_df_map[2])[10];

// Size of the map.
int venc_df_map_size;

// Initialize deblocking filter dump.
void venc_df_init_dump(f265_enc_thread *t)
{
    f265_gen_data *gd = &t->enc->gd;

    static int init_flag = 0;
    if (init_flag) return;
    init_flag = 1;

    for (int i = 0; i < 2; i++) venc_df_b4_dim[i] = gd->pix_dim[i]>>2;

    venc_df_map_size = 2*venc_df_b4_dim[0]*venc_df_b4_dim[1]*10;

    uint8_t *map_buf = (uint8_t*)malloc(venc_df_map_size);
    venc_f265_df_map[0] = (int8_t(*)[10])map_buf;
    venc_f265_df_map[1] = (int8_t(*)[10])(map_buf + (venc_df_map_size>>1));

    map_buf = (uint8_t*)malloc(venc_df_map_size);
    venc_hm_df_map[0] = (int8_t(*)[10])map_buf;
    venc_hm_df_map[1] = (int8_t(*)[10])(map_buf + (venc_df_map_size>>1));

    void venc_open_dump_file(FILE **f, const char *path, const char *mode);
    venc_open_dump_file(&venc_df_dump_file, VAN_DUMP_DEBLOCK_FILTER_PATH, "rb");
}

// Read from the HM dump file.
static void venc_df_read(void *dst, int size, FILE *file)
{
    if (fread(dst, size, 1, file) != 1)
    {
        printf("HM deblocking filter end-of-file.\n");
        exit(1);
    }
}

// Initialize frame dump.
void venc_df_init_frame_dump()
{
    memset(venc_f265_df_map[0], 0, venc_df_map_size);

    venc_df_read(venc_hm_df_map[0], venc_df_map_size, venc_df_dump_file);
    uint32_t sanity, expected = 0x42dec042;
    venc_df_read(&sanity, 4, venc_df_dump_file);
    if (sanity != expected)
    {
        printf("DF sanity mismatch, expected %x got %x.\n", expected, sanity);
        exit(1);
    }
}

// Finish frame dump.
void venc_df_finish_frame_dump(f265_enc_thread *t)
{
    for (int dir = 1; dir >= 0; dir--)
    {
        for (int y = 0; y < venc_df_b4_dim[1]; y++)
        {
            for (int x = 0; x < venc_df_b4_dim[0]; x++)
            {
                int nb_cmp = 10;

                // Do not consider the horizontal dE* values while the vertical
                // filtering isn't fully functional.
                #if 0
                if (!dir) nb_cmp = 7;
                #endif

                int b4_idx = y*venc_df_b4_dim[0] + x;
                int8_t *va = venc_f265_df_map[dir][b4_idx];
                int8_t *ha = venc_hm_df_map[dir][b4_idx];

                // Another kludge. Skip the edge test because HM's data is
                // unreliable.
                #if 1
                if (memcmp(va+1, ha+1, nb_cmp-1))
                #else
                if (memcmp(va, ha, nb_cmp))
                #endif
                {
                    printf("Mismatch at frame %d dir %d pix (%d,%d) 4x4 (%d, %d) b4_idx %d b8_idx %d.\n",
                           (int)t->src_frame->abs_poc, dir, x<<2, y<<2, x, y, b4_idx,
                           (y>>1)*(venc_df_b4_dim[0]>>1) + (x>>1));
                    int8_t *arrays[2] = { va, ha };
                    const char *labels[2] = { "f265:",  "hm:  " };
                    for (int i = 0; i < 2; i++)
                    {
                        printf("%s ", labels[i]);
                        for (int j = 0; j < nb_cmp; j++) printf("%d ", arrays[i][j]);
                        printf("\n");
                    }
                    exit(1);
                }
            }
        }
    }
}
#endif

// Save the deblocking filter information contained in the CB transform tree.
void venc_save_deblock_tt(f265_enc_thread *t, f265_cb *cb, int lg_bs, int cb_ox, int cb_oy)
{
    uint8_t tn_val = *t->tt.tn++;
    int split_tt_flag = !!(tn_val&8);
    int b4_stride = t->enc->gd.pix_dim[0]>>2;
    int b8_stride = t->enc->gd.pix_dim[0]>>3;

    // HM compatibility kludge. It turns out that HM *sometimes* treats a 64x64
    // inter CB with no residual as a 64x64 TB instead of 4 32x32 TBs. This
    // works because the boundary strength is 0 for the non-PB edges. When
    // debugging, we verify that we flag the same edges as HM so we emulate this
    // behavior here. Remove this eventually.
    //
    // Update after further debugging: it looks like that the behavior described
    // above isn't consistent. I'm removing the HM edge validation test because
    // the rest of the code seems to work. This condition still needs need to be
    // there to prevent non-edge mismatches with HM.
    int hm_null_inter_64_flag = !(cb->flags&F265_CB_INTRA) && lg_bs == 6 && !(tn_val&7);

    // Identify the TB edges, excluding 4x4 TBs.
    if (lg_bs == 3 || lg_bs > 3 && !split_tt_flag || hm_null_inter_64_flag)
    {
        int b8_size = 1<<(lg_bs-3);
        uint16_t *b8_map = t->src_frame->b8_map + venc_get_block_map_offset(t, cb, cb_ox, cb_oy, 3, b8_stride);
        for (int i = 0; i < b8_size; i++)
        {
            b8_map[i] |= 5<<3;
            b8_map[i*b8_stride] |= 5<<4;
        }
    }

    // Split the transform tree.
    if (split_tt_flag && !hm_null_inter_64_flag)
    {
        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
            venc_save_deblock_tt(t, cb, lg_bs-1, cb_ox + (i&1)*sbs, cb_oy + (i>>1)*sbs);
        return;
    }

    // Save the TB luma non-zero flag. Eventually, replace this code by code
    // that splatters with uint64_t (worst case is 8 stores).
    int nz_flag = tn_val&1;
    int b4_size = 1<<(lg_bs-2);
    uint8_t *b4_map = t->src_frame->b4_map + venc_get_block_map_offset(t, cb, cb_ox, cb_oy, 2, b4_stride);
    for (int i = 0; i < b4_size; i++)
        for (int j = 0; j < b4_size; j++)
            F265_SET_FLAG(b4_map[i*b4_stride + j], 1<<6, nz_flag);
}

// Save the deblocking filter information contained in the CB partition tree.
void venc_save_deblock_pb(f265_enc_thread *t, f265_cb *cb)
{
    // Intra partitions and the first inter partition are already marked via the
    // transform tree. We assume inter HV isn't allowed. Thus, only the second
    // inter partition needs to be marked. Since the deblocking filter only
    // considers edges on the 8x8 grid, we skip partitions unaligned on 8x8.
    if ((cb->flags&F265_CB_INTRA) || cb->inter_part == F265_PART_UN) return;

    int part_size[2], part_off[2];
    venc_get_part_loc(part_size, part_off, cb->inter_part, 1, 1<<cb->lg_bs);

    // Edge direction. 0xc4: 1100 0100.
    int dir = (0xc4>>cb->inter_part)&1;

    // Unaligned edge.
    if (part_size[!dir]&4) return;

    // Mark the edges.
    int b8_size = part_size[dir]>>3;
    int b8_stride = t->enc->gd.pix_dim[0]>>3;
    uint16_t *b8_map = t->src_frame->b8_map + venc_get_block_map_offset(t, cb, part_off[0], part_off[1], 3, b8_stride);
    if (dir) for (int i = 0; i < b8_size; i++) b8_map[i*b8_stride] |= 1<<4;
    else for (int i = 0; i < b8_size; i++) b8_map[i] |= 1<<3;
}

// Build the deblocking frame identification map.
void venc_build_deblock_frame_id_map(f265_enc_thread *t)
{
    f265_frame *f = t->src_frame;

    // Map the unused reference index (-1) to the last position (16).
    t->deblock_frame_id_map[0][0] = t->deblock_frame_id_map[1][0] = 16;

    // Map the reference indices naively (this is fast enough).
    for (int list = 0; list < t->nb_lists; list++)
        for (int ref_idx = 0; ref_idx < t->nb_ref_idx[list]; ref_idx++)
            for (int frame_idx = 0; frame_idx < f->dpb_size; frame_idx++)
                if (f->dpb[frame_idx] == t->ref_ctx[list][ref_idx].frame)
                {
                    t->deblock_frame_id_map[list][1+ref_idx] = frame_idx;
                    break;
                }
}

// Return true if the specified pair of motion vectors is distant.
static int venc_is_distant_mv_pair(f265_mv a, f265_mv b)
{
    for (int i = 0; i < 2; i++) if (F265_ABS(a.v[i]-b.v[i]) >= 4) return 1;
    return 0;
}

// Return the boundary strength associated to the current edge.
static int venc_get_boundary_strength(f265_enc_thread *t, int b4_idx[2], int b8_idx[2], int dir)
{
    f265_frame *f = t->src_frame;
    int b8_vals[2], b4_vals[2];
    for (int i = 0; i < 2; i++)
    {
        b4_vals[i] = f->b4_map[b4_idx[i]];
        b8_vals[i] = f->b8_map[b8_idx[i]];
    }

    // Intra CB, BS=2.
    if ((b8_vals[0]|b8_vals[1])&1) return 2;

    // Transform block boundary and a luma TB is non-zero, BS=1.
    if (b8_vals[1]&(1<<(5+dir)) && (b4_vals[0]|b4_vals[1])&(1<<6)) return 1;

    // Both blocks are inter. Map their reference indices to frame identifiers
    // and fetch their motion vectors.
    int frame_ids[2][2];
    f265_mv mvs[2][2];
    for (int i = 0; i < 2; i++)
        for (int list = 0; list < 2; list++)
        {
            frame_ids[i][list] = t->deblock_frame_id_map[list][1+f->ref_idx[b4_idx[i]][list]];
            mvs[i][list] = f->mv[b4_idx[i]][list];
        }

    // Sort the frame identifiers and motion vectors.
    for (int i = 0; i < 2; i++)
        if (frame_ids[i][0] > frame_ids[i][1])
        {
            F265_SWAP(int, frame_ids[i][0], frame_ids[i][1]);
            F265_SWAP(f265_mv, mvs[i][0], mvs[i][1]);
        }

    // If the frame identifiers mismatch, BS=1.
    if (memcmp(frame_ids[0], frame_ids[1], 2*sizeof(int))) return 1;

    // If the motion vectors are distant, then BS=1, otherwise BS=0.

    // Bi-prediction on one reference frame.
    if (frame_ids[0][0] == frame_ids[0][1])
    {
        return (venc_is_distant_mv_pair(mvs[0][0], mvs[1][0]) || venc_is_distant_mv_pair(mvs[0][1], mvs[1][1])) &&
               (venc_is_distant_mv_pair(mvs[0][0], mvs[1][1]) || venc_is_distant_mv_pair(mvs[0][1], mvs[1][0]));
    }

    // Single prediction or bi-prediction on different frames.
    return venc_is_distant_mv_pair(mvs[0][0], mvs[1][0]) ||
           (frame_ids[0][1] != 16 && venc_is_distant_mv_pair(mvs[0][1], mvs[1][1]));
}

// Reverse the 4x4 block vertically in place.
//      [0 1 2 3]    [3 2 1 0]
//      [4 5 6 7] => [7 6 5 4]
//      [8 9 A B]    [B A 9 8]
//      [C D E F]    [F E D C]
static void venc_reverse_4x4_block(f265_pix pix[4*4])
{
    f265_pix tmp[4*4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            tmp[i*4 + j] = pix[i*4 + 3-j];
    memcpy(pix, tmp, sizeof(tmp));
}

// Load the 4x4 source pixels on each side of the edge. For simplicity, we
// tranpose the horizontal case to the vertical case, and we reverse the
// pixels vertically in the 'P' block so that indexing is the same as the
// 'Q' block.
static void venc_load_edge_block(f265_pix pix[2][4*4], int b_pos[2][2], f265_pix *rec_plane, int src_stride, int dir)
{
    for (int i = 0; i < 2; i++)
    {
        f265_pix *src = rec_plane + b_pos[i][1]*src_stride + b_pos[i][0];
        if (dir) venc_copy_block(pix[i], 4, src, src_stride, 4, 4);
        else venc_transpose_block(pix[i], 4, src, src_stride, 4, 4);
    }
    venc_reverse_4x4_block(pix[0]);
}

// Store the deblocked pixels on each side of the edge.
static void venc_store_edge_block(f265_pix pix[2][4*4], int b_pos[2][2], f265_pix *rec_plane, int src_stride, int dir)
{
    venc_reverse_4x4_block(pix[0]);
    for (int i = 0; i < 2; i++)
    {
        f265_pix *src = rec_plane + b_pos[i][1]*src_stride + b_pos[i][0];
        if (dir) venc_copy_block(src, src_stride, pix[i], 4, 4, 4);
        else venc_transpose_block(src, src_stride, pix[i], 4, 4, 4);
    }
}

// Compute the luma edge filtering parameters (in order: dE, dEp, dEq).
static void venc_compute_luma_edge_params(f265_enc_thread *t, f265_pix pix[2][4*4], int beta, int tc, int params[3])
{
    // Straight from the specification, but the X,Y indices are reversed.
    f265_pix (*p)[4] = (f265_pix(*)[4])pix[0], (*q)[4] = (f265_pix(*)[4])pix[1];
    int dp0 = F265_ABS(p[0][2] - 2*p[0][1] + p[0][0]);
    int dp3 = F265_ABS(p[3][2] - 2*p[3][1] + p[3][0]);
    int dq0 = F265_ABS(q[0][2] - 2*q[0][1] + q[0][0]);
    int dq3 = F265_ABS(q[3][2] - 2*q[3][1] + q[3][0]);
    int dpq0 = dp0 + dq0;
    int dpq3 = dp3 + dq3;
    int dp = dp0 + dp3;
    int dq = dq0 + dq3;
    int d = dpq0 + dpq3;

    if (d >= beta)
    {
        params[0] = params[1] = params[2] = 0;
        return;
    }

    // Get dE.
    params[0] = 2;
    int dpq_array[2] = { dpq0<<1, dpq3<<1 };
    for (int i = 0; i < 2; i++)
    {
        int dpq = dpq_array[i];
        int r = i*3;
        int decision = dpq < (beta>>2) &&
                       F265_ABS(p[r][3] - p[r][0]) + F265_ABS(q[r][0] - q[r][3]) < (beta>>3) &&
                       F265_ABS(p[r][0] - q[r][0]) < ((5*tc + 1)>>1);
        if (!decision)
        {
            params[0] = 1;
            break;
        }
    }

    // Get dEp/dEq.
    int d_array[2] = { dp, dq };
    for (int i = 0; i < 2; i++) params[1+i] = d_array[i] < ((beta + (beta>>1)) >> 3);
}

// Deblock a line of pixels of a luma edge in place.
static void venc_deblock_luma_line(f265_pix *lp, f265_pix *lq, int params[3], int tc, int uf)
{
    // Fix the bit depth issues eventually.
    int max_pix = 255;

    // Compute the value of Delta1 (cannot be unified below due to the shift).
    int delta1 = (9*(lq[0] - lp[0]) - 3*(lq[1] - lp[1]) + 8) >> 4;

    // Get the filtered pixels. For simplicity, we translate the 'Q' case into
    // the 'P' case.
    f265_pix *lines[2] = { lp, lq };
    f265_pix filtered[2][3];
    int nb_filtered[2] = { 0, 0 };
    for (int i = 0; i < 2; i++)
    {
        // Skip the unfiltered edge.
        if ((uf>>i)&1) continue;

        f265_pix *p = lines[i], *q = lines[!i], *f = filtered[i];
        int *n = nb_filtered + i;

        // Strong filtering.
        if (params[0] == 2)
        {
            f[0] = F265_CLAMP((p[2] + 2*p[1] + 2*p[0] + 2*q[0] + q[1] + 4) >> 3, p[0]-2*tc, p[0]+2*tc);
            f[1] = F265_CLAMP((p[2] + p[1] + p[0] + q[0] + 2) >> 2, p[1]-2*tc, p[1]+2*tc);
            f[2] = F265_CLAMP((2*p[3] + 3*p[2] + p[1] + p[0] + q[0] + 4) >> 3, p[2]-2*tc, p[2]+2*tc);
            *n = 3;
        }

        // Weak filtering.
        else
        {
            // Skipped.
            if (F265_ABS(delta1) >= 10*tc) continue;

            // 1 pixel filtered.
            int delta2 = F265_CLAMP(delta1, -tc, tc);
            if (i) delta2 = -delta2;
            f[0] = F265_CLAMP(p[0] + delta2, 0, max_pix);
            *n = 1;

            // 2 pixels filtered.
            if (params[1+i])
            {
                int dp = F265_CLAMP((((p[2]+p[0]+1) >> 1) - p[1] + delta2) >> 1, -(tc>>1), tc>>1);
                f[1] = F265_CLAMP(p[1] + dp, 0, max_pix);
                *n = 2;
            }
        }
    }

    // Replace the filtered pixels.
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < nb_filtered[i]; j++)
            lines[i][j] = filtered[i][j];
}

// Deblock a luma edge.
static void venc_deblock_luma_edge(f265_enc_thread *t, int dir, int b_pos[2][2], int beta, int tc, int uf)
{
    // Load the block.
    int src_stride = t->me.ref_stride;
    f265_pix *rec_plane = t->src_frame->rec_planes[0];
    f265_pix pix[2][4*4];
    venc_load_edge_block(pix, b_pos, rec_plane, src_stride, dir);

    // Get the edge parameters. Bail out if dE is null.
    int params[3];
    venc_compute_luma_edge_params(t, pix, beta, tc, params);
    if (!params[0]) return;
    #ifdef VAN_DUMP_DEBLOCK_FILTER
    {
        int b4_idx = (b_pos[1][1]>>2)*venc_df_b4_dim[0] + (b_pos[1][0]>>2);
        venc_f265_df_map[dir][b4_idx][7] = params[0];
        venc_f265_df_map[dir][b4_idx][8] = params[1];
        venc_f265_df_map[dir][b4_idx][9] = params[2];
    }
    #endif

    // Filter each line.
    for (int i = 0; i < 4; i++) venc_deblock_luma_line(pix[0] + 4*i, pix[1] + 4*i, params, tc, uf);

    // Store the block.
    venc_store_edge_block(pix, b_pos, rec_plane, src_stride, dir);
}

// Deblock a chroma edge.
static void venc_deblock_chroma_edge(f265_enc_thread *t, int dir, int b_pos[2][2], int tc, int uf)
{
    // Fix the bit depth issues eventually.
    int max_pix = 255;
    int src_stride = t->me.ref_stride;

    // Process each component.
    for (int comp = 1; comp < 3; comp++)
    {
        // Load the block.
        f265_pix *rec_plane = t->src_frame->rec_planes[3+comp];
        f265_pix pix[2][4*4];
        venc_load_edge_block(pix, b_pos, rec_plane, src_stride, dir);

        // Filter each line.
        for (int i = 0; i < 4; i++)
        {
            f265_pix *p = pix[0] + 4*i, *q = pix[1] + 4*i;
            int delta = F265_CLAMP((((q[0] - p[0]) << 2) + p[1] - q[1] + 4) >> 3, -tc, tc);
            if (!(uf&1)) p[0] = F265_CLAMP(p[0] + delta, 0, max_pix);
            if (!(uf&2)) q[0] = F265_CLAMP(q[0] - delta, 0, max_pix);
        }

        // Store the block.
        venc_store_edge_block(pix, b_pos, rec_plane, src_stride, dir);
    }
}

// Deblock the edge at the specified position.
static void venc_deblock_edge(f265_enc_thread *t, int dir, int px, int py)
{
    // Compute the position and the index on the 4x4/8x8 grid of the two blocks
    // on each side of the edge. The second entry in the arrays refers to the
    // block associated to the edge (the right/bottom block, referred to as 'q'
    // in the specification).
    f265_gen_data *gd = &t->enc->gd;
    f265_frame *f = t->src_frame;
    int b4_stride = gd->pix_dim[0]>>2;
    int b8_stride = gd->pix_dim[0]>>3;
    int b_pos[2][2] = { { px - (dir<<2), py - (!dir<<2) }, { px, py } };
    int b_pos_c[2][2] = { { (px>>1) - (dir<<2), (py>>1) - (!dir<<2) }, { px>>1, py>>1 } };
    int b4_idx[2], b8_idx[2];
    for (int i = 0; i < 2; i++)
    {
        b4_idx[i] = (b_pos[i][1]>>2)*b4_stride + (b_pos[i][0]>>2);
        b8_idx[i] = (b_pos[i][1]>>3)*b8_stride + (b_pos[i][0]>>3);
    }

    // Bail out if this is not a deblocking edge.
    if (!(f->b8_map[b8_idx[1]]&(1<<(3+dir)))) return;
    #ifdef VAN_DUMP_DEBLOCK_FILTER
    venc_f265_df_map[dir][b4_idx[1]][0] = 1;
    #endif

    // Compute the boundary strength. Bail out if zero.
    int bs = venc_get_boundary_strength(t, b4_idx, b8_idx, dir);
    if (!bs) return;

    // Compute whether this is a deblocked chroma edge.
    int chroma_edge_flag = dir ? !(px&8) && !(py&4) : !(py&8) && !(px&4);

    // Extract the unfiltered flags.
    int uf = 0;
    for (int i = 0; i < 2; i++) uf |= ((f->b8_map[b8_idx[i]]>>7)&1)<<i;

    // Compute the average luma QP.
    int qpl = (f->qmap[b8_idx[0]] + f->qmap[b8_idx[1]] + 1)>>1;

    // Fix the bit depth issues below eventually.

    // Compute Beta. We avoid the table lookup with range tests.
    int beta_q = F265_CLAMP(qpl + gd->deblock_off[0], 0, 51);
    int beta = beta_q<16 ? 0 : beta_q<29 ? beta_q-10 : (beta_q<<1)-38;

    // Compute luma tC. The table is too irregular for an arithmetic solution,
    // so table lookups seem unavoidable. Assembly: use VGATHERPS so we can at
    // least do 8 lookups somewhat in parallel.
    int luma_tc_q = F265_CLAMP(qpl + ((bs-1)<<1) + gd->deblock_off[1], 0, 53);
    int luma_tc = f265_tc_table[luma_tc_q];
    #ifdef VAN_DUMP_DEBLOCK_FILTER
    venc_f265_df_map[dir][b4_idx[1]][1] = bs;
    venc_f265_df_map[dir][b4_idx[1]][2] = qpl;
    venc_f265_df_map[dir][b4_idx[1]][4] = beta;
    venc_f265_df_map[dir][b4_idx[1]][5] = luma_tc;
    #endif

    // Process the luma edge.
    venc_deblock_luma_edge(t, dir, b_pos, beta, luma_tc, uf);

    // Process the chroma edge.
    if (chroma_edge_flag && bs == 2)
    {
        // Compute chroma tC. FIXME: analyze if the QP table index can overflow
        // here.
        int qpc = gd->chroma_qp_table[qpl + gd->chroma_qp_idx_off];
        int chroma_tc_q = F265_CLAMP(qpc + ((bs-1)<<1) + gd->deblock_off[1], 0, 53);
        int chroma_tc = f265_tc_table[chroma_tc_q];
        #ifdef VAN_DUMP_DEBLOCK_FILTER
        venc_f265_df_map[dir][b4_idx[1]][3] = qpc;
        venc_f265_df_map[dir][b4_idx[1]][6] = chroma_tc;
        #endif

        venc_deblock_chroma_edge(t, dir, b_pos_c, chroma_tc, uf);
    }
}

// Deblock the whole frame in one direction.
static void venc_deblock_frame_dir(f265_enc_thread *t, int dir)
{
    // Deblock every line on the 8x8 grid except the first and the last, and
    // every 4-pixel edge along the current line.
    int32_t pix_dim[2] = { t->enc->gd.pix_dim[0], t->enc->gd.pix_dim[1] };
    for (int i = 8; i < pix_dim[!dir]; i += 8)
        for (int j = 0; j < pix_dim[dir]; j += 4)
            venc_deblock_edge(t, dir, dir ? i : j, dir ? j : i);
}

// Deblock the whole frame.
void venc_deblock_frame(f265_enc_thread *t)
{
    // Build the deblocking frame identification map.
    venc_build_deblock_frame_id_map(t);

    // Deblock the vertical/horizontal edges.
    for (int dir = 1; dir >= 0; dir--) venc_deblock_frame_dir(t, dir);
    #ifdef VAN_DUMP_DEBLOCK_FILTER
    venc_df_finish_frame_dump(t);
    #endif
}

