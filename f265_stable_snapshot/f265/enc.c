// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Miscellaneous stuff that doesn't fit elsewhere:
// - Debugging functions.
// - Initialization.
// - High-level logic.

#include "f265/enc.h"

f265_large_tables f265_lt;

static int venc_encode_frame_mt_none(f265_enc *e, f265_frame *f);
static int venc_encode_section(f265_enc_thread *t, int section_idx);
static void venc_encode_ctb(f265_enc_thread *t);
static f265_enc_thread* venc_set_enc_thread(f265_enc *e, f265_frame *f);
static int venc_check_enc_result(f265_enc_thread *t, f265_frame *next_frame);

// Dump the YUV frame specified. Silent failures.
#ifdef F265_HAVE_STDIO
void venc_dump_yuv_file(FILE *file, f265_pix *planes[3], int64_t frame_poc, int stride, int width, int height,
                        int sx, int sy)
{
    int64_t nb_pix = width*height;
    int64_t frame_size = nb_pix + ((nb_pix<<1)>>(sx+sy));
    int64_t seek_pos = frame_poc*frame_size;

    if (fseek(file, seek_pos, SEEK_SET)) return;
    for (int i = 0; i < 3; i++)
    {
        int w = i ? width>>sx : width;
        int h = i ? height>>sy : height;
        f265_pix *p = planes[i];
        for (int y = 0; y < h; y++, p += stride) fwrite(p, w*F265_PSIZE, 1, file);
    }
}
#endif


// Temporary object to hold the encoder initialization data.
typedef struct f265_enc_mem_data
{
    // Chroma scale factor in the X,Y directions (0 or 1).
    int64_t csf[2];

    // Frame width and height in terms of pixels (multiple of CB size).
    int64_t pix_dim[2];

    // CTB block counts.
    int64_t ctb_size;
    int64_t ctb_dim[2];
    int64_t nb_ctb;

    // True if subpel motion estimation is used.
    int64_t subpel_me_flag;

    // Number of sections in the frame.
    int64_t nb_sections;

    // Total lookahead delay.
    int64_t la_delay;

    // Number of encoding and lookahead thread objects.
    int64_t nb_enc_threads;
    int64_t nb_la_threads;

    // Number of frame threads. 0 if frame multithreading is not used.
    int64_t nb_frame_threads;

    // Number of frames encoded in parallel.
    int64_t nb_enc_frames;

    // YUV planes stride, in pixels.
    int64_t stride;

    // Size of the padded luma/chroma plane buffers.
    int64_t plane_size[2];

    // Offset of the first non-padded pixel of the plane from the luma/chroma
    // plane buffers (YUV).
    int64_t plane_off[3];

    // Number of reconstructed luma planes.
    int64_t nb_rec_luma;

    // Size allocated per CTB in the chunk buffer.
    int64_t ctb_chunk_size;

    // Size of the external bitstream.
    int64_t bs_mem;

    // Object counts.
    int64_t nb_frame_objs;
    int64_t nb_src_objs;
    int64_t nb_la_objs;
    int64_t nb_enc_objs;
    int64_t nb_dyn_objs; // Frame, source, lookahead, encoding.
    int64_t nb_deblock_objs;
    int64_t nb_chunk_objs;

    // Object sizes.
    int64_t src_obj_size;
    int64_t la_obj_size;
    int64_t enc_obj_size;
    int64_t deblock_obj_size;
    int64_t chunk_obj_size;
    int64_t fmap_obj_size;

    // Object pointers.
    f265_enc *enc;
    f265_enc_thread *enc_threads;
    f265_frame *frame_objs;
    uint8_t *src_objs;
    uint8_t *la_objs;
    uint8_t *enc_objs;
    uint8_t *deblock_objs;
    uint8_t *chunk_objs;
    uint8_t *me_buf_objs;
    f265_frame_map *fmap;
    f265_enc_thread** md_enc_threads;
    f265_frame** md_enc_frames;
    uint8_t** md_unused_objs;
    f265_frame** la_display;
    f265_frame** la_coded;
    f265_frame** la_queues;

} f265_enc_mem_data;

// Analyze the encoder memory requirements.
//
// The size of an object has the same alignment as the greatest alignment of its
// fields. Thus, the size is aligned on a cache line boundary if the object
// contains a field aligned on a cache line boundary.
//
// For simplicity, we allocate large objects contiguously on cache line
// boundaries to avoid false sharing and alignment issues.
static void venc_analyze_enc_mem(f265_enc_params *p, f265_enc_mem_data *d, uint8_t *buf)
{
    // Set subpel motion estimation tentatively. FIXME.
    d->subpel_me_flag = 1;

    // Maximum number of slice headers in the frame. FIXME.
    int64_t nb_slice_headers = 1;

    // number of sections in the frame. FIXME.
    d->nb_sections = 1;


    // Chroma scaling factors (X shift, Y shift, YUV size numerator, YUV size
    // denominator shift).
    int csf_table[4][4] = { {0,0,1,0}, {1,1,3,1}, {1,0,2,0}, {0,0,3,0} };
    d->csf[0] = csf_table[p->chroma_format][0];
    d->csf[1] = csf_table[p->chroma_format][1];
    int csf_num = csf_table[p->chroma_format][2];
    int csf_shift = csf_table[p->chroma_format][3];


    // Frame and CTB dimensions.
    int min_cb_lg = p->cb_range[0];
    int max_cb_lg = p->cb_range[1];
    d->ctb_size = 1<<max_cb_lg;
    for (int i = 0; i < 2; i++)
    {
        d->pix_dim[i] = ((p->clip_dim[i] + (1<<min_cb_lg) - 1) >> min_cb_lg) << min_cb_lg;
        d->ctb_dim[i] = (d->pix_dim[i] + d->ctb_size - 1) >> max_cb_lg;
    }
    d->nb_ctb = d->ctb_dim[0]*d->ctb_dim[1];

    // Frame 4x4/8x8 block dimensions.
    int64_t b4_dim[2] = { d->pix_dim[0]>>2, d->pix_dim[1]>>2 };
    int64_t nb_b4 = b4_dim[0]*b4_dim[1];
    int64_t nb_b8 = nb_b4>>2;


    // Total lookahead delay.
    d->la_delay = p->la_decision_delay + p->la_process_delay;

    // Threading.
    d->nb_enc_threads = p->nb_workers[0] + (p->mt_mode != F265_MT_ENC_FRAME);
    // There is one lookahead thread object even if the lookahead is not used.
    d->nb_la_threads = F265_MAX(p->nb_workers[1], 1);
    d->nb_frame_threads = (p->mt_mode == F265_MT_ENC_FRAME) ? d->nb_enc_threads : 0;
    d->nb_enc_frames = (p->mt_mode == F265_MT_ENC_FRAME) ? d->nb_enc_threads : 1;


    // Object counts and sizes.

    // The next frame to encode uses one object. There are la_delay frames in
    // the lookahead after analysis, plus one unification frame if it's used.
    // Each reference frame uses one object. Each frame thread uses one object.
    // The last encoded frame uses one object.
    //
    // Be careful with assumptions. The last encoded frame is not necessarily
    // used for reference. The lookahead unification frame is not necessarily
    // the last encoded frame or used for reference. Since we explicitly keep an
    // object for the last encoded frame, we cover the case where the DPB
    // contains more frames than the number of reference frames before the RPS
    // is applied.
    d->nb_frame_objs = 1 + d->la_delay + p->la_flag + p->nb_refs + d->nb_frame_threads + 1;

    // The next frame to encode uses one object. There are la_delay frames in
    // the lookahead after analysis, plus one unification frame if it's used.
    // Each frame thread uses one object.
    d->nb_src_objs = 1 + d->la_delay + p->la_flag + d->nb_frame_threads;
    d->nb_la_objs = d->nb_src_objs;

    // FIXME: the first frame thread does NOT consume an encode object. Revisit
    // code below when implementing frame multithreading.

    // The next frame to encode uses one object. Each reference frame uses one
    // object. Each frame thread uses one object.
    d->nb_enc_objs = 1 + p->nb_refs + d->nb_frame_threads;
    d->nb_dyn_objs = d->nb_frame_objs + d->nb_src_objs + d->nb_la_objs + d->nb_enc_objs;

    // Each frame encoded in parallel uses one object.
    d->nb_chunk_objs = d->nb_enc_frames;
    d->nb_deblock_objs = d->nb_enc_frames;

    // YUV planes stride, aligned on a cache line boundary.
    d->stride = F265_ALIGN_VAL((d->pix_dim[0] + 2*F265_LUMA_PLANE_PADDING)*F265_PSIZE, 64) >> (F265_PSIZE-1);

    // Luma and chroma plane size. The extra cache line is needed to align the
    // planes.
    d->plane_size[0] = d->stride*(d->pix_dim[1] + 2*F265_LUMA_PLANE_PADDING)*F265_PSIZE + 64;
    d->plane_size[1] = p->chroma_format ? (d->plane_size[0]<<1)>>(d->csf[0]+d->csf[1]) : 0;

    // Offset of the first non-padded pixel in each plane.
    d->plane_off[0] = F265_ALIGN_VAL((d->stride+1)*F265_LUMA_PLANE_PADDING*F265_PSIZE, 64);
    if (p->chroma_format == 3)
    {
        d->plane_off[1] = d->plane_off[0];
        d->plane_off[2] = d->plane_size[0] + d->plane_off[0];
    }
    else
    {
        if (p->chroma_format == 1) d->plane_off[1] = d->plane_off[0]>>1;
        else d->plane_off[1] = F265_ALIGN_VAL((d->stride*F265_LUMA_PLANE_PADDING +
                                               (F265_LUMA_PLANE_PADDING>>1))*F265_PSIZE, 32);
        d->plane_off[2] = d->plane_off[1] + ((d->stride*F265_PSIZE)>>1);
    }

    // Source object size. Luma + chroma.
    d->src_obj_size = d->plane_size[0] + d->plane_size[1];

    // Lookahead object size.
    d->la_obj_size = 0;

    // Number of reconstructed luma planes.
    d->nb_rec_luma = 1 + 3*d->subpel_me_flag;

    // Encoding object size.
    d->enc_obj_size = F265_ALIGN_VAL(// Cache line aligned.
                                     d->nb_rec_luma*d->plane_size[0] +  // Luma planes.
                                     d->plane_size[1] +                 // Chroma planes.
                                     // 32-bit aligned.
                                     nb_b4*4*2 +                        // mv.
                                     // 16-bit aligned.
                                     nb_b8*2 +                          // b8_map.
                                     // Byte aligned.
                                     d->nb_ctb*sizeof(f265_ctb_sao) +   // sao.
                                     nb_b8*1 +                          // qmap.
                                     nb_b4*(2+1),                       // ref_idx, b4_map.
                                     64);

    // Deblocking context object size.
    d->deblock_obj_size = 0;

    // Chunk buffer and external bitstream size.
    {
        // Miscellaneous frame NAL overhead (AU delimiter, VPS, SPS, PPS). The
        // scaling list data takes the most space (~6K). We assume 1 NAL unit of
        // each kind. 16K ought to be plenty.
        int64_t frame_nal_size = 16384;

        // NAL overhead per CTB.
        // - Assuming every CTB has its own segment (without a slice header).
        // - NAL formatting: 6 bytes.
        // - Dependent segment: 15 bytes.
        //   - Header: address 4 bytes, flags and PPS 4 bytes.
        //   - Flush: 2 bytes.
        //   - Total escaped: (4+4+2)*3/2 = 15.
        int64_t ctb_nal_size = 21;

        // Slice header overhead.
        // - Address and entry points: shadowed by the dependent segments.
        // - Miscellaneous stuff: 50 bytes.
        // - 16 reference frames:
        //   - RPS: 10 bytes for the POC and the flags.
        //   - Weighting: 3*(2+2+2) bytes for the weights, offsets and flags.
        // - Total escaped: (50 + 16*(10 + 3*(2+2+2))) * 3/2 = 747.
        int64_t slice_header_size = 747;

        // Number of 4x4/8x8 blocks per CTB.
        int64_t b4_per_ctb = 1<<((max_cb_lg-2)<<1);
        int64_t b8_per_ctb = 1<<((max_cb_lg-3)<<1);

        // Number of 4x4 transform blocks per CTB (all image components).
        int64_t tb_per_ctb = (b4_per_ctb*csf_num)>>csf_shift;

        // Count the maximum number of context and bypass bins per CTB.
        int64_t context_bins = 0;
        int64_t bypass_bins = 0;

        // SAO.
        context_bins += 5;
        bypass_bins += 126;

        // We can use 8x8 PCM CBs without PBs and QPs in the worst case. We add
        // 24 bypass bins per 8x8 CB for CABAC synchronization and mode
        // encoding.
        if (p->pcm_range[0] != -1)
        {
            bypass_bins += b8_per_ctb*24 + tb_per_ctb*16*p->bit_depth[2];
        }

        // Same as the CTB CABAC buffer, with the appropriate object counts and
        // without the correction bits.
        else
        {
            context_bins += b8_per_ctb*11 + b4_per_ctb*23 + tb_per_ctb*34;
            bypass_bins  += b8_per_ctb*14 + b4_per_ctb*34 + tb_per_ctb*16*35;
        }

        // Account for escaping and the worst case for the context bins.
        int64_t data_per_ctb = ((context_bins*3/2 + bypass_bins)*3/2 + 7)/8;

        // Account for the segmentation and the chunk data.
        d->ctb_chunk_size = 5*4 + data_per_ctb;

        // Chunk buffer size.
        d->chunk_obj_size = F265_ALIGN_VAL(d->nb_ctb*d->ctb_chunk_size, 64);

        // External bitstream size.
        d->bs_mem = frame_nal_size + nb_slice_headers*slice_header_size + d->nb_ctb*(ctb_nal_size + data_per_ctb);
    }

    // Frame map.
    d->fmap_obj_size = sizeof(f265_frame_map) + d->nb_sections*4*2 + d->nb_ctb*sizeof(f265_frame_ctb);


    // Memory reservation/allocation.

    // Align the buffer, if any.
    uint8_t *b = (uint8_t*)F265_ALIGN_VAL((uint64_t)buf, 64);

    // Track the current allocation offset.
    int64_t b_off = 0;

    // Set the pointer P of type T to a memory area of size S aligned on A.
    // Update the offset counter with S. S needs not be aligned on A.
    #define SET_OFF(T, P, S, A) b_off = F265_ALIGN_VAL(b_off, (A)); d->P = (T*)(b + b_off); b_off += (S)

    // Flush the current cache line.
    #define FLUSH_CL() b_off = F265_ALIGN_VAL(b_off, 64)

    // Track the offset and size of each object group.
    SET_OFF(f265_enc, enc, sizeof(f265_enc), 64);
    SET_OFF(f265_enc_thread, enc_threads, d->nb_enc_threads*sizeof(f265_enc_thread), 64);
    SET_OFF(f265_frame, frame_objs, d->nb_frame_objs*sizeof(f265_frame), 64);
    SET_OFF(uint8_t, src_objs, d->nb_src_objs*d->src_obj_size, 64);
    SET_OFF(uint8_t, enc_objs, d->nb_enc_objs*d->enc_obj_size, 64);
    SET_OFF(uint8_t, deblock_objs, d->nb_deblock_objs*d->deblock_obj_size, 64);
    SET_OFF(uint8_t, chunk_objs, d->nb_chunk_objs*d->chunk_obj_size, 64);
    SET_OFF(f265_frame_map, fmap, d->fmap_obj_size, 64);
    FLUSH_CL();
    SET_OFF(f265_enc_thread*, md_enc_threads, d->nb_enc_threads*sizeof(void*), sizeof(void*));
    SET_OFF(f265_frame*, md_enc_frames, (d->nb_enc_frames+1)*sizeof(void*), sizeof(void*));
    SET_OFF(uint8_t*, md_unused_objs, d->nb_dyn_objs*sizeof(void*), sizeof(void*));
    FLUSH_CL();
    SET_OFF(f265_frame*, la_display, (1+p->la_decision_delay+1)*sizeof(void*), sizeof(void*));
    SET_OFF(f265_frame*, la_coded, (p->la_decision_delay+1)*sizeof(void*), sizeof(void*));
    SET_OFF(f265_frame*, la_queues, 2*(1+p->la_process_delay)*sizeof(void*), sizeof(void*));
    FLUSH_CL();
    #undef SET_OFF
    #undef FLUSH_CL

    // Total memory required. The extra cache line is used to align the buffer.
    int64_t enc_mem = b_off + 64;

    // Set the memory requirements.
    if (!buf)
    {
        p->enc_mem = enc_mem;
        p->bs_mem = d->bs_mem;
    }
}

// Compute the intra availability map.
static void f265_compute_intra_avail_map(f265_gen_data *gd)
{
    // For each block, we compute the size of the top/left blocks and the
    // position of the current block with respect to those blocks.
    //
    // The width of the left block is given by the horizontal alignment of the
    // current block, e.g. there's a 16x16 block left of a block at X=16 or
    // X=48. The vertical position of the current block with respect to the left
    // block is given by the Y offset modulo the left block width.
    //
    // The width of the top block is twice the vertical alignment of the current
    // block, e.g. theres'a 32x32 block above a block at Y=16 or Y=48. The
    // horizontal position of the current block with respect to the top block is
    // given by the X offset modulo the top block width.
    //
    // The CTB size determines the top/left block width at offset 0.

    // Pass each 4x4 block of the CTB.
    int b4_stride = gd->ctb_size>>2;
    for (int by = 0; by < b4_stride; by++)
    {
        for (int bx = 0; bx < b4_stride; bx++)
        {
            // Offsets of the current block.
            int boff[2] = { bx, by };

            // Process each direction (top/left).
            int map_val = 0;
            for (int i = 0; i < 2; i++)
            {
                // Get the log of the neighbour block width (bit scan forward).
                int off = boff[!i]|b4_stride;
                int nl = 0;
                while (!(off&(1<<nl))) nl++;

                // Twice the alignment for the top block.
                nl += !i;

                // Find the position of the current block relative to the
                // neighbour block.
                int pos = boff[i] & (0xff >> (8-nl));

                // Number of 4x4 blocks available in the neighbour.
                int nb_avail = (1<<nl) - pos;

                // Get the log of the number of available 4x4 blocks (bit scan
                // reverse). We don't care about fractions.
                int nb_avail_log = 5;
                while (!(nb_avail&(1<<nb_avail_log))) nb_avail_log--;

                // Store the number of available pixels.
                map_val |= (nb_avail_log + 2)<<(i*4);
            }

            gd->intra_avail_map[by*b4_stride + bx] = map_val;
        }
    }
}

// Set the CTB segmentation flags.
static void venc_set_ctb_seg_flags(f265_enc *enc, f265_frame_ctb *ctb, int ctb_x, int ctb_y, int ctb_xy)
{
    f265_gen_data *gd = &enc->gd;
    ctb->seg_flags = 0;

    // Frame start/end.
    if (ctb_xy == 0) ctb->seg_flags |= F265_SC_CABAC_INIT;
    if (ctb_xy == gd->nb_ctb - 1) ctb->seg_flags |= F265_SC_SEGMENT_END|F265_SC_CHUNK_END;

    // WPP.
    if (gd->eflags&F265_PF_WPP)
    {
        // Start of a row and not the first row. If the second CTB of the
        // previous row is fully available, load its CABAC contexts. Here we
        // assume the previous row is available. Otherwise, initialize the
        // CABAC contexts.
        if (ctb_x == 0 && ctb_y != 0)
        {
            if (gd->pix_dim[0] >= 2*gd->ctb_size) ctb->seg_flags |= F265_SC_CABAC_LOAD;
            else ctb->seg_flags |= F265_SC_CABAC_INIT;
        }

        // Second CTB of the row. Save the CABAC contexts.
        if (ctb_x == 1) ctb->seg_flags |= F265_SC_CABAC_SAVE;

        // Last CTB of the row. End the current chunk.
        if (ctb_x == gd->ctb_dim[0] - 1) ctb->seg_flags |= F265_SC_CHUNK_END;
    }
}

// Initialize the encoder memory.
static void venc_init_enc_mem(f265_enc_params *p, f265_enc_mem_data *d, char **error_handle)
{
    // General data.
    f265_enc *enc = d->enc;
    f265_gen_data *gd = &enc->gd;
    gd->hbd_flag = p->hbd_flag;
    for (int i = 0; i < 2; i++) gd->csf[i] = d->csf[i];
    gd->ctb_size = d->ctb_size;
    gd->stride = d->stride;
    for (int i = 0; i < 2; i++) gd->plane_size[i] = d->plane_size[i];
    for (int i = 0; i < 3; i++) gd->plane_off[i] = d->plane_off[i];
    gd->nb_ctb = d->nb_ctb;
    gd->ctb_chunk_size = d->ctb_chunk_size;
    for (int i = 0; i < 2; i++)
    {
        gd->ctb_dim[i] = d->ctb_dim[i];
        gd->pix_dim[i] = d->pix_dim[i];
        gd->clip_dim[i] = p->clip_dim[i];
        gd->nb_workers[i] = p->nb_workers[i];
    }
    gd->mt_mode = p->mt_mode;
    gd->chroma_format = p->chroma_format;
    for (int i = 0; i < 4; i++) gd->bit_depth[i] = p->bit_depth[i];
    for (int i = 0; i < 2; i++)
    {
        gd->cb_range[i] = p->cb_range[i];
        gd->pcm_range[i] = p->pcm_range[i];
        gd->tb_range[i] = p->tb_range[i];
        gd->tb_depth[i] = p->tb_depth[i];
    }
    gd->qg_log = p->qg_log;
    gd->nb_refs = p->nb_refs;
    gd->nb_b_frames = p->nb_b_frames;
    gd->profile_idc = p->profile_idc;
    gd->level_idc = p->level_idc;
    gd->chroma_qp_idx_off = p->chroma_qp_idx_off;
    for (int i = 0; i < 2; i++) gd->deblock_off[i] = p->deblock_off[i];
    gd->merge_cand = p->merge_cand;
    gd->parallel_merge_level = p->parallel_merge_level;
    gd->poc_bits = 8; // FIXME, hardcoded for HM compatibility.
    gd->default_nb_ref_idx[0] = F265_MAX(p->nb_refs, 1);
    gd->default_nb_ref_idx[1] = 1;
    gd->nb_reordered_frames = !!p->nb_b_frames;
    gd->frame_rate_num = p->frame_rate_num;
    gd->frame_rate_den = p->frame_rate_den;
    gd->algo = p->algo;

    f265_me_search_func sf[4] = {venc_me_dia_search, venc_me_xdia_search, venc_me_hex_search, venc_me_square_search};
    for (int i = 0; i < 3; i++)
    {
        gd->early_me_params[0].dist_func_ids[i] = (p->me_dist[i]<<1) + (i==2 && p->chroma_me_flag);
        gd->early_me_params[0].nb_iters[i] = p->me_iter[i];
        gd->early_me_params[0].search_funcs[i] = sf[p->me_algo[i]];
    }

    // Encoder flags.
    gd->eflags = 0;
    #define SET_PF(flag, value) F265_SET_FLAG(gd->eflags, flag, value)
    SET_PF(F265_PF_DEBLOCK, p->deblock_flag);
    SET_PF(F265_PF_DEBLOCK_ALL, 1);
    SET_PF(F265_PF_SAO, p->sao_flag);
    SET_PF(F265_PF_SCALING, p->scaling_flag);
    SET_PF(F265_PF_PCM, p->pcm_flag);
    SET_PF(F265_PF_TRANSQUANT_BYPASS, p->transquant_bypass_flag);
    SET_PF(F265_PF_SIGN_HIDING, p->sign_hiding_flag);
    SET_PF(F265_PF_TRANSFORM_SKIP, p->transform_skip_flag);
    SET_PF(F265_PF_SMOOTH_INTRA, p->smooth_intra_flag);
    SET_PF(F265_PF_AMP, p->amp_flag);
    SET_PF(F265_PF_TMV, p->tmv_flag);
    SET_PF(F265_PF_CHROMA_ME, p->chroma_me_flag);
    SET_PF(F265_PF_SUBPEL_ME, d->subpel_me_flag);
    SET_PF(F265_PF_LA, p->la_flag);
    SET_PF(F265_PF_ALLOW_FRAME_RECODE, p->allow_recode_flag[0]);
    SET_PF(F265_PF_ALLOW_CTB_RECODE, p->allow_recode_flag[1]);
    SET_PF(F265_PF_RDOQ, p->rdoq_flag);
    SET_PF(F265_PF_WPP, p->wpp_flag);
    #undef SET_PF

    // Chroma QP table. Taken from the spec. Assumes BD 8.
    for (int qpy = 0; qpy < 52; qpy++)
    {
        const uint8_t table[52] =
        {
            0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
            26, 27, 28, 29, 29, 30, 31, 32, 33, 33, 34, 34, 35,
            35, 36, 36, 37, 37, 38, 39, 40, 41, 42, 43, 44, 45
        };
        gd->chroma_qp_table[qpy] = table[F265_CLAMP(qpy + gd->chroma_qp_idx_off, 0, 51)];
    }

    // Intra availability map.
    f265_compute_intra_avail_map(gd);

    // Frame map.
    {
        f265_frame_map *fmap = d->fmap;
        uint8_t *o = (uint8_t*)fmap; o += sizeof(f265_frame_map);
        fmap->section_map = (int(*)[2])o; o += d->nb_sections*4*2;
        fmap->ctb_map = (f265_frame_ctb*)o; o += d->nb_ctb*sizeof(f265_frame_ctb);
        int pix_dim[2] = { gd->pix_dim[0], gd->pix_dim[1] };
        int ctb_size = gd->ctb_size;

        // FIXME: we assume one tile/slice per thread.
        fmap->nb_sections = d->nb_sections;
        fmap->section_map[0][0] = 0;
        fmap->section_map[0][1] = gd->nb_ctb;
        for (int y = 0; y < gd->ctb_dim[1]; y++)
        {
            for (int x = 0; x < gd->ctb_dim[0]; x++)
            {
                int ctb_xy = y*gd->ctb_dim[0] + x;
                f265_frame_ctb *ctb = fmap->ctb_map + ctb_xy;

                // Position.
                ctb->pos[0] = x;
                ctb->pos[1] = y;

                // Neighbours (D, B, C, A).
                ctb->neighbours = ((y>0&&x>0)<<0) | ((y>0)<<1) | ((y>0&&x<gd->ctb_dim[0]-1)<<2) | ((x>0)<<3);

                // Set the CTB segmentation flags.
                venc_set_ctb_seg_flags(enc, ctb, x, y, ctb_xy);

                // The quarterpel filter refers to 3 pixels left/above of the
                // current pixel and 4 pixels right/below.
                //
                // Example: frame with 64+64+32=160 pixels, ctb_size=64.
                //          Let F265_LUMA_PLANE_PADDING=72, F265_SEARCH_OOB=0
                //          so that the MC and ME bounds are symmetrical.
                //          Padded pixels: left [-72..-1], right [160..231].
                //
                //          CTB 1 has 64 pixels left, 160-64=96 pixels right.
                //          ME left:  64 - (64+72-4)    = -68. Pixels [-68..-4].
                //          ME right: 64 + (96+72-4-64) = 164. Pixels [164..227].
                //
                //            ME left     *           ME right
                //          | |<==>| |    *         | |<==>| |
                //          | <------+----*---------+->      |
                //          |4  64  4| 64 * 64 | 32 |4  64  4|
                //          |        |    |<==>|    |        |
                //          ^        ^    Block     ^        ^
                //         -72       0             160      131

                // Number of pixels in the frame around the CTB.
                int l = x*ctb_size, t = y*ctb_size, r = pix_dim[0]-l, b = pix_dim[1]-t;

                // Motion estimation bounds.
                int range = F265_MAX_MV_FPEL - F265_SEARCH_OOB;
                int min_border = F265_LUMA_PLANE_PADDING - 4 - F265_SEARCH_OOB;
                int max_border = min_border - 64;
                ctb->me_bounds[0] = -(F265_MIN(range, l + min_border)<<2);
                ctb->me_bounds[1] = -(F265_MIN(range, t + min_border)<<2);
                ctb->me_bounds[2] =  (F265_MIN(range, r + max_border)<<2);
                ctb->me_bounds[3] =  (F265_MIN(range, b + max_border)<<2);

                // Motion compensation bounds.
                ctb->mc_bounds[0] = -(l + 64 + 3)<<2;
                ctb->mc_bounds[1] = -(t + 64 + 3)<<2;
                ctb->mc_bounds[2] =  (r + 2)<<2;
                ctb->mc_bounds[3] =  (b + 2)<<2;
            }
        }
    }

    // Main data.
    f265_main_data *md = &enc->md;
    md->enc_threads = d->md_enc_threads;
    for (int i = 0; i < d->nb_enc_threads; i++)
    {
        md->enc_threads[i] = d->enc_threads + i;
        venc_me_set_params(enc, &md->enc_threads[i]->me);
    }
    md->enc_frames = d->md_enc_frames;
    md->enc_frames[0] = NULL;
    md->la_ref_frame = NULL;
    md->abs_poc_counter = 0;
    md->nb_past_frames = 0;
    md->nb_enc_threads = d->nb_enc_threads;
    md->nb_enc_frames = 0;
    md->nb_la_frames = 0;
    md->la_process_delay = p->la_process_delay;
    md->enc_exit_flag = 0;
    md->unused_objs[F265_UNUSED_FRM_OBJ] = d->md_unused_objs;
    md->unused_objs[F265_UNUSED_SRC_OBJ] = md->unused_objs[F265_UNUSED_FRM_OBJ] + d->nb_frame_objs;
    md->unused_objs[F265_UNUSED_LA_OBJ]  = md->unused_objs[F265_UNUSED_SRC_OBJ] + d->nb_src_objs;
    md->unused_objs[F265_UNUSED_ENC_OBJ] = md->unused_objs[F265_UNUSED_LA_OBJ]  + d->nb_la_objs;
    md->nb_unused_objs[F265_UNUSED_FRM_OBJ] = d->nb_frame_objs;
    md->nb_unused_objs[F265_UNUSED_SRC_OBJ] = d->nb_src_objs;
    md->nb_unused_objs[F265_UNUSED_LA_OBJ] = d->nb_la_objs;
    md->nb_unused_objs[F265_UNUSED_ENC_OBJ] = d->nb_enc_objs;
    for (int i = 0; i < d->nb_frame_objs; i++) md->unused_objs[F265_UNUSED_FRM_OBJ][i] = (uint8_t*)(d->frame_objs + i);
    for (int i = 0; i < d->nb_src_objs; i++) md->unused_objs[F265_UNUSED_SRC_OBJ][i] = d->src_objs + i*d->src_obj_size;
    for (int i = 0; i < d->nb_enc_objs; i++) md->unused_objs[F265_UNUSED_ENC_OBJ][i] = d->enc_objs + i*d->enc_obj_size;

    // Stream-level rate control initialization.
    venc_rc_init_stream(enc, p);

    // Lookahead.
    f265_lookahead *la = &enc->la;
    la->display = d->la_display;
    la->display[0] = NULL;
    la->coded = d->la_coded;

    // The output queue is initialized with process_delay null frames for
    // initial buffering.
    la->in_queue = d->la_queues;
    la->out_queue = d->la_queues + (1+p->la_process_delay);
    la->in_queue_len = 0;
    la->out_queue_len = p->la_process_delay;
    memset(la->out_queue, 0, la->out_queue_len*sizeof(void*));

    la->key_interval = p->key_interval;
    la->key_type = p->key_type;
    la->key_countdown = 1;
    la->key_idr_flag = 1;
    la->process_delay = p->la_process_delay;
    la->decision_delay = p->la_decision_delay;
    la->nb_undecided = 0;
    la->nb_committed = 0;
    la->nb_b_frames = p->nb_b_frames;

    // Encoding thread objects.
    for (int i = 0; i < d->nb_enc_threads; i++)
    {
        f265_enc_thread *t = d->enc_threads + i;
        t->plane_stride = d->stride;
        t->enc = enc;
        t->me.t = t;
        t->me.ref_stride = d->stride;
        t->chunk_buf = (p->mt_mode == F265_MT_ENC_FRAME) ? d->chunk_objs + i*d->chunk_obj_size : d->chunk_objs;
        t->cbs.mode = F265_CABAC_RAW;
        t->cbs.encode_context_bin = venc_encode_context_bin_raw;
        t->cbs.encode_bypass_bin = venc_encode_bypass_bin_raw;
        t->cbs.encode_bypass_bins = venc_encode_bypass_bins_raw;
        t->cbs.encode_term0_bin = venc_encode_term0_bin_raw;
        t->store = t->store_buf;

        // Set the quality settings.
        t->an.rdm_flag = !p->rdo_level;
        t->an.hm_me_flag = p->hm_me_flag;
        t->an.all_intra_flag = p->all_intra_flag;
        t->an.nullify_inter_tb_flag = p->nullify_inter_tb_flag;

        // Initialize the prediction map.
        for (int j = 0; j < 100; j++) t->pmap[j] = F265_UNAVAIL_CB_IDX;

        // Initialize the unavailable CB. We set its size to 64x64 so that it is
        // considered unsplit.
        t->cb[F265_UNAVAIL_CB_IDX].lg_bs = 6;
        t->cb[F265_UNAVAIL_CB_IDX].enc_idx = 127;
    }

    // Frame objects.
    for (int i = 0; i < d->nb_frame_objs; i++)
    {
        f265_frame *f = d->frame_objs + i;
        f->fmap = d->fmap;
    }

    // YUV dump.
    if (p->yuv_dump_path)
    {
        gd->yuv_dump_file = fopen(p->yuv_dump_path, "wb");
        if (!gd->yuv_dump_file)
        {
            *error_handle = "cannot open YUV dump file";
            return;
        }
    }

    // HM GOP file specified.
    if (p->hm_gop_path)
    {
        venc_parse_hm_gop_file(enc, p->hm_gop_path);

        // Intra-only compatibility.
        if (!gd->hm_gop_size)
            enc->gd.default_nb_ref_idx[0] = enc->gd.default_nb_ref_idx[1] = 4;
    }

    #ifdef VAN_USE_HM_ME
    // Memory allocation for BD conversion to 16s. Needed for bi-dir
    // inter search HM compatible.
    {
        md->hm_last_poc = -1;
        int stride = enc->gd.stride;
        int *pix_dim = enc->gd.pix_dim;

        int src_size = stride * pix_dim[1];
        md->hm_src_plane_16s = (f265_pix*)malloc(src_size * 2);
        md->hm_new_hm_src_plane_16s = (f265_pix*)malloc(src_size * 2);

        int ref_size = stride * (pix_dim[1] + 2 * F265_LUMA_PLANE_PADDING);
        for (int l = 0; l < 2; l++)
            for (int r = 0; r < 16; r++)
            {
                md->hm_ref_plane_16s[l][r] = (f265_pix*)malloc(ref_size * 2);
            }
        md->hm_ref_dual_plane_16s = (f265_pix*)malloc(src_size * 2);
    }
    #endif
}

// Analyze the parameters and initialize the encoder if a memory buffer is
// provided.
static void venc_alloc_enc_mem(f265_enc_params *p, uint8_t *buf, f265_enc **enc_handle, char **error_handle)
{
    f265_enc_mem_data d;
    venc_analyze_enc_mem(p, &d, buf);
    if (buf)
    {
        *enc_handle = d.enc;
        *error_handle = NULL;
        venc_init_enc_mem(p, &d, error_handle);
        if (*error_handle) venc_deinit_enc(d.enc);
    }
}

// Note: the initialization code also calls this function to clean up on
// failure, so fields must be zeroed properly.
void venc_deinit_enc(f265_enc *enc)
{
    if (!enc) return;
    f265_gen_data *gd = &enc->gd;

    #ifdef F265_HAVE_STDIO
    if (gd->yuv_dump_file) { fclose(gd->yuv_dump_file); gd->yuv_dump_file = NULL; }
    #endif

    #ifdef VAN_USE_HM_ME
    // Free the memory allocated for BD conversion.
    f265_main_data *md = &enc->md;

    if (md->hm_src_plane_16s) free(md->hm_src_plane_16s);
    if (md->hm_new_hm_src_plane_16s) free(md->hm_new_hm_src_plane_16s);

    for (int l = 0; l < 2; l++)
        for (int r = 0; r < 16; r++)
            if (md->hm_ref_plane_16s[l][r]) free(md->hm_ref_plane_16s[l][r]);
    if (md->hm_ref_dual_plane_16s) free(md->hm_ref_dual_plane_16s);
    #endif

    #ifdef VAN_LOAD_CTB_ANALYSIS
    extern FILE *van_an_dump_file;
    // To prevent valgrind errors.
    if (van_an_dump_file) fclose(van_an_dump_file);
    #endif
}

void venc_analyze_params(f265_enc_params *params)
{
    venc_alloc_enc_mem(params, NULL, NULL, NULL);
}

void venc_init_enc(f265_enc_params *params, uint8_t *buf, f265_enc **enc_handle, char **error_handle)
{
    venc_alloc_enc_mem(params, buf, enc_handle, error_handle);
}

// Get an unused object of the category specified.
uint8_t* venc_get_unused_obj(f265_enc *e, int cat)
{
    f265_main_data *md = &e->md;
    assert(md->nb_unused_objs[cat]);
    return md->unused_objs[cat][--md->nb_unused_objs[cat]];
}

// Link an unused source object to the frame specified.
void venc_link_src_obj(f265_enc *e, f265_frame *f)
{
    f265_gen_data *gd = &e->gd;
    uint8_t *o = venc_get_unused_obj(e, F265_UNUSED_SRC_OBJ);
    f->mem_buf[1] = o;
    f->src_planes[0] = (f265_pix*)(o + gd->plane_off[0]); o += gd->plane_size[0];
    for (int i = 0; i < 2; i++) f->src_planes[1+i] = (f265_pix*)(o + gd->plane_off[1+i]);
    o += gd->plane_size[1];
}

// Link an unused encoding object to the frame specified.
void venc_link_enc_obj(f265_enc *e, f265_frame *f)
{
    f265_gen_data *gd = &e->gd;
    int nb_luma_rec = 1 + F265_GET_FLAG(gd->eflags, F265_PF_SUBPEL_ME)*3;
    int64_t b4_dim[2] = { gd->pix_dim[0]>>2, gd->pix_dim[1]>>2 };
    int64_t nb_b4 = b4_dim[0]*b4_dim[1];
    int64_t nb_b8 = nb_b4>>2;

    uint8_t *o = venc_get_unused_obj(e, F265_UNUSED_ENC_OBJ);
    f->mem_buf[3] = o;
    for (int i = 0; i < nb_luma_rec; i++, o += gd->plane_size[0]) f->rec_planes[i] = (f265_pix*)(o + gd->plane_off[0]);
    for (int i = 0; i < 2; i++) f->rec_planes[4+i] = (f265_pix*)(o + gd->plane_off[1+i]);
    o += gd->plane_size[1];
    f->mv = (f265_mv(*)[2])o; o += nb_b4*4*2;
    f->b8_map = (uint16_t*)o; o += nb_b8*2;
    f->sao = (f265_ctb_sao*)o; o += gd->nb_ctb*sizeof(f265_ctb_sao);
    f->qmap = (int8_t*)o; o += nb_b8;
    f->b4_map = (uint8_t*)o; o += nb_b4;
    f->ref_idx = (int8_t(*)[2])o; o += nb_b4*2;
}

// Release the memory no longer needed. This includes the frame object itself.
// The function assumes the frame needs not be encoded.
void venc_release_frame_mem(f265_enc *e, f265_frame *f)
{
    #define REL(CAT, FLAGS) \
        if (unlikely(!(f->main_flags&(FLAGS)) && f->mem_buf[CAT])) \
        { \
            e->md.unused_objs[CAT][e->md.nb_unused_objs[CAT]++] = f->mem_buf[CAT]; \
            f->mem_buf[CAT] = NULL; \
        }
    REL(F265_UNUSED_FRM_OBJ, F265_FF_LA_REF|F265_FF_DPB_REF|F265_FF_LAST_ENC);
    REL(F265_UNUSED_SRC_OBJ, F265_FF_LA_REF);
    REL(F265_UNUSED_LA_OBJ,  F265_FF_LA_REF);
    REL(F265_UNUSED_ENC_OBJ, F265_FF_DPB_REF);
    #undef REL
}

int venc_process_enc_req(f265_enc *e, f265_enc_req *req)
{
    // Set the request in the encoder.
    int res;
    f265_main_data *md = &e->md;
    md->req = req;

    do
    {
        // Encode a new frame. Bail out if the input frame is null and the encoder
        // has been fully flushed.
        if (req->input == NULL && md->nb_la_frames + md->nb_enc_frames == 0)
        {
            res = F265_RET_EMPTY;
            break;
        }

        // Prepare the frame for insertion in the lookahead.
        f265_frame *in_frame = req->input ? venc_la_make_frame(e) : NULL;

        // Add the frame to the lookahead if the frame is non-null or the lookahead
        // is being flushed.
        f265_frame *next_frame = NULL;
        if (in_frame || md->nb_la_frames)
        {
            // Retrieve the next frame to encode. If the input frame is null, we
            // must keep flushing until the lookahead has produced a frame.
            while (1)
            {
                next_frame = venc_la_add_frame(e, in_frame);
                if (in_frame || next_frame) break;
            }

            // The lookahead has not produced a frame, so it is still buffering.
            if (!next_frame)
            {
                res = F265_RET_EMPTY;
                break;
            }
        }

        // Encode the frame, if any.
        res = venc_encode_frame_mt_none(e, next_frame);

    } while (0);

    return res;
}

// Set the DPB of the current frame. Improve the logic eventually.
static void venc_set_frame_dpb(f265_enc *enc, f265_frame *f)
{
    // Set the DPB reference flag of the current frame.
    F265_SET_FLAG(f->main_flags, F265_FF_DPB_REF, F265_GET_FLAG(f->gen_flags, F265_FF_REF));

    // Get the list of frames currently in the DPB, i.e. the frames in the DPB
    // of the previous frame and the previous frame itself if it is used for
    // reference.
    f265_frame *frames[17];
    int nb_frames = 0;
    if (f->abs_poc)
    {
        f265_frame *pf = enc->md.enc_frames[enc->md.nb_enc_frames];
        int add_pf_flag = F265_GET_FLAG(pf->gen_flags, F265_FF_REF);
        for (int i = 0; i < pf->dpb_size; i++)
        {
            f265_frame *ref = pf->dpb[i];
            if (add_pf_flag && ref->abs_poc > pf->abs_poc) { frames[nb_frames++] = pf; add_pf_flag = 0; }
            frames[nb_frames++] = ref;
        }
        if (add_pf_flag) frames[nb_frames++] = pf;
    }

    // Flag indicating whether a frame is kept in the DPB.
    uint8_t kept[17];
    // Bogus 'if' because gcc is buggy (__warn_memset_zero_len).
    if (nb_frames) memset(kept, 1, nb_frames);

    // Track whether we need to vacate a frame.
    int max_nb_refs = enc->gd.nb_refs;
    int vacate_flag = nb_frames > max_nb_refs;

    // Vacate all frames.
    if (f->gen_flags&F265_FF_IDR)
    {
        vacate_flag = 0;
        memset(kept, 0, nb_frames);
    }

    // Vacate all but the last CRA frame.
    else if (nb_frames && (frames[nb_frames-1]->gen_flags&F265_FF_CRA) && f->abs_poc > frames[nb_frames-1]->abs_poc)
    {
        vacate_flag = 0;
        memset(kept, 0, nb_frames-1);
    }

    // Vacate frames using the HM logic.
    if (enc->gd.hm_gop_size)
    {
        for (int i = 0; i < nb_frames; i++)
        {
            if (!kept[i]) continue;
            f265_frame *ref = frames[i];
            f265_hm_gop_entry *entry = f->gop_entry;
            int found_flag = 0;
            for (int j = 0; j < entry->nb_refs && !found_flag; j++)
                found_flag = f->abs_poc + entry->refs[j] == ref->abs_poc;
            kept[i] = found_flag;
        }
    }

    // Vacate the oldest frame.
    else if (vacate_flag) kept[0] = 0;

    // Set the frames used for reference in the DPB. Free the others if they are
    // not being encoded.
    f->dpb_neg = f->dpb_pos = f->dpb_size = 0;
    for (int i = 0; i < nb_frames; i++)
    {
        f265_frame *ref = frames[i];
        if (kept[i])
        {
            f->dpb[f->dpb_size++] = ref;
            if (ref->abs_poc < f->abs_poc) f->dpb_neg++;
            else f->dpb_pos++;
        }
        else
        {
            F265_SET_FLAG(ref->main_flags, F265_FF_DPB_REF, 0);
            if (ref->main_flags&F265_FF_ENC) venc_release_frame_mem(enc, ref);
        }
    }
}

// Build the default reference lists. Each list contains each negative and
// positive frame once, independently of the frame type.
static void venc_get_default_ref_lists(f265_frame *f, f265_frame *def_lists[2][16])
{
    int neg_pos = f->dpb_neg-1;
    for (int i = 0; i < f->dpb_neg; i++)
    {
        f265_frame *ref = f->dpb[neg_pos-i];
        def_lists[0][i] = ref;
        def_lists[1][f->dpb_pos+i] = ref;
    }
    for (int i = 0; i < f->dpb_pos; i++)
    {
        f265_frame *ref = f->dpb[neg_pos+1+i];
        def_lists[0][f->dpb_neg+i] = ref;
        def_lists[1][i] = ref;
    }
}

// Build the typical reference lists.
static void venc_get_typical_ref_lists(f265_enc_thread *t, f265_frame *f, f265_frame *act_lists[2][16],
                                       f265_frame *def_lists[2][16])
{
    // HM logic.
    if (t->enc->gd.hm_gop_size)
    {
        int nb_refs = F265_MIN(f->dpb_size, f->gop_entry->nb_active);
        for (int i = 0; i < 2; i++)
        {
            memcpy(act_lists[i], def_lists[i], nb_refs*sizeof(void*));
            t->nb_ref_idx[i] = nb_refs;
        }
        if (f->frame_type == F265_FRAME_P) t->nb_ref_idx[1] = 0;
        return;
    }

    // P frame.
    if (f->frame_type == F265_FRAME_P)
    {
        t->nb_ref_idx[0] = f->dpb_size;
        t->nb_ref_idx[1] = 0;
        memcpy(act_lists[0], def_lists[0], t->nb_ref_idx[0]*sizeof(void*));
    }

    // B frame.
    else
    {
        t->nb_ref_idx[0] = F265_MAX(f->dpb_neg, 1);
        t->nb_ref_idx[1] = F265_MAX(f->dpb_pos, 1);
        memcpy(act_lists[0], def_lists[0], t->nb_ref_idx[0]*sizeof(void*));
        memcpy(act_lists[1], def_lists[1], t->nb_ref_idx[1]*sizeof(void*));
    }
}

// Build the inter configuration.
static void venc_set_frame_inter(f265_enc_thread *t, f265_frame *f)
{
    // I frame.
    if (f->frame_type == F265_FRAME_I)
    {
        t->nb_lists = t->nb_ref_idx[0] = t->nb_ref_idx[1] = 0;
    }

    // P/B frame.
    else
    {
        // Build the default and the actual reference lists.
        f265_frame *def_lists[2][16], *act_lists[2][16];
        venc_get_default_ref_lists(f, def_lists);
        venc_get_typical_ref_lists(t, f, act_lists, def_lists);

        // Set the number of lists.
        t->nb_lists = 1 + (f->frame_type == F265_FRAME_B);

        // Set the reference contexts.
        for (int list = 0; list < t->nb_lists; list++)
            for (int idx = 0; idx < t->nb_ref_idx[list]; idx++)
            {
                f265_ref_ctx *rc = t->ref_ctx[list] + idx;
                f265_frame *ref = act_lists[list][idx];
                rc->frame = ref;
                for (int i = 0; i < 6; i++) rc->planes[i] = ref->rec_planes[i];
                rc->weights = f265_lt.wp_default;
            }
    }

    // Cache the information in the frame.
    f->nb_lists = t->nb_lists;
    for (int i = 0; i < 2; i++) f->nb_ref_idx[i] = t->nb_ref_idx[i];

    // Set up the distance scale factors.
    venc_frame_set_up_inter_pred(t);
}

// Encode the current frame with the main thread. The input frame cannot be
// null. Return VALID, ABORT, ERROR.
static int venc_encode_frame_mt_none(f265_enc *e, f265_frame *f)
{
    f265_enc_thread *t = venc_set_enc_thread(e, f);

    // Kludge to load HM DF map before encoding to normalize it during encoding.
    #ifdef VAN_DUMP_DEBLOCK_FILTER
    void venc_df_init_dump(f265_enc_thread *t);
    void venc_df_init_frame_dump();
    venc_df_init_dump(t);
    venc_df_init_frame_dump();
    #endif

    // Trace progression.
    #ifdef VAN_TRACE_FRAME_ENCODE
    char fr_type[3] = {'I', 'P', 'B'};
    printf("Encoding frame POC %d type %c.\n", (int)t->src_frame->abs_poc, fr_type[t->src_frame->frame_type]);
    #endif

    // Trace reference lists.
    #ifdef VAN_TRACE_REF_LIST
    printf("L0: ");
    for (int r = 0; r < t->nb_ref_idx[0]; r++)
        printf("%ld ", t->ref_ctx[0][r].frame->abs_poc);
    printf("\n");
    printf("L1: ");
    for (int r = 0; r < t->nb_ref_idx[1]; r++)
        printf("%ld ", t->ref_ctx[1][r].frame->abs_poc);
    printf("\n");
    #endif

    // Frame re-encoding loop.
    while (1)
    {
        // FIXME: assuming one section.
        int res = venc_encode_section(t, 0);
        if (res == F265_RET_ABORT) return res;
        res = venc_check_enc_result(t, f);
        if (res != F265_RET_RETRY) return res;
    }
}

// Encode the current section. Return VALID or ABORT.
static int venc_encode_section(f265_enc_thread *t, int section_idx)
{
    // Import the section data.
    f265_frame *f = t->src_frame;
    f265_frame_map *fmap = f->fmap;
    t->ctb_idx = fmap->section_map[section_idx][0];
    int nb_ctbs = fmap->section_map[section_idx][1];

    // Get the chunk buffer offsets.
    uint8_t *chunks;
    venc_get_chunk_buf_offsets(t->enc, t->chunk_buf, &t->nb_segs, &t->seg, &chunks, t->ctb_idx, nb_ctbs);

    // Set the segment count to 0 and remember we need to initialize a segment.
    *t->nb_segs = 0;
    t->init_segment_flag = 1;

    // Initialize the current chunk and the CABAC engine.
    venc_init_seg_chunk(&t->chunk, chunks);
    venc_init_cabac_engine(&t->cbs, t->cabac_buf);

    // Encode each CTB.
    for (int i = 0; i < nb_ctbs; i++, t->ctb_idx++) venc_encode_ctb(t);

    return F265_RET_VALID;
}

// Return the plane offset corresponding to the specified block offset in the
// CTB.
int venc_get_ctb_block_plane_off(f265_enc_thread *t, int comp, int ct_ox, int ct_oy)
{
    int csf = !!comp;
    return ((t->ctb_off[1]>>csf) + ct_oy)*t->me.ref_stride + (t->ctb_off[0]>>csf) + ct_ox;
}

// Return the plane offset corresponding to the specified block offset in the
// CB.
int venc_get_cb_block_plane_off(f265_enc_thread *t, f265_cb *cb, int comp, int cb_ox, int cb_oy)
{
    int csf = !!comp;
    return (((t->ctb_off[1] + cb->cb_off[1])>>csf) + cb_oy)*t->me.ref_stride +
            ((t->ctb_off[0] + cb->cb_off[0])>>csf) + cb_ox;
}

// Return the depth scan index of the block (x,y) given the number of blocks on
// a side. For example, the output is 9 for x=1, y=2, size=4.
//         0145
//         2367
//         89CD
//         ABEF
// Remove this eventually.
int venc_get_depth_scan_idx(int x, int y, int size)
{
    if (size == 1) return 0;

    int depth = 0;

    // Half the size.
    int h = size>>1;

    if (x >= h)
    {
        x -= h;
        depth += h*h;
    }

    if (y >= h)
    {
        y -= h;
        depth += 2*h*h;
    }

    return depth + venc_get_depth_scan_idx(x, y, h);
}

// Return the partitioning mode of the coding block. Optimize this eventually.
int venc_get_cb_part_mode(f265_cb *cb)
{
    if (cb->flags&F265_CB_INTRA) return (cb->intra_luma_mode[1] == -1) ? F265_PART_UN : F265_PART_HV;
    return cb->inter_part;
}

// Get the CB offsets in 4x4 blocks from the start of the CTB. Optimize this
// eventually.
void venc_get_cb_loc4(int cb_loc[2], f265_cb *cb)
{
    cb_loc[0] = cb->cb_off[0]>>2;
    cb_loc[1] = cb->cb_off[1]>>2;
}

// Return the partition sizes and offsets in 4x4 blocks. Optimize this eventually.
void venc_get_part_loc4(int part_size[2], int part_off[2], int part_mode, int part_idx, int cb_size)
{
    int tmp = f265_part_table[part_mode][part_idx][0];
    part_size[0] = ((tmp&15)*cb_size)>>4;
    part_size[1] = ((tmp>>4)*cb_size)>>4;
    tmp = f265_part_table[part_mode][part_idx][1];
    part_off[0] = ((tmp&15)*cb_size)>>4;
    part_off[1] = ((tmp>>4)*cb_size)>>4;
}

// Return the partition sizes and offsets in pixels. Optimize this eventually.
void venc_get_part_loc(int part_size[2], int part_off[2], int part_mode, int part_idx, int cb_size)
{
    int tmp = f265_part_table[part_mode][part_idx][0];
    part_size[0] = ((tmp&15)*cb_size)>>2;
    part_size[1] = ((tmp>>4)*cb_size)>>2;
    tmp = f265_part_table[part_mode][part_idx][1];
    part_off[0] = ((tmp&15)*cb_size)>>2;
    part_off[1] = ((tmp>>4)*cb_size)>>2;
}

// Update the prediction map entry at base + offset. The offsets are in 4x4
// blocks.
static void venc_update_pmap_entry(uint16_t *base, int b4_ox, int b4_oy, int cb_idx, int part_idx)
{
    // Prediction map index.
    int pmap_idx = (b4_oy>>1)*10 + (b4_ox>>1);

    // 4x4 block mask shift.
    int shift = 8 + ((((b4_oy&1)<<1) + (b4_ox&1))<<1);

    // Update the entry.
    uint16_t *e = base + pmap_idx;
    *e = (*e&0xff00)|cb_idx;
    *e = (*e&~(3<<shift))|(part_idx<<shift);
}

// Update the prediction map for the unsplit CB specified. Optimize this
// eventually.
void venc_update_pmap_unsplit_cb(f265_enc_thread *t, f265_cb *cb)
{
    int cb_idx = cb - t->cb;
    int part_mode = venc_get_cb_part_mode(cb);
    int nb_parts = f265_nb_parts[part_mode];
    int cb_size = 1<<cb->lg_bs;
    int cb_loc[2];
    venc_get_cb_loc4(cb_loc, cb);

    // CB entry in the prediction map.
    uint16_t *cb_entry = t->pmap + 11 + (cb_loc[1]>>1)*10 + (cb_loc[0]>>1);

    // Pass each partition.
    for (int part_idx = 0; part_idx < nb_parts; part_idx++)
    {
        // Get the partition sizes and offsets in 4x4 blocks.
        int part_size[2], part_off[2];
        venc_get_part_loc4(part_size, part_off, part_mode, part_idx, cb_size);

        // Pass each 4x4 block and update the entry.
        for (int by = 0; by < part_size[1]; by++)
            for (int bx = 0; bx < part_size[0]; bx++)
                venc_update_pmap_entry(cb_entry, part_off[0] + bx, part_off[1] + by, cb_idx, part_idx);
    }
}

// Update the prediction map recursively for the CB specified.
void venc_update_pmap_cb(f265_enc_thread *t, f265_cb *cb)
{
    if (cb->flags&F265_CB_SPLIT)
        for (int i = 0; i < 4; i++)
            venc_update_pmap_cb(t, t->cb + cb->child_idx + i);
    else
        venc_update_pmap_unsplit_cb(t, cb);
}

// Return the CB index and the partition index given the 4x4 block position from
// the CB origin.
void venc_lookup_pmap(f265_enc_thread *t, int x, int y, int *cb_idx, int *part_idx)
{
    uint16_t entry = t->pmap[11 + (y>>1)*10 + (x>>1)];
    int shift = 8 + ((((y&1)<<1) + (x&1))<<1);
    *cb_idx = entry&255;
    *part_idx = (entry>>shift)&3;
}

// Load the neighbour partition data.
static void venc_load_neighbour_part(f265_cb *cb, f265_frame *f, int b4_off, int part_idx, int intra_flag)
{
    if (intra_flag)
    {
        cb->intra_luma_mode[part_idx] = f->b4_map[b4_off]&63;
    }

    else
    {
        memcpy(cb->ref_idx[part_idx], f->ref_idx[b4_off], 2);
        memcpy(cb->mv[part_idx], f->mv[b4_off], 2*sizeof(f265_mv));
    }
}

// Return true if the frame data doesn't match the first partition data.
static int venc_is_ctb_neighbour_cb_split(f265_cb *cb, f265_frame *f, int b4_off, int intra_flag)
{
    if (intra_flag)
    {
        return cb->intra_luma_mode[0] != (f->b4_map[b4_off]&63);
    }

    else
    {
        return memcmp(cb->ref_idx[0], f->ref_idx[b4_off], 2) ||
               memcmp(cb->mv[0], f->mv[b4_off], 2*sizeof(f265_mv));
    }
}

// Load the neighbour CB data and update the prediction map. loc is the location
// of the neighbour (0: top, 1: top-right, 2: left, 3: top-left). n_off is the
// neighbour offset from the first 8x8 block in the CTB. The function returns
// the number of 8x8 blocks in the neighbour.
int venc_load_ctb_neighbour_cb(f265_enc_thread *t, int *next_cb_idx, int n_off, int loc)
{
    f265_gen_data *gd = &t->enc->gd;
    f265_frame *f = t->src_frame;
    int b4_stride = gd->pix_dim[0]>>2;
    int b8_stride = gd->pix_dim[0]>>3;

    // Allocate the CB entry.
    int cb_idx = (*next_cb_idx)++;
    f265_cb *cb = t->cb + cb_idx;

    // Offsets of the neighbour 4x4/8x8 block from the first 4x4/8x8 block of
    // the CTB.
    int b4_x = (loc >= 2) ? -1 : n_off<<1;
    int b4_y = (loc != 2) ? -1 : n_off<<1;
    int b8_x = b4_x>>1, b8_y = b4_y>>1;

    // Offset of the neighbour 8x8/4x4 block in the frame cache.
    int b4_off = ((t->ctb_off[1]>>2) + b4_y)*b4_stride + (t->ctb_off[0]>>2) + b4_x;
    int b8_off = ((t->ctb_off[1]>>3) + b8_y)*b8_stride + (t->ctb_off[0]>>3) + b8_x;

    // Set the neighbour information.
    int b8_val = f->b8_map[b8_off];
    int intra_flag = b8_val&1;
    cb->lg_bs = ((b8_val>>1)&3) + 3;
    cb->enc_idx = -1;
    int cb_size = 1<<cb->lg_bs;
    int skip_flag = F265_GET_FLAG(b8_val, 1<<8);

    cb->flags = F265_CB_PRESENT;
    if (intra_flag) cb->flags |= F265_CB_INTRA;
    if (skip_flag) cb->flags |= F265_CB_SKIP;

    // Number of 4x4 blocks to update in the prediction map.
    int cb_nb_b4 = (loc&1) ? 1 : cb_size>>2;

    // Number of 4x4 blocks in partition 0.
    int part0_nb_b4 = cb_nb_b4;

    // Increment to get to the next 4x4 block in the frame cache.
    int b4_next_inc = (!loc) ? 1 : b4_stride;

    // Increment to get to the next 4x4 block in the prediction map.
    int b4_next_x = !loc, b4_next_y = !!loc;

    // Import the first partition data.
    venc_load_neighbour_part(cb, f, b4_off, 0, intra_flag);

    // Check if we use the second partition by probing all possible splits. We
    // do not care about the actual partitioning mode.
    if (!(loc&1))
    {
        for (int frac = 1; frac < 4; frac++)
        {
            // Skip assymetrical splits if the CB is 8x8 or intra.
            if (frac != 2 && (cb_size == 8 || intra_flag)) continue;

            // Tentative partition 0 size in 4x4 blocks.
            int ps = (cb_size*frac)>>4;

            // Offset of the second partition 4x4 block in the frame cache.
            int part1_b4_off = b4_off + ps*b4_next_inc;

            // Split detected. Copy the second partition data and update the
            // first partition size.
            if (venc_is_ctb_neighbour_cb_split(cb, f, part1_b4_off, intra_flag))
            {
                venc_load_neighbour_part(cb, f, part1_b4_off, 1, intra_flag);
                part0_nb_b4 = ps;
                break;
            }
        }
    }

    // Update the prediction map.
    for (int i = 0; i < cb_nb_b4; i++)
        venc_update_pmap_entry(t->pmap + 11, b4_x + b4_next_x*i, b4_y + b4_next_y*i, cb_idx, i >= part0_nb_b4);

    return cb_size>>3;
}

// Load the CTB neighbour coding blocks.
void venc_load_ctb_neighbours(f265_enc_thread *t)
{
    f265_gen_data *gd = &t->enc->gd;
    uint16_t *pmap = t->pmap;
    int ctb_d_avail_flag = !!(t->ctb_neighbours & 1);
    int ctb_b_avail_flag = !!(t->ctb_neighbours & 2);
    int ctb_c_avail_flag = !!(t->ctb_neighbours & 4);
    int ctb_a_avail_flag = !!(t->ctb_neighbours & 8);

    // Index of the next CB object allocated as a neighbour.
    int next_cb_idx = 85;

    // Number of 8x8 blocks on the side of the CTB.
    int b8_width = gd->ctb_size>>3;

    // Number of 8x8 blocks to import on the borders.
    int nb_import_top  = ctb_b_avail_flag*(F265_MIN(gd->pix_dim[0]-t->ctb_off[0], gd->ctb_size)>>3);
    int nb_import_left = ctb_a_avail_flag*(F265_MIN(gd->pix_dim[1]-t->ctb_off[1], gd->ctb_size)>>3);

    // Top-left.
    if (ctb_d_avail_flag) venc_load_ctb_neighbour_cb(t, &next_cb_idx, 0, 3);
    else pmap[0] = F265_UNAVAIL_CB_IDX;

    // Top.
    for (int i = 0; i < b8_width;)
    {
        if (i < nb_import_top) i += venc_load_ctb_neighbour_cb(t, &next_cb_idx, i, 0);
        else pmap[1+i++] = F265_UNAVAIL_CB_IDX;
    }

    // Top-right.
    if (ctb_c_avail_flag) venc_load_ctb_neighbour_cb(t, &next_cb_idx, b8_width, 1);
    else pmap[1+b8_width] = F265_UNAVAIL_CB_IDX;

    // Left.
    for (int i = 0; i < b8_width;)
    {
        if (i < nb_import_left) i += venc_load_ctb_neighbour_cb(t, &next_cb_idx, i, 2);
        else pmap[(1+i++)*10] = F265_UNAVAIL_CB_IDX;
    }
}

// Save the CB data to the frame. Optimize this eventually.
void venc_save_cb(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return;

    // Handle split CB.
    if (cb->flags&F265_CB_SPLIT)
    {
        for (int i = 0; i < 4; i++) venc_save_cb(t, t->cb + cb->child_idx + i);
        return;
    }

    // Handle unsplit CB.
    f265_gen_data *gd = &t->enc->gd;
    f265_frame *f = t->src_frame;
    int intra_flag = F265_GET_FLAG(cb->flags, F265_CB_INTRA);
    int skip_flag = F265_GET_FLAG(cb->flags, F265_CB_SKIP);
    int part_mode = venc_get_cb_part_mode(cb);
    int nb_parts = f265_nb_parts[part_mode];
    int cb_size = 1<<cb->lg_bs;
    int b8_val = (intra_flag<<0) | ((cb->lg_bs-3)<<1) | (skip_flag<<8);
    int b4_stride = gd->pix_dim[0]>>2;
    int b8_stride = gd->pix_dim[0]>>3;
    int cb_loc[2];
    venc_get_cb_loc4(cb_loc, cb);

    // Pass each partition.
    for (int part_idx = 0; part_idx < nb_parts; part_idx++)
    {
        // Get the partition sizes and offsets in 4x4 blocks.
        int part_size[2], part_off[2];
        venc_get_part_loc4(part_size, part_off, part_mode, part_idx, cb_size);

        // Get the partition data.
        uint8_t intra_mode = cb->intra_luma_mode[part_idx];
        f265_mv *mv = cb->mv[part_idx];
        int8_t *ref_idx = cb->ref_idx[part_idx];

        // Pass each 4x4 block.
        for (int by = 0; by < part_size[1]; by++)
        {
            for (int bx = 0; bx < part_size[0]; bx++)
            {
                // Block offsets from the start of the frame.
                int b4_x = (t->ctb_off[0]>>2) + cb_loc[0] + part_off[0] + bx;
                int b4_y = (t->ctb_off[1]>>2) + cb_loc[1] + part_off[1] + by;
                int b4_off = b4_y*b4_stride + b4_x;
                int b8_off = (b4_y>>1)*b8_stride + (b4_x>>1);

                // Update the frame cache.
                f->b8_map[b8_off] = b8_val;

                if (intra_flag)
                {
                    f->b4_map[b4_off] = intra_mode;
                    memset(f->ref_idx[b4_off], -1, 2);
                }
                else
                {
                    memcpy(f->ref_idx[b4_off], ref_idx, 2);
                    memcpy(f->mv[b4_off], mv, 2*sizeof(f265_mv));
                }
            }
        }
    }

    // Save the deblocking parameters.
    if (f->gen_flags&F265_FF_DEBLOCK)
    {
        venc_save_deblock_tt(t, cb, cb->lg_bs, 0, 0);
        venc_save_deblock_pb(t, cb);
    }
}

// Load the frame CTB data in the thread.
static void venc_load_ctb(f265_enc_thread *t)
{
    f265_gen_data *gd = &t->enc->gd;
    f265_me_ctx *me = &t->me;
    f265_frame *f = t->src_frame;
    f265_frame_ctb *ctb = f->fmap->ctb_map + t->ctb_idx;
    f265_cabac_bs *cbs = &t->cbs;
    int sx = gd->csf[0], sy = gd->csf[1];
    int stride = t->me.ref_stride;

    // Initialize.
    t->ctb_x = ctb->pos[0];
    t->ctb_y = ctb->pos[1];
    t->ctb_xy = t->ctb_y*gd->ctb_dim[0] + t->ctb_x;
    t->ctb_off[0] = t->ctb_x*gd->ctb_size;
    t->ctb_off[1] = t->ctb_y*gd->ctb_size;
    t->mc_bounds64 = ctb->mc_bounds64;
    t->seg_flags = ctb->seg_flags;
    t->ctb_neighbours = ctb->neighbours;
    me->me_bounds64 = ctb->me_bounds64;
    me->me_bounds_packs[0] = venc_pack_mv_range(-me->me_bounds[0], -me->me_bounds[1]);
    me->me_bounds_packs[1] = venc_pack_mv_range(me->me_bounds[2], me->me_bounds[3]);

    // Set the pointers to the reconstructed planes.
    for (int i = 0; i < 4; i++)
        t->rec_planes[i] = f->rec_planes[i] + t->ctb_off[1]*stride + t->ctb_off[0];
    for (int i = 0; i < 2; i++)
        t->rec_planes[4+i] = f->rec_planes[4+i] + (t->ctb_off[1]*stride>>sy) + (t->ctb_off[0]>>sx);

    // Set the intra availability.
    t->intra_avail_cutoff[0] = F265_MIN(gd->pix_dim[0] - t->ctb_off[0], gd->ctb_size<<1);
    t->intra_avail_cutoff[1] = F265_MIN(gd->pix_dim[1] - t->ctb_off[1], gd->ctb_size);

    // Mark CBs unavailable by default. Optimize this eventually.
    for (int i = 0; i < 100; i++)
        t->pmap[i] = F265_UNAVAIL_CB_IDX;

    // Load the CTB neighbour coding blocks.
    venc_load_ctb_neighbours(t);

    // Load the temporal motion vectors.
    venc_ctb_set_up_inter_pred(t);

    // Initialize the coding blocks inside the CTB. Optimize this eventually. We
    // could skip the initialization if the current and last CTBs encoded were
    // fully inside the frame.
    for (int range_idx = gd->cb_range[1], range_cb_idx = 0; range_idx >= gd->cb_range[0]; range_idx--)
    {
        int lg_bs = range_idx;
        int bs = 1<<lg_bs;
        int b8_width = 1<<(gd->cb_range[1]-3);
        int sb_width = 1<<(gd->cb_range[1]-range_idx);
        int last_level_flag = range_idx == gd->cb_range[0];
        for (int sb_y = 0; sb_y < sb_width; sb_y++)
        {
            for (int sb_x = 0; sb_x < sb_width; sb_x++)
            {
                int depth_idx = venc_get_depth_scan_idx(sb_x, sb_y, sb_width);
                int cb_idx = range_cb_idx + depth_idx;
                int child_idx = range_cb_idx + sb_width*sb_width + (depth_idx<<2);
                int px = t->ctb_off[0] + sb_x*bs;
                int py = t->ctb_off[1] + sb_y*bs;
                int present_flag = px<gd->pix_dim[0] && py<gd->pix_dim[1];
                int forbidden_flag = present_flag && !last_level_flag && (px+bs>gd->pix_dim[0] || py+bs>gd->pix_dim[1]);

                f265_cb *cb = t->cb + cb_idx;
                cb->lg_bs = lg_bs;
                cb->child_idx = child_idx;
                cb->cb_off[0] = sb_x*bs;
                cb->cb_off[1] = sb_y*bs;
                cb->enc_idx = venc_get_depth_scan_idx(cb->cb_off[0]>>3, cb->cb_off[1]>>3, b8_width);
                cb->flags = 0;
                F265_SET_FLAG(cb->flags, F265_CB_PRESENT, present_flag);
                F265_SET_FLAG(cb->flags, F265_CB_FORBIDDEN|F265_CB_SPLIT, forbidden_flag);
            }
        }

        range_cb_idx += sb_width*sb_width;
    }

    // Execute the segmentation operations.
    if (unlikely(t->seg_flags&F265_SC_CABAC_INIT)) venc_init_cabac_contexts(cbs, f->frame_type, f->qp);
    if (unlikely(t->seg_flags&F265_SC_CABAC_LOAD)) venc_load_cabac_contexts(cbs, &t->wpp_snap);
    if (unlikely(t->init_segment_flag)) venc_init_segment(t, F265_GET_FLAG(t->seg_flags, F265_SC_DEPENDENT));
    if (likely(t->sflags&F265_PF_ALLOW_CTB_RECODE)) venc_save_cabac(cbs, &t->ctb_snap);

    // Set the rate control data.
    //venc_rc_load_ctb(t);
}

// Save the thread CTB data to the frame.
static void venc_save_ctb(f265_enc_thread *t)
{
    // Apply the deblocking filter on full rows. FIXME.
    #if 0
    if (t->ctb_x == ctb_width - 1 && t->deblock_ctx.in_line)
        venc_deblock_ctb_row(&t->deblock_ctx, t->ctb_y);
    #endif

    // Reset the transform node pointer to save the deblocking information.
    t->tt.tn = t->tmap;

    // Save the CBs recursively.
    venc_save_cb(t, t->cb);

    // Save the QP information. Frame QP for now.
    {
        int qp = t->src_frame->qp;
        int b8_size_x = F265_MIN(t->enc->gd.ctb_size, t->enc->gd.pix_dim[0] - t->ctb_off[0])>>3;
        int b8_size_y = F265_MIN(t->enc->gd.ctb_size, t->enc->gd.pix_dim[1] - t->ctb_off[1])>>3;
        int b8_stride = t->enc->gd.pix_dim[0]>>3;
        int8_t *qmap = t->src_frame->qmap + (t->ctb_off[1]>>3)*b8_stride + (t->ctb_off[0]>>3);
        for (int i = 0; i < b8_size_y; i++)
            for (int j = 0; j < b8_size_x; j++)
                qmap[i*b8_stride + j] = qp;
    }
}

// Write and commit the CTB as required. Return F265_RET_RETRY if the analysis
// must be redone. Return F265_RET_VALID on commit.
static int venc_check_ctb(f265_enc_thread *t)
{
    f265_cabac_bs *cbs = &t->cbs;
    int segment_end_flag = F265_GET_FLAG(t->seg_flags, F265_SC_SEGMENT_END);
    int chunk_end_flag = F265_GET_FLAG(t->seg_flags, F265_SC_CHUNK_END);

    // Maximum size of a segment (example: MTU size - NAL packaging overhead). 0
    // if unlimited.
    int max_seg_size = 0;

    // Tentative chunk data.
    f265_seg_chunk chunk = t->chunk;

    // Tentatively write the CTB and the end-of-stream flags in the chunk.
    venc_write_ctb_chunk(t, &chunk, segment_end_flag, chunk_end_flag);

    // The rate control rejects the CTB. Restore the full CABAC state.
    #if 0
    if (unlikely(venc_rc_check_ctb(t) == F265_RET_RETRY))
    #else
    if (0)
    #endif
    {
        venc_load_cabac(cbs, &t->ctb_snap);
        return F265_RET_RETRY;
    }

    // Check if we're busting the maximum segment size. We don't care if we bust
    // and we have only one chunk.
    if (unlikely(max_seg_size && !chunk.start_flag))
    {
        // FIXME: decide whether we want to write the segment header
        // (dependent/independent) at the beginning of the frame and compute the
        // size, adding some bytes for the address, or do that when we start a
        // new segment.
        int hdr_size = 0x1234;

        // Size of the header and the chunk bytestream.
        int seg_size = hdr_size + chunk.cur - chunk.base;

        // If we're not ending the chunk, factor in the outstanding bytes and
        // the bytes affected by a carry operation. Those bytes may become null
        // and require escaping. Add 3 more bytes to end the stream.
        if (!chunk_end_flag) seg_size += 6*(cbs->raw.chunks_outstanding + 1) + 3;

        // Busted!
        if (seg_size > max_seg_size)
        {
            // Restore the full CABAC state and the chunk state.
            venc_load_cabac(cbs, &t->ctb_snap);
            chunk = t->chunk;

            // Write the end-of-stream flags of the last CTB of the chunk and
            // commit. Create a new dependent segment, update the engine
            // snapshot and write the current CTB again to get back on the
            // standard path.
            venc_write_stream_end(cbs, 1, 1);
            venc_write_seg_chunk(cbs, &chunk, 1);
            venc_commit_seg_chunk(t, &chunk, 1, 1);
            venc_init_segment(t, 1);
            venc_save_cabac_engine(cbs, &t->ctb_snap);
            venc_write_ctb_chunk(t, &chunk, segment_end_flag, chunk_end_flag);
        }
    }

    // Commit the chunk.
    venc_commit_seg_chunk(t, &chunk, segment_end_flag, chunk_end_flag);
    t->chunk = chunk;

    // Save the CABAC contexts if required.
    if (unlikely(t->seg_flags&F265_SC_CABAC_SAVE)) venc_save_cabac_contexts(cbs, &t->wpp_snap);

    return F265_RET_VALID;
}

#ifdef VAN_VERIFY_RDO_ANALYSIS
f265_pix venc_verify_analysis_rec[2][3][64*64];
uint8_t venc_verify_analysis_ctx[2][F265_NB_CABAC_CTX];

// Save the reconstruction and CABAC contexts at the location specified.
void venc_verify_analysis_save(f265_enc_thread *t, int loc)
{
    int ref_stride = t->me.ref_stride;
    for (int comp = 0; comp < 3; comp++)
    {
        int csf = !!comp;
        int bs = t->enc->gd.ctb_size>>csf;
        f265_pix *plane = t->src_frame->rec_planes[comp ? 3+comp : 0] + venc_get_ctb_block_plane_off(t, comp, 0, 0);
        venc_copy_block(venc_verify_analysis_rec[loc][comp], bs, plane, ref_stride, bs, bs);
    }

    memcpy(venc_verify_analysis_ctx[loc], t->cbs.contexts, F265_NB_CABAC_CTX);
}

// Verify if the analysis matches the encoding.
void venc_verify_analysis_check(f265_enc_thread *t)
{
    for (int comp = 0; comp < 3; comp++)
    {
        int csf = !!comp;
        int bs = t->enc->gd.ctb_size>>csf;

        for (int y = 0; y < bs; y++)
        {
            for (int x = 0; x < bs; x++)
            {
                int off = y*bs + x;
                int p0 = venc_verify_analysis_rec[0][comp][off];
                int p1 = venc_verify_analysis_rec[1][comp][off];
                if (p0 != p1)
                {
                    printf("Frame %d CTB %d comp %d pixel (%d,%d) mismatch.\n",
                           (int)t->src_frame->abs_poc, t->ctb_xy, comp, x, y);
                    printf("Analysis:\n");
                    print_mb_pix(venc_verify_analysis_rec[0][comp], bs, bs);
                    printf("\n");
                    printf("Encode:\n");
                    print_mb_pix(venc_verify_analysis_rec[1][comp], bs, bs);
                    printf("\n");
                    exit(1);
                }
            }
        }
    }

    for (int i = 0; i < F265_NB_CABAC_CTX; i++)
    {
        int c0 = venc_verify_analysis_ctx[0][i];
        int c1 = venc_verify_analysis_ctx[1][i];
        if (c0 != c1)
        {
            printf("Frame %d CTB %d context 0x%x mismatch (expected %d, analysis %d).\n",
                   (int)t->src_frame->abs_poc, t->ctb_xy, i, c1, c0);
            exit(1);
        }
    }
}
#endif

// Encode the current CTB. Failure is not an option.
static void venc_encode_ctb(f265_enc_thread *t)
{
    // Load the frame CTB data in the thread.
    venc_load_ctb(t);

    // Loop until the CTB is successfully encoded.
    while (1)
    {
        // Choose the encoding mode.
        #ifdef VAN_LOAD_CTB_ANALYSIS
        venc_load_hm_ctb_analysis(t);
        int rec_flag = 1;
        #elif defined(VAN_DUMP_INTRA)
        venc_hm_set_intra(t);
        int rec_flag = 1;
        #else
        int rec_flag = venc_analyze_ctb(t);
        //int rec_flag = venc_hm_an_ctb(t);
        #endif

        #if 0
        // Test HM bit-exact motion estimation and inter search.
        // We keep this code until we completly drop HM bit-exact mode.
        if (t->src_frame->frame_type != F265_FRAME_I) venc_hm_analyze_inter_ctb(t);
        #endif

        // Reconstruct the CTB if this was not done during the analysis.
        if (rec_flag) venc_rec_ctb(t);

        // Check if the CTB is acceptable. If not, redo the analysis. Otherwise,
        // commit.
        if (venc_check_ctb(t) != F265_RET_RETRY) break;
    }

    #ifdef VAN_VERIFY_RDO_ANALYSIS
    venc_verify_analysis_save(t, 1);
    venc_verify_analysis_check(t);
    #endif

    // Save the thread CTB data to the frame.
    venc_save_ctb(t);
}

// Function to determine the lambdas used for RDO.
// The QP factor modifies the lambda (but not the QP) of the current frame.
// The QP offset modifies the QP (and the lambda) of the current frame.
static double venc_calc_lambda(f265_gen_data *gd, f265_frame *f, int8_t *qp)
{
    int gop = gd->nb_b_frames + 1;
    int poc = f->abs_poc;
    double qp_temp = (double)*qp;

    // Calculate depth of frame in the GOP.
    int idx = poc % gop;
    int depth = 0;
    if (idx > 0)
        for (int step = gop, i = gop>>1; i >= 1; i >>= 1, step >>= 1, depth++)
            for (int j = i; j < gop; j += step) if (j == poc) { i=0; break; }

    // Compute QP factor.
    double qp_fac = 1.0;
    if (f->frame_type == F265_FRAME_I)
    {
        // Compute scale.
        double scale = 0.05 * (double)gd->nb_b_frames;
        scale = 1.0 - F265_CLAMP(scale, 0.0, 0.5);
        qp_fac = 0.57 * scale;
    }
    else if (gd->hm_gop_size)
    {
        qp_fac = f->gop_entry->qp_factor;

        // Add QP offset.
        qp_temp += f->gop_entry->qp_offset;
    }

    // Update QP.
    *qp = F265_CLAMP((int8_t)floor(qp_temp + 0.5), 0, 51);

    // Compute lambda.
    qp_temp -= 12;
    double lambda = pow(2.0, qp_temp/3.0) * qp_fac;
    if (depth > 0) lambda *= F265_CLAMP(qp_temp / 6.0, 2.0, 4.0);

    // FIXME This hard-coding should be replaced by the value of the HADME flag obtained from the config file.
    uint8_t isHADME = 1;
    if (isHADME == 0 && f->frame_type != F265_FRAME_I) lambda *= 0.95;

    return lambda;
}

// Function to determine the lambdas used for chroma RDO.
static double venc_calc_lambda_chroma(int chroma_qp, int luma_qp, double luma_lambda)
{
    int delta_qp = chroma_qp - luma_qp;
    double weight = pow(2.0, delta_qp/3.0);
    return luma_lambda / weight;
}

// Assign the frame to the next idle encoding thread and prepare encoding.
// Return the thread.
static f265_enc_thread* venc_set_enc_thread(f265_enc *e, f265_frame *f)
{
    f265_main_data *md = &e->md;

    // Set the DPB.
    venc_set_frame_dpb(e, f);

    // Link the encoding object.
    venc_link_enc_obj(e, f);

    // Get the previous frame.
    f265_frame *prev = md->enc_frames[0];

    // We have one more encoding thread active.
    f265_enc_thread *t = md->enc_threads[md->nb_enc_frames++];

    // Assign the frame.
    t->src_frame = f;
    f->frame_type = f->la_frame_type;

    // Assign the flags tentatively.
    t->mflags = t->sflags = e->gd.eflags;

    // Set the frame deblocking flag.
    int deblock_flag = t->mflags&F265_PF_DEBLOCK && (t->mflags&F265_PF_DEBLOCK_ALL || f->gen_flags&F265_FF_REF);
    F265_SET_FLAG(f->gen_flags, F265_FF_DEBLOCK, deblock_flag);

    // Set the frame QP.
    f->qp = venc_rc_frame_start(t, prev);

    // Set the QP.
    t->qp[0] = f->qp;
    t->qp[1] = e->gd.chroma_qp_table[f->qp];

    // Set the RDO values.
    // The frame QP is modified if the QP offset is set.
    t->hm_lambda[0] = venc_calc_lambda(&e->gd, f, &f->qp);
    t->hm_lambda[1] = venc_calc_lambda_chroma(t->qp[0],t->qp[1], t->hm_lambda[0]);

    // Set the chroma weight.
    t->hm_wcd = pow(2.0, (t->qp[0] - t->qp[1]) / 3.0);

    // Build the inter configuration.
    venc_set_frame_inter(t, f);

    return t;
}

// Check the encoding result and take the appropriate action. On failure, dipose
// of the next frame to encode, if any. Return VALID, RETRY, ABORT, ERROR.
static int venc_check_enc_result(f265_enc_thread *t, f265_frame *next_frame)
{
    f265_enc *e = t->enc;
    f265_gen_data *gd = &e->gd;
    f265_main_data *md = &e->md;
    int res = F265_RET_VALID;

    do
    {
        // Abort if requested.
        if (md->enc_exit_flag)
        {
            res = F265_RET_ABORT;
            break;
        }

        // Encode the NALs.
        venc_encode_frame_nals(e, t->src_frame, t->chunk_buf);

        // We have a valid result. Ask the rate control about it.
        res = venc_rc_check_frame(t);

        // The rate control chocked or wants to retry.
        if (res != F265_RET_VALID) break;

        // The rate control wants to commit.
        // Release the encoding thread.
        md->nb_enc_frames--;

        // Clear the last encoded frame and release memory.
        {
            f265_frame *f = md->enc_frames[0];
            if (f)
            {
                F265_SET_FLAG(f->main_flags, F265_FF_LAST_ENC, 0);
                venc_release_frame_mem(e, f);
            }
        }

        // Set the last encoded frame and release memory.
        f265_frame *f = t->src_frame;
        F265_SET_FLAG(f->main_flags, F265_FF_ENC|F265_FF_LAST_ENC, 1);
        md->enc_frames[0] = f;
        venc_release_frame_mem(e, f);

        // Deblock the frame.
        if (f->gen_flags&F265_FF_DEBLOCK) venc_deblock_frame(t);

        #ifdef F265_HAVE_STDIO
        f265_pix *rec_planes[3] = { f->rec_planes[0], f->rec_planes[4], f->rec_planes[5] };
        if (gd->yuv_dump_file) venc_dump_yuv_file(gd->yuv_dump_file, rec_planes, f->abs_poc, gd->stride,
                                                  gd->clip_dim[0], gd->clip_dim[1], gd->csf[0], gd->csf[1]);
        #endif

        // Handle padding and the halfpel interpolation.
        if (gd->mt_mode != F265_MT_ENC_FRAME)
        {
            // Pad and interpolate.
            if (f->gen_flags&F265_FF_REF)
            {
                // Pad the YUV planes. This must be done prior to the halfpel
                // interpolation.
                venc_pad_rec_yuv_planes(t);

                // Perform the halfpel interpolation.
                if (t->sflags&F265_PF_SUBPEL_ME) venc_compute_hpel_planes(t);
            }
        }

    } while (0);

    // Release the next frame if we won't encode it.
    if (next_frame && (res == F265_RET_ABORT || res == F265_RET_ERROR))
        venc_release_frame_mem(e, next_frame);

    return res;
}

