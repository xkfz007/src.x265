// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Block analysis prior to encoding.

// FIXME:
// Multiply the distortion itself so you can elimiate a shift for the lambda
// multiplication and the chroma weighting factor.

// Compile the encoding functions in RDO mode.
#define F265_RDO
#include "entropy.c"


///////////////////////////////////////////////////////////////////////////////
// Temporary debugging code.

// Print the analysis costs.
//#define VAN_TRACE_ANALYSIS

// Print the transform block pixels during the analysis.
//#define VAN_TRACE_TB_ANALYSIS

// Print the stash operations.
//#define VAN_TRACE_STASH

// Verify the CB inter cost in RDM mode. No effect in RDO mode.
//#define VAN_VERIFY_CB_INTER_COST

// Perform the chroma analysis. If not defined, use the luma mode for chroma.
#define VAN_CHROMA_ANALYSIS

// Add more precision to the bit costs. Turn off to ease debugging.
#define VAN_SCALE_UP_PRECISION

#ifdef VAN_TRACE_ANALYSIS
// Set to 0 to disable tracing.
int venc_trace_analysis_flag = 0;
#endif

void venc_an_cb_loc(f265_enc_thread *t, f265_cb *cb, int comp)
{
    int csf = !!comp;
    printf("CB (%d,%d) size %d pix (%d,%d): ",
           cb->cb_off[0]>>csf, cb->cb_off[1]>>csf, (1<<cb->lg_bs)>>csf,
           (t->ctb_off[0]+cb->cb_off[0])>>csf, (t->ctb_off[1]+cb->cb_off[1])>>csf);
}

void venc_an_tt_loc(f265_enc_thread *t, f265_cb *cb, int comp, int cb_ox, int cb_oy, int lg_bs)
{
    int csf = !!comp;
    printf("CB (%d,%d) TT (%d,%d) size %d pix (%d,%d): ",
           cb->cb_off[0]>>csf, cb->cb_off[1]>>csf, cb_ox, cb_oy, (1<<lg_bs),
           ((t->ctb_off[0]+cb->cb_off[0])>>csf)+cb_ox, ((t->ctb_off[1]+cb->cb_off[1])>>csf)+cb_oy);
}

void venc_an_print_tn(f265_enc_thread *t, int spaces)
{
    int split_flag = !!((*t->tt.tn++)&8);
    for (int i = 0; i < spaces; i++) printf(" ");
    printf("%d\n", split_flag);
    if (split_flag)
        for (int i = 0; i < 4; i++)
            venc_an_print_tn(t, spaces+1);
}

void venc_an_print_tt(f265_enc_thread *t)
{
    uint8_t *saved_tn = t->tt.tn;
    t->tt.tn = t->tmap;
    printf("==== Transform tree ====\n");
    venc_an_print_tn(t, 0);
    printf("Last node at offset %d, current node at offset %d.\n",
           (int)(saved_tn - t->tmap), (int)(t->tt.tn - t->tmap));

    printf("========================\n");
    t->tt.tn = saved_tn;
}

#ifdef VAN_TRACE_STASH
void venc_stash_trace(f265_enc_thread *t, const char *what)
{
    f265_stash *s = &t->stash;
    printf("* %s, offsets (%d,%d,%d) tn_off %d block (%d,%d,%d,%d).\n",
           what, s->stash_off[0], s->stash_off[1], s->stash_off[2], s->tn_off, s->px, s->py, s->width, s->height);
}
#endif


///////////////////////////////////////////////////////////////////////////////


// Number of bytes to copy in a stash frame.
#define F265_STASH_FRAME_SIZE (8*4)

// Save or load the CABAC contexts from the stash. We might want to refine this
// eventually, e.g. we'll never need the SAO contexts.
static inline void venc_stash_copy_cabac(f265_enc_thread *t, uint8_t **wp, int save_flag)
{
    if (save_flag) memcpy(*wp, t->cbs.contexts, F265_NB_CABAC_CTX);
    else memcpy(t->cbs.contexts, *wp, F265_NB_CABAC_CTX);
    *wp += F265_NB_CABAC_CTX;
}

// Save or load the transform nodes from the stash.
static inline void venc_stash_copy_tn(f265_enc_thread *t, uint8_t **wp, int save_flag)
{
    f265_stash *s = &t->stash;

    if (save_flag)
    {
        uint8_t *tn = t->tt.tn;
        int nb_tn = (tn - t->tmap) - s->tn_off;
        memcpy(*wp, &nb_tn, 4); *wp += 4;
        memcpy(*wp, tn - nb_tn, nb_tn); *wp += nb_tn;
    }

    else
    {
        uint8_t *tn = t->tmap + s->tn_off;
        int nb_tn;
        memcpy(&nb_tn, *wp, 4); *wp += 4;
        memcpy(tn, *wp, nb_tn); *wp += nb_tn;
        t->tt.tn = tn + nb_tn;
    }
}

// Save or load the reconstructed pixels from the stash.
static inline void venc_stash_copy_rec(f265_enc_thread *t, uint8_t **wp, int save_flag, int flags)
{
    f265_stash *s = &t->stash;
    for (int comp = 0; comp < 3; comp++)
    {
        if (!(flags&(1<<comp))) continue;

        int csf = !!comp;
        int px = s->px>>csf;
        int py = s->py>>csf;
        int width = s->width>>csf;
        int height = s->height>>csf;
        int ref_stride = t->me.ref_stride;
        f265_pix *rec_plane = t->src_frame->rec_planes[comp ? 3+comp : 0] + py*ref_stride + px;
        f265_pix *stash_rec = (f265_pix*)*wp;

        if (save_flag) venc_copy_block(stash_rec, width, rec_plane, ref_stride, width, height);
        else venc_copy_block(rec_plane, ref_stride, stash_rec, width, width, height);
        *wp += width*height*F265_PSIZE;
    }
}

// Push a new stash frame.
void venc_stash_push(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "pushing");
    #endif
    f265_stash *s = &t->stash;
    memcpy(s->buf + s->stash_off[2], &s->stash_off[0], F265_STASH_FRAME_SIZE);
    s->stash_off[0] = s->stash_off[1] = s->stash_off[2] = s->stash_off[2] + F265_STASH_FRAME_SIZE;
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "pushed");
    #endif
}

// Pop the current stash frame.
void venc_stash_pop(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "popping");
    #endif
    f265_stash *s = &t->stash;
    memcpy(&s->stash_off[0], s->buf + s->stash_off[0] - F265_STASH_FRAME_SIZE, F265_STASH_FRAME_SIZE);
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "popped");
    #endif
}

// Save the initial state of the block. The block location is specified in
// terms of luma pixels.
void venc_stash_init(f265_enc_thread *t, int px, int py, int width, int height)
{
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "initializing");
    #endif
    f265_stash *s = &t->stash;

    // Set the block location.
    s->px = px;
    s->py = py;
    s->width = width;
    s->height = height;

    // Save the transform node offset.
    s->tn_off = t->tt.tn - t->tmap;

    // Save the CABAC contexts.
    if (s->flags&16)
    {
        uint8_t *wp = s->buf + s->stash_off[0];
        venc_stash_copy_cabac(t, &wp, 1);
        s->stash_off[1] = s->stash_off[2] = wp - s->buf;
    }
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "initialized");
    #endif
}

// Frontend for venc_stash_init(). Set the luma block location given the CB and
// the luma block offsets from the CB origin.
void venc_stash_init_cb_off(f265_enc_thread *t, f265_cb *cb, int cb_ox, int cb_oy, int width, int height)
{
    venc_stash_init(t,
                    t->ctb_off[0] + cb->cb_off[0] + cb_ox,
                    t->ctb_off[1] + cb->cb_off[1] + cb_oy,
                    width, height);
}

// Save as above without the offsets.
void venc_stash_init_cb(f265_enc_thread *t, f265_cb *cb)
{
    venc_stash_init_cb_off(t, cb, 0, 0, 1<<cb->lg_bs, 1<<cb->lg_bs);
}

// Rollback to the initial state.
void venc_stash_reset(f265_enc_thread *t)
{
    f265_stash *s = &t->stash;
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "reset");
    #endif

    // Load the transform node pointer.
    t->tt.tn = t->tmap + s->tn_off;

    // Load the CABAC contexts.
    if (s->flags&16)
    {
        uint8_t *wp = s->buf + s->stash_off[0];
        venc_stash_copy_cabac(t, &wp, 0);
    }
}

// Save the analyzed state of the block. The previously saved state is
// discarded, if any.
void venc_stash_save(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "saving");
    #endif
    f265_stash *s = &t->stash;
    uint32_t flags = s->flags;
    uint8_t *wp = s->buf + s->stash_off[1];

    if (flags&8) venc_stash_copy_tn(t, &wp, 1);
    if (flags&7) venc_stash_copy_rec(t, &wp, 1, flags);
    if (flags&16) venc_stash_copy_cabac(t, &wp, 1);

    s->stash_off[2] = wp - s->buf;
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "saved");
    #endif
}

// Rollback to the saved state.
void venc_stash_restore(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_STASH
    venc_stash_trace(t, "restored");
    #endif
    f265_stash *s = &t->stash;
    uint32_t flags = s->flags;
    uint8_t *wp = s->buf + s->stash_off[1];

    if (flags&8) venc_stash_copy_tn(t, &wp, 0);
    if (flags&7) venc_stash_copy_rec(t, &wp, 0, flags);
    if (flags&16) venc_stash_copy_cabac(t, &wp, 0);
}

// Convenience function for save and reset.
void venc_stash_save_reset(f265_enc_thread *t)
{
    venc_stash_save(t);
    venc_stash_reset(t);
}

// Set the transform block pointers to the beginning of the arrays to discard
// the data after the encoding. We might want to optimize this eventually.
void venc_set_tmp_tb(f265_enc_thread *t)
{
    f265_tt_enc *tt = &t->tt;
    tt->tb = t->tb;
    tt->sb = t->sb;
    tt->levels = t->levels;
}

// Kludge to compute the cost of a syntax element. Fix this eventually.
#define F265_RDO_COST(var, eval) t->cbs.rdo.bits = 0; eval; var = venc_rdo_bit_cost(t, t->cbs.rdo.bits)

// Compute the RDO cost of the entropy bits specified (using 32768 bit
// fractions for now).
static int64_t venc_rdo_bit_cost(f265_enc_thread *t, int rdo_bits)
{
    // 64-bit arithmetic to avoid overflows before shift.
    int64_t rdo_lambda = t->an.rdo_lambda;
    #ifdef VAN_SCALE_UP_PRECISION
    return (rdo_bits*rdo_lambda)>>(15+8-5);
    #else
    return (rdo_bits*rdo_lambda)>>(15+8);
    #endif
}

// Return the cost of encoding the context bin specified without encoding it.
static finline int64_t venc_get_context_bin_cost(f265_enc_thread *t, int ctx_idx, int bin)
{
    return venc_rdo_bit_cost(t, f265_cabac_entropy_table[t->cbs.contexts[ctx_idx]^bin]);
}

// Get the minimum and maximum tested depths for an intra partition transform
// tree. All transform blocks within those bounds are tested.
static void venc_get_intra_part_depth_range(f265_enc_thread *t, int lg_bs, int depth_range[2])
{
    int min_lg_bs = t->an.tb_range[0], max_lg_bs = t->an.tb_range[1], max_depth = t->an.tb_depth[0];
    depth_range[0] = F265_MAX(0, lg_bs-max_lg_bs);
    depth_range[1] = F265_MAX(depth_range[0], F265_MIN(lg_bs-min_lg_bs, max_depth));
}

// Same as above for inter.
static void venc_get_inter_part_depth_range(f265_enc_thread *t, int inter_part, int lg_bs, int depth_range[2])
{
    int min_lg_bs = t->an.tb_range[0], max_lg_bs = t->an.tb_range[1], max_depth = t->an.tb_depth[1];
    depth_range[0] = F265_MAX(0, lg_bs-max_lg_bs);
    depth_range[1] = F265_MAX(depth_range[0], F265_MIN(lg_bs-min_lg_bs, max_depth));

    // If the partitioning mode is not UN and the maximum depth is 0, force at
    // least one split level.
    if (inter_part != F265_PART_UN && !max_depth)
    {
        depth_range[0] = F265_MAX(1, depth_range[0]);
        depth_range[1] = F265_MAX(1, depth_range[1]);
    }
}

// Return the SSD of the block specified.
static int64_t venc_an_block_ssd(f265_enc_thread *t, f265_pix *src0, int src0_stride, f265_pix *src1, int src1_stride,
                                 int bs, int comp)
{
    // Compute the reconstruction distortion.
    int64_t dist = venc_ssd(src0, src0_stride, src1, src1_stride, bs, bs, 0, 8);
    #ifdef VAN_SCALE_UP_PRECISION
    dist <<= 5;
    #endif

    // Scale the chroma distortion.
    if (comp) dist *= t->hm_wcd;

    return dist;
}

//#define VAN_INTRA_TT_EXPERIMENT
#ifdef VAN_INTRA_TT_EXPERIMENT
f265_pix pred_buf[4][64*64];
#endif

// Return the cost of the intra transform block specified. Set the non-zero
// flag.
static int64_t venc_analyze_intra_tb(f265_enc_thread *t, f265_cb *cb, int *nz_flag, int mode, int comp,
                                     int lg_bs, int cb_ox, int cb_oy)
{
    f265_pix nbuf[2][129], pred[32*32];
    int csf = !!comp;
    int ct_ox = (cb->cb_off[0]>>csf) + cb_ox;
    int ct_oy = (cb->cb_off[1]>>csf) + cb_oy;
    int bs = 1<<lg_bs;
    int ref_stride = t->me.ref_stride;
    int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox, ct_oy);
    int filter_edge_flag, filter_neighbour_flag, neighbour_bilinear_flag;
    int64_t cost;
    f265_pix *neighbours;

    int fast_neighbour_flag = t->an.rdm_flag;

    // Predict the neighbours.
    if (fast_neighbour_flag)
    {
        // Do not filter anything when using the source pixels. Depth 0
        // optimization here eventually (compute those neighbours once).
        venc_predict_intra_neighbours(t, nbuf[0], 0, comp, bs, ct_ox, ct_oy);
        filter_edge_flag = 0;
        neighbours = nbuf[0];
    }

    else
    {
        int smooth_intra_flag = F265_GET_FLAG(t->enc->gd.eflags, F265_PF_SMOOTH_INTRA);
        venc_get_intra_filter_flags(&filter_edge_flag, &filter_neighbour_flag, &neighbour_bilinear_flag,
                                    comp, lg_bs, mode, smooth_intra_flag);
        venc_predict_intra_neighbours(t, nbuf[0], 1, comp, bs, ct_ox, ct_oy);
        if (filter_neighbour_flag) venc_filter_intra_neighbours(nbuf[1], nbuf[0], bs, 8, neighbour_bilinear_flag);
        neighbours = nbuf[filter_neighbour_flag];
    }

    // Predict the block pixels.
    venc_predict_intra_mode(pred, neighbours, lg_bs, 8, mode, filter_edge_flag);

    // Compute the cost.
    if (t->an.rdm_flag)
    {
        // Compute the prediction distortion.
        f265_pix *src = t->src_frame->src_planes[comp] + plane_off;
        cost = t->me.dist[t->an.mode_metric](src, ref_stride, pred, bs, bs, bs, 8);

        // Avoid Valgrind errors. FIXME.
        *nz_flag = 1;
    }

    else
    {

#ifdef VAN_INTRA_TT_EXPERIMENT
if (!comp)
{
for (int i = 0; i < bs; i++)
    for (int j = 0; j < bs; j++)
        pred_buf[lg_bs-2][(cb_oy+i)*64+cb_ox+j] = pred[i*bs+j];

// Experimental early termination.
if (lg_bs > 2)
{
    // Where the prediction samples will end up.
    int offset = (bs>>1)-1;

    int sad = 0;

    // Row SAD.
    for (int i = 0; i < bs; i++)
        sad += F265_ABS(pred_buf[lg_bs-3][(cb_oy+offset)*64+cb_ox+i] - pred_buf[lg_bs-2][(cb_oy+offset)*64+cb_ox+i]);

    // Column SAD.
    for (int i = 0; i < bs; i++)
        sad += F265_ABS(pred_buf[lg_bs-3][(cb_oy+i)*64+cb_ox+offset] - pred_buf[lg_bs-2][(cb_oy+i)*64+cb_ox+offset]);

    // The idea is that the greater the SAD, the less interesting
    // the bigger prediction signal should be.
    int threshold = 4*bs;
printf("sad=%d, threshold=%d\n", sad, threshold);
    if (sad > threshold) printf("split should win: ");
    else printf("unsplit should win: ");
}

}
#endif
        // Compute the reconstruction.
        int dst_flag, order;
        venc_get_intra_encode_flags(&dst_flag, &order, comp, lg_bs, mode);
        venc_set_tmp_tb(t);
        *nz_flag = venc_rec_tb(t, pred, 1<<lg_bs, comp, lg_bs, dst_flag, order, 0, ct_ox, ct_oy, 0, 1, 0);

        // Compute the encoding cost.
        venc_set_tmp_tb(t);
        int64_t coeff_cost = 0;
        if (*nz_flag) { F265_RDO_COST(coeff_cost, venc_write_tb(&t->cbs, &t->tt, !!comp)); }

        // Compute the reconstruction distortion.
        f265_pix *src = t->src_frame->src_planes[comp] + plane_off;
        f265_pix *rec = t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off;
        int64_t dist = venc_an_block_ssd(t, rec, ref_stride, src, ref_stride, bs, comp);

        cost = coeff_cost + dist;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, comp, cb_ox, cb_oy, lg_bs);
            printf("encoded TB cost %d (coeffs %d, dist %d), nz %d.\n",
                   (int)cost, (int)coeff_cost, (int)dist, *nz_flag);
        }
        #endif
    }

    return cost;
}

// Return the cost of the intra chroma transform tree. Set the two chroma
// non-zero flags.
static int64_t venc_analyze_intra_chroma_tt(f265_enc_thread *t, f265_cb *cb, int *chroma_nz, int mode, int depth,
                                            int lg_bs, int cb_ox, int cb_oy)
{
    int64_t cost;
    int tt_nz = 0;

    // Get the split flag. Skip the 4x4 subblocks if present.
    int split_flag = !!((*t->tt.tn++)&8);
    if (split_flag && lg_bs == 2)
    {
        split_flag = 0;
        t->tt.tn += 4;
    }

    if (split_flag)
    {
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 1, cb_ox, cb_oy, lg_bs);
            printf("analyzing split TT depth %d.\n", depth);
        }
        #endif

        // Packed subtree non-zero flags, in order (VU VU VU VU).
        int packed_nz_flags = 0;

        // Analyze the subtree.
        int64_t subtree_cost = 0;
        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
        {
            int nz_flags;
            subtree_cost += venc_analyze_intra_chroma_tt(t, cb, &nz_flags, mode, depth+1, lg_bs-1,
                                                         cb_ox + (i&1)*sbs, cb_oy + (i>>1)*sbs);
            packed_nz_flags |= nz_flags<<(i<<1);
        }

        // Encode the subtree CBFs when at least one component is non-zero. The
        // CBFs at different depths use different contexts, so in RDO mode they
        // can be encoded in any order. The UV flags are interleaved, so we
        // process each component following this order.
        int64_t cbf_cost = 0;
        if (!t->an.rdm_flag && packed_nz_flags)
        {
            // Determine the value of the chroma non-zero flags.
            int comp_nz[2];
            for (int i = 0; i < 2; i++)
            {
                comp_nz[i] = !!((packed_nz_flags>>i)&0x55);
                tt_nz |= comp_nz[i]<<i;
            }

            // Encode the CBFs in order.
            t->cbs.rdo.bits = 0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 2; j++)
                {
                    // Skip the component when it's all zero.
                    if (!comp_nz[j]) continue;

                    int cbf = (packed_nz_flags>>((i<<1)+j))&1;
                    venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+depth+1, cbf);
                }
            cbf_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
        }

        cost = cbf_cost + subtree_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 1, cb_ox, cb_oy, lg_bs);
            printf("analyzed split TT depth %d cost %d (cbf %d, subtree %d).\n",
                   depth, (int)cost, (int)cbf_cost, (int)subtree_cost);
        }
        #endif
    }

    else
    {
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 1, cb_ox, cb_oy, lg_bs);
            printf("analyzing unsplit TT depth %d.\n", depth);
        }
        #endif

        // Encode each transform block.
        cost = 0;
        for (int i = 0; i < 2; i++)
        {
            int nz_flag;
            cost += venc_analyze_intra_tb(t, cb, &nz_flag, mode, 1+i, lg_bs, cb_ox, cb_oy);
            tt_nz |= nz_flag<<i;
        }

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 1, cb_ox, cb_oy, lg_bs);
            printf("analyzed unsplit TT depth %d TB %d.\n", depth, (int)cost);
        }
        #endif
    }

    *chroma_nz = tt_nz;
    return cost;
}

// Return the cost of the chroma intra block. The function depends on the luma
// transform nodes being set correctly. The transform tree pointer is unaffected
// by this function.
static int64_t venc_analyze_intra_cb_chroma(f265_enc_thread *t, f265_cb *cb, uint8_t *init_tn)
{
    #ifndef VAN_CHROMA_ANALYSIS
    cb->intra_chroma_mode = cb->intra_luma_mode[0];
    return 0;
    #endif

    int luma_mode = cb->intra_luma_mode[0];
    int best_mode;
    int64_t best_cost = F265_MAX_SSD;
    int stash_flag = !!t->stash.flags;

    if (stash_flag) venc_stash_init_cb(t, cb);

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 1);
        printf("analyzing chroma modes.\n");
    }
    #endif

    // Test every mode.
    for (int i = 0; i < 5; i++)
    {
        int mode;
        int64_t mode_cost;

        // Test the luma mode first since it is the most probable.
        if (!i)
        {
            mode = luma_mode;

            if (t->an.rdm_flag) mode_cost = t->an.se_costs[F265_SE_INTRA_CHROMA_MODE+0];
            else { F265_RDO_COST(mode_cost, venc_encode_context_bin(&t->cbs, F265_CO_INTRA_CHROMA_PRED, 0)); }
        }

        // Test the standard chroma modes.
        else
        {
            mode = f265_chroma_mode_table[i-1];
            if (mode == luma_mode) mode = 34;

            if (t->an.rdm_flag) mode_cost = t->an.se_costs[F265_SE_INTRA_CHROMA_MODE+1];
            else
            {
                t->cbs.rdo.bits = 0;
                venc_encode_context_bin(&t->cbs, F265_CO_INTRA_CHROMA_PRED, 1);
                venc_encode_bypass_bins(&t->cbs, 0, 2);
                mode_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
            }
        }

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 1);
            printf("testing mode %d.\n", mode);
        }
        #endif

        // Set the current transform node.
        t->tt.tn = init_tn;

        // Analyze the transform tree.
        int chroma_nz;
        int64_t tt_cost = venc_analyze_intra_chroma_tt(t, cb, &chroma_nz, mode, 0, cb->lg_bs-1, 0, 0);

        // Encode the two chroma CBFs.
        int64_t cbf_cost = 0;
        if (!t->an.rdm_flag)
        {
            t->cbs.rdo.bits = 0;
            for (int i = 0; i < 2; i++)
                venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+0, (chroma_nz>>i)&1);
            cbf_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
        }

        int64_t cost = mode_cost + cbf_cost + tt_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 1);
            printf("analyzed mode %d cost %d (mode %d, CBF %d, subtree %d).\n",
                   mode, (int)cost, (int)mode_cost, (int)cbf_cost, (int)tt_cost);
        }
        #endif

        if (cost < best_cost)
        {
            if (stash_flag) venc_stash_save(t);
            best_mode = mode;
            best_cost = cost;
        }

        if (stash_flag) venc_stash_reset(t);
    }

    if (stash_flag) venc_stash_restore(t);
    cb->intra_chroma_mode = best_mode;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 1);
        printf("kept chroma mode %d cost %d.\n", best_mode, (int)best_cost);
    }
    #endif

    return best_cost;
}

// Return the cost of the intra luma transform tree. This function currently
// stores the transform nodes unconditionally.
static int64_t venc_analyze_intra_luma_tt(f265_enc_thread *t, f265_cb *cb, int mode, int depth_range[2], int depth,
                                          int lg_bs, int cb_ox, int cb_oy)
{
    // FIXME. For the prediction distortion case, when we change the TT layout
    // parameters, we might want to account for the skewed flag costs.
    int split_flag = depth < depth_range[1];
    int unsplit_flag = depth >= depth_range[0];
    int split_present_flag = split_flag && unsplit_flag;
    int stash_flag = t->stash.flags && split_present_flag;
    int64_t split_cost = F265_MAX_SSD, unsplit_cost = F265_MAX_SSD, best_cost;

    #ifdef VAN_INTRA_TT_EXPERIMENT
    printf("venc_analyze_intra_luma_tt lg_bs=%d\n", lg_bs);
    #endif

    if (stash_flag) venc_stash_init_cb_off(t, cb, cb_ox, cb_oy, 1<<lg_bs, 1<<lg_bs);

    // Get the split flag costs.
    uint16_t split_flag_costs[2] = { 0, 0 };
    #if 0
    if (split_present_flag)
    #else
    if (t->an.rdm_flag && split_present_flag)
    #endif
    {
        int off = (5-lg_bs)<<1;
        for (int i = 0; i < 2; i++) split_flag_costs[i] = t->an.se_costs[F265_SE_SPLIT_TRANSFORM + off + i];
    }

    // FIXME: consider extracting this code in a function.
    if (split_flag)
    {
        // Add a split transform node.
        *t->tt.tn++ = 15;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzing split TT depth %d.\n", depth);
        }
        #endif

        int64_t split_flag_cost = 0;
        if (t->an.rdm_flag) split_flag_cost = split_flag_costs[1];
        else if (split_present_flag)
        {
            F265_RDO_COST(split_flag_cost,
                          venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 1));
        }

        venc_stash_push(t);
        int64_t subtree_cost = 0;
        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
            subtree_cost += venc_analyze_intra_luma_tt(t, cb, mode, depth_range, depth+1, lg_bs-1,
                                                       cb_ox + (i&1)*sbs, cb_oy + (i>>1)*sbs);

#ifdef VAN_INTRA_TT_EXPERIMENT
//{
//        int block_size = 1 << lg_bs;
//        for (int i = 0; i < block_size; i++)
//            for (int j = 0; j < block_size; j++)
//                printf("%3d%s", pred_buf[lg_bs-3][(cb_oy+i)*64+cb_ox+j], j==(block_size-1) ? "\n" : " ");
//}
#endif
        venc_stash_pop(t);

        split_cost = split_flag_cost + subtree_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzed split TT depth %d cost %d (split_flag %d, subtree %d).\n",
                   depth, (int)split_cost, (int)split_flag_cost, (int)subtree_cost);
            #ifdef VAN_TRACE_TB_ANALYSIS
            {
                printf("======> Split part result.\n");
                int plane_off = venc_get_cb_block_plane_off(t, cb, 0, cb_ox, cb_oy);
                f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
                print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
            }
            #endif
        }
        #endif
    }

    if (unsplit_flag)
    {
        if (stash_flag) venc_stash_save_reset(t);

        // Add an unsplit transform node.
        *t->tt.tn++ = 7;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzing unsplit TT depth %d.\n", depth);
        }
        #endif

        // For the prediction distortion case, we have to estimate the cost of
        // splitting the transform tree. The prediction distortion alone won't
        // account for the overhead of syntax elements like last_sig_coeff.
        //
        // For now, use 2 bypass bins per TB.
        int64_t unsplit_flag_cost = 0;
        int64_t fudge = 0;
        if (t->an.rdm_flag)
        {
            unsplit_flag_cost = split_flag_costs[0];
            fudge = 2*t->an.se_costs[F265_SE_BYPASS];
        }
        else if (split_present_flag)
        {
            F265_RDO_COST(unsplit_flag_cost,
                          venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 0));
        }

        int nz_flag;
        int64_t tb_cost = venc_analyze_intra_tb(t, cb, &nz_flag, mode, 0, lg_bs, cb_ox, cb_oy);
#ifdef VAN_INTRA_TT_EXPERIMENT
//{
//        int block_size = 1 << lg_bs;
//        for (int i = 0; i < block_size; i++)
//            for (int j = 0; j < block_size; j++)
//                printf("%3d%s", pred_buf[lg_bs-2][(cb_oy+i)*64+cb_ox+j], j==(block_size-1) ? "\n" : " ");
//}
#endif

        int64_t cbf_cost = 0;
        if (!t->an.rdm_flag)
        {
            F265_RDO_COST(cbf_cost, venc_encode_context_bin(&t->cbs, F265_CO_CBF_LUMA + (cb->lg_bs==lg_bs), nz_flag));
        }

        unsplit_cost = unsplit_flag_cost + cbf_cost + fudge + tb_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzed unsplit TT depth %d cost %d (split_flag %d, cbf %d, fudge %d, TB %d).\n",
                   depth, (int)unsplit_cost, (int)unsplit_flag_cost, (int)cbf_cost, (int)fudge, (int)tb_cost);
            #ifdef VAN_TRACE_TB_ANALYSIS
            {
                printf("======> Unsplit part result.\n");
                int plane_off = venc_get_cb_block_plane_off(t, cb, 0, cb_ox, cb_oy);
                f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
                print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
            }
            #endif
        }
        #endif
    }

    if (split_cost < unsplit_cost)
    {
        if (stash_flag) venc_stash_restore(t);
        best_cost = split_cost;
#ifdef VAN_INTRA_TT_EXPERIMENT
printf("split won\n");
{
        // Transfer the prediction signal one level up.
        int block_size = 1 << lg_bs;
        for (int i = 0; i < block_size; i++)
            for (int j = 0; j < block_size; j++)
                pred_buf[lg_bs-2][(cb_oy+i)*64+cb_ox+j] = pred_buf[lg_bs-3][(cb_oy+i)*64+cb_ox+j];
}
#endif
    }

    else
    {
        best_cost = unsplit_cost;
#ifdef VAN_INTRA_TT_EXPERIMENT
printf("unsplit won\n");
#endif
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("analyzed TT depth %d best %d (split %d, unsplit %d).\n",
               depth, (int)best_cost, (int)split_cost, (int)unsplit_cost);
    }
    #endif

    return best_cost;
}

// Compute the rough cost of using a specific intra prediction mode.
int64_t venc_analyze_intra_rough_cost(f265_enc_thread *t, f265_cb *cb, int part_idx, int lg_bs,
                                      int cb_ox, int cb_oy, int mode)
{
    int nz_flag = 0;
    int64_t cost = 0;
    if (lg_bs == 6)
    {
        for (int y = 0; y < 64; y += 32)
            for (int x = 0; x < 64; x += 32)
                cost += venc_analyze_intra_tb(t, cb, &nz_flag, mode, 0, 5, cb_ox+x, cb_oy+y);
    }

    else if (lg_bs > 2)
        cost = venc_analyze_intra_tb(t, cb, &nz_flag, mode, 0, lg_bs, cb_ox, cb_oy);

    else
    {
        int y = part_idx > 1 ? 2 : 0;
        int x = part_idx & 1 ? 2 : 0;
        cost = venc_analyze_intra_tb(t, cb, &nz_flag, mode, 0, lg_bs, cb_ox+x, cb_oy+y);
    }

    return cost;
}

// Compute the cost of using a specific intra prediction mode.
int64_t venc_analyze_intra_mode_cost(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx, int lg_bs,
                                     int cb_ox, int cb_oy, int mode, uint16_t mode_costs[35], int rmd_flag)
{
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("testing mode %d.\n", mode);
    }
    #endif

    int64_t mode_cost = 0;
    if (t->an.rdm_flag) mode_cost = mode_costs[mode];
    else
    {
        // Optimize this eventually (and extract code).
        t->cbs.rdo.bits = 0;
        cb->intra_luma_mode[part_idx] = mode;
        int val = venc_write_intra_luma_pred(t, cb, part_idx);
        if (val < 3) venc_write_intra_mpm_idx(t, val);
        else venc_write_intra_rem_mode(t, val-3);
        mode_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
    }

    int64_t tt_cost = 0;
    if (rmd_flag) tt_cost = venc_analyze_intra_rough_cost(t, cb, part_idx, lg_bs, cb_ox, cb_oy, mode);

    else
    {
        venc_stash_push(t);
        tt_cost = venc_analyze_intra_luma_tt(t, cb, mode, depth_range, 0, lg_bs, cb_ox, cb_oy);
        venc_stash_pop(t);
    }

    int64_t cost = mode_cost + tt_cost;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("analyzed mode %d cost %d (mode %d, subtree %d).\n", mode, (int)cost, (int)mode_cost, (int)tt_cost);
        #ifdef VAN_TRACE_TB_ANALYSIS
        {
            printf("======> Part result.\n");
            int plane_off = venc_get_cb_block_plane_off(t, cb, 0, cb_ox, cb_oy);
            f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
            print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
        }
        #endif
    }
    #endif

    return cost;
}

// Refine the intra prediction mode to find a better one.
// FIXME. Rough mode decision not supported.
void venc_analyze_intra_mode_refinement(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx, int lg_bs,
                                        int cb_ox, int cb_oy, uint16_t mode_costs[35], int modes[2],
                                        int64_t costs[2], int iterations)
{
    int stash_flag = !!t->stash.flags;

    for (; iterations > 0; iterations--)
    {
        if (stash_flag) venc_stash_reset(t);

        int mode = (modes[0] + modes[1]) >> 1;
        int64_t cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                     cb_ox, cb_oy, mode, mode_costs, 0);

        // The middle cost is better than the current best.
        if (cost < costs[0])
        {
            costs[1] = costs[0];
            modes[1] = modes[0];
            costs[0] = cost;
            modes[0] = mode;
            if (stash_flag) venc_stash_save(t);
        }

        // This middle mode is better than the 2nd best mode.
        else if (cost < costs[1])
        {
            costs[1] = cost;
            modes[1] = mode;
        }

        // The middle mode performs worse. Assume current best is global min.
        else
            break;
    }
}

// Refine the intra prediction mode to find a better one.
// Force a whole-depth binary search.
// Use the best mode.
// FIXME. Rough mode decision not supported.
void venc_analyze_intra_mode_refinement2(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx, int lg_bs,
                                         int cb_ox, int cb_oy, uint16_t mode_costs[35], int modes[2],
                                         int64_t costs[2], int iterations)
{
    int stash_flag = !!t->stash.flags;

    // Split the interval of 7 modes into equal parts. We are going to test
    // three modes, then an additional mode after analyzing the results.
    int quarter_modes[5];
    quarter_modes[0] = modes[0];
    quarter_modes[4] = modes[1];
    quarter_modes[2] = (quarter_modes[0] + quarter_modes[4]) >> 1;
    quarter_modes[1] = (quarter_modes[0] + quarter_modes[2]) >> 1;
    quarter_modes[3] = (quarter_modes[2] + quarter_modes[4]) >> 1;

    // Test the three modes that split the initial interval equally.
    // Keep track of the best mode while doing this.
    int64_t quarter_costs[5] = {costs[0], F265_MAX_SSD, F265_MAX_SSD, F265_MAX_SSD, costs[1]};

    int best_mode_idx = 0;
    int64_t best_cost = quarter_costs[0];
    for (int i = 1; i < 4; i++)
    {
        quarter_costs[i] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                        cb_ox, cb_oy, quarter_modes[i], mode_costs, 0);
        if (quarter_costs[i] < best_cost)
        {
            best_mode_idx = i;
            best_cost = quarter_costs[i];
            if (stash_flag) venc_stash_save(t);
        }

        if (stash_flag) venc_stash_reset(t);
    }

    // Determine if we are checking to the left or the right of the best mode.
    int delta = 1;
    if (best_mode_idx > 0 && quarter_costs[best_mode_idx-1] < quarter_costs[best_mode_idx+1]) delta = -1;

    // Derive the last mode to be tested.
    int last_mode = (quarter_modes[best_mode_idx] + quarter_modes[best_mode_idx+delta]) >> 1;

    // Fetch the cost of that mode.
    int64_t last_cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                     cb_ox, cb_oy, last_mode, mode_costs, 0);

    // Either use the last mode or the best mode prior to the final check.
    if (last_cost < best_cost)
    {
        costs[0] = last_cost;
        modes[0] = last_mode;
        if (stash_flag) venc_stash_save_reset(t);
    }
    else
    {
        costs[0] = best_cost;
        modes[0] = quarter_modes[best_mode_idx];
    }
}

// Opportunistic search, take 2!
// Check "1-away". Stop if current is better.
// Check "4-away".
// Check "3-away" if "4-away" is better than "1-away".
// Check "2-away" if "1-away" is better than "4-away".
// Select the best option.
// FIXME. Rough mode decision not supported.
void venc_analyze_intra_mode_refinement3(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx, int lg_bs,
                                         int cb_ox, int cb_oy, uint16_t mode_costs[35], int modes[2], int64_t costs[2])
{
    int stash_flag = !!t->stash.flags;
    int inc = modes[0] < modes[1] ? 1 : -1;

    int near_modes[5] = { modes[0] };
    for (int i = 1; i < 5; i++) near_modes[i] = near_modes[i-1] + inc;

    int64_t near_costs[5] = { costs[0], F265_MAX_SSD, F265_MAX_SSD, F265_MAX_SSD, F265_MAX_SSD };

    near_costs[1] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                 cb_ox, cb_oy, near_modes[1], mode_costs, 0);
    if (near_costs[0] <= near_costs[1]) return;

    if (stash_flag) venc_stash_save_reset(t);

    near_costs[4] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                 cb_ox, cb_oy, near_modes[4], mode_costs, 0);

    int best_mode_idx = 1;
    if (near_costs[4] < near_costs[1])
    {
        best_mode_idx = 4;
        if (stash_flag) venc_stash_save_reset(t);
        near_costs[3] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                 cb_ox, cb_oy, near_modes[3], mode_costs, 0);
        if (near_costs[3] < near_costs[4])
        {
            best_mode_idx = 3;
            if (stash_flag) venc_stash_save(t);
        }
    }
    else
    {
        if (stash_flag) venc_stash_reset(t);
        near_costs[2] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                 cb_ox, cb_oy, near_modes[2], mode_costs, 0);
        if (near_costs[2] < near_costs[1])
        {
            best_mode_idx = 2;
            if (stash_flag) venc_stash_save(t);
        }
    }

    modes[0] = near_modes[best_mode_idx];
    costs[0] = near_costs[best_mode_idx];
}

// Refine the intra prediction mode to find a better one.
void venc_analyze_intra_breadth_refinement(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx,
                                           int lg_bs, int cb_ox, int cb_oy, uint16_t mode_costs[35],
                                           int *best_mode, int64_t *best_cost, int iterations)
{
    int stash_flag = !!t->stash.flags;
    int algo_flags = t->enc->gd.algo;
    int rough_flag = (algo_flags&4)>>2;

    // Check both modes adjacent to the current best.
    int64_t left_cost = LLONG_MAX, right_cost = LLONG_MAX;

    // Temporary buffers for rough mode decision.
    int temp_modes[6] = {best_mode[0], best_mode[1], best_mode[2], best_mode[0], best_mode[0], -1};
    int64_t temp_costs[6] = {best_cost[0], best_cost[1], best_cost[2], LLONG_MAX, LLONG_MAX, LLONG_MAX};

    // Analyze the angular mode before the best mode.
    if (*best_mode > 2)
    {
        if (rough_flag)
        {
            temp_modes[3] -= 1;
            temp_costs[3] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs, cb_ox, cb_oy,
                                                         temp_modes[3], mode_costs, rough_flag);
        }

        else
        {
            left_cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs, cb_ox, cb_oy,
                                                     *best_mode-1, mode_costs, rough_flag);
            if (stash_flag)
            {
                if (left_cost < *best_cost) venc_stash_save(t);
                venc_stash_reset(t);
            }
        }
    }

    // Analyze the angular mode after the best mode.
    if (*best_mode < 34)
    {
        if (rough_flag)
        {
            temp_modes[4] += 1;
            temp_costs[4] = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs, cb_ox, cb_oy,
                                                         temp_modes[4], mode_costs, rough_flag);
        }

        else
        {
            right_cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs, cb_ox, cb_oy,
                                                      *best_mode+1, mode_costs, rough_flag);
            if (stash_flag)
            {
                if (right_cost < *best_cost && right_cost < left_cost) venc_stash_save(t);
                venc_stash_reset(t);
            }
        }
    }

    // Determine the search direction.
    int inc = -1;
    if (!rough_flag)
    {
        // If current best is a local minimum, stop.
        if (left_cost >= *best_cost && *best_cost <= right_cost) return;

        // Determine which interval to continue exploring.
        if (left_cost < *best_cost && left_cost <= right_cost)
        {
            *best_cost = left_cost;
        }
        else
        {
            *best_cost = right_cost;
            inc = 1;
        }
        *best_mode += inc;
    }

    else
    {
        if (temp_costs[4] < temp_costs[3])
        {
            inc = 1;
            F265_SWAP(int64_t, temp_costs[3], temp_costs[4]);
            F265_SWAP(int, temp_modes[3], temp_modes[4]);
        }

        // Merge sort. Keep the best 3 modes.
        int64_t *in_costs = &temp_costs[0], *ref_costs = &temp_costs[3];
        int *in_modes = &temp_modes[0], *ref_modes = &temp_modes[3];
        for (int n = 0; n < 3; n++)
        {
            if (*in_costs <= *ref_costs)
            {
                best_mode[n] = *in_modes; in_modes++;
                best_cost[n] = *in_costs; in_costs++;
            }
            else
            {
                best_mode[n] = *ref_modes; ref_modes++;
                best_cost[n] = *ref_costs; ref_costs++;
            }
        }

        // If current best is a local minimum, stop.
        if (best_mode[0] == temp_modes[0]) return;
    }

    // Search for a local minimum.
    for (; iterations > 0; iterations--)
    {
        int mode = *best_mode + inc;
        if (mode < 2 || mode > 34) break;
        int64_t cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs, cb_ox, cb_oy,
                                                    mode, mode_costs, rough_flag);
        if (cost >= *best_cost)
        {
            if (rough_flag)
            {
                if (cost < best_cost[1])
                {
                    best_mode[2] = best_mode[1];
                    best_mode[1] = mode;
                    best_cost[2] = best_cost[1];
                    best_cost[1] = cost;
                }
                else if (cost < best_cost[2])
                {
                    best_cost[2] = cost;
                    best_mode[2] = mode;
                }
            }
            break;
        }

        // Prepare the next iteration.
        if (!rough_flag)
        {
            *best_cost = cost;
            *best_mode = mode;
            if (stash_flag) venc_stash_save_reset(t);
        }

        else
        {
            best_cost[2] = best_cost[1];
            best_cost[1] = best_cost[0];
            best_cost[0] = cost;
            best_mode[2] = best_mode[1];
            best_mode[1] = best_mode[0];
            best_mode[0] = mode;
        }
    }
}

// Return the cost of the luma intra partition.
static int64_t venc_analyze_intra_part_luma(f265_enc_thread *t, f265_cb *cb, int depth_range[2], int part_idx,
                                            int lg_bs, int cb_ox, int cb_oy)
{
    int mode_list[35];
    int nb_modes;
    int algo_flags = t->enc->gd.algo;
    int rough_flag = (algo_flags&4)>>2;
    int bsearch_flag = !((algo_flags&2)>>1);
    int refine_flag = (algo_flags&32)>>5;

    // Best rough costs kept for actual RDO.
    int best_rdm_modes[3] = { 0, 0, 0 };
    int64_t best_rdm_costs[3] = { LLONG_MAX, LLONG_MAX, LLONG_MAX };

    // Storage when binary-like search is used.
    int64_t intra_cost[35];

    if (t->an.all_intra_flag)
    {
        for (int i = 0; i < 35; i++) mode_list[i] = i;
        nb_modes = 35;
    }

    else
    {
        // Use 7 "base" modes.
        if (algo_flags&1)
        {
            mode_list[0] = 0;
            mode_list[1] = 1;
            mode_list[2] = 2;
            mode_list[3] = 10;
            mode_list[4] = 18;
            mode_list[5] = 26;
            mode_list[6] = 34;
            nb_modes = 7;
        }

        // Use 5 "base" modes.
        else
        {
            mode_list[0] = 0;
            mode_list[1] = 1;
            mode_list[2] = 2;
            mode_list[3] = 18;
            mode_list[4] = 34;
            nb_modes = 5;
        }
    }

    int best_mode = 0;
    int64_t best_cost = F265_MAX_SSD;
    int stash_flag = !!t->stash.flags;

    if (stash_flag) venc_stash_init_cb_off(t, cb, cb_ox, cb_oy, 1<<lg_bs, 1<<lg_bs);

    // Predict the MPMs and set up the mode cost table.
    int mpm_list[3];
    uint16_t mode_costs[35];
    #if 1
    if (t->an.rdm_flag)
    {
    #endif
    venc_get_intra_pred_mode(t, cb, part_idx, mpm_list);
    for (int i = 0; i < 35; i++) mode_costs[i] = t->an.se_costs[F265_SE_INTRA_LUMA_MODE+3];
    for (int i = 0; i < 3; i++) mode_costs[mpm_list[i]] = t->an.se_costs[F265_SE_INTRA_LUMA_MODE+i];
    #if 1
    }
    #endif

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("analyzing luma part %d, MPM list (%d,%d,%d) est (%d,%d,%d|%d).\n",
               part_idx, mpm_list[0], mpm_list[1], mpm_list[2],
               t->an.se_costs[F265_SE_INTRA_LUMA_MODE+0], t->an.se_costs[F265_SE_INTRA_LUMA_MODE+1],
               t->an.se_costs[F265_SE_INTRA_LUMA_MODE+2], t->an.se_costs[F265_SE_INTRA_LUMA_MODE+3]);
    }
    #endif

    // Test the initial modes.
    // All modes are tested when t->an.all_intra_flag is raised.
    for (int i = 0; i < nb_modes; i++)
    {
        int mode = mode_list[i];
        int64_t cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                    cb_ox, cb_oy, mode, mode_costs, rough_flag);

        if (bsearch_flag) intra_cost[i] = cost;

        if (cost < best_cost)
        {
            // Eventually, this could be optimized by not resetting&saving for
            // the last mode.
            if (stash_flag && !rough_flag) venc_stash_save(t);
            best_mode = mode;
            best_cost = cost;
        }

        // Reset contexts to test the next mode.
        if (stash_flag && !rough_flag) venc_stash_reset(t);

        // Keep track of the N best modes.
        else if (rough_flag)
        {
            if (cost < best_rdm_costs[0])
            {
                best_rdm_costs[2] = best_rdm_costs[1];
                best_rdm_costs[1] = best_rdm_costs[0];
                best_rdm_costs[0] = cost;

                best_rdm_modes[2] = best_rdm_modes[1];
                best_rdm_modes[1] = best_rdm_modes[0];
                best_rdm_modes[0] = mode;
            }

            else if (cost < best_rdm_costs[1])
            {
                best_rdm_costs[2] = best_rdm_costs[1];
                best_rdm_costs[1] = cost;
                best_rdm_modes[2] = best_rdm_modes[1];
                best_rdm_modes[1] = mode;
            }

            else if (cost < best_rdm_costs[2])
            {
                best_rdm_costs[2] = cost;
                best_rdm_modes[2] = mode;
            }
        }

    }

    // The assumption is that only angular modes can be refined.
    if (!t->an.all_intra_flag && refine_flag && bsearch_flag && best_mode > 1)
    {
        int modes[2] = {best_mode, 18};
        int64_t costs[2] = {best_cost, intra_cost[3]};

        // Rules for the 7 "base" modes.
        if (algo_flags&1)
        {
            // These rules are a little more complex since we have to choose one of
            // four intervals.
            if (best_mode == 2)
            {
                modes[1] = 10;
                costs[1] = intra_cost[3];
            }
            else if (best_mode == 10)
            {
                if (intra_cost[2] <= intra_cost[4])
                {
                    modes[1] = 2;
                    costs[1] = intra_cost[2];
                }
                else
                {
                    modes[1] = 18;
                    costs[1] = intra_cost[4];
                }
            }
            else if (best_mode == 18)
            {
                if (intra_cost[3] <= intra_cost[5])
                {
                    modes[1] = 10;
                    costs[1] = intra_cost[3];
                }
                else
                {
                    modes[1] = 26;
                    costs[1] = intra_cost[5];
                }
            }
            else if (best_mode == 26)
            {
                if (intra_cost[4] <= intra_cost[6])
                {
                    modes[1] = 18;
                    costs[1] = intra_cost[4];
                }
                else
                {
                    modes[1] = 34;
                    costs[1] = intra_cost[6];
                }
            }
            else
            {
                modes[1] = 26;
                costs[1] = intra_cost[5];
            }

            // Opportunistic 2: whole-depth binary-search.
            if (algo_flags&8)
            {
                // Force a log2(n) search of the selected interval.
                venc_analyze_intra_mode_refinement2(t, cb, depth_range, part_idx, lg_bs,
                                                    cb_ox, cb_oy, mode_costs, modes, costs, 3);
            }

            // Opportunistic 3: assume uniform mode distribution and bell-shaped
            //                  curve that models costs.
            else if (algo_flags&16)
            {
                // Force a log2(n) search of the selected interval.
                venc_analyze_intra_mode_refinement3(t, cb, depth_range, part_idx, lg_bs,
                                                    cb_ox, cb_oy, mode_costs, modes, costs);
            }

            // Actual binary-like search.
            else
            {
                // Run up to 3 iterations to find a better mode than the current best.
                venc_analyze_intra_mode_refinement(t, cb, depth_range, part_idx, lg_bs,
                                                   cb_ox, cb_oy, mode_costs, modes, costs, 3);
            }
        }

        // Rules for the 5 "base" modes.
        else
        {
            // The rules here only work for the specified starting values. Changing
            // these values requires updating the rules. The rules are simple here
            // because mode 18 *has* to be used. We only have to update the values
            // when mode 18 is the best one to use.
            if (best_mode == 18)
            {
                if (intra_cost[2] <= intra_cost[4])
                {
                    modes[1] = 2;
                    costs[1] = intra_cost[2];
                }
                else
                {
                    modes[1] = 34;
                    costs[1] = intra_cost[4];
                }
            }

            // Actual binary-like search.
            // Run up to 4 iterations to find a better mode than the current best.
            venc_analyze_intra_mode_refinement(t, cb, depth_range, part_idx, lg_bs,
                                               cb_ox, cb_oy, mode_costs, modes, costs, 4);
        }

        // Transfer back the results.
        best_mode = modes[0];
        best_cost = costs[0];
    }

    // The assumption is that only angular modes can be refined.
    if (!t->an.all_intra_flag && refine_flag && !bsearch_flag && best_mode > 1)
    {
        int iterations = algo_flags&1 ? 15 : 7;

        // Use straight up RDO.
        if (!rough_flag)
        {
            venc_analyze_intra_breadth_refinement(t, cb, depth_range, part_idx, lg_bs,
                                                  cb_ox, cb_oy, mode_costs, &best_mode,
                                                  &best_cost, iterations);
        }

        // Run rough mode decision.
        else
        {
            venc_analyze_intra_breadth_refinement(t, cb, depth_range, part_idx, lg_bs,
                                                  cb_ox, cb_oy, mode_costs, best_rdm_modes,
                                                  best_rdm_costs, iterations);
        }
    }

    // Analyze the transform tree (RDO or RDM).
    if (rough_flag)
    {
        best_mode = best_rdm_modes[0];
        best_cost = LLONG_MAX;

        // Just analyze the first mode when delaying transform tree
        // computations.
        int nb_tt_modes = (t->enc->gd.algo & (1<<11)) ? 1 : 3;

        for (int n = 0; n < nb_tt_modes; n++)
        {
            int64_t cost = venc_analyze_intra_mode_cost(t, cb, depth_range, part_idx, lg_bs,
                                                        cb_ox, cb_oy, best_rdm_modes[n], mode_costs, 0);
            if (cost < best_cost)
            {
                if (stash_flag) venc_stash_save(t);
                best_mode = best_rdm_modes[n];
                best_cost = cost;
            }

            if (stash_flag) venc_stash_reset(t);
        }
    }

    if (stash_flag) venc_stash_restore(t);

    // Store the best mode.
    cb->intra_luma_mode[part_idx] = best_mode;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("kept luma mode %d cost %d.\n", best_mode, (int)best_cost);
    }
    #endif

    return best_cost;
}

// Return the cost of encoding the CB with intra prediction.
static int64_t venc_analyze_intra_cb(f265_enc_thread *t, f265_cb *cb)
{
    int lg_bs = cb->lg_bs;
    int split_flag = lg_bs == t->enc->gd.cb_range[0];
    int64_t split_cost = F265_MAX_SSD, unsplit_cost, best_cost;
    int depth_range[2];

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing intra CB.\n");
    }
    #endif

    // Set the stash flags to analyze luma and chroma.
    //
    // For luma, save the transform nodes since we currently need them for
    // chroma. Do not save the chroma reconstruction.
    //
    // For chroma, do not save the transform nodes since they're set with luma.
    // Do not save the luma reconstruction.
    //
    // FIXME: the final encode still needs a valid transform tree when we don't
    // save the transform nodes.
    int saved_stash_flags = t->stash.flags;
    int luma_stash_flags = (saved_stash_flags | 8) & ~6;
    int chroma_stash_flags = saved_stash_flags & ~9;
    int stash_flag = split_flag;

    // Remember the initial transform tree position.
    uint8_t *init_tn = t->tt.tn;

    // Set the CB mode to intra.
    F265_SET_FLAG(cb->flags, F265_CB_INTRA, 1);

    // Add the CB mode cost.
    int64_t cb_mode_cost;
    if (t->an.rdm_flag)
    {
        int skip_idx = F265_SE_CU_SKIP + (venc_get_cu_skip_ctx_off(t, cb)<<1);
        cb_mode_cost = t->an.se_costs[skip_idx+0] + t->an.se_costs[F265_SE_PRED_MODE+1];
    }

    else
    {
        t->cbs.rdo.bits = 0;
        if (t->src_frame->frame_type != F265_FRAME_I) venc_write_cu_skip_flag(t, cb, 0);
        venc_write_pred_mode_flag(t, cb, 1);
        cb_mode_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
    }

    if (stash_flag) venc_stash_init_cb(t, cb);

    // Test HV. Consider putting HV last, there are less modes to save/restore
    // if worse than UN.
    uint8_t hv_luma_mode[4], hv_chroma_mode;
    if (split_flag)
    {
        // Add a split transform node.
        *t->tt.tn++ = 15;

        // Set the partitioning mode to HV and update the prediction map.
        // This needs to be optimized.
        cb->intra_luma_mode[1] = 0;
        venc_update_pmap_unsplit_cb(t, cb);

        // Set the analyzed depth range.
        venc_get_intra_part_depth_range(t, lg_bs-1, depth_range);

        int64_t split_flag_cost;
        if (t->an.rdm_flag) split_flag_cost = t->an.se_costs[F265_SE_INTRA_PART+0];
        else { F265_RDO_COST(split_flag_cost, venc_write_intra_part_mode(t, cb)); }

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzing HV (depth range %d,%d).\n", depth_range[0], depth_range[1]);
        }
        #endif

        venc_stash_push(t);

        t->stash.flags = luma_stash_flags;
        int64_t luma_cost = 0;
        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
            luma_cost += venc_analyze_intra_part_luma(t, cb, depth_range, i, lg_bs-1, (i&1)*sbs, (i>>1)*sbs);

        t->stash.flags = chroma_stash_flags;
        int64_t chroma_cost = 0;
        if (!t->an.rdm_flag || t->an.chroma_me_flag)
            chroma_cost = venc_analyze_intra_cb_chroma(t, cb, init_tn);

        venc_stash_pop(t);

        split_cost = split_flag_cost + luma_cost + chroma_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzed HV cost %d (split_flag %d, luma %d, chroma %d).\n",
                   (int)split_cost, (int)split_flag_cost, (int)luma_cost, (int)chroma_cost);
            #ifdef VAN_TRACE_TB_ANALYSIS
            {
                printf("======> HV result.\n");
                int plane_off = venc_get_cb_block_plane_off(t, cb, 0, 0, 0);
                f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
                print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
            }
            #endif
        }
        #endif

        for (int i = 0; i < 4; i++) hv_luma_mode[i] = cb->intra_luma_mode[i];
        hv_chroma_mode = cb->intra_chroma_mode;
    }

    // Test UN.
    if (stash_flag)
    {
        t->stash.flags = saved_stash_flags;
        venc_stash_save_reset(t);
    }
    cb->intra_luma_mode[1] = -1;
    venc_get_intra_part_depth_range(t, lg_bs, depth_range);

    int64_t unsplit_flag_cost;
    if (t->an.rdm_flag) unsplit_flag_cost = t->an.se_costs[F265_SE_INTRA_PART+1];
    else { F265_RDO_COST(unsplit_flag_cost, venc_write_intra_part_mode(t, cb)); }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing UN (depth range %d,%d).\n", depth_range[0], depth_range[1]);
    }
    #endif

    venc_stash_push(t);

    t->stash.flags = luma_stash_flags;
    int64_t luma_cost = venc_analyze_intra_part_luma(t, cb, depth_range, 0, lg_bs, 0, 0);

    t->stash.flags = chroma_stash_flags;
    int64_t chroma_cost = 0;
    if (!t->an.rdm_flag || t->an.chroma_me_flag)
        chroma_cost = venc_analyze_intra_cb_chroma(t, cb, init_tn);

    venc_stash_pop(t);

    unsplit_cost = unsplit_flag_cost + luma_cost + chroma_cost;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed UN cost %d (split_flag %d, luma %d, chroma %d).\n",
               (int)unsplit_cost, (int)unsplit_flag_cost, (int)luma_cost, (int)chroma_cost);
        #ifdef VAN_TRACE_TB_ANALYSIS
        {
            printf("======> UN result.\n");
            int plane_off = venc_get_cb_block_plane_off(t, cb, 0, 0, 0);
            f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
            print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
        }
        #endif
    }
    #endif

    t->stash.flags = saved_stash_flags;

    if (split_cost < unsplit_cost)
    {
        if (stash_flag) venc_stash_restore(t);
        for (int i = 0; i < 4; i++) cb->intra_luma_mode[i] = hv_luma_mode[i];
        cb->intra_chroma_mode = hv_chroma_mode;
        best_cost = split_cost;
    }

    else
    {
        best_cost = unsplit_cost;
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("finished CB intra best %d (split %d, unsplit %d).\n",
               (int)best_cost, (int)split_cost, (int)unsplit_cost);
        #ifdef VAN_TRACE_TB_ANALYSIS
        {
            printf("======> Intra result.\n");
            int plane_off = venc_get_cb_block_plane_off(t, cb, 0, 0, 0);
            f265_pix *rec_plane = t->src_frame->rec_planes[0] + plane_off;
            print_mb_pix(rec_plane, t->me.ref_stride, 1<<lg_bs);
        }
        #endif
    }
    #endif

    return best_cost + cb_mode_cost;
}

// Do the motion compensation of the current inter mode in the reconstructed
// planes.
static void venc_inter_mode_mc(f265_enc_thread *t, f265_cb *cb)
{
    for (int comp = 0; comp < 3; comp++)
    {
        int plane_off = venc_get_ctb_block_plane_off(t, comp, cb->cb_off[0]>>!!comp, cb->cb_off[1]>>!!comp);
        f265_pix *rec = t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off;
        venc_mc_cb(t, cb, rec, t->me.ref_stride, comp);
    }
}

// Return the cost of the inter transform block specified. Set the non-zero
// flag. If transform_flag is false, the transform is skipped.
static int64_t venc_analyze_inter_tb(f265_enc_thread *t, f265_cb *cb, int *nz_flag, int transform_flag, int comp,
                                     int lg_bs, int cb_ox, int cb_oy)
{
    int csf = !!comp;
    int ct_ox = (cb->cb_off[0]>>csf) + cb_ox;
    int ct_oy = (cb->cb_off[1]>>csf) + cb_oy;
    int bs = 1<<lg_bs;
    int ref_stride = t->me.ref_stride;
    int dist_rec_stride = ref_stride;
    int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox, ct_oy);
    f265_pix *src = t->src_frame->src_planes[comp] + plane_off;
    f265_pix *rec = t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off;
    f265_pix *dist_rec = rec;
    F265_ALIGN64 f265_pix tmp_rec[32*32];

    // Do the transform.
    int64_t coeff_cost = 0;
    if (transform_flag)
    {
        // Compute the reconstruction. The sooner this code is gone, the better.
        venc_set_tmp_tb(t);
        {
            dist_rec = tmp_rec;
            dist_rec_stride = bs;

            int qp = t->qp[!!comp];
            f265_rec_params rp;
            rp.t = t;
            rp.src = src;
            rp.src_stride = ref_stride;
            rp.pred = rec;
            rp.pred_stride = ref_stride;
            rp.rec = dist_rec;
            rp.rec_stride = dist_rec_stride;
            rp.tt = &t->tt;
            rp.cabac_contexts = t->cbs.contexts;
            rp.comp = comp;
            rp.lg_bs = lg_bs;
            rp.tb_depth = 0;
            rp.dst_flag = 0;
            rp.order = 0;
            rp.zero_flag = 0;
            rp.sign_hiding_flag = F265_GET_FLAG(t->enc->gd.eflags, F265_PF_SIGN_HIDING);
            rp.rdoq_flag = 0;
            rp.qp = qp;
            rp.bd = 8;
            rp.iframe_flag = 0;
            rp.dump_flag = 0;
            rp.intra_flag = 0;
            rp.lambda = t->hm_lambda[!!comp];
            rp.compute_ssd_from_idct2 = 0;
            *nz_flag = venc_rec_block(&rp);
        }

        // Compute the encoding cost.
        venc_set_tmp_tb(t);
        if (*nz_flag) { F265_RDO_COST(coeff_cost, venc_write_tb(&t->cbs, &t->tt, !!comp)); }
    }

    else
    {
        // Only for tracing.
        *nz_flag = 0;
    }

    // Compute the reconstruction distortion.
    int64_t dist = venc_an_block_ssd(t, dist_rec, dist_rec_stride, src, ref_stride, bs, comp);

    int64_t cost = coeff_cost + dist;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, comp, cb_ox, cb_oy, lg_bs);
        printf("encoded TB cost %d (coeffs %d, dist %d), nz %d.\n", (int)cost, (int)coeff_cost, (int)dist, *nz_flag);
    }
    #endif

    return cost;
}

// Return the cost of the inter transform block specified. Set the non-zero
// flag. This function tries both values of the CBF when requested.
static int64_t venc_analyze_inter_tb_full(f265_enc_thread *t, f265_cb *cb, int *nz_flag, int comp, int depth,
                                          int lg_bs, int cb_ox, int cb_oy)
{
    if (!t->an.nullify_inter_tb_flag)
        return venc_analyze_inter_tb(t, cb, nz_flag, 1, comp, lg_bs, cb_ox, cb_oy);

    int stash_flag = !!t->stash.flags;
    venc_stash_push(t);

    // The function venc_stash_init_cb_off() expects luma offsets. In the case
    // of chroma, convert the chroma offsets into luma offsets.
    if (stash_flag) venc_stash_init_cb_off(t, cb, cb_ox<<!!comp, cb_oy<<!!comp, 1<<(lg_bs+!!comp), 1<<(lg_bs+!!comp));

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, comp, cb_ox, cb_oy, lg_bs);
        printf("analyzing non-zero TB comp %d.\n", comp);
    }
    #endif

    // Try the non-zero case.
    int64_t nz_cost = venc_analyze_inter_tb(t, cb, nz_flag, 1, comp, lg_bs, cb_ox, cb_oy);
    int64_t ret_cost = nz_cost;

    // Try the zero case.
    if (*nz_flag)
    {
        if (stash_flag) venc_stash_save_reset(t);

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, comp, cb_ox, cb_oy, lg_bs);
            printf("analyzing zero TB comp %d.\n", comp);
        }
        #endif

        int64_t z_cost = venc_analyze_inter_tb(t, cb, nz_flag, 0, comp, lg_bs, cb_ox, cb_oy);

        // Add the CBF cost. This is an approximation.
        int ctx_idx = comp ? F265_CO_CBF_CHROMA+depth : F265_CO_CBF_LUMA+!depth;
        int64_t nz_cbf_cost = venc_get_context_bin_cost(t, ctx_idx, 1);
        int64_t z_cbf_cost =  venc_get_context_bin_cost(t, ctx_idx, 0);
        int64_t est_nz_cost = nz_cost + nz_cbf_cost;
        int64_t est_z_cost = z_cost + z_cbf_cost;

        *nz_flag = est_nz_cost < est_z_cost;
        if (*nz_flag)
        {
            if (stash_flag) venc_stash_restore(t);
        }

        else ret_cost = z_cost;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, comp, cb_ox, cb_oy, lg_bs);
            printf("analyzed TB comp %d cost %d (nz_cbf %d + nz_tb %d = %d, z_cbf %d + z_tb %d = %d).\n",
                   comp, (int)ret_cost, (int)nz_cbf_cost, (int)nz_cost, (int)est_nz_cost, (int)z_cbf_cost,
                   (int)z_cost, (int)est_z_cost);
        }
        #endif
    }

    venc_stash_pop(t);

    return ret_cost;
}

// Return the cost of the inter transform tree. This function currently stores
// the transform nodes unconditionally.
static int64_t venc_analyze_inter_sub_tt(f265_enc_thread *t, f265_cb *cb, int *yuv_flags, int depth_range[2],
                                         int depth, int block_idx, int lg_bs, int cb_ox, int cb_oy)
{
    int split_flag = depth < depth_range[1];
    int unsplit_flag = depth >= depth_range[0];
    int split_present_flag = split_flag && unsplit_flag;
    int stash_flag = t->stash.flags && split_present_flag;
    int64_t split_cost = F265_MAX_SSD, unsplit_cost = F265_MAX_SSD, best_cost;
    int split_yuv_flags = 0, unsplit_yuv_flags = 0;

    if (stash_flag) venc_stash_init_cb_off(t, cb, cb_ox, cb_oy, 1<<lg_bs, 1<<lg_bs);

    if (split_flag)
    {
        // Add a split transform node.
        uint8_t *tn = t->tt.tn++;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzing split TT depth %d.\n", depth);
        }
        #endif

        // Encode the split flag.
        int64_t split_flag_cost = 0;
        if (split_present_flag)
        {
            F265_RDO_COST(split_flag_cost,
                          venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 1));
        }

        // Analyze the subtree.
        venc_stash_push(t);

        // Packed chroma subtree non-zero flags, in order (VU VU VU VU).
        int packed_chroma_nz = 0;

        int64_t subtree_cost = 0;
        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
        {
            int nz_flags;
            subtree_cost += venc_analyze_inter_sub_tt(t, cb, &nz_flags, depth_range, depth+1, i, lg_bs-1,
                                                      cb_ox + (i&1)*sbs, cb_oy + (i>>1)*sbs);
            split_yuv_flags |= nz_flags&1;
            packed_chroma_nz |= (nz_flags>>1)<<(i<<1);
        }

        venc_stash_pop(t);

        // Encode the chroma CBFs of the subtree.
        int64_t chroma_cbf_cost = 0;
        if (lg_bs != 3)
        {
            if (packed_chroma_nz)
            {
                // Determine the value of the chroma non-zero flags.
                int comp_nz[2];
                for (int i = 0; i < 2; i++)
                {
                    comp_nz[i] = !!((packed_chroma_nz>>i)&0x55);
                    split_yuv_flags |= comp_nz[i]<<(1+i);
                }

                // Encode the CBFs in order.
                t->cbs.rdo.bits = 0;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 2; j++)
                    {
                        // Skip the component when it's all zero.
                        if (!comp_nz[j]) continue;

                        int cbf = (packed_chroma_nz>>((i<<1)+j))&1;
                        venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+depth+1, cbf);
                    }
                chroma_cbf_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
            }
        }

        // Get the chroma flags from the last block.
        else
        {
            split_yuv_flags |= (packed_chroma_nz>>6)<<1;
        }

        split_cost = split_flag_cost + chroma_cbf_cost + subtree_cost;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 1, cb_ox, cb_oy, lg_bs);
            printf("analyzed split TT depth %d cost %d (split_flag %d, chroma_cbf %d, subtree %d).\n",
                   depth, (int)split_cost, (int)split_flag_cost, (int)chroma_cbf_cost, (int)subtree_cost);
        }
        #endif

        *tn = 8|split_yuv_flags;
    }

    if (unsplit_flag)
    {
        if (stash_flag) venc_stash_save_reset(t);

        // Add an unsplit transform node.
        uint8_t *tn = t->tt.tn++;

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzing unsplit TT depth %d.\n", depth);
        }
        #endif

        int64_t unsplit_flag_cost = 0;
        if (split_present_flag)
        {
            F265_RDO_COST(unsplit_flag_cost,
                          venc_encode_context_bin(&t->cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 0));
        }

        // Encode the luma transform block.
        int nz_flag;
        int64_t tb_cost = venc_analyze_inter_tb_full(t, cb, &nz_flag, 0, depth, lg_bs, cb_ox, cb_oy);
        unsplit_yuv_flags |= nz_flag;

        // Encode the chroma transform blocks.
        // FIXME: we should consider the tentative chroma flags for the
        // split/unsplit comparison (without passing up the tentative costs).
        if (lg_bs > 2 || block_idx == 3)
        {
            int chroma_lg_bs = F265_MAX(lg_bs-1, 2);
            int chroma_off[2] = { (cb_ox>>1)&~3, (cb_oy>>1)&~3 };

            for (int comp = 1; comp < 3; comp++)
            {
                tb_cost += venc_analyze_inter_tb_full(t, cb, &nz_flag, comp, depth, chroma_lg_bs,
                                                      chroma_off[0], chroma_off[1]);
                unsplit_yuv_flags |= nz_flag<<comp;
            }
        }

        // Encode the luma CBF unless it is inferred by root_cbf.
        int64_t luma_cbf_cost = 0;
        if (depth || unsplit_yuv_flags&6)
        {
            F265_RDO_COST(luma_cbf_cost,
                          venc_encode_context_bin(&t->cbs, F265_CO_CBF_LUMA+!depth, unsplit_yuv_flags&1));
        }

        unsplit_cost = unsplit_flag_cost + luma_cbf_cost + tb_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
            printf("analyzed unsplit TT depth %d cost %d (split_flag %d, luma_cbf %d, TBs %d).\n",
                   depth, (int)unsplit_cost, (int)unsplit_flag_cost, (int)luma_cbf_cost, (int)tb_cost);
        }
        #endif
        *tn = unsplit_yuv_flags;
    }

    if (split_cost < unsplit_cost)
    {
        if (stash_flag) venc_stash_restore(t);
        best_cost = split_cost;
        *yuv_flags = split_yuv_flags;
    }

    else
    {
        best_cost = unsplit_cost;
        *yuv_flags = unsplit_yuv_flags;
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_tt_loc(t, cb, 0, cb_ox, cb_oy, lg_bs);
        printf("analyzed TT depth %d best %d (split %d, unsplit %d).\n",
               depth, (int)best_cost, (int)split_cost, (int)unsplit_cost);
    }
    #endif

    return best_cost;
}

// Return the cost of the inter transform tree with a residual (including
// root_cbf if needed). The cost is maximal if it turns out there is no
// residual.
static int64_t venc_analyze_inter_tt_residual(f265_enc_thread *t, f265_cb *cb, int root_cbf_present_flag)
{
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing residual cost.\n");
    }
    #endif

    // Get the depth range.
    int depth_range[2];
    venc_get_inter_part_depth_range(t, cb->inter_part, cb->lg_bs, depth_range);

    // Analyze the subtree.
    int nz_flags;
    int64_t subtree_cost = venc_analyze_inter_sub_tt(t, cb, &nz_flags, depth_range, 0, 0, cb->lg_bs, 0, 0);

    // When no residual is present, complications arise.
    //
    // If we haven't split the transform tree, then we encoded something like
    // root_cbf=1, cbf_cb=0, cbf_cr=0. This is not a valid encoding.
    //
    // If we have split the transform tree, then we encoded something like
    // root_cbf=1, cbf_cb=0, cbf_cr=0, cbf_luma=0 x4. This is a valid encoding
    // as far as the specification is concerned. It could even be RD-optimal
    // compared to encoding root_cbf=0, although it seems fairly unlikely. The
    // final encoding would currently convert this case to root_cbf=0.
    // Preventing this is tricky, and undesired if we run RDOQ in the final
    // encoding but not in the analysis.
    //
    // The whole situation is a can of worms. We keep it simple and fallback to
    // root_cbf=0 when all YUV components are zero.
    if (!nz_flags)
    {
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("residual is zero, discarding residual encoding.\n");
        }
        #endif

        return F265_MAX_SSD;
    }

    // Encode root_cbf.
    int64_t root_cbf_cost = 0;
    if (root_cbf_present_flag)
    {
        F265_RDO_COST(root_cbf_cost, venc_encode_context_bin(&t->cbs, F265_CO_ROOT_CBF, 1));
    }

    // Encode the two chroma CBFs.
    t->cbs.rdo.bits = 0;
    for (int i = 0; i < 2; i++)
        venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+0, (nz_flags>>(1+i))&1);
    int64_t chroma_cbf_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);

    int64_t residual_cost = root_cbf_cost + chroma_cbf_cost + subtree_cost;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed residual cost %d (root_cbf %d, chroma CBF %d, subtree %d).\n",
               (int)residual_cost, (int)root_cbf_cost, (int)chroma_cbf_cost, (int)subtree_cost);
    }
    #endif

    return residual_cost;
}

// Return the cost of the inter transform tree without a residual (excluding
// root_cbf).
static int64_t venc_analyze_inter_tt_no_residual(f265_enc_thread *t, f265_cb *cb)
{
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing no-residual cost.\n");
    }
    #endif

    // Add an unsplit node.
    if (cb->lg_bs == 6)
    {
        *t->tt.tn++ = 8;
        for (int i = 0; i < 4; i++) *t->tt.tn++ = 0;
    }

    else *t->tt.tn++ = 0;

    // Compute the distortion.
    int64_t dist_cost = 0;
    for (int comp = 0, dummy; comp < 3; comp++)
        dist_cost += venc_analyze_inter_tb(t, cb, &dummy, 0, comp, cb->lg_bs - !!comp, 0, 0);

    return dist_cost;
}

// Return the cost of the inter transform tree for the current inter mode.
static int64_t venc_analyze_inter_tt(f265_enc_thread *t, f265_cb *cb)
{
    // Using a variable just in case we change this code.
    int stash_flag = 1;
    int64_t residual_cost, no_residual_cost, best_cost;

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing inter TT.\n");
    }
    #endif

    if (stash_flag) venc_stash_init_cb(t, cb);

    // Encode with a residual. Consider swapping the cases to avoid a restore
    // statistically.
    {
        venc_stash_push(t);
        residual_cost = venc_analyze_inter_tt_residual(t, cb, 1);
        venc_stash_pop(t);

        if (stash_flag)
        {
            if (residual_cost != F265_MAX_SSD) venc_stash_save_reset(t);
            else venc_stash_reset(t);
        }
    }

    // Encode without a residual.
    {
        // Get the distortion.
        int64_t dist_cost = venc_analyze_inter_tt_no_residual(t, cb);

        // Encode root_cbf.
        F265_RDO_COST(int64_t root_cbf_cost, venc_encode_context_bin(&t->cbs, F265_CO_ROOT_CBF, 0));

        no_residual_cost = root_cbf_cost + dist_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzed no-residual cost %d (root_cbf %d, dist %d).\n",
                   (int)no_residual_cost, (int)root_cbf_cost, (int)dist_cost);
        }
        #endif
    }

    // Residual case is better.
    if (residual_cost < no_residual_cost)
    {
        if (stash_flag) venc_stash_restore(t);
        best_cost = residual_cost;
    }

    // No-residual case is better.
    else
    {
        best_cost = no_residual_cost;
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed inter TT cost %d (residual %d, no-residual %d).\n",
               (int)best_cost, (int)residual_cost, (int)no_residual_cost);
    }
    #endif

    return best_cost;
}

// Return the cost of encoding the partition with unidirectional prediction
// using the reference index specified.
static int venc_analyze_inter_uni(f265_enc_thread *t, f265_cb *cb, int part_idx, int list, int ref_idx)
{
    f265_me_ctx *me = &t->me;
    f265_inter_block *ib = &t->ib;
    f265_mv *pmv = ib->part_mv;
    f265_ref_ctx *ref_ctx = t->ref_ctx[list] + ref_idx;

    // Set the reference data.
    venc_me_set_ref(t, ref_ctx, 0);

    // Predict the PMVs and choose the best using the SAD metric.
    venc_get_pmv(t, part_idx, ib->neighbour_mv, ref_idx, list, pmv);
    int pmv_cost;
    int init_pmv_idx = venc_me_test_pmv(pmv, me, 0, &pmv_cost);
    venc_me_set_pmv(me, pmv[init_pmv_idx], t->qp[0], 0);

    // Set the best motion vector and its cost.
    me->ref[0].best_mv = pmv[init_pmv_idx];
    me->best_cost = pmv_cost + venc_me_mv_cost_test(me, me->ref[0].best_mv, 0);

    // Regular motion estimation.
    if (!t->an.hm_me_flag)
    {
        // Do the search. WARNING: works currently because the list won't be
        // modified by venc_early_me().
        venc_early_me(me, pmv + init_pmv_idx, 1, NULL);
    }

    // HM-compatible motion estimation.
    else
    {
        // Back up the motion estimation context.
        f265_me_ctx saved_me = *me;

        // Set up ME.
        f265_frac_ref_block rb;
        f265_weight wp;
        f265_ref_ctx hm_ctx[2];
        hm_ctx[0] = *ref_ctx;

        rb.bipred_flag = 0;
        for (int i = 0; i < 2; i++)
        {
            rb.ctx[i] = hm_ctx + i;
            rb.pos[i] = me->pos[i];
            rb.size[i] = me->dim[i];
        }
        rb.comp = 0;

        venc_hm_setup_me(t, cb, me->dim[1], 0, list, ref_idx, &wp);
        me->dist_func_id = t->an.mode_metric;

        // Perform ME.
        int hm_mv_bits, hm_bits = 0, hm_cost;
        venc_hm_me_search(t, cb, &rb, 0, &me->ref[0].pmv, &hm_mv_bits, &hm_bits, &hm_cost);
        f265_mv best_mv = me->ref[0].hm_best_mv;

        // Restore the motion estimation context.
        *me = saved_me;
        me->ref[0].best_mv = best_mv;
        me->best_cost = hm_cost;
    }

    // Add the best PMV cost.
    // Eventually consider refining if we change the PMV.
    ib->part_pmv_idx[0] = venc_me_find_nearest_pmv(pmv, t->an.se_costs + F265_SE_MVP, me);

    return me->best_cost;
}

// Stash the data of the current inter mode if it is better than the previous
// inter mode and rollback to the initial state.
//
// Eventually use a single function to save best unsplit CB mode (intra, PCM,
// bypass, inter) conditionally. We're saving too much.
static void venc_keep_best_inter_mode(f265_enc_thread *t, f265_cb *cb, int64_t cost, int inter_mode)
{
    f265_inter_block *ib = &t->ib;
    f265_stash *s = &t->stash;
    int stash_flag = !!s->flags;

    if (cost < ib->cb_best_cost)
    {
        ib->cb_best_cost = cost;
        ib->cb_inter_mode = inter_mode;
        memcpy(ib->cb_mv, cb->mv, sizeof(cb->mv));
        memcpy(ib->cb_ref_idx, cb->ref_idx, sizeof(cb->ref_idx));
        memcpy(ib->cb_inter_bf, cb->inter_bf, sizeof(cb->inter_bf));
        if (stash_flag) venc_stash_save(t);
    }

    if (stash_flag) venc_stash_reset(t);
}

// Restore the best coding block inter mode.
static void venc_restore_best_inter_mode(f265_enc_thread *t, f265_cb *cb)
{
    f265_inter_block *ib = &t->ib;

    memcpy(cb->mv, ib->cb_mv, sizeof(cb->mv));
    memcpy(cb->ref_idx, ib->cb_ref_idx, sizeof(cb->ref_idx));
    memcpy(cb->inter_bf, ib->cb_inter_bf, sizeof(cb->inter_bf));

    if (ib->cb_inter_mode == -1)
    {
        cb->inter_part = F265_PART_UN;
        F265_SET_FLAG(cb->flags, F265_CB_SKIP, 1);
    }

    else
    {
        cb->inter_part = ib->cb_inter_mode;
        F265_SET_FLAG(cb->flags, F265_CB_SKIP, 0);
    }
}

// Analyze the UN skip and merge modes.
static void venc_analyze_inter_skip(f265_enc_thread *t, f265_cb *cb, int64_t parent_base_costs[2])
{
    int stash_flag = !!t->stash.flags;

    // Base syntax element costs. Index 0 corresponds to the skip case, index 1
    // corresponds to the merge case (root_cbf inferred to be 1).
    //   base_costs[0]: cu_skip_flag=1
    //   base_costs[1]: cu_skip_flag=0, pred_mode_flag=0, part_mode=1,
    //                  merge_flag=1.
    int64_t base_costs[2];

    // CABAC contexts to restore if merge is selected (part_mode, merge_flag).
    uint8_t merge_cabac_ctx[2];

    // Indices of the contexts above.
    uint8_t merge_cabac_ctx_idx[2] = { F265_CO_PART_MODE, F265_CO_MERGE_FLAG };

    // Set the skip cost.
    base_costs[0] = parent_base_costs[0];

    // Exclude the merge flag in RDM mode, this is handled below.
    if (t->an.rdm_flag)
    {
        base_costs[1] = parent_base_costs[1] + t->an.se_costs[F265_SE_INTER_PART+0];
    }

    // Encode both cases and assume skip will be selected.
    else
    {
        t->cbs.rdo.bits = 0;

        // Save, encode, cache and restore each context.
        for (int i = 0; i < 2; i++)
        {
            int ctx_idx = merge_cabac_ctx_idx[i];
            int orig_ctx = t->cbs.contexts[ctx_idx];
            venc_encode_context_bin(&t->cbs, ctx_idx, 1);
            merge_cabac_ctx[i] = t->cbs.contexts[ctx_idx];
            t->cbs.contexts[ctx_idx] = orig_ctx;
        }

        base_costs[1] = parent_base_costs[1] + venc_rdo_bit_cost(t, t->cbs.rdo.bits);
    }

    // Set the partitioning mode to UN.
    cb->inter_part = F265_PART_UN;

    // Set the motion estimation data. We don't need to set up the partition for
    // RDO since we do the motion compensation.
    if (t->an.rdm_flag)
    {
        venc_me_set_partition(t, cb, 0);
        t->me.dist_func_id = t->an.mode_metric;
    }

    else
        venc_get_neighbour_mvs(t, cb, 0, t->ib.neighbour_mv);

    // Get the merge candidates. The memset() is required for memcmp() (MVs
    // aren't zeroed currently).
    f265_inter_neighbour_mv cands[5];
    memset(cands, 0, sizeof(cands));
    venc_get_merge_candidate(t, cb, 0, t->ib.neighbour_mv, cands);

    // Test each merge candidate.
    for (int merge_idx = 0; merge_idx < t->enc->gd.merge_cand; merge_idx++)
    {
        f265_inter_neighbour_mv *cand = cands + merge_idx;

        // Skip the merge candidate if it's the same as the previous one (the
        // merge index has a negligible effect on the RD cost). Consider using a
        // faster algorithm here. We could also remember the previous cost and
        // just update the cost difference as needed.
        if (merge_idx && !memcmp(cands + merge_idx - 1, cand, sizeof(f265_inter_neighbour_mv)))
            continue;

        // Set the merge candidate data.
        cb->inter_bf[0] = merge_idx;
        for (int i = 0; i < 2; i++)
        {
            cb->ref_idx[0][i] = cand->ref_idx[i];
            cb->mv[0][i] = cand->mv[i];
        }

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("Testing UN merge candidate %d, ref %d (%d,%d), ref %d (%d,%d).\n",
                   merge_idx,
                   cb->ref_idx[0][0], cb->mv[0][0].x, cb->mv[0][0].y,
                   cb->ref_idx[0][1], cb->mv[0][1].x, cb->mv[0][1].y);
        }
        #endif

        // In RDM, do not test for skip explicitly. The final encode will
        // transform the merge case into the skip case as needed.
        if (t->an.rdm_flag)
        {
            int merge_idx_cost = t->an.se_costs[F265_SE_MERGE_IDX+merge_idx];
            int dist = venc_me_merge_cand_dist(t, *cand);
            int merge_cost = base_costs[1] + merge_idx_cost + dist;

            #ifdef VAN_TRACE_ANALYSIS
            if (venc_trace_analysis_flag)
            {
                venc_an_cb_loc(t, cb, 0);
                printf("analyzed merge cost %d (base %d, merge_flag_and_idx %d, dist %d).\n",
                       (int)merge_cost, (int)base_costs[1], (int)merge_idx_cost, (int)dist);
            }
            #endif
            venc_keep_best_inter_mode(t, cb, merge_cost, F265_PART_UN);
        }

        else
        {
            // Do the motion compensation.
            venc_inter_mode_mc(t, cb);

            // Encode the merge index and cache its context.
            F265_RDO_COST(int64_t merge_idx_cost, venc_write_merge_idx(t, cb, merge_idx));
            int merge_idx_ctx = t->cbs.contexts[F265_CO_MERGE_IDX];

            // Test merge.
            {
                // Encode the transform tree.
                venc_stash_push(t);
                int64_t tt_cost = venc_analyze_inter_tt_residual(t, cb, 0);
                venc_stash_pop(t);

                // There is a residual.
                if (tt_cost != F265_MAX_SSD)
                {
                    int64_t merge_cost = base_costs[1] + merge_idx_cost + tt_cost;
                    #ifdef VAN_TRACE_ANALYSIS
                    if (venc_trace_analysis_flag)
                    {
                        venc_an_cb_loc(t, cb, 0);
                        printf("analyzed merge cost %d (base %d, merge_idx %d, TT %d).\n",
                               (int)merge_cost, (int)base_costs[1], (int)merge_idx_cost, (int)tt_cost);
                    }
                    #endif

                    // Replace the merge contexts. Make this conditional
                    // eventually (based on best cost). The merge index context
                    // is already up-to-date.
                    for (int i = 0; i < 2; i++) t->cbs.contexts[merge_cabac_ctx_idx[i]] = merge_cabac_ctx[i];

                    venc_keep_best_inter_mode(t, cb, merge_cost, F265_PART_UN);
                }

                // There is no residual.
                else
                    if (stash_flag) venc_stash_reset(t);
            }

            // Test skip.
            {
                int64_t dist_cost = venc_analyze_inter_tt_no_residual(t, cb);
                int64_t skip_cost = base_costs[0] + merge_idx_cost + dist_cost;
                #ifdef VAN_TRACE_ANALYSIS
                if (venc_trace_analysis_flag)
                {
                    venc_an_cb_loc(t, cb, 0);
                    printf("analyzed skip cost %d (base %d, merge_idx %d, dist %d).\n",
                           (int)skip_cost, (int)base_costs[0], (int)merge_idx_cost, (int)dist_cost);
                }
                #endif

                // Update the merge index context.
                t->cbs.contexts[F265_CO_MERGE_IDX] = merge_idx_ctx;

                venc_keep_best_inter_mode(t, cb, skip_cost, -1);
            }
        }
    }
}

// Return the cost of the inter partition specified.
//
// During the partition analysis, we approximate the cost of the transform tree
// using the prediction distortion. Theoretically it might be possible to gain
// some quality by computing the reconstruction distortion on a virtual
// non-square transform tree centered on the partition, but the cost would be
// prohibitive and the implementation awkward. It seems more productive to
// refine the motion vectors in a second pass after each partition has been
// analyzed.
//
// Since the transform tree cost is imprecise, there is no point in using
// extremely precise costs for the other syntax elements. Hence, we approximate
// all costs and we do not touch the CABAC contexts during the partition
// analysis.
//
// This function is not split into smaller functions so that the compiler can
// optimize the use of local variables.
static int venc_analyze_inter_part(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    f265_me_ctx *me = &t->me;
    f265_inter_block *ib = &t->ib;

    // Set the partition data.
    venc_me_set_partition(t, cb, part_idx);

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing part %d.\n", part_idx);
    }
    #endif

    // Find the best motion vector in each reference list. Store the data in the
    // CB and cache the PMVs in the IB.
    int uni_cost[2] = { F265_MAX_SAD, F265_MAX_SAD };
    for (int list = 0; list < t->nb_lists; list++)
    {
        for (int ref_idx = 0; ref_idx < t->nb_ref_idx[list]; ref_idx++)
        {
            int cost = venc_analyze_inter_uni(t, cb, part_idx, list, ref_idx);
            cost += t->an.se_costs[F265_SE_REF_IDX + list*16 + ref_idx];

            if (cost < uni_cost[list])
            {
                uni_cost[list] = cost;
                cb->mv[part_idx][list] = me->ref[0].best_mv;
                cb->ref_idx[part_idx][list] = ref_idx;
                for (int i = 0; i < 2; i++) ib->uni_pmv[list][i] = ib->part_mv[i];
                ib->uni_pmv_idx[list] = ib->part_pmv_idx[0];
            }

            #ifdef VAN_TRACE_ANALYSIS
            if (venc_trace_analysis_flag)
            {
                venc_an_cb_loc(t, cb, 0);
                printf("analyzed part %d list %d ref_idx %d cost %d MV (%d,%d).\n",
                       part_idx, list, ref_idx, cost, me->ref[0].best_mv.x, me->ref[0].best_mv.y);
            }
            #endif
        }
    }

    int best_uni_list = uni_cost[1] < uni_cost[0];
    int best_uni_cost = uni_cost[best_uni_list];
    int store_uni_flag = 1;
    int part_cost = best_uni_cost;

    // Test bidirectional prediction.
    // FIXME. Derive best_uni_list inside this branch (without the 4x8 test), by
    // adding the second bin context value of inter_pred_idc. If not 4x8, add
    // the cost of the first bin context. Keep the best cost overall. The P case
    // must fast track.
    if (t->nb_lists == 2 && (cb->lg_bs > 3 || cb->inter_part == F265_PART_UN))
    {
        // Use the mode comparison metric by default.
        me->dist_func_id = t->an.mode_metric;

        // Assuming the reference indices and motion vectors stay the same as
        // the unidirectional prediction.
        for (int list = 0; list < 2; list++)
        {
            ib->part_mv[list] = cb->mv[part_idx][list];
            ib->part_pmv_idx[list] = ib->uni_pmv_idx[list];
        }

        // Set the references and PMVs.
        for (int list = 0; list < 2; list++)
        {
            venc_me_set_ref(t, t->ref_ctx[list] + cb->ref_idx[part_idx][list], list);
            venc_me_set_pmv(me, ib->uni_pmv[list][ib->uni_pmv_idx[list]], t->qp[0], list);
        }

        // Get the cost.
        int bi_cost = venc_me_mv_total_cost_bi(ib->part_mv, me);
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzed part %d bi cost %d.\n", part_idx, bi_cost);
        }
        #endif

        // Do the cost comparison here because the compiler doesn't know that
        // the unidirectional cost is lower when biprediction isn't tested.
        if (bi_cost < best_uni_cost)
        {
            store_uni_flag = 0;
            part_cost = bi_cost;
        }
    }

    // Keep the data of the best prediction mode so far.
    if (store_uni_flag)
    {
        int pmv_idx = ib->uni_pmv_idx[best_uni_list];
        ib->mode_pmv[part_idx][best_uni_list] = ib->uni_pmv[best_uni_list][pmv_idx];
        cb->ref_idx[part_idx][!best_uni_list] = -1;
        cb->inter_bf[part_idx] = (pmv_idx<<(3+best_uni_list)) | 5;
    }

    else
    {
        for (int list = 0; list < 2; list++)
        {
            cb->mv[part_idx][list] = ib->part_mv[list];
            ib->mode_pmv[part_idx][list] = ib->uni_pmv[list][ib->part_pmv_idx[list]];
        }
        cb->inter_bf[part_idx] = (ib->part_pmv_idx[1]<<4) | (ib->part_pmv_idx[0]<<3) | 5;
    }

    // Add the non-merge flag cost.
    part_cost += t->an.se_costs[F265_SE_MERGE_IDX+5];

    // Test the merge candidates unless the partitioning mode is UN.
    if (cb->inter_part != F265_PART_UN)
    {
        // Use the mode comparison metric.
        me->dist_func_id = t->an.mode_metric;

        // Get the merge candidates.
        f265_inter_neighbour_mv cands[5];
        memset(cands, 0, sizeof(cands));
        venc_get_merge_candidate(t, cb, part_idx, t->ib.neighbour_mv, cands);

        for (int merge_idx = 0; merge_idx < t->enc->gd.merge_cand; merge_idx++)
        {
            f265_inter_neighbour_mv *cand = cands + merge_idx;

            // Skip the merge candidate if it's the same as the previous one.
            if (merge_idx && !memcmp(cands + merge_idx - 1, cand, sizeof(f265_inter_neighbour_mv)))
                continue;

            // Compute the cost.
            int dist = venc_me_merge_cand_dist(t, *cand);
            int cost = t->an.se_costs[F265_SE_MERGE_IDX+merge_idx] + dist;
            #ifdef VAN_TRACE_ANALYSIS
            if (venc_trace_analysis_flag)
            {
                venc_an_cb_loc(t, cb, 0);
                printf("analyzed merge index %d cost %d.\n", merge_idx, cost);
            }
            #endif

            // Keep the merge candidate.
            if (cost < part_cost)
            {
                part_cost = cost;

                cb->inter_bf[part_idx] = merge_idx;
                for (int i = 0; i < 2; i++)
                {
                    cb->ref_idx[part_idx][i] = cand->ref_idx[i];
                    cb->mv[part_idx][i] = cand->mv[i];
                }
            }
        }
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed part %d, estimated cost %d.\n", part_idx, part_cost);
    }
    #endif

    return part_cost;
}

// Return the RDO cost of the partitioning mode specified.
static int64_t venc_rdo_inter_part_mode(f265_enc_thread *t, f265_cb *cb, int part_mode, int64_t base_cost)
{
    // Do the motion compensation.
    venc_inter_mode_mc(t, cb);

    // Analyze the transform tree.
    venc_stash_push(t);
    int64_t tt_cost = venc_analyze_inter_tt(t, cb);
    venc_stash_pop(t);

    // Compute the prediction block cost.
    t->cbs.rdo.bits = 0;
    venc_write_inter_part_mode(t, cb);
    int nb_parts = f265_nb_parts[part_mode];
    for (int i = 0; i < nb_parts; i++) venc_write_inter_part(t, cb, i);
    int64_t pb_cost = venc_rdo_bit_cost(t, t->cbs.rdo.bits);

    int64_t total_cost = base_cost + pb_cost + tt_cost;
    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed inter mode %d cost %d (base %d, PB %d, TT %d).\n",
               part_mode, (int)total_cost, (int)base_cost, (int)pb_cost, (int)tt_cost);
    }
    #endif

    return total_cost;
}

// Refine the UN mode by testing different motion vectors.
static void venc_refine_inter_un_mode(f265_enc_thread *t, f265_cb *cb, int64_t base_cost)
{
    f265_inter_block *ib = &t->ib;

    // Set the partitioning mode and update the prediction map. This needs to be
    // optimized.
    cb->inter_part = F265_PART_UN;
    venc_update_pmap_unsplit_cb(t, cb);

    // Set the partition data.
    venc_me_set_partition(t, cb, 0);

    // Keep the original inter data.
    int uni_list = cb->ref_idx[0][0] == -1;
    f265_mv orig_mv = ib->cb_mv[0][uni_list];

    // Restore the rest of the UN data.
    for (int list = 0; list < 2; list++)
    {
        cb->ref_idx[0][list] = ib->cb_ref_idx[0][list];
        cb->inter_bf[0] = ib->cb_inter_bf[0];
    }

    // Test the diamond positions.
    const int8_t dia_pos[4][2] = { {0, -1}, {0, 1}, { -1, 0}, {1, 0} };
    for (int pos = 0; pos < 4; pos++)
    {
        const int8_t *mv_off = dia_pos[pos];

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("Refining UN mode with MV offset (%d,%d).\n", mv_off[0], mv_off[1]);
        }
        #endif

        // Offset the motion vector.
        cb->mv[0][uni_list] = orig_mv;
        for (int i = 0; i < 2; i++) cb->mv[0][uni_list].v[i] += mv_off[i];

        int64_t total_cost = venc_rdo_inter_part_mode(t, cb, 0, base_cost);
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("Refined cost %d, prev best was %d.\n", (int)total_cost, (int)ib->cb_best_cost);
        }
        #endif

        venc_keep_best_inter_mode(t, cb, total_cost, F265_PART_UN);
    }

    // FIXME: add pmv test.
}

// Analyze the inter partitioning mode specified.
static void venc_analyze_inter_part_mode(f265_enc_thread *t, f265_cb *cb, int part_mode, int64_t base_cost)
{
    int64_t total_cost = base_cost;

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing inter mode %d.\n", part_mode);
    }
    #endif

    // Set the partitioning mode and update the prediction map. This needs to be
    // optimized.
    cb->inter_part = part_mode;
    venc_update_pmap_unsplit_cb(t, cb);

    // Analyze the partitions.
    int nb_parts = f265_nb_parts[part_mode];
    for (int part_idx = 0; part_idx < nb_parts; part_idx++)
        total_cost += venc_analyze_inter_part(t, cb, part_idx);

    // Add the estimated partitioning mode cost in RDM mode.
    if (t->an.rdm_flag)
    {
        int min_bs_flag = cb->lg_bs == t->enc->gd.cb_range[0];
        total_cost += t->an.se_costs[F265_SE_INTER_PART + 3*!min_bs_flag + part_mode];
    }

    // Compute the real costs in RDO mode. We keep the motion estimation data we
    // have found but ignore the imprecise RDM costs.
    else
        total_cost = venc_rdo_inter_part_mode(t, cb, part_mode, base_cost);

    venc_keep_best_inter_mode(t, cb, total_cost, part_mode);
}

// Sanity test: verify that the actual cost of the inter CB matches the RDM
// analysis cost.
#ifdef VAN_VERIFY_CB_INTER_COST
static void venc_verify_inter_cb_cost(f265_enc_thread *t, f265_cb *cb, int orig_cost)
{
    f265_me_ctx *me = &t->me;
    f265_inter_block *ib = &t->ib;
    int min_bs_flag = cb->lg_bs == t->enc->gd.cb_range[0];
    int skip_idx = F265_SE_CU_SKIP + (venc_get_cu_skip_ctx_off(t, cb)<<1);
    int part_mode = cb->inter_part;
    int nb_parts = f265_nb_parts[part_mode];
    int cb_cost = 0;
    int debug = 0;
    #ifdef VAN_TRACE_ANALYSIS
    debug = venc_trace_analysis_flag;
    #endif

    // Update the prediction map.
    venc_update_pmap_unsplit_cb(t, cb);

    // Do the motion compensation.
    venc_inter_mode_mc(t, cb);

    if (debug) printf("Actual inter mode %d.\n", part_mode);

    // Compute the distortion. We have to follow the partition layout since the
    // SATD depends on the block size.
    {
        int dist_func_id = t->an.mode_metric;
        int nb_comp = 1 + (t->an.chroma_me_flag<<1);

        for (int part_idx = 0; part_idx < nb_parts; part_idx++)
        {
            int part_size[2], part_off[2];
            venc_get_part_loc(part_size, part_off, part_mode, part_idx, 1<<cb->lg_bs);

            for (int comp = 0; comp < nb_comp; comp++)
            {
                int csf = !!comp;
                int ct_ox = (cb->cb_off[0] + part_off[0])>>csf;
                int ct_oy = (cb->cb_off[1] + part_off[1])>>csf;
                int bs_x = part_size[0]>>csf;
                int bs_y = part_size[1]>>csf;
                int ref_stride = t->me.ref_stride;
                int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox, ct_oy);
                f265_pix *src = t->src_frame->src_planes[comp] + plane_off;
                f265_pix *rec = t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off;
                int dist = (me->dist[dist_func_id])(src, ref_stride, rec, ref_stride, bs_x, bs_y, me->bit_depth[0]);
                cb_cost += dist;
                if (debug) printf("Actual dist %d.\n", dist);

                // Dump the pixels for debugging.
                #if 0
                if (debug)
                {
                    printf("Comp %d src:\n", comp);
                    print_block_pix(src, ref_stride, bs_x, bs_y);
                    printf("Comp %d ref:\n", comp);
                    print_block_pix(rec, ref_stride, bs_x, bs_y);
                }
                #endif
            }
        }
    }

    // Skip.
    if (cb->flags&F265_CB_SKIP)
    {
        int skip_cost = t->an.se_costs[skip_idx+1];
        if (debug) printf("Actual skip cost %d.\n", skip_cost);
        cb_cost += skip_cost;
    }

    // Non-skip.
    else
    {
        int skip_cost = t->an.se_costs[skip_idx+0];
        int pred_cost = t->an.se_costs[F265_SE_PRED_MODE+0];
        int part_cost = t->an.se_costs[F265_SE_INTER_PART + 3*!min_bs_flag + part_mode];
        if (debug) printf("Actual skip/pred/part cost %d/%d/%d.\n", skip_cost, pred_cost, part_cost);
        cb_cost += skip_cost + pred_cost + part_cost;

        for (int part_idx = 0; part_idx < nb_parts; part_idx++)
        {
            venc_me_set_partition(t, cb, part_idx);
            int inter_bf = cb->inter_bf[part_idx];
            int merge_idx = inter_bf&7;
            int merge_idx_cost = t->an.se_costs[F265_SE_MERGE_IDX+merge_idx];
            if (debug) printf("Actual merge idx %d cost %d.\n", merge_idx, merge_idx_cost);
            cb_cost += merge_idx_cost;

            // Non-merge.
            if (merge_idx == 5)
            {
                for (int list = 0; list < 2; list++)
                {
                    int ref_idx = cb->ref_idx[part_idx][list];
                    if (ref_idx == -1) continue;

                    int pmv_idx = (inter_bf>>(3+list))&1;
                    f265_mv pmv_array[2];
                    venc_get_pmv(t, part_idx, ib->neighbour_mv, ref_idx, list, pmv_array);
                    f265_mv pmv = pmv_array[pmv_idx];
                    venc_me_set_pmv(me, pmv, t->qp[0], 0);
                    f265_mv mv = cb->mv[part_idx][list];

                    int ref_idx_cost = t->an.se_costs[F265_SE_REF_IDX + list*16 + ref_idx];
                    int mv_cost = venc_me_mv_cost_test(me, mv, 0);
                    int pmv_idx_cost = t->an.se_costs[F265_SE_MVP + pmv_idx];
                    if (debug) printf("Actual ref_idx/mv/pmv %d/(%d,%d)/%d:(%d,%d) cost %d/%d/%d.\n",
                                      ref_idx, mv.x, mv.y, pmv_idx, pmv.x, pmv.y,
                                      ref_idx_cost, mv_cost, pmv_idx_cost);
                    cb_cost += ref_idx_cost + mv_cost + pmv_idx_cost;
                }
            }

            // Merge.
            else
            {
                for (int list = 0; list < 2; list++)
                {
                    int ref_idx = cb->ref_idx[part_idx][list];
                    if (ref_idx == -1) continue;

                    f265_mv mv = cb->mv[part_idx][list];
                    if (debug) printf("Actual ref_idx/mv %d/(%d,%d).\n", ref_idx, mv.x, mv.y);
                }
            }
        }
    }

    if (orig_cost != cb_cost || debug)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("frame %d ctb %d cost verification: expected %d, actual %d.\n",
               (int)t->src_frame->abs_poc, t->ctb_xy, orig_cost, cb_cost);
        if (orig_cost != cb_cost) exit(1);
    }
}
#endif

// Return the cost of encoding the CB with inter prediction.
static int64_t venc_analyze_inter_cb(f265_enc_thread *t, f265_cb *cb)
{
    f265_inter_block *ib = &t->ib;
    f265_stash *s = &t->stash;
    int stash_flag = !!s->flags;

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzing inter CB.\n");
    }
    #endif

    // Set the CB mode to inter.
    F265_SET_FLAG(cb->flags, F265_CB_INTRA, 0);

    // Precompute the skip and non-skip base costs.

    // Base syntax element costs.
    //   base_costs[0]: cu_skip_flag=1
    //   base_costs[1]: cu_skip_flag=0, pred_mode_flag=0
    int64_t base_costs[2];

    // CABAC contexts to restore if skip is selected in RDO mode.
    int skip_ctx_idx;
    uint8_t skip_cabac_ctx[2];

    // Use the approximated costs directly.
    if (t->an.rdm_flag)
    {
        int skip_idx = F265_SE_CU_SKIP + (venc_get_cu_skip_ctx_off(t, cb)<<1);
        base_costs[0] = t->an.se_costs[skip_idx+1];
        base_costs[1] = t->an.se_costs[skip_idx+0] + t->an.se_costs[F265_SE_PRED_MODE+0];
    }

    // Encode both cases and assume the non-skip case will be selected.
    else
    {
        // Get the context index of the skip flag.
        skip_ctx_idx = F265_CO_CU_SKIP + venc_get_cu_skip_ctx_off(t, cb);

        // Preserve the original skip context.
        uint8_t orig_skip_ctx = t->cbs.contexts[skip_ctx_idx];

        // Encode the skip case and save its CABAC contexts.
        F265_RDO_COST(base_costs[0], venc_encode_context_bin(&t->cbs, skip_ctx_idx, 1));
        skip_cabac_ctx[0] = t->cbs.contexts[skip_ctx_idx];
        skip_cabac_ctx[1] = t->cbs.contexts[F265_CO_PRED_MODE];

        // Restore the original skip context and encode the non-skip case.
        t->cbs.contexts[skip_ctx_idx] = orig_skip_ctx;
        t->cbs.rdo.bits = 0;
        venc_encode_context_bin(&t->cbs, skip_ctx_idx, 0);
        venc_encode_context_bin(&t->cbs, F265_CO_PRED_MODE, 0);
        base_costs[1] = venc_rdo_bit_cost(t, t->cbs.rdo.bits);
    }

    // Snapshot the state after the skip/non-skip precomputation.
    if (stash_flag) venc_stash_init_cb(t, cb);

    // Initialize the best inter cost.
    ib->cb_best_cost = F265_MAX_SSD;

    // Test the skip/merge modes.
    venc_analyze_inter_skip(t, cb, base_costs);

    // Test UN/H2/V2.
    for (int part_mode = 0; part_mode < 3; part_mode++)
        venc_analyze_inter_part_mode(t, cb, part_mode, base_costs[1]);

    // Test the AMP modes.
    if (cb->lg_bs > t->enc->gd.cb_range[0] && t->enc->gd.eflags&F265_PF_AMP)
        for (int part_mode = 4; part_mode < 8; part_mode++)
            venc_analyze_inter_part_mode(t, cb, part_mode, base_costs[1]);

    // Refine the non-merge, non-bi-predicted UN mode if it is the best.
    // Disabled for now, this is merely a proof-of-concept.
    if (0 &&
        ib->cb_inter_mode == F265_PART_UN &&
        (ib->cb_inter_bf[0]&7)==5 &&
        (ib->cb_ref_idx[0][0] == -1 || ib->cb_ref_idx[0][1] == -1))
    {
        venc_refine_inter_un_mode(t, cb, base_costs[1]);
    }

    // Restore the best mode data.
    if (stash_flag) venc_stash_restore(t);
    venc_restore_best_inter_mode(t, cb);

    // Restore the skip CABAC contexts.
    if (ib->cb_inter_mode == -1 && !t->an.rdm_flag)
    {
        t->cbs.contexts[skip_ctx_idx] = skip_cabac_ctx[0];
        t->cbs.contexts[F265_CO_PRED_MODE] = skip_cabac_ctx[1];
    }

    // Consider adding a redundancy prevention test eventually (two partitions
    // with the same inter data).

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("analyzed inter CB best mode %d cost %d.\n", ib->cb_inter_mode, (int)ib->cb_best_cost);
    }
    #endif

    return ib->cb_best_cost;
}

// Return the cost of an unsplit CB.
static int64_t venc_analyze_cb(f265_enc_thread *t, f265_cb *cb)
{
    int intra_flag = 1;
    int inter_flag = t->src_frame->frame_type != F265_FRAME_I;
    int stash_flag = intra_flag && inter_flag;
    int64_t best_cost, intra_cost = F265_MAX_SSD, inter_cost = F265_MAX_SSD;
    uint8_t *init_tn = t->tt.tn;

    if (stash_flag) venc_stash_init_cb(t, cb);

    if (inter_flag)
    {
        // Do not save the inter reconstruction. During the transform tree
        // exploration, the motion compensation is saved into the reconstructed
        // planes but the IDCT result is discarded. Hence the prediction is
        // preserved during the search. We reconstruct the CB below if required.
        int saved_stash_flags = t->stash.flags;
        t->stash.flags &= ~7;

        venc_stash_push(t);
        inter_cost = venc_analyze_inter_cb(t, cb);
        venc_stash_pop(t);

        t->stash.flags = saved_stash_flags;

        // Optional verification.
        #ifdef VAN_VERIFY_CB_INTER_COST
        if (t->an.rdm_flag) venc_verify_inter_cb_cost(t, cb, inter_cost);
        #endif
    }

    if (intra_flag)
    {
        if (stash_flag) venc_stash_save_reset(t);

        // Optimize those pushes eventually.
        venc_stash_push(t);
        intra_cost = venc_analyze_intra_cb(t, cb);
        venc_stash_pop(t);
    }

    if (inter_cost < intra_cost)
    {
        if (stash_flag)
        {
            venc_stash_restore(t);
            F265_SET_FLAG(cb->flags, F265_CB_INTRA, 0);
        }

        // Reconstruct the CB with the best inter mode.
        t->tt.tn = init_tn;

        // Don't reconstruct if delaying the transform tree exploration.
        if (t->enc->gd.algo & (1<<11))
        {
        }

        // Fast transform tree exploration.
        else if (t->an.rdm_flag)
        {
            // Do the motion compensation.
            venc_inter_mode_mc(t, cb);

            // We need to obtain a transform tree in RDM mode. Using the RDO
            // logic for now.
            f265_cabac_bs saved_cbs = t->cbs;
            int saved_stash_flags = t->stash.flags;
            t->stash.flags = 16+8;

            int depth_range[2];
            venc_get_inter_part_depth_range(t, cb->inter_part, cb->lg_bs, depth_range);

            venc_stash_push(t);
            int nz_flags;
            venc_analyze_inter_sub_tt(t, cb, &nz_flags, depth_range, 0, 0, cb->lg_bs, 0, 0);
            venc_stash_pop(t);

            t->stash.flags = saved_stash_flags;
            t->cbs = saved_cbs;
        }

        // Just do the reconstruction.
        else
        {
            // FIXME: we can skip the reconstruction altogether for the CBs are not
            // used for the intra prediction of other CBs. The final encode will
            // take care of reconstructing correctly.
            //
            // FIXME WARNING: this code is setting the YUV flags in the transform
            // nodes to their actual non-zero values. This is not what we want for
            // the final RDOQ. Also, we're writing the coefficients uselessly.
            venc_set_tmp_tb(t);
            venc_rec_inter_cb(t, cb);
        }

        best_cost = inter_cost;
    }

    else
    {
        F265_SET_FLAG(cb->flags, F265_CB_SKIP, 0);
        best_cost = intra_cost;

        // If we haven't analyzed intra chroma, do it now to get a valid mode.
        if (t->an.rdm_flag && !t->an.chroma_me_flag && !(t->enc->gd.algo & (1<<11)))
            venc_analyze_intra_cb_chroma(t, cb, init_tn);
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("finished CB mode analysis best %d (intra %d, inter %d).\n",
               (int)best_cost, (int)intra_cost, (int)inter_cost);
    }
    #endif

    return best_cost;
}

// Return the cost of the CB (split or unsplit).
static int64_t venc_visit_cb_bottom_up(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return 0;

    int split_flag = cb->lg_bs > t->enc->gd.cb_range[0];
    int unsplit_flag = !F265_GET_FLAG(cb->flags, F265_CB_FORBIDDEN);
    int stash_flag = t->stash.flags && unsplit_flag && split_flag;
    int64_t split_cost = F265_MAX_SSD, unsplit_cost = F265_MAX_SSD, best_cost;

    // Get the split flag costs.
    uint16_t split_flag_costs[2] = { 0, 0 };
    #if 0
    if (split_flag && unsplit_flag)
    #else
    if (t->an.rdm_flag && split_flag && unsplit_flag)
    #endif
    {
        int off = venc_get_split_cu_ctx_off(t, cb)<<1;
        for (int i = 0; i < 2; i++) split_flag_costs[i] = t->an.se_costs[F265_SE_SPLIT_CU + off + i];
    }

    if (stash_flag) venc_stash_init_cb(t, cb);

    // Analyze split CB.
    if (split_flag)
    {
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("@@ analyzing split CB (size %d).\n", 1<<cb->lg_bs);
        }
        #endif

        int64_t split_flag_cost;
        if (t->an.rdm_flag) split_flag_cost = split_flag_costs[1];
        else { F265_RDO_COST(split_flag_cost, venc_write_split_cu_flag(t, cb, 1)); }

        venc_stash_push(t);
        int64_t subtree_cost = 0;
        for (int i = 0; i < 4; i++) subtree_cost += venc_visit_cb_bottom_up(t, t->cb + cb->child_idx + i);
        venc_stash_pop(t);

        split_cost = split_flag_cost + subtree_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzed split CB cost %d (split_flag %d, subtree %d).\n",
                   (int)split_cost, (int)split_flag_cost, (int)subtree_cost);
        }
        #endif
    }

    // Analyze unsplit CB.
    if (unsplit_flag)
    {
        if (stash_flag) venc_stash_save_reset(t);

        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("@@ analyzing unsplit CB (size %d).\n", 1<<cb->lg_bs);
        }
        #endif

        int64_t unsplit_flag_cost;
        if (t->an.rdm_flag) unsplit_flag_cost = split_flag_costs[0];
        else { F265_RDO_COST(unsplit_flag_cost, venc_write_split_cu_flag(t, cb, 0)); }

        venc_stash_push(t);
        int64_t cb_cost = venc_analyze_cb(t, cb);
        venc_stash_pop(t);

        unsplit_cost = unsplit_flag_cost + cb_cost;
        #ifdef VAN_TRACE_ANALYSIS
        if (venc_trace_analysis_flag)
        {
            venc_an_cb_loc(t, cb, 0);
            printf("analyzed unsplit CB cost %d (split_flag %d, CB %d).\n",
                   (int)unsplit_cost, (int)unsplit_flag_cost, (int)cb_cost);
        }
        #endif
    }

    if (split_cost < unsplit_cost)
    {
        if (stash_flag) venc_stash_restore(t);
        best_cost = split_cost;
        F265_SET_FLAG(cb->flags, F265_CB_SPLIT, 1);
        // FIXME: we might want to back this up instead of doing recursive
        // calls? 6 32-byte stores for the whole CTB minus the first row would
        // do the job.
        venc_update_pmap_cb(t, cb);
    }

    else
    {
        best_cost = unsplit_cost;
        F265_SET_FLAG(cb->flags, F265_CB_SPLIT, 0);
        venc_update_pmap_unsplit_cb(t, cb);
    }

    #ifdef VAN_TRACE_ANALYSIS
    if (venc_trace_analysis_flag)
    {
        venc_an_cb_loc(t, cb, 0);
        printf("finished CB analysis best %d (split %d, unsplit %d).\n",
               (int)best_cost, (int)split_cost, (int)unsplit_cost);
    }
    #endif

    return best_cost;
}

// Obtain valid transform trees for the chosen CTB layout.
static void venc_analyze_ctb_tt(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return;

    // Split CB.
    if (cb->flags&F265_CB_SPLIT)
    {
        for (int i = 0; i < 4; i++) venc_analyze_ctb_tt(t, t->cb + cb->child_idx + i);
        return;
    }

    // Unsplit CB.
    int lg_bs = cb->lg_bs;
    uint8_t *init_tn = t->tt.tn;

    // Intra.
    if (cb->flags&F265_CB_INTRA)
    {
        // HV.
        if (cb->intra_luma_mode[1] != -1)
        {
            // Add a split transform node.
            *t->tt.tn++ = 15;

            int depth_range[2];
            venc_get_intra_part_depth_range(t, lg_bs-1, depth_range);

            for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
                venc_analyze_intra_luma_tt(t, cb, cb->intra_luma_mode[i], depth_range, 0,
                                           lg_bs-1, (i&1)*sbs, (i>>1)*sbs);
        }

        // UN.
        else
        {
            int depth_range[2];
            venc_get_intra_part_depth_range(t, lg_bs, depth_range);
            venc_analyze_intra_luma_tt(t, cb, cb->intra_luma_mode[0], depth_range, 0, lg_bs, 0, 0);
        }

        // Set the intra chroma mode and transform trees.
        venc_analyze_intra_cb_chroma(t, cb, init_tn);
    }

    // Inter.
    else
    {
        // Do the motion compensation.
        venc_inter_mode_mc(t, cb);

        // Analyze the transform tree.
        uint8_t *init_tn = t->tt.tn;
        venc_analyze_inter_tt(t, cb);

        // Reconstruct.
        t->tt.tn = init_tn;
        venc_set_tmp_tb(t);
        venc_rec_inter_cb(t, cb);
    }
}

// Validation pass.
static void venc_analyze_ctb_single_pass(f265_enc_thread *t)
{
    f265_gen_data *gd = &t->enc->gd;
    f265_an_block *an = &t->an;

    // Set the analyzed TB range to the actual parameter values.
    for (int i = 0; i < 2; i++)
    {
        an->tb_range[i] = gd->tb_range[i];
        an->tb_depth[i] = gd->tb_depth[i];
    }

    // Do the analysis.
    venc_visit_cb_bottom_up(t, t->cb);

    // Obtain valid transform trees for the chosen CTB layout.
    if (t->enc->gd.algo & (1<<11))
    {
        // Switch to RDO.
        int saved_stash_flags = t->stash.flags;
        t->stash.flags = 31;
        t->an.rdm_flag = 0;

        // Reset the transform tree.
        t->tt.tn = t->tmap;

        // Analyze.
        venc_analyze_ctb_tt(t, t->cb);

        // Switch to RDM.
        t->stash.flags = saved_stash_flags;
        t->an.rdm_flag = 1;
    }
}

// Analyze the current CTB. Return true if the CTB needs to be reconstructed
// (probably a design error -- fix eventually).
int venc_analyze_ctb(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_SYNTAX
    venc_trace_syntax_flag = 0;
    #endif

    #ifdef VAN_TRACE_ANALYSIS
    // FIXME. Used for debug.
    venc_trace_analysis_flag = t->src_frame->abs_poc == 1 && t->ctb_xy == 0;
    if (venc_trace_analysis_flag)
    {
        printf("\nFrame %d CTB %d (%d,%d) pix (%d,%d).\n",
               (int)t->src_frame->abs_poc, t->ctb_xy, t->ctb_x, t->ctb_y, t->ctb_off[0], t->ctb_off[1]);
    }
    #endif

    // Back up the raw CABAC object.
    f265_cabac_bs saved_cbs = t->cbs;

    // Remember whether we consider the chroma cost in the RDM analysis.
    t->an.chroma_me_flag = t->me.early_me_params->dist_func_ids[2]&1;

    // Set the mode comparison metric. Currently using the qpel distortion
    // metric.
    t->an.mode_metric = t->me.early_me_params->dist_func_ids[2]>>1;

    // Set the lambdas.
    t->an.rdo_lambda = t->hm_lambda[0]*256.0f + 0.5f;
    t->me.lambda = t->an.rdm_lambda = (sqrt(t->hm_lambda[0])*256.0f) + 0.5f;

    // Update the costs based on the fixed/initial/current CABAC context values.
    // Experimental code.
    if (t->ctb_xy == 0 || 0) venc_update_se_costs(t);

    // Set the stash flags.
    t->stash.flags = t->an.rdm_flag ? 8 : 31;

    // Initialize the transform tree.
    venc_init_transform_tree(t);

    // Do a single pass.
    venc_analyze_ctb_single_pass(t);

    // Optionally verify the RDO analysis.
    #ifdef VAN_VERIFY_RDO_ANALYSIS
    venc_verify_analysis_save(t, 0);
    #endif

    // Restore the CABAC object.
    t->cbs = saved_cbs;

    return 1;
}

// Return the entropy associated to the context and bin value specified in
// fractions of 2048.
static inline int venc_get_cabac_entropy(f265_enc_thread *t, int ctx_idx, int bin)
{
    // Experiments.

    #if 0
    int ctx = t->cbs.contexts[ctx_idx];
    return f265_cabac_entropy_table[ctx^bin]>>4;
    #endif

    #if 0
    int ctx = t->cbs.contexts[ctx_idx];
    return f265_cabac_unskewed_entropy_table[ctx^bin]>>4;
    #endif

    #if 1
    return 32768>>4;
    #endif
}

// Update the syntax element costs based on the current values of the CABAC
// contexts.
//
// The analysis of the worst case for maximizing the precision is as follow.
// The maximum entropy loss is less than 5.8 bits for a context bin. The
// worst-case binarization is for "inter_part", which uses 3 context bins and 1
// bypass bin. The maximum bit cost is 19 (3*5.8 + 1). The maximum lambda value
// is less than 100.
void venc_update_se_costs(f265_enc_thread *t)
{
    uint16_t *se_costs = t->an.se_costs;
    int i_flag = t->src_frame->frame_type == F265_FRAME_I;
    int b_flag = t->src_frame->frame_type == F265_FRAME_B;
    int bypass_val = 32768>>4;
    int min_entropy_val = f265_cabac_entropy_table[127]>>4;

    // Intra.
    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            se_costs[F265_SE_SPLIT_CU+2*i+j] = venc_get_cabac_entropy(t, F265_CO_SPLIT_CU+i, j);
            se_costs[F265_SE_SPLIT_TRANSFORM+2*i+j] = venc_get_cabac_entropy(t, F265_CO_SPLIT_TRANSFORM+i, j);
        }

        se_costs[F265_SE_PRED_MODE+j] = i_flag*venc_get_cabac_entropy(t, F265_CO_PRED_MODE, j);
        se_costs[F265_SE_INTRA_PART+j] = venc_get_cabac_entropy(t, F265_CO_PART_MODE, j);
    }

    {
        int luma_pred_0 = venc_get_cabac_entropy(t, F265_CO_INTRA_LUMA_PRED, 0);
        int luma_pred_1 = venc_get_cabac_entropy(t, F265_CO_INTRA_LUMA_PRED, 1);
        se_costs[F265_SE_INTRA_LUMA_MODE+0] = luma_pred_1 + bypass_val;
        se_costs[F265_SE_INTRA_LUMA_MODE+1] = se_costs[F265_SE_INTRA_LUMA_MODE+2] = luma_pred_1 + 2*bypass_val;
        se_costs[F265_SE_INTRA_LUMA_MODE+3] = luma_pred_0 + 5*bypass_val;

        int chroma_pred_0 = venc_get_cabac_entropy(t, F265_CO_INTRA_CHROMA_PRED, 0);
        int chroma_pred_1 = venc_get_cabac_entropy(t, F265_CO_INTRA_CHROMA_PRED, 1);
        se_costs[F265_SE_INTRA_CHROMA_MODE+0] = chroma_pred_0;
        se_costs[F265_SE_INTRA_CHROMA_MODE+1] = chroma_pred_1 + 2*bypass_val;
    }
    se_costs[F265_SE_BYPASS] = bypass_val;

    // Inter.
    if (!i_flag)
    {
        // CB skip.
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 3; i++)
                se_costs[F265_SE_CU_SKIP+2*i+j] = venc_get_cabac_entropy(t, F265_CO_CU_SKIP+i, j);

        // Partitioning mode.
        {
            uint16_t ctx_costs[3][2];
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 2; i++)
                    ctx_costs[j][i] = venc_get_cabac_entropy(t, F265_CO_PART_MODE+j, i);

            // UN.
            se_costs[F265_SE_INTER_PART+0] = se_costs[F265_SE_INTER_PART+3] = ctx_costs[0][1];

            // Minimum block size, H2 and V2.
            se_costs[F265_SE_INTER_PART+1] = ctx_costs[0][0] + ctx_costs[1][1];
            se_costs[F265_SE_INTER_PART+2] = ctx_costs[0][0] + ctx_costs[1][0];

            // Since we know we don't use inter HV, we set the minimum cost.
            if (t->enc->gd.cb_range[0] != 3) se_costs[F265_SE_INTER_PART+2] += min_entropy_val;

            // Non-minimum block size, AMP enabled.
            if (t->enc->gd.eflags&F265_PF_AMP)
            {
                // Exclude UN and ignore HV.
                for (int mode = 1; mode < 8; mode++)
                {
                    uint32_t ctx_str = 0x00220460;
                    int cost = mode > 3 ? bypass_val : 0;
                    for (int i = 0; i < 3; i++) cost += ctx_costs[i][(ctx_str>>((mode<<2)+i))&1];
                    se_costs[F265_SE_INTER_PART+3+mode] = cost;
                }
            }

            // Non-minimum block size, AMP disabled.
            else
            {
                se_costs[F265_SE_INTER_PART+3+1] = ctx_costs[0][0] + ctx_costs[1][1];
                se_costs[F265_SE_INTER_PART+3+2] = ctx_costs[0][0] + ctx_costs[1][0];
            }
        }

        // Merge.
        {
            int nb_cand = t->enc->gd.merge_cand;

            int flag_ctx[2], idx_ctx[2];
            for (int i = 0; i < 2; i++)
            {
                flag_ctx[i] = venc_get_cabac_entropy(t, F265_CO_MERGE_FLAG, i);
                idx_ctx[i] = venc_get_cabac_entropy(t, F265_CO_MERGE_IDX, i);
            }

            se_costs[F265_SE_MERGE_IDX+5] = flag_ctx[0];

            if (nb_cand == 1)
                se_costs[F265_SE_MERGE_IDX+0] = flag_ctx[1];

            else
            {
                int gt0_cost = flag_ctx[1] + idx_ctx[1];
                se_costs[F265_SE_MERGE_IDX+0] = flag_ctx[1] + idx_ctx[0];
                for (int i = 1; i < nb_cand - 1; i++) se_costs[F265_SE_MERGE_IDX+i] = gt0_cost + i*bypass_val;
                se_costs[F265_SE_MERGE_IDX+nb_cand-1] = gt0_cost + (nb_cand-2)*bypass_val;
            }
        }

        // MVP.
        for (int i = 0; i < 2; i++)
            se_costs[F265_SE_MVP+i] = venc_get_cabac_entropy(t, F265_CO_MVP, i);

        // Reference indices.
        {
            uint16_t ctx_costs[2][2];
            for (int j = 0; j < 2; j++)
                for (int i = 0; i < 2; i++)
                    ctx_costs[j][i] = venc_get_cabac_entropy(t, F265_CO_REF_IDX+j, i);

            for (int list = 0; list < t->nb_lists; list++)
            {
                uint16_t *list_costs = se_costs + F265_SE_REF_IDX + list*16;
                int nb_ref_idx = t->nb_ref_idx[list];

                if (nb_ref_idx == 1)
                    list_costs[0] = 0;

                if (nb_ref_idx > 1)
                {
                    list_costs[0] = ctx_costs[0][0];
                    list_costs[1] = ctx_costs[0][1];
                }

                if (nb_ref_idx > 2)
                {
                    list_costs[1] += ctx_costs[1][0];

                    int gt1_cost = ctx_costs[0][1] + ctx_costs[1][1];
                    for (int i = 2; i < nb_ref_idx - 1; i++) list_costs[i] = gt1_cost + (i-1)*bypass_val;
                    list_costs[nb_ref_idx-1] = gt1_cost + (nb_ref_idx-3)*bypass_val;
                }
            }
        }

        // Inter prediction indicator.
        if (b_flag)
        {
            // FIXME. Store the flag values for every context (5*2 entries).
            // Null on P frames.
        }
    }

    int shift = 11+8;
    for (int i = 0; i < F265_SE_SIZE; i++)
        se_costs[i] = ((int64_t)se_costs[i]*t->an.rdm_lambda + (1<<(shift-1)))>>shift;
}

// HM RDO analysis.

#define F265_MAX_RD_COST 10000000000.0

static int32_t venc_hm_an_cb(f265_enc_thread *t, f265_cb *cb);
static int32_t venc_hm_an_intra_tt(f265_enc_thread *t, f265_cb *cb, uint8_t lg_bs, int mode,
                                   int ox, int oy, int depth);

#ifdef VAN_TRACE_ANALYSIS
int venc_print_hm_rdo_trace;

static finline void venc_trace_rdo_bits(const char *where, uint64_t count)
{
    if (venc_print_hm_rdo_trace) printf("%s bits: %llu\n", where, (unsigned long long)count);
}

static finline void venc_trace_rdo_cost(const char *where, double cost)
{
    if (venc_print_hm_rdo_trace) printf("%s cost: %.4f\n", where, cost);
}

static finline void venc_trace_rdo_dist(const char *where, int dist)
{
    if (venc_print_hm_rdo_trace) printf("%s dist: %d\n", where, dist);
}

static finline void venc_trace_rdo_path_cb(int bs)
{
    if (venc_print_hm_rdo_trace) printf("CB %dx%d\n", bs, bs);
}

static finline void venc_trace_rdo_path_tb(const char *where, int bs, int mode, f265_enc_thread *t,
                                           int ox, int oy, int chroma_flag)
{
    if (venc_print_hm_rdo_trace)
        printf("%s TB %dx%d (mode=%d)\n(ox,oy)=(%d,%d)\n", where, bs, bs, mode,
               (t->ctb_off[0]>>chroma_flag)+ox, (t->ctb_off[1]>>chroma_flag)+oy);
}

static finline void venc_trace_rdo_luma_modes(f265_cb *cb)
{
    if (!venc_print_hm_rdo_trace) return;

    if (cb->intra_luma_mode[1] == -1)
        printf("Best mode = %d\n", (int)cb->intra_luma_mode[0]);

    else
        printf("Best mode = { %d, %d, %d, %d }\n", (int)cb->intra_luma_mode[0],
                                                   (int)cb->intra_luma_mode[1],
                                                   (int)cb->intra_luma_mode[2],
                                                   (int)cb->intra_luma_mode[3]);
}

static finline void venc_trace_rdo_chroma_mode(f265_cb *cb)
{
    if (venc_print_hm_rdo_trace) printf("Best chroma mode = %d\n", (int)cb->intra_chroma_mode);
}

#else
static finline void venc_trace_rdo_bits(const char *where, uint64_t count) {}
static finline void venc_trace_rdo_cost(const char *where, double cost) {}
static finline void venc_trace_rdo_dist(const char *where, int dist) {}
static finline void venc_trace_rdo_path_cb(int bs) {}
static finline void venc_trace_rdo_path_tb(const char *where, int bs, int mode, f265_enc_thread *t,
                                           int ox, int oy, int chroma_flag)
{}
static finline void venc_trace_rdo_luma_modes(f265_cb *cb) {}
static finline void venc_trace_rdo_chroma_mode(f265_cb *cb) {}
#endif

// Calculate rate-distortion cost.
static finline double venc_hm_an_cost(f265_enc_thread *t, f265_cabac_bs *cbs, int dist)
{
    // The "0.5" rounding factor is to mimic HM.
    return (double)dist + floor(t->hm_lambda[0] * (double)(cbs->rdo.bits>>15) + 0.5);
}

// Calculate sum of squared error.
static finline int32_t venc_hm_an_dist(f265_enc_thread *t, int lg_bs, int ct_ox, int ct_oy, int comp)
{
    int bs = 1<<lg_bs;
    int bd = t->enc->gd.bit_depth[0];
    int stride = t->me.ref_stride;
    int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox, ct_oy);
    int32_t dist = venc_ssd(t->src_frame->src_planes[comp] + plane_off, stride,
                            t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off, stride,
                            bs, bs, 0, bd);
    return !comp ? dist : (int32_t)(dist * t->hm_wcd);
}

// Signal the luma prediction mode.
static finline void venc_hm_an_luma_pred_cost(f265_enc_thread *t, f265_cb *cb, int part)
{
    int mode_val = venc_write_intra_luma_pred(t, cb, part);
    if (mode_val < 3) venc_write_intra_mpm_idx(t, mode_val);
    else venc_write_intra_rem_mode(t, mode_val-3);
}

// Recursively signal the chroma CBF flags from the current transform tree
// location.
static void venc_signal_intra_chroma_cbf(f265_enc_thread *t, int lg_bs, int depth, int cbf_mask)
{
    uint8_t tn = *t->tt.tn++;

    // Signal flags for current TB.
    if (lg_bs >= 2)
    {
        // FIXME. There should be an finline function in entropy.c to signal
        //        the cbf_cb and cbf_cr flags.
        for (int c = 1; c < 3; c++)
            if (depth == 0 || cbf_mask&(1<<c))
            {
                venc_trace_syntax_element(c == 1 ? "cbf_cb" : "cbf_cr");
                venc_encode_context_bin(&t->cbs, F265_CO_CBF_CHROMA+depth, (tn>>c)&1);
            }
    }

    // Signal flags for child TBs.
    if (tn&8)
    {
        int sub_lg_bs = lg_bs-1, next_depth = depth+1;
        for (int part = 0; part < 4; part++) venc_signal_intra_chroma_cbf(t, sub_lg_bs, next_depth, tn);
    }
}

// Reconstruct either Cb, Cr or both components following the luma transform
// tree. Recursively update the CBF flags found in the tree.
// NOTE The HM computes the distortion at each TB instead of computing it at
//      the CB level. Here, we mimic that process because of the weighted chroma
//      distortion. w*A + w*B + w*C <= w*(A+B+C) because the distortion is not
//      kept with floating point precision. Now because of this, the distortion
//      is computed too often in less than optimal fashion. When bit-exactness
//      is no longer required, factor out the distortion from this function and
//      have the caller compute it afterwards.
static int32_t venc_hm_an_intra_chroma_tt(f265_enc_thread *t, int mode, int lg_bs,
                                          int ox, int oy, int uv_mask, int *dist, int depth)
{
    uint8_t *tn = t->tt.tn++;
    int split_flag = *tn&8;

    // Clear Cb, Cr or both flags.
    *tn &= ~(uv_mask<<1);

    // Process unsplit TBs or split luma 8x8 TBs.
    if ((!split_flag && lg_bs > 1) || (split_flag && lg_bs == 2))
    {
        if (uv_mask == 1) venc_trace_rdo_path_tb("Intra Chroma", 1<<lg_bs, mode, t, ox, oy, 1);

        f265_cabac_bs *cbs = &t->cbs;
        for (int comp = 1; comp < 3; comp++)
            if (!!(comp&uv_mask))
            {
                venc_set_tmp_tb(t);
                int nz_flag = venc_rec_intra_tb(t, comp, lg_bs, mode, 0, ox, oy, depth);
                if (nz_flag)
                {
                    venc_set_tmp_tb(t);
                    venc_write_tb(cbs, &t->tt, 1);
                    *tn |= 1<<comp;
                }
                *dist += venc_hm_an_dist(t, lg_bs, ox, oy, comp);
            }
    }

    // Recursively process child TBs.
    if (split_flag)
    {
        int nz_flag = 0;
        for (int c = 0, sub_lg_bs = lg_bs-1; c < 4; c++)
            nz_flag |= venc_hm_an_intra_chroma_tt(t, mode, sub_lg_bs, ox+((c&1)<<sub_lg_bs), oy+((c>>1)<<sub_lg_bs),
                                                  uv_mask, dist, depth+1);
        *tn |= nz_flag;
    }

    return *tn&7;
}

// Intra chroma analysis. Follows the transform tree selected with the luma
// component. Try all 5 modes.
static int32_t venc_hm_an_intra_chroma_mode(f265_enc_thread *t, f265_cb *cb)
{
    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    // Save input state.
    f265_cabac_bs *cbs = &t->cbs;
    f265_cabac_snapshot cs;
    venc_save_cabac_contexts(cbs, &cs);
    uint8_t *tn = t->tt.tn;
    uint64_t bits = cbs->rdo.bits;

    // Try all available chroma angles.
    int lg_bs = cb->lg_bs-1;
    int cb_ox = cb->cb_off[0]>>1, cb_oy = cb->cb_off[1]>>1;
    uint8_t best_mode = 0;
    const uint8_t chroma_modes[5] = { 0, 26, 10, 1, 36 };
    for (int mode = 0; mode < 5; mode++)
    {
        // Convert mode index to actual value.
        int chroma_mode = chroma_modes[mode];
        if (chroma_mode == cb->intra_luma_mode[0]) chroma_mode=34;
        else if (chroma_mode == 36) chroma_mode = cb->intra_luma_mode[0];

        // Signal intra chroma pred mode.
        cb->intra_chroma_mode = chroma_mode;
        venc_write_intra_chroma_mode(t, cb);

        // Optimize chroma components independently.
        dist[1] = 0;
        for (int comp = 1; comp < 3; comp++)
        {
            venc_hm_an_intra_chroma_tt(t, cb->intra_chroma_mode, lg_bs, cb_ox, cb_oy, comp, &dist[1], 0);
            t->tt.tn = tn;
        }
        venc_signal_intra_chroma_cbf(t, lg_bs, 0, 0);
        cost[1] = venc_hm_an_cost(t, cbs, dist[1]);

        venc_trace_rdo_dist("Chroma", dist[1]);
        venc_trace_rdo_bits("Chroma", cbs->rdo.bits);
        venc_trace_rdo_cost("Chroma", cost[1]);

        // Track best mode.
        if (cost[1] < cost[0])
        {
            cost[0] = cost[1];
            dist[0] = dist[1];
            best_mode = chroma_mode;
        }

        // Restore input state.
        venc_load_cabac_contexts(cbs, &cs);
        t->tt.tn = tn;
        cbs->rdo.bits = bits;
    }

    venc_trace_rdo_dist("Best chroma", dist[0]);
    venc_trace_rdo_cost("Best chroma", cost[0]);

    // Encode chroma components correctly using the best selected mode.
    cb->intra_chroma_mode = best_mode;
    venc_write_intra_chroma_mode(t, cb);
    venc_hm_an_intra_chroma_tt(t, cb->intra_chroma_mode, lg_bs, cb_ox, cb_oy, 3, &dist[1], 0);
    t->tt.tn = tn;
    venc_signal_intra_chroma_cbf(t, lg_bs, 0, 0);
    cost[0] = venc_hm_an_cost(t, cbs, dist[0]);

    return dist[0];
}

// Split the current TB into four equal parts. Derive cost. Return distortion.
static finline void venc_hm_an_intra_tt_split(f265_enc_thread *t, f265_cb *cb, uint8_t lg_bs, int mode,
                                              int cb_ox, int cb_oy, int32_t *dist, double *cost, int depth)
{
    f265_cabac_bs *cbs = &t->cbs;

    // Signal split TT.
    *t->tt.tn &= 0xf0;
    *t->tt.tn++ |= 8;

    // Process child TBs.
    uint64_t bits = 0;
    for (int part = 0, sub_lg_bs = lg_bs-1; part < 4; part++)
    {
        *dist += venc_hm_an_intra_tt(t, cb, sub_lg_bs, mode, cb_ox+((part&1)<<sub_lg_bs), cb_oy+((part>>1)<<sub_lg_bs),
                                     depth);
         bits += cbs->rdo.bits>>15;

        // Keep "unoutputted" engine state. (This mimics the HM.)
        cbs->rdo.bits &= 32767;
    }
    cbs->rdo.bits += bits<<15;

    // FIXME. There should be an finline function in entropy.c to signal
    //        the split_transform_flag.
    if (lg_bs <= t->enc->gd.tb_range[1])
    {
        venc_trace_syntax_element("split_transform_flag");
        venc_encode_context_bin(cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 1);
    }

    // Compute cost of 4 optimized children.
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Split TT", *dist);
    venc_trace_rdo_bits("Split TT", cbs->rdo.bits);
    venc_trace_rdo_cost("Split TT", *cost);
}

// Analyze distortion and cost of current TB.
static finline void venc_hm_an_intra_tt_unsplit(f265_enc_thread *t, f265_cb *cb, uint8_t lg_bs, int mode,
                                                int cb_ox, int cb_oy, int32_t *dist, double *cost, int depth)
{
    f265_cabac_bs *cbs = &t->cbs;
    uint8_t *tn = t->tt.tn++;

    // Signal unsplit TT.
    *tn &= 0xf0;

    // Signal split_transform_flag.
    // FIXME. There should be an finline function in entropy.c to signal
    //        the split_transform_flag.
    f265_gen_data *gd = &t->enc->gd;
    if (lg_bs > gd->tb_range[0] && cb->lg_bs - lg_bs < gd->tb_depth[0])
    {
        venc_trace_syntax_element("split_transform_flag");
        venc_encode_context_bin(cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, 0);
    }

    // DCT/Quant/Dequant/IDCT.
    venc_set_tmp_tb(t);
    int ct_ox = cb->cb_off[0] + cb_ox, ct_oy = cb->cb_off[1] + cb_oy;
    int nz_flag = venc_rec_intra_tb(t, 0, lg_bs, mode, 0, ct_ox, ct_oy, depth);

    // Write TB.
    // FIXME. There should be an finline function in entropy.c to signal
    //        the cbf_luma.
    int ctx_off = cb->lg_bs == lg_bs;
    venc_trace_syntax_element("cbf_luma");
    venc_encode_context_bin(cbs, F265_CO_CBF_LUMA+ctx_off, nz_flag);
    if (nz_flag)
    {
        venc_set_tmp_tb(t);
        venc_write_tb(cbs, &t->tt, 0);
        *tn |= 1;
    }

    *dist = venc_hm_an_dist(t, lg_bs, ct_ox, ct_oy, 0);
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Unsplit TT", *dist);
    venc_trace_rdo_bits("Unsplit TT", cbs->rdo.bits);
    venc_trace_rdo_cost("Unsplit TT", *cost);
}

// Recursively analyze the transform tree using the selected luma prediction
// mode. Keep the best layout.
static int32_t venc_hm_an_intra_tt(f265_enc_thread *t, f265_cb *cb, uint8_t lg_bs,
                                   int mode, int cb_ox, int cb_oy, int depth)
{
    int intra_hv = cb->intra_luma_mode[1] != -1;
    f265_gen_data *gd = &t->enc->gd;
    int unsplit_tt_flag = lg_bs <= gd->tb_range[1];
    int split_tt_flag = lg_bs > gd->tb_range[0] && cb->lg_bs - lg_bs  < gd->tb_depth[0] + intra_hv;

    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    venc_trace_rdo_path_tb("Intra Luma", 1<<lg_bs, mode, t, cb->cb_off[0]+cb_ox, cb->cb_off[1]+cb_oy, 0);

    // Remember current transform tree offset and CABAC contexts.
    venc_stash_init_cb_off(t, cb, cb_ox, cb_oy, 1<<lg_bs, 1<<lg_bs);

    // Remember input bit count.
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;

    // Unsplit TT.
    if (unsplit_tt_flag)
    {
        venc_stash_push(t);
        venc_hm_an_intra_tt_unsplit(t, cb, lg_bs, mode, cb_ox, cb_oy, &dist[0], &cost[0], depth);
        venc_stash_pop(t);
    }

    // Split TT.
    if (split_tt_flag)
    {
        // Stash unsplit TT context (if available).
        if (unsplit_tt_flag)
        {
            // Save transform tree, luma samples and CABAC contexts.
            // Revert to input transform tree position and CABAC contexts.
            venc_stash_save_reset(t);

            // Remember current bit count. Reload input bit count.
            F265_SWAP(uint64_t, cbs->rdo.bits, bits);
        }

        venc_stash_push(t);
        venc_hm_an_intra_tt_split(t, cb, lg_bs, mode, cb_ox, cb_oy, &dist[1], &cost[1], depth+1);
        venc_stash_pop(t);

        // Restore unsplit TT results if split TT doesn't perform better.
        if (cost[0] <= cost[1])
        {
            // Restore unplit transform tree, luma samples and CABAC contexts.
            venc_stash_restore(t);

            // Reload bit count following unsplit TT exploration.
            cbs->rdo.bits = bits;
        }

        // Split distortion is better.
        else
            return dist[1];
    }

    // Return best distortion.
    return dist[0];
}

// Intra luma prediction mode analysis.
// Try all 35 modes: planar, DC and 33 angular modes.
static int32_t venc_hm_an_intra_luma_mode(f265_enc_thread *t, f265_cb *cb, int split_flag, int part_idx)
{
    int lg_bs = cb->lg_bs-split_flag;
    int cb_ox = (part_idx&0x1)<<lg_bs, cb_oy = (part_idx>>1)<<lg_bs;

    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    // Save current state (transform tree pointer offset, CABAC contexts).
    venc_stash_init_cb_off(t, cb, cb_ox, cb_oy, 1<<lg_bs, 1<<lg_bs);

    // Remember input bit count.
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits, best_bits;

    int best_mode = 0;
    for (int angle = 0; angle < 35; angle ++)
    {
        // Signal intra prediction.
        cb->intra_luma_mode[part_idx] = angle;
        venc_hm_an_luma_pred_cost(t, cb, part_idx);

        // Transform tree analysis.
        venc_stash_push(t);
        dist[1] = venc_hm_an_intra_tt(t, cb, lg_bs, angle, cb_ox, cb_oy, 0);
        venc_stash_pop(t);
        cost[1] = venc_hm_an_cost(t, &t->cbs, dist[1]);

        // Keep best mode.
        if (cost[1] < cost[0])
        {
            // Save transform tree, luma samples, CABAC contexts.
            venc_stash_save(t);

            dist[0] = dist[1];
            cost[0] = cost[1];

            // Remember current bit count and current mode.
            best_bits = cbs->rdo.bits;
            best_mode = angle;
        }

        // Reload input state (transform tree offset and CABAC contexts).
        venc_stash_reset(t);

        // Restore input bit count.
        cbs->rdo.bits = bits;
    }

    // Restore best mode's transform tree, luma samples and CABAC contexts.
    venc_stash_restore(t);

    // Remember best mode.
    cb->intra_luma_mode[part_idx] = best_mode;

    // Restore bit count.
    cbs->rdo.bits = best_bits;

    return dist[0];
}

// Reencode the luma residual following the known transform tree.
static void venc_hm_an_intra_hv_bits(f265_enc_thread *t, f265_cb *cb, int lg_bs, int mode, int cb_ox, int cb_oy,
                                     int depth)
{
    f265_cabac_bs *cbs = &t->cbs;
    uint8_t tn_val = *t->tt.tn++;
    int split_flag = !!(tn_val&8);

    // Signal the split_transform_flag (when present).
    // Note that the split_transform_flag at depth 0 is
    // derived because the CB uses 4 PBs. The flags here
    // are for depths >= 1.
    // FIXME. There should be an finline function in entropy.c to signal
    //        the split_transform_flag.
    f265_gen_data *gd = &t->enc->gd;
    if (lg_bs <= gd->tb_range[1] &&
        lg_bs > gd->tb_range[0] &&
        cb->lg_bs - lg_bs < gd->tb_depth[0] + 1)
    {
        venc_trace_syntax_element("split_transform_flag");
        venc_encode_context_bin(cbs, F265_CO_SPLIT_TRANSFORM+5-lg_bs, split_flag);
    }

    if (split_flag)
        for (int part = 0, sub_lg_bs = lg_bs-1; part < 4; part++)
            venc_hm_an_intra_hv_bits(t, cb, sub_lg_bs, mode, cb_ox+((part&1)<<sub_lg_bs), cb_oy+((part>>1)<<sub_lg_bs),
                                     depth+1);

    else
    {
        int nz_flag = tn_val&1;

        // FIXME. There should be an finline function in entropy.c to signal
        //        the cbf_luma syntax element.
        int ctx_off = cb->lg_bs == lg_bs;
        venc_trace_syntax_element("cbf_luma");
        venc_encode_context_bin(cbs, F265_CO_CBF_LUMA+ctx_off, nz_flag);

        // DCT/Quant/Dequant/IDCT + Write TB.
        if (nz_flag)
        {
            int ct_ox = cb->cb_off[0] + cb_ox, ct_oy = cb->cb_off[1] + cb_oy;
            venc_set_tmp_tb(t);
            venc_rec_intra_tb(t, 0, lg_bs, mode, 0, ct_ox, ct_oy, depth);
            venc_set_tmp_tb(t);
            venc_write_tb(cbs, &t->tt, 0);
        }
    }
}

// Perform intra HV analysis.
static finline void venc_hm_an_intra_hv(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;
    f265_cabac_snapshot cs;
    venc_save_cabac_contexts(cbs, &cs);

    // Signal split TT.
    *t->tt.tn++ = 8;

    // Signal part_mode is NxN.
    // NOTE The bits to signal this info are only part of the 1st NxN PB
    //      bit count as the loop below returns to the bit count prior
    //      to signaling this info (mimic HM). We could wait until after
    //      the independent analysis to signal this.
    cb->intra_luma_mode[1] = 0;
    venc_write_intra_part_mode(t, cb);

    // Run analysis on each HV part independently. Reset the CABAC contexts
    // in between each part. Mimic HM.
    uint8_t *tn = t->tt.tn;
    for (int part = 0; part < 4; part++)
    {
        *dist += venc_hm_an_intra_luma_mode(t, cb, 1, part);

        // Reload contexts between subparts (mimic HM).
        venc_load_cabac_contexts(cbs, &cs);

        // Reset bits for next iteration.
        cbs->rdo.bits = bits;

        // Update prediction map.
        venc_update_pmap_cb(t, cb);
    }

    // Revert TT.
    t->tt.tn = tn;

    // Signal part_mode and luma prediction modes.
    venc_write_intra_part_mode(t, cb);
    for (int part = 0; part < 4; part++) venc_hm_an_luma_pred_cost(t, cb, part);

    // Encode each part using the selected best mode, correctly updating the
    // CABAC contexts along the way.
    for (int part = 0, lg_bs = cb->lg_bs-1; part < 4; part++)
        venc_hm_an_intra_hv_bits(t, cb, lg_bs, cb->intra_luma_mode[part], (part&1)<<lg_bs, (part>>1)<<lg_bs, 0);

    // Optimize chroma prediction mode.
    // NOTE The HM resets the bits for chroma exploration.
    //      I'm not sure this really changes anything as it only
    //      affects the base (reference). All modes would start
    //      with the same count.
    uint64_t luma_bits = cbs->rdo.bits;
    cbs->rdo.bits = bits;
    t->tt.tn = tn-1;
    *dist += venc_hm_an_intra_chroma_mode(t, cb);
    cbs->rdo.bits = (cbs->rdo.bits-bits) + luma_bits;

    // Cost for best luma and best chroma.
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Best Intra HV", *dist);
    venc_trace_rdo_bits("Best Intra HV", cbs->rdo.bits);
    venc_trace_rdo_cost("Best Intra HV", *cost);
    venc_trace_rdo_luma_modes(cb);
}

// Perform intra UN analysis.
static finline void venc_hm_an_intra_un(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    // Save input state.
    uint8_t *tn = t->tt.tn;
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;

    // FIXME Signal PCM is not used (when PCM is enabled).

    // Optimize intra angular mode.
    cb->intra_luma_mode[1] = -1;
    venc_write_intra_part_mode(t, cb);
    *dist = venc_hm_an_intra_luma_mode(t, cb, 0, 0);

    // Optimize chroma prediction mode.
    // NOTE The HM resets the bits for chroma exploration.
    //      I'm not sure this really changes anything as it only
    //      affects the base (reference). All modes would start
    //      with the same count.
    uint64_t luma_bits = cbs->rdo.bits;
    cbs->rdo.bits = bits;
    t->tt.tn = tn;
    *dist += venc_hm_an_intra_chroma_mode(t, cb);
    cbs->rdo.bits = cbs->rdo.bits - bits + luma_bits;

    // Cost for best luma and best chroma.
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Best Intra UN", *dist);
    venc_trace_rdo_bits("Best Intra UN", cbs->rdo.bits);
    venc_trace_rdo_cost("Best Intra UN", *cost);
    venc_trace_rdo_luma_modes(cb);
}

// Perform intra RDO analysis.
static int32_t venc_hm_an_intra(f265_enc_thread *t, f265_cb *cb)
{
    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    // Save transform tree pointer offset and CABAC contexts.
    venc_stash_init_cb(t, cb);

    // Remember input bit count.
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;

    // Intra UN.
    venc_stash_push(t);
    venc_hm_an_intra_un(t, cb, &dist[0], &cost[0]);
    venc_stash_pop(t);

    // Intra HV.
    if (cb->lg_bs == t->enc->gd.cb_range[0])
    {
        // Save current state so we may revert.
        int un_luma_mode = cb->intra_luma_mode[0], un_chroma_mode = cb->intra_chroma_mode;
        venc_stash_save_reset(t);

        // Remember intra UN bit count. Reload input bit count.
        F265_SWAP(uint64_t, cbs->rdo.bits, bits);

        venc_stash_push(t);
        venc_hm_an_intra_hv(t, cb, &dist[1], &cost[1]);
        venc_stash_pop(t);

        // Reload intra UN results if intra HV doesn't perform better.
        if (cost[0] <= cost[1])
        {
            // Restore transform tree, reconstructed samples, CABAC contexts.
            venc_stash_restore(t);

            // Reinstate bit count after intra UN exploration.
            cbs->rdo.bits = bits;

            // Reload best modes.
            cb->intra_luma_mode[0] = un_luma_mode;
            cb->intra_luma_mode[1] = -1;
            cb->intra_chroma_mode = un_chroma_mode;
        }
        else
        {
            dist[0] = dist[1];
            cost[0] = cost[1];
        }
    }

    // FIXME Intra PCM.

    return dist[0];
}

// Setup intra RDO analysis. Signal prediction mode, fetch distortion from
// intra RDO analysis and compute RDO cost.
static finline void venc_hm_an_pred_intra(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;
    F265_SET_FLAG(cb->flags, F265_CB_INTRA, 1);
    venc_write_pred_mode_flag(t, cb, 1);
    *dist = venc_hm_an_intra(t, cb);
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Best Intra", *dist);
    venc_trace_rdo_bits("Best Intra", cbs->rdo.bits);
    venc_trace_rdo_cost("Best Intra", *cost);
}

// Recursively analyze the inter transform tree.
// Figure out the best partitioning.
// FIXME: set static when called.
int32_t venc_hm_an_inter_tt(f265_enc_thread *t, f265_cb *cb)
{
    // TODO.
    return 0;
}

// Run RDO analysis for the inter mode set in the CB.
static finline void venc_hm_an_inter_mode(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
#ifdef VAN_TRACE_ANALYSIS
    char *str;
    int inter_part = cb->inter_part;
    if (inter_part == F265_PART_UN) str = "2Nx2N";
    else if (inter_part == F265_PART_H2) str = "2NxN";
    else if (inter_part == F265_PART_V2) str = "Nx2N";
    else if (inter_part == F265_PART_H1) str = "2NxnU";
    else if (inter_part == F265_PART_H3) str = "2NxnD";
    else if (inter_part == F265_PART_V1) str = "nLx2N";
    else if (inter_part == F265_PART_V3) str = "nRx2N";

    f265_cabac_bs *cbs = &t->cbs;
    venc_trace_rdo_dist(str, *dist);
    venc_trace_rdo_bits(str, cbs->rdo.bits);
    venc_trace_rdo_cost(str, *cost);
#endif
}

// Inter skip RDO analysis.
static finline void venc_hm_an_inter_skip(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;
    venc_trace_rdo_dist("Best skip", *dist);
    venc_trace_rdo_bits("Best skip", cbs->rdo.bits);
    venc_trace_rdo_cost("Best skip", *cost);
}

// FIXME.
// Perform inter RDO analysis.
static int32_t venc_hm_an_inter(f265_enc_thread *t, f265_cb *cb)
{
    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    uint8_t *tn = t->tt.tn;
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;
    f265_cabac_snapshot cs;
    venc_save_cabac_contexts(cbs, &cs);

    int best_part = F265_PART_UN;
    venc_hm_an_inter_skip(t, cb, &dist[0], &cost[0]);
    for (int mode = F265_PART_UN; mode <= F265_PART_V3; mode++)
    {
        t->tt.tn = tn;
        venc_load_cabac_contexts(cbs, &cs);
        cbs->rdo.bits = bits;
        cb->inter_part = mode;
        venc_hm_an_inter_mode(t, cb, &dist[1], &cost[1]);
        if (cost[1] < cost[0])
        {
            cost[0] = cost[1];
            dist[0] = dist[1];
            best_part = mode;

            // FIXME Save the info associated to the best mode.
        }
    }

    if (cb->inter_part != best_part)
    {
        cb->inter_part = best_part;

        // FIXME Reload the info associated to the best mode.
    }

    return dist[0];
}

// Call inter RDO analysis.
static finline void venc_hm_an_pred_inter(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;
    F265_SET_FLAG(cb->flags, F265_CB_INTRA, 0);
    venc_write_pred_mode_flag(t, cb, 0);
    *dist = venc_hm_an_inter(t, cb);
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Best Inter", *dist);
    venc_trace_rdo_bits("Best Inter", cbs->rdo.bits);
    venc_trace_rdo_cost("Best Inter", *cost);
}

// Choose between intra and inter prediction based on RDO analysis and
// available modes.
static int32_t venc_hm_an_pred(f265_enc_thread *t, f265_cb *cb)
{
    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;

    // FIXME Signal transquant bypass (when the option is enabled).

    // FIXME Signal skip flag (when t->src_frame->frame_type != F265_FRAME_I).

    // Save input state (transform tree pointer offset, CABAC contexts).
    venc_stash_init_cb(t, cb);

    f265_frame *f = t->src_frame;
    if (f->frame_type != F265_FRAME_I)
    {
        // Inter analysis.
        venc_stash_push(t);
        venc_hm_an_pred_inter(t, cb, &dist[0], &cost[0]);
        venc_stash_pop(t);

        // Save inter transform tree, reconstructed samples, CABAC contexts.
        // Revert to input state for intra exploration.
        venc_stash_save_reset(t);

        // Remember bit count after inter exploration. Reload input bit count.
        F265_SWAP(uint64_t, cbs->rdo.bits, bits);
    }

    // Intra analysis.
    venc_stash_push(t);
    venc_hm_an_pred_intra(t, cb, &dist[1], &cost[1]);
    venc_stash_pop(t);

    // Reload inter results if intra doesn't perform better.
    int inter_flag = cost[0] <= cost[1];
    if (inter_flag)
    {
        // Restore transform tree, reconstructed samples, CABAC contexts from
        // inter exploration.
        venc_stash_restore(t);

        // Restore bit count following inter exploration.
        cbs->rdo.bits = bits;

        // Reinstate inter flag.
        F265_SET_FLAG(cb->flags, F265_CB_INTRA, 0);
    }

    return dist[!inter_flag];
}

// Perform RDO analysis on the children of current the CB.
static finline void venc_hm_an_split_cb(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;

    // Explore children. CABAC state is continued from one child to the next.
    uint64_t bits = 0;
    F265_SET_FLAG(cb->flags, F265_CB_SPLIT, 1);
    f265_cb *child_cb = t->cb + cb->child_idx;
    for (int part = 0; part < 4; part++)
    {
        *dist += venc_hm_an_cb(t, child_cb++);

        // Reset bits in between each child (mimic HM).
        bits += cbs->rdo.bits>>15;
        cbs->rdo.bits &= 32767;
    }
    cbs->rdo.bits += bits<<15;

    venc_write_split_cu_flag(t, cb, 1);
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Split", *dist);
    venc_trace_rdo_bits("Split", cbs->rdo.bits&32767);
    venc_trace_rdo_cost("Split", *cost);
}

// Perform RDO analysis on current the CB.
static finline void venc_hm_an_unsplit_cb(f265_enc_thread *t, f265_cb *cb, int32_t *dist, double *cost)
{
    f265_cabac_bs *cbs = &t->cbs;
    F265_SET_FLAG(cb->flags, F265_CB_SPLIT, 0);
    *dist = venc_hm_an_pred(t, cb);
    venc_write_split_cu_flag(t, cb, 0);
    *cost = venc_hm_an_cost(t, cbs, *dist);

    venc_trace_rdo_dist("Unsplit", *dist);
    venc_trace_rdo_bits("Unsplit", cbs->rdo.bits&32767);
    venc_trace_rdo_cost("Unsplit", *cost);
}

// Recursively run analysis on CB and its children.
// Don't analyze missing CBs from partial CTBs.
static int32_t venc_hm_an_cb(f265_enc_thread *t, f265_cb *cb)
{
    if (!F265_GET_FLAG(cb->flags, F265_CB_PRESENT)) return 0;

    venc_trace_rdo_path_cb(1<<cb->lg_bs);

    double cost[2] = { F265_MAX_RD_COST, F265_MAX_RD_COST };
    int32_t dist[2] = { 0, 0 };

    // Save the block position and size.
    // Save the position of the tmap pointer.
    // Save the CABAC contexts.
    venc_stash_init_cb(t, cb);

    // Remember the current bit count.
    f265_cabac_bs *cbs = &t->cbs;
    uint64_t bits = cbs->rdo.bits;

    // RDO analysis on unsplit CB.
    f265_gen_data *gd = &t->enc->gd;
    if (!F265_GET_FLAG(cb->flags, F265_CB_FORBIDDEN))
    {
        venc_stash_push(t);
        venc_hm_an_unsplit_cb(t, cb, &dist[0], &cost[0]);
        venc_stash_pop(t);
    }

    // RDO analysis on split CB.
    if (cb->lg_bs > gd->cb_range[0])
    {
        // Save unsplit results (when present).
        if (cost[0] != F265_MAX_RD_COST)
        {
            // Save the transform tree state.
            // Copy the reconstructred samples (Y, Cb and Cr) to the stash.
            // Save the CABAC contexts following unsplit exploration.
            // Reposition the tmap pointer at the location prior to unsplit
            // exploration.
            // Reload the CABAC contexts to the state prior to unsplit
            // exploration.
            venc_stash_save_reset(t);

            // Remember bit count after unsplit exploration.
            // Reload bit count prior to unsplit exploration.
            F265_SWAP(uint64_t, cbs->rdo.bits, bits);
        }

        venc_stash_push(t);
        venc_hm_an_split_cb(t, cb, &dist[1], &cost[1]);
        venc_stash_pop(t);

        // Reload unsplit results if splitting doesn't perform better.
        if (cost[0] <= cost[1])
        {
            // Reload the transform tree, the reconstructed samples and the
            // CABAC contexts.
            venc_stash_restore(t);

            // Reload bit count following unsplit exploration.
            cbs->rdo.bits = bits;
        }
    }

    // Update CB hierarchy according to best option.
    int split_flag = cost[1] < cost[0];
    F265_SET_FLAG(cb->flags, F265_CB_SPLIT, split_flag);
    if (!split_flag) venc_update_pmap_unsplit_cb(t, cb);

    venc_trace_rdo_dist("Best CB", dist[split_flag]);
    venc_trace_rdo_bits("Best CB", t->cbs.rdo.bits&32767);
    venc_trace_rdo_cost("Best CB", cost[split_flag]);

    return dist[split_flag];
}

// Perform RDO analysis on the current CTB, updating the CB flags and the
// transform tree. CABAC engine is restored to calling state prior to return.
int venc_hm_an_ctb(f265_enc_thread *t)
{
    #ifdef VAN_TRACE_ANALYSIS
    venc_print_hm_rdo_trace = 1;
    #endif

    // Remember pre-analysis CABAC state.
    f265_cabac_bs *cbs = &t->cbs;
    f265_cabac_bs snapshot = *cbs;

    // Run analysis.
    cbs->rdo.bits = 0;
    t->stash.flags = 31;
    venc_init_transform_tree(t);
    venc_hm_an_cb(t, t->cb);

    // Restore pre-analysis CABAC state.
    t->cbs = snapshot;

    #ifdef VAN_TRACE_ANALYSIS
    venc_print_hm_rdo_trace = 0;
    #endif

    return 1;
}
