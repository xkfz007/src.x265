// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include "enc.h"

// TODO: unify the chroma scaling factors correctly.

// Define to trace motion vector costs during analysis.
//#define VAN_TRACE_MV_COST
#ifdef VAN_TRACE_MV_COST
extern int venc_trace_analysis_flag;
#endif


void venc_me_interpol_plane(int16_t *dst, int dst_stride, f265_mv mv, f265_me_ctx *me, int csf[2],
                            f265_pix *ref_plane, int comp, int width, int height)
{
    // Get the integer part of the MV.
    int xmv = mv.x >> (2 + csf[0]);
    int ymv = mv.y >> (2 + csf[1]);

    // Get the fractional part of the MV.
    int xfrac = mv.x & ((1<<(csf[0]+2))-1);
    int yfrac = mv.y & ((1<<(csf[1]+2))-1);

    // Get a pointer to the top left pixel of the block.
    f265_pix *ref = ref_plane + xmv + ymv * me->ref_stride;

    if (comp)
        venc_interpol_chroma(ref, me->ref_stride, dst, dst_stride, width, height, xfrac, yfrac, me->bit_depth[1]);
    else
        venc_interpol_luma(ref, me->ref_stride, dst, dst_stride, width, height, xfrac, yfrac, me->bit_depth[0]);
}

// Interpolate the reference pixels at the MV position.
void venc_me_interpol(f265_pix *dst, int dst_stride, f265_mv mv, f265_me_ctx *me, int comp)
{
    int is_chroma = !!comp;

    // Find the chroma scale factor.
    int csf[2];
    for (int i = 0; i < 2; i++)
        csf[i] = me->csf[i] * is_chroma;

    // Find the block actual dimension and stride.
    int width  = me->dim[0] >> csf[0];
    int height = me->dim[1] >> csf[1];

    int16_t buf[64*64];
    venc_me_interpol_plane(buf, 64, mv, me, csf, me->ref[0].ref_planes[comp], comp, width, height);

    // Weight the prediction.
    f265_weight *weights = me->ref[0].weights + comp;
    if (weights->used_flag) venc_weight_pix_uni(buf, 64, dst, dst_stride, width, height, weights, me->bit_depth[is_chroma]);
    else venc_scale_pix_uni(buf, 64, dst, dst_stride, width, height, me->bit_depth[is_chroma]);
}

// Interpolate the reference pixels at the MV position. Each MV and me_ctx represents a reference.
void venc_me_interpol_bi(f265_pix *dst, int dst_stride, f265_mv mv[2], f265_me_ctx *me, int comp)
{
    int is_chroma = !!comp;

    // Find the chroma scale factor.
    int csf[2];
    for (int i = 0; i < 2; i++)
        csf[i] = me->csf[i] * is_chroma;

    // Find the block actual dimension and stride.
    int width  = me->dim[0] >> csf[0];
    int height = me->dim[1] >> csf[1];

    int16_t buf[2][64*64];
    for (int i = 0; i < 2; i++)
        venc_me_interpol_plane(buf[i], 64, mv[i], me, csf, me->ref[i].ref_planes[comp], comp, width, height);

    // Weight the prediction.
    f265_weight *weights_l0 = me->ref[0].weights + comp;
    f265_weight *weights_l1 = me->ref[1].weights + comp;
    if (weights_l0->used_flag || weights_l1->used_flag)
        venc_weight_pix_bi(buf[0], buf[1], 64, dst, dst_stride, width, height, weights_l0, weights_l1, me->bit_depth[is_chroma]);
    else
        venc_scale_pix_bi(buf[0], buf[1], 64, dst, dst_stride, width, height, me->bit_depth[is_chroma]);
}

// Return a pointer to the fake luma reference pixels of the current
// uni-predicted block. The function sets "out_stride" to the stride of the
// reference block. "tmp_buf" is used to store the interpolation result as
// needed.
static f265_pix* venc_fake_luma_ref_p(f265_enc_thread *t, int *out_stride, f265_pix *tmp_buf, int mx, int my)
{
    f265_me_ctx *me = &t->me;
    f265_pix **ref_planes = me->ref[0].planes;
    int weight_flag = me->ref[0].weights[0].used_flag;
    int plane_stride = t->plane_stride;
    int packed_dims = me->packed_dims[0];
    int width = (packed_dims>>8)&0xff, awi = packed_dims>>24;
    f265_pix *out = tmp_buf;

    // Get the reference plane information.
    int ref_off = (my>>2)*plane_stride + (mx>>2);
    int hpel_idx = ((my&3)<<2) + (mx&3);

    // Get the first halfpel plane.
    f265_pix *ref0 = ref_planes[f265_hpel_src0[hpel_idx]] + ref_off + (((my&3) == 3) ? plane_stride : 0);

    // Average the two closest halfpel reference planes.
    if (likely(hpel_idx&5))
    {
        *out_stride = width;
        f265_pix *ref1 = ref_planes[f265_hpel_src1[hpel_idx]] + ref_off + ((mx&3) == 3);
        venc_avg_pix[awi](out, ref0, plane_stride, ref1, plane_stride, packed_dims);

        // Weight in place. FIXME.
        if (weight_flag)
        {
        }
    }

    // Use only the first halfpel plane.
    else
    {
        // Weight. FIXME.
        if (weight_flag)
        {
        }

        // Use the halfpel plane directly.
        else
        {
            out = ref0;
            *out_stride = plane_stride;
        }
    }

    return out;
}

// Return the fake luma distortion cost for the specified motion vector.
static int venc_get_fake_luma_block_dist(f265_enc_thread *t, f265_mv mv)
{
    f265_me_ctx *me = &t->me;
    int packed_dims = me->packed_dims[0];
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff, awi = packed_dims>>24;
    int aligned_block_size = F265_ALIGN_VAL(width*height, 64);
    int alloc_size = aligned_block_size*sizeof(f265_pix);
    f265_pix *tmp_buf = (f265_pix*)t->store;

    t->store += alloc_size;
    int ref_stride;
    f265_pix *ref = venc_fake_luma_ref_p(t, &ref_stride, tmp_buf, mv.x, mv.y);
    int dist = venc_fsad[awi](me->src_planes[0], t->plane_stride, ref, ref_stride, packed_dims);
    t->store -= alloc_size;

    return dist;
}

// Return the distortion using the current distortion metric.
int venc_me_get_dist(f265_me_ctx *me, f265_pix *src0, int32_t stride0, f265_pix *src1,
                     int32_t stride1, int32_t width, int32_t height, int32_t bitdepth)
{
    return (me->dist[me->dist_func_id])(src0, stride0, src1, stride1, width, height, bitdepth);
}

// Return the distortion cost of a MV on a reference luma plane.
int venc_me_luma_cost(f265_mv mv, f265_me_ctx *me)
{
    // General case.
    if (mv.x&3 || mv.y&3 || me->ref[0].weights->used_flag)
    {
        f265_pix buf[64*64];
        venc_me_interpol(buf, 64, mv, me, 0);
        return venc_me_get_dist(me, me->src_planes[0], me->ref_stride, buf, 64, me->dim[0], me->dim[1],
                                me->bit_depth[0]);
    }

    // Unweighted fullpel case.
    else
    {
        f265_pix *ref = me->ref[0].ref_planes[0] + (mv.y>>2)*me->ref_stride + (mv.x>>2);
        return venc_me_get_dist(me, me->src_planes[0], me->ref_stride, ref, me->ref_stride, me->dim[0], me->dim[1],
                                me->bit_depth[0]);
    }
}

// Return the distortion cost of a MV on a reference chroma plane.
int venc_me_chroma_cost(f265_mv mv, f265_me_ctx *me, int comp)
{
    f265_pix buf[64*64];
    venc_me_interpol(buf, 64, mv, me, comp);
    return venc_me_get_dist(me, me->src_planes[comp], me->ref_stride, buf, 64,
                            me->dim[0]>>me->csf[0], me->dim[1]>>me->csf[1], me->bit_depth[1]);
}

// Return the distortion of the merge candidate specified.
int venc_me_merge_cand_dist(f265_enc_thread *t, f265_inter_neighbour_mv cand)
{
    f265_me_ctx *me = &t->me;
    int bi_flag = cand.ref_idx[0] != -1 && cand.ref_idx[1] != -1;
    int uni_list = cand.ref_idx[0] == -1;
    int uni_ref_idx = cand.ref_idx[uni_list];
    int dist = 0;

    // Set the references and clip the MVs in place.
    if (bi_flag)
        for (int list = 0; list < 2; list++)
        {
            venc_me_set_ref(t, t->ref_ctx[list] + cand.ref_idx[list], list);
            f265_clip_mv(cand.mv[list], t->mc_bounds64);
        }

    else
    {
        venc_me_set_ref(t, t->ref_ctx[uni_list] + uni_ref_idx, 0);
        f265_clip_mv(cand.mv[uni_list], t->mc_bounds64);
    }

    // Fast track.
    if ((t->enc->gd.algo&0x400) && !bi_flag)
        return venc_get_fake_luma_block_dist(t, cand.mv[uni_list]);

    // Pass each component.
    int nb_comp = 1 + (me->chroma_flag<<1);
    for (int comp = 0; comp < nb_comp; comp++)
    {
        // Do the interpolation.
        f265_pix buf[64*64];
        if (bi_flag) venc_me_interpol_bi(buf, 64, cand.mv, me, comp);
        else venc_me_interpol(buf, 64, cand.mv[uni_list], me, comp);

        // Compute the distortion.
        int is_chroma = !!comp;
        int scale_x = me->csf[0] * is_chroma;
        int scale_y = me->csf[1] * is_chroma;
        dist += (me->dist[me->dist_func_id])(me->src_planes[comp], me->ref_stride, buf, 64,
                                             me->dim[0]>>scale_x, me->dim[1]>>scale_y, me->bit_depth[0]);
    }

    return dist;
}

// Return the cost of the MV (classic algorithm).
int venc_me_mv_cost(f265_mv mv, uint16_t *cost_table[2])
{
    return cost_table[0][mv.x] + cost_table[1][mv.y];
}

// Compute the MV length like HM.
static int32_t venc_me_mv_hm_len(int32_t val)
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

// Compute the MV length exactly like the spec.
static int32_t venc_me_mv_spec_len(int val)
{
    // Adapted from venc_egk_bin().
    int k = 1;
    int o = val + (1<<k);
    int n = 31 - __builtin_clz(o);
    return 1 + (n<<1) - k;
}

// Experimental MV cost computations. FIXME. Leaving debug for now.
//
// HM computes the motion vector costs as follow. The probabilities for the
// greater0/1 flags are assumed to be 50%. Thus each flag costs 1 bit
// independently of its value. The sign bin and the exponential golomb order 1
// bins cost 1 bit each (bypass). Example table:
//
//   mvd 0 => 1 (gt0).
//   mvd 1 => 3 (gt0, gt1, sign).
//   mvd 2 => 5 (gt0, gt1, sign, eg1:00).
//   mvd 3 => 5 (gt0, gt1, sign, eg1:01).
//   mvd 4 => 7 (gt0, gt1, sign, eg1:1000).
//   mvd 5 => 7 (gt0, gt1, sign, eg1:1001).
//   mvd 6 => 7 (gt0, gt1, sign, eg1:1010).
//   mvd 7 => 7 (gt0, gt1, sign, eg1:1011).
//
// It might be possible to improve on that algorithm by using better
// approximations for the greater0/1 flag costs.
//
// Experimentally, using the current CABAC contexts for computing the motion
// vector costs in RDO mode is not always productive. Perhaps this nullifies
// the negative feedback loop that steers the CABAC contexts toward the
// statistical model average which tend to be globally optimal. Consequently
// the CABAC contexts get stuck in a local minimum and the quality suffers.
//
// Using the initial CABAC contexts is another alternative. Theorically they
// reflect the statistical average.
//
// In RDM mode, both approaches seems to help a little but there are no
// observed benefits in RDO mode. To be revisited.

int venc_me_mv_cost_test(f265_me_ctx *me, f265_mv mv, int ref_id)
{
    int algo = 1;
    f265_enc_thread *t = me->t;
    t = t;

    // H.264 style for CABAC.
    if (algo == 0)
        return venc_me_mv_cost(mv, me->ref[ref_id].mv_costs);

    // HM-style.
    else if (algo == 1)
    {
        f265_mv pmv = me->ref[ref_id].pmv;
        int sum = 0;

        #if 0
        if (venc_trace_analysis_flag) printf("===> MV (%d,%d) PMV (%d,%d) lambda %d.\n",
                                             mv.x, mv.y, pmv.x, pmv.y, t->me.lambda);
        #endif

        for (int comp = 0; comp < 2; comp++)
        {
            int bits = venc_me_mv_hm_len(mv.v[comp] - pmv.v[comp]);
            sum += bits;

            #if 0
            int mvd = F265_ABS(mv.v[comp] - pmv.v[comp]);
            if (venc_trace_analysis_flag) printf("Comp %d, mvd %d, cost %d.\n", comp, mvd, bits);
            #endif
        }

        int cost = ((int64_t)sum*me->lambda)>>8;
        #if 0
        if (venc_trace_analysis_flag) printf("Total cost %d.\n", cost);
        #endif

        return cost;
    }

    // CABAC-aware.
    else
    {
        // Get the context costs. The first index is the context offset, the second
        // index is the bin value.
        int ctx_costs[2][2];

        // Dynamic CABAC contexts.
        #if 1
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                int ctx_idx = F265_CO_MVD_GREATER0+i;
                int ctx_val = t->cbs.contexts[ctx_idx];
                ctx_costs[i][j] = f265_cabac_entropy_table[ctx_val^j];
            }
        #endif

        // HM-like CABAC contexts.
        #if 0
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                ctx_costs[i][j] = 32768;
        #endif

        // Chosen CABAC contexts. Hardcoded for QP 30 P frame.
        #if 0
        ctx_costs[0][0] = 32768*1.45114;
        ctx_costs[0][1] = 32768*0.65683;
        ctx_costs[1][0] = 32768*0.80499;
        ctx_costs[1][1] = 32768*1.22556;
        #endif

        // Get the costs for each component.
        f265_mv pmv = me->ref[ref_id].pmv;
        int sum_ctx_cost = 0, sum_bypass_count = 0;
        #if 0
        if (venc_trace_analysis_flag) printf("===> MV (%d,%d) PMV (%d,%d) lambda %d ctx %.5f/%.5f | %.5f/%.5f.\n",
                                             mv.x, mv.y, pmv.x, pmv.y, t->me.lambda,
                                             ctx_costs[0][0]/32768.0, ctx_costs[0][1]/32768.0,
                                             ctx_costs[1][0]/32768.0, ctx_costs[1][1]/32768.0);
        #endif
        for (int comp = 0; comp < 2; comp++)
        {
            int mvd = F265_ABS(mv.v[comp] - pmv.v[comp]);

            // Get the context cost and the bypass bin count.
            int ctx_cost, bypass_count;
            if (mvd == 0)
            {
                ctx_cost = ctx_costs[0][0];
                bypass_count = 0;
            }

            else if (mvd == 1)
            {
                ctx_cost = ctx_costs[0][1] + ctx_costs[1][0];
                bypass_count = 1;
            }

            else
            {
                ctx_cost = ctx_costs[0][1] + ctx_costs[1][1];
                bypass_count = venc_me_mv_spec_len(mvd-2) + 1;
            }

            #if 0
            if (venc_trace_analysis_flag) printf("Comp %d, mvd %d, ctx_cost %.2f, bypass %d.\n",
                                                 comp, mvd, ctx_cost/32768.0, bypass_count);
            #endif

            sum_ctx_cost += ctx_cost;
            sum_bypass_count += bypass_count;
        }

        // Scale and multiply by lambda.
        int shift = 15+8;
        int cost = ((int64_t)(sum_ctx_cost + (sum_bypass_count<<15))*me->lambda + (1<<(shift-1)))>>shift;

        #if 0
        if (venc_trace_analysis_flag) printf("Total cost %d.\n", cost);
        #endif

        return cost;
    }
}

// Return the total cost of the MV. Include all required plane distortions, as set with the me->chroma_flag.
int venc_me_mv_total_cost(f265_mv mv, f265_me_ctx *me)
{
    int mv_cost = venc_me_mv_cost_test(me, mv, 0);
    int dist = venc_me_luma_cost(mv, me);
    if (me->chroma_flag)
        for (int comp = 1; comp < 3; comp++)
            dist += venc_me_chroma_cost(mv, me, comp);
    #ifdef VAN_TRACE_MV_COST
    if (venc_trace_analysis_flag)
        printf("Search (%d,%d): cost %d, MV %d, dist %d.\n", mv.x, mv.y, mv_cost + dist, mv_cost, dist);
    #endif
    return mv_cost + dist;
}

// Return the total cost of a bi-predicted MV. Include all required plane distortions, as set with the me->chroma_flag.
int venc_me_mv_total_cost_bi(f265_mv mv[2], f265_me_ctx *me)
{
    int cost = 0;

    // Get the cost of encoding the MV.
    for (int i = 0; i < 2; i++)
        cost += venc_me_mv_cost_test(me, mv[i], i);

    f265_pix dst[64*64];

    // Get the distortion of each plane.
    for (int i = 0; i < 1+2*me->chroma_flag; i++)
    {
        // Get the interpolated plane.
        venc_me_interpol_bi(dst, 64, mv, me, i);

        // Get the distortion.
        int is_chroma = !!i;
        int scale_x = me->csf[0] * is_chroma;
        int scale_y = me->csf[1] * is_chroma;
        cost += (*me->dist[me->dist_func_id])(me->src_planes[i], me->ref_stride,
                                              dst, 64,
                                              me->dim[0]>>scale_x, me->dim[1]>>scale_y,
                                              me->bit_depth[is_chroma]);
    }

    return cost;
}

// Return the index of the PMV with the lowest distortion. We ignore the PMV
// index cost. This function does not save the best PMV in the ME context. The
// distortion cost is set.
int venc_me_test_pmv(f265_mv pmv[2], f265_me_ctx *me, int dist, int *cost)
{
    f265_enc_thread *t = me->t;

    // Set the required distortion function.
    me->dist_func_id = dist;

    // Find the cost of both candidates. We clip to the motion estimation
    // bounds, not the motion compensation bounds here.
    int costs[2];

    if (t->enc->gd.algo&0x400)
    {
        for (int i = 0; i < 2; i++)
            costs[i] = venc_get_fake_luma_block_dist(t, f265_clip_mv(pmv[i], me->me_bounds64));
    }

    else
    {
        for (int i = 0; i < 2; i++)
            costs[i] = venc_me_luma_cost(f265_clip_mv(pmv[i], me->me_bounds64), me);
    }

    // Return the index of the best PMV. Choose PMV 0 in priority.
    int idx = costs[0] > costs[1];
    *cost = costs[idx];
    return idx;
}

// Return the index of the nearest PMV based on the MV cost and the PMV index
// cost. Update the best cost in the MV context, but do not save the best PMV.
// Note that the best cost now includes the base PMV cost.
int venc_me_find_nearest_pmv(f265_mv pmv[2], uint16_t base_costs[2], f265_me_ctx *me)
{
    // Get the distortion.
    f265_me_ref *ref = me->ref;
    f265_mv mv = ref->best_mv;

    // FIXME. Do not use the table for the cost computations until we've
    // finished experimenting with the MV costs.
    #if 0
    int dist = me->best_cost - ref->mv_costs[0][mv.x] - ref->mv_costs[1][mv.y];

    // Get the costs for each PMV.
    uint16_t *cost_table = me->mv_cost_table + (F265_NB_MV_COSTS >> 1);
    int pmv_costs[2];
    for (int i = 0; i < 2; i++)
    {
        int mv_cost = cost_table[mv.x - pmv[i].x] + cost_table[mv.y - pmv[i].y];
        pmv_costs[i] = base_costs[i] + mv_cost;
        #ifdef VAN_TRACE_MV_COST
        if (venc_trace_analysis_flag)
            printf("PMV index %d PMV (%d,%d) MV (%d,%d) base/MV/dist/total %d/%d/%d/%d.\n",
                   i, pmv[i].x, pmv[i].y, mv.x, mv.y, base_costs[i], mv_cost, dist, pmv_costs[i] + dist);
        #endif
    }
    #else
    int dist = me->best_cost - venc_me_mv_cost_test(me, mv, 0);
    int pmv_costs[2];
    for (int i = 0; i < 2; i++)
    {
        venc_me_set_pmv(me, pmv[i], me->t->qp[0], 0);
        int mv_cost = venc_me_mv_cost_test(me, mv, 0);
        pmv_costs[i] = base_costs[i] + mv_cost;
        #ifdef VAN_TRACE_MV_COST
        if (venc_trace_analysis_flag)
            printf("PMV index %d PMV (%d,%d) MV (%d,%d) base/MV/dist/total %d/%d/%d/%d.\n",
                   i, pmv[i].x, pmv[i].y, mv.x, mv.y, base_costs[i], mv_cost, dist, pmv_costs[i] + dist);
        #endif
    }
    #endif

    // Choose PMV 0 in priority.
    int pmv_idx = pmv_costs[0] > pmv_costs[1];

    // Update the best cost.
    me->best_cost = dist + pmv_costs[pmv_idx];

    return pmv_idx;
}

// Test and save the MV if it has a better cost.
void venc_me_keep_best_mv(f265_mv cand_mv, f265_me_ctx *me)
{
    int cand_cost = venc_me_mv_total_cost(cand_mv, me);
    if (cand_cost < me->best_cost)
    {
        #ifdef VAN_TRACE_MV_COST
        if (venc_trace_analysis_flag)
            printf("Kept best MV (%d,%d) cost %d.\n", cand_mv.x, cand_mv.y, cand_cost);
        #endif
        me->best_cost = cand_cost;
        me->ref[0].best_mv = cand_mv;
    }
}

// Test all offsets centered on the MV.
void venc_me_test_all_positions(f265_mv base_mv, f265_me_ctx *me, const int8_t offset[][2], int nb_offset)
{
    for (int pos = 0; pos < nb_offset; pos++)
    {
        f265_mv cand_mv = base_mv;
        cand_mv.x += offset[pos][0];
        cand_mv.y += offset[pos][1];
        venc_me_keep_best_mv(cand_mv, me);
    }
}

// Dia search.
void venc_me_dia_search(int32_t nb_iter, int32_t scale, f265_me_ctx *me)
{
    const int8_t dia_pos[3][4][2] =
    {
        { {0, -1}, {0, 1}, { -1, 0}, {1, 0} },
        { {0, -2}, {0, 2}, { -2, 0}, {2, 0} },
        { {0, -4}, {0, 4}, { -4, 0}, {4, 0} }
    };

    for (int iter = 0; iter < nb_iter; iter++)
    {
        f265_mv ini_mv = me->ref[0].best_mv;
        venc_me_test_all_positions(ini_mv, me, dia_pos[scale], 4);

        if (ini_mv.v[0] == me->ref[0].best_mv.v[0] && ini_mv.v[1] == me->ref[0].best_mv.v[1]) break;
        if (venc_mv_out_of_range(me->ref[0].best_mv, me->me_bounds_packs)) break;
    }
}

// X-dia search.
void venc_me_xdia_search(int32_t nb_iter, int32_t scale, f265_me_ctx *me)
{
    const int8_t x_pos[3][4][2] =
    {
        { { -1, 1}, {1, 1}, { -1, -1}, {1, -1} },
        { { -2, 2}, {2, 2}, { -2, -2}, {2, -2} },
        { { -4, 4}, {4, 4}, { -4, -4}, {4, -4} }
    };

    const int8_t dia_pos[3][4][2] =
    {
        { {0, -1}, {0, 1}, { -1, 0}, {1, 0} },
        { {0, -2}, {0, 2}, { -2, 0}, {2, 0} },
        { {0, -4}, {0, 4}, { -4, 0}, {4, 0} }
    };

    f265_mv ini_mv = me->ref[0].best_mv;
    venc_me_test_all_positions(ini_mv, me, x_pos[scale], 4);

    ini_mv = me->ref[0].best_mv;
    venc_me_test_all_positions(ini_mv, me, dia_pos[scale], 4);
}

// Hex search.
void venc_me_hex_search(int32_t nb_iter, int32_t scale, f265_me_ctx *me)
{
    // The positions on each side of the array are used for unification.
    const int8_t hex_pos[3][8][2] =
    {
        { { -1, -2}, { -2, 0}, { -1, 2}, {1, 2}, {2, 0}, {1, -2}, { -1, -2}, { -2, 0} },
        { { -2, -4}, { -4, 0}, { -2, 4}, {2, 4}, {4, 0}, {2, -4}, { -2, -4}, { -4, 0} },
        { { -4, -8}, { -8, 0}, { -4, 8}, {4, 8}, {8, 0}, {4, -8}, { -4, -8}, { -8, 0} }
    };
    // Notice that the next iteration of the hexagon search overlaps with three
    // positions of the previous iteration. Those positions must be skipped.
    // Label the hexagon directions as follow: W:0, NW:1, NE:2, E:3, SE:4, SW:5.
    // Let 'D' be the direction taken by the previous iteration. The three
    // directions to search are (D-1)%6, D, (D+1)%6.
//    const uint8_t next_hex_dir[8] = { 5, 0, 1, 2, 3, 4, 5, 0 };
    const int8_t square_pos[3][8][2] =
    {
        { {0, -1}, {0, 1}, { -1, 0}, {1, 0}, { -1, -1}, { -1, 1}, {1, -1}, {1, 1} },
        { {0, -2}, {0, 2}, { -2, 0}, {2, 0}, { -2, -2}, { -2, 2}, {2, -2}, {2, 2} },
        { {0, -4}, {0, 4}, { -4, 0}, {4, 0}, { -4, -4}, { -4, 4}, {4, -4}, {4, 4} }
    };

    for (int iter = 0; iter < nb_iter; iter++)
    {
        f265_mv ini_mv = me->ref[0].best_mv;
        venc_me_test_all_positions(ini_mv, me, hex_pos[scale] + 1, 6);

        if (ini_mv.v[0] == me->ref[0].best_mv.v[0] && ini_mv.v[1] == me->ref[0].best_mv.v[1]) break;
        if (venc_mv_out_of_range(me->ref[0].best_mv, me->me_bounds_packs)) break;
        //TODO: skip the three known positions.
    }

    f265_mv ini_mv = me->ref[0].best_mv;
    venc_me_test_all_positions(ini_mv, me, square_pos[scale], 8);
}

// Get the halfpel luma reference with the extra pixels specified.
static void venc_me_get_luma_hpel_ref(f265_me_ctx *me, f265_pix *dst, f265_pix *base_ref, f265_mv mv,
                                      int extra_x, int extra_y)
{
    F265_ALIGN64 int16_t tmp[66*66];
    int dst_stride = 66;
    int width = me->dim[0]+extra_x, height = me->dim[1]+extra_y;
    int csf[2] = { 0, 0 };

    venc_me_interpol_plane(tmp, dst_stride, mv, me, csf, base_ref, 0, width, height);

    f265_weight *weights = me->ref[0].weights;
    if (weights->used_flag) venc_weight_pix_uni(tmp, dst_stride, dst, dst_stride, width, height, weights,
                                                me->bit_depth[0]);
    else venc_scale_pix_uni(tmp, dst_stride, dst, dst_stride, width, height, me->bit_depth[0]);
}

// Get the quarterpel luma reference approximation.
static void venc_me_get_luma_qpel_ref(f265_me_ctx *me, f265_pix *dst, f265_pix *base_ref, f265_mv mv)
{
    // Get the two halfpel positions used for the qpel approximation. For a
    // diagonal approximation, we always use the fullpel position as one of the
    // two positions, since it is more precise than a halfpel position.
    //
    //      FqHqF                         F H F H F
    //      qqqqq   F=fullpel
    //      HqHqH   H=halfpel             H H H H H
    //      qqqqq   q=quarterpel
    //      FqHqF                         F H F H F
    //
    //                                    H H H H H
    //
    //                                    F H F H F
    //
    // The first index is the halfpel position. The second index is the motion
    // vector component.
    const int8_t hpel_pos_table[2][7] =
    {
        { -4, -2,  0, 0, 0, 2, 4 },
        { -2, -2, -2, 0, 2, 2, 2 }
    };

    f265_mv hpel_mv[2];
    for (int pos = 0; pos < 2; pos++)
        for (int comp = 0; comp < 2; comp++)
        {
            int v = mv.v[comp];
            int sign = v < 0 ? -1 : 1;
            int vi = sign*(F265_ABS(v)&~3);
            int vf = sign*(F265_ABS(v)&3);
            hpel_mv[pos].v[comp] = vi + hpel_pos_table[pos][3+vf];
        }

    // Interpolate between the two halfpel locations (to be replaced by
    // pre-computed halfpel planes).
    F265_ALIGN64 f265_pix tmp[2][66*66];
    int dst_stride = 66;
    int width = me->dim[0], height = me->dim[1];
    for (int pos = 0; pos < 2; pos++)
        venc_me_get_luma_hpel_ref(me, tmp[pos], base_ref, hpel_mv[pos], 0, 0);

    // Get the average. This is H.264's algorithm.
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            int off = y*dst_stride + x;
            dst[off] = (tmp[0][off] + tmp[1][off] + 1)>>1;
        }
}

// Helper function for venc_me_square_search(). Test the luma reference
// specified and update the best cost as needed.
static void venc_me_get_luma_ref_cost(f265_me_ctx *me, f265_pix *ref, int ref_stride, f265_mv mv)
{
    int mv_cost = venc_me_mv_cost_test(me, mv, 0);
    int dist = venc_me_get_dist(me, me->src_planes[0], me->ref_stride, ref, ref_stride,
                                me->dim[0], me->dim[1], me->bit_depth[0]);
    int cost = mv_cost + dist;
    if (cost < me->best_cost)
    {
        me->best_cost = cost;
        me->ref[0].best_mv = mv;
    }
}

// 1-iteration square search for halfpel/quarterpel.
void venc_me_square_search(int32_t nb_iter, int32_t scale, f265_me_ctx *me)
{
    // Order: N, S, W, E, NW, SW, NE, SE.
    const int8_t square_pos[2][8][2] =
    {
        { {0, -1}, {0, 1}, { -1, 0}, {1, 0}, { -1, -1}, { -1, 1}, {1, -1}, {1, 1} },
        { {0, -2}, {0, 2}, { -2, 0}, {2, 0}, { -2, -2}, { -2, 2}, {2, -2}, {2, 2} },
    };

    assert(scale == 0 || scale == 1);

    // Center the MV on the nearest fpel/hpel.
    f265_mv orig_mv = me->ref[0].best_mv;
    f265_mv center_mv = orig_mv;
    for (int comp = 0; comp < 2; comp++)
    {
        if (scale == 0) center_mv.v[comp] &= ~1;
        else center_mv.v[comp] = (center_mv.v[comp] + 2) & ~3;
    }

    // Test the center MV if we moved.
    if (center_mv.p != orig_mv.p)
        venc_me_keep_best_mv(center_mv, me);

    // Generic case. Disabling qpel until we have precomputed halfpel planes.
    if (scale == 0 || me->chroma_flag)
    {
        venc_me_test_all_positions(center_mv, me, square_pos[scale], 8);
        return;
    }

    // Luma reference plane centered on the block.
    f265_pix *base_ref = me->ref[0].ref_planes[0];

    // Buffer big enough to cover one extra pixel in each direction.
    F265_ALIGN64 f265_pix ref[66*66];
    int ref_stride=66;

    // Compute the surrounding MVs, in order.
    f265_mv mvs[8];
    for (int i = 0; i < 8; i++)
        for (int comp = 0; comp < 2; comp++)
            mvs[i].v[comp] = center_mv.v[comp] + square_pos[scale][i][comp];

    // Halfpel.
    if (scale == 1)
    {
        // Vertical. Interpolate at (0,-2).
        {
            venc_me_get_luma_hpel_ref(me, ref, base_ref, mvs[0], 0, 1);
            venc_me_get_luma_ref_cost(me, ref, ref_stride, mvs[0]);
            venc_me_get_luma_ref_cost(me, ref+ref_stride, ref_stride, mvs[1]);
        }

        // Horizontal. Interpolate at (-2,0).
        {
            venc_me_get_luma_hpel_ref(me, ref, base_ref, mvs[2], 1, 0);
            venc_me_get_luma_ref_cost(me, ref, ref_stride, mvs[2]);
            venc_me_get_luma_ref_cost(me, ref+1, ref_stride, mvs[3]);
        }

        // Diagonal. Interpolate at (-2,-2).
        {
            venc_me_get_luma_hpel_ref(me, ref, base_ref, mvs[4], 1, 1);
            venc_me_get_luma_ref_cost(me, ref, ref_stride, mvs[4]);
            venc_me_get_luma_ref_cost(me, ref+ref_stride, ref_stride, mvs[5]);
            venc_me_get_luma_ref_cost(me, ref+1, ref_stride, mvs[6]);
            venc_me_get_luma_ref_cost(me, ref+ref_stride+1, ref_stride, mvs[7]);
        }
    }

    // Quarterpel.
    else
    {
        for (int i = 0; i < 8; i++)
        {
            venc_me_get_luma_qpel_ref(me, ref, base_ref, mvs[i]);
            venc_me_get_luma_ref_cost(me, ref, ref_stride, mvs[i]);
        }

        // Optionally, correct the final cost. It seems not very useful in RDO,
        // but pretty useful in RDM.
        if (1)
        {
            me->best_cost = F265_MAX_SAD;
            venc_me_keep_best_mv(me->ref[0].best_mv, me);
        }
    }
}

// During the search, we pack the cost of a position with its ID inside a single
// integer ((cost<<shift|id). That way, we can update the best cost and the best
// position with a single conditional move operation.
//
// The following functions implement that scheme. We use macros to force the
// compiler to generate the code we want.
#define COST(idx) ((costs[idx]<<shift)+idx+1)
#define CMP(idx) do { if (COST(idx)<best_cost) best_cost = COST(idx); } while (0)
static inline int venc_update_best_cost3(int best_cost, int *costs)
{ int shift = 3; CMP(0); CMP(1); CMP(2); return best_cost; }
static inline int venc_update_best_cost4(int best_cost, int *costs)
{ int shift = 3; CMP(0); CMP(1); CMP(2); CMP(3); return best_cost; }
static inline int venc_update_best_cost6(int best_cost, int *costs)
{ int shift = 3; CMP(0); CMP(1); CMP(2); CMP(3); CMP(4); CMP(5); return best_cost; }
static inline int venc_update_best_cost8(int best_cost, int *costs)
{ int shift = 4; CMP(0); CMP(1); CMP(2); CMP(3); CMP(4); CMP(5); CMP(6); CMP(7); return best_cost; }
#undef COST
#undef CMP

// Experimental code.
void venc_fast_fullpel_search(f265_enc_thread *t)
{
    f265_me_ctx *me = &t->me;
    int packed_dims = me->packed_dims[0];
    int awi = packed_dims>>24;
    venc_sad3_func sad3 = venc_sad3[awi];
    venc_sad4_func sad4 = venc_sad4[awi];
    f265_pix *src = me->src_planes[0];
    f265_pix *fref = me->ref[0].planes[0];
    f265_pix *refs[8];
    int costs[8];
    int stride = me->ref_stride;
    int bmx = me->ref[0].best_mv.x>>2;
    int bmy = me->ref[0].best_mv.y>>2;
    int bc = me->best_cost;
    int nb_iter = 16;
    int idx;
    uint32_t mvbp[2] = { me->me_bounds_packs[0], me->me_bounds_packs[1] };

    // Call the SADx function.
    #define SAD3(off) sad3(costs + off, src, stride, refs + off, stride, packed_dims)
    #define SAD4(off) sad4(costs + off, src, stride, refs + off, stride, packed_dims)

    // Set the (ox, oy) reference from the best motion vector at the array
    // offset specified.
    #define REF(off, ox, oy)\
    do\
    {\
        refs[off] = fref + bmy*stride + bmx + (oy)*stride + (ox);\
    } while (0)

    // Add the (ox, oy) reference motion vector cost at the array offset
    // specified.
    #define MVC(off, ox, oy)\
    do\
    {\
        f265_mv tmv = {{ (bmx+(ox))<<2, (bmy+(oy))<<2 }};\
        costs[off] += venc_me_mv_cost_test(me, tmv, 0);\
    } while (0)

    // Compute the costs of 3 or 4 references. Do not use a loop for this. The
    // "merge SAD costs" operation has a long latency in assembly, so it's
    // preferable to compute the motion vector costs simultaneously.
    #define REF3(off, ox0, oy0, ox1, oy1, ox2, oy2)\
    do\
    {\
        REF(off+0, ox0, oy0); REF(off+1, ox1, oy1); REF(off+2, ox2, oy2); SAD3(off);\
        MVC(off+0, ox0, oy0); MVC(off+1, ox1, oy1); MVC(off+2, ox2, oy2);\
    } while(0)

    #define REF4(off, ox0, oy0, ox1, oy1, ox2, oy2, ox3, oy3)\
    do\
    {\
        REF(off+0, ox0, oy0); REF(off+1, ox1, oy1); REF(off+2, ox2, oy2); REF(off+3, ox3, oy3); SAD4(off);\
        MVC(off+0, ox0, oy0); MVC(off+1, ox1, oy1); MVC(off+2, ox2, oy2); MVC(off+3, ox3, oy3);\
    } while(0)


    // Hexagon search:
    //               Search the following 6 positions around the center 'X':
    //     NW NE       W:  (-2,  0)  E:  (2,  0).
    //      sss        NW: (-1,  2)  NE: (1,  2).
    //     WsXsE       SW: (-1, -2)  NW: (1, -2).
    //      sss      When 'X' matches better than all the hexagon pixels:
    //     SW SE       Search the 8 pixels surrounding X, i.e. square search.
    if (1)
    {
        // The positions on each side of the array are used for unification.
        const int8_t hex_pos[8][2] = { {-1, -2}, {-2, 0}, {-1, 2}, {1, 2}, {2, 0}, {1, -2}, {-1, -2}, {-2, 0} };

        // Notice that the next iteration of the hexagon search overlaps with
        // three positions of the previous iteration. Those positions must be
        // skipped. Label the hexagon directions as follow: W:0, NW:1, NE:2,
        // E:3, SE:4, SW:5. Let 'D' be the direction taken by the previous
        // iteration. The three directions to search are (D-1)%6, D, (D+1)%6.
        // The unstable test order slightly lowers quality.
        const uint8_t next_hex_dir[8] = { 5, 0, 1, 2, 3, 4, 5, 0 };

        // The first position is used for unification.
        const int8_t square_pos[9][2] =
            { {0, 0}, {0, -1}, {0, 1}, {-1, 0}, {1, 0}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1} };

        // Hexagon search.
        bc <<= 3;
        do
        {
            // First hexagon iteration.
            REF3(0, -2,0, -1,2, 1,2); REF3(3, 2,0, 1,-2, -1,-2);
            bc = venc_update_best_cost6(bc, costs); idx = bc&7; bc &= ~7;
            if (!idx) break;
            bmx += hex_pos[idx][0]; bmy += hex_pos[idx][1];
            int dir = next_hex_dir[idx];

            // Next hexagon iterations.
            for (int iter = 1; iter < nb_iter && !venc_mv_out_of_range2(bmx<<2, bmy<<2, mvbp); iter++)
            {
                REF3(0, hex_pos[dir+0][0], hex_pos[dir+0][1],
                        hex_pos[dir+1][0], hex_pos[dir+1][1],
                        hex_pos[dir+2][0], hex_pos[dir+2][1]);
                bc = venc_update_best_cost3(bc, costs); idx = bc&7; bc &= ~7;
                if (!idx) break;
                bmx += hex_pos[dir + idx - 1][0]; bmy += hex_pos[dir + idx - 1][1];
                dir = next_hex_dir[dir + idx - 1];
            }

        } while (0);

        // Square search.
        REF4(0, 0,-1, 0,1, -1,0, 1,0); REF4(4, -1,-1, -1,1, 1,-1, 1,1);
        bc = venc_update_best_cost8(bc<<1, costs); idx = bc&15; bc >>= 4;
        bmx += square_pos[idx][0]; bmy += square_pos[idx][1];
    }

    me->ref[0].best_mv.x = bmx<<2;
    me->ref[0].best_mv.y = bmy<<2;
    me->best_cost = bc;

    #undef SAD3
    #undef SAD4
    #undef REF
    #undef MVC
    #undef REF3
    #undef REF4
}

// Return the halfpel plane corresponding to the motion vector specified.
static f265_pix* venc_get_halfpel_plane_from_mv(f265_pix **ref_planes, int mx, int my, int stride)
{
    int ref_off = (my>>2)*stride + (mx>>2);
    int plane_off = (my&2) + ((mx&2)>>1);
    return ref_planes[plane_off] + ref_off;
}

// Experimental code.
static void venc_fast_halfpel_search(f265_enc_thread *t)
{
    // The first position is used for unification.
    const int8_t x_pos[5][2] = { {0, 0}, {-2, 2}, {2, 2}, {-2, -2}, {2, -2} };
    const int8_t dia_pos[5][2] = { {0, 0}, {0, -2}, {0, 2}, { -2, 0}, {2, 0} };

    f265_me_ctx *me = &t->me;
    int packed_dims = me->packed_dims[0];
    int awi = packed_dims>>24;
    venc_fsad_func sad = venc_fsad[awi];
    venc_sad4_func sad4 = venc_sad4[awi];
    f265_pix *src = me->src_planes[0];
    f265_pix **ref_planes = me->ref[0].planes;
    f265_pix *refs[4];
    int costs[4];
    int stride = me->ref_stride;
    int bc = me->best_cost;
    int idx;

    #if 0
    // Center the MV on the nearest fullpel.
    //
    // We probably already tested the center position, but we lost its cost.
    // Pretend it is the best position temporarily for the search. The
    // understimated center position cost does not seem to affect quality much.
    //
    // This is only done for the halfpel search, to avoid searching quarterpel
    // positions. The quarterpel search is faster if the search center is a
    // quarterpel position since less interpolation needs to be done, and
    // quality improves as well.
    //
    // FIXME: validate those results with more tests.
    f265_mv orig_mv = me->ref[0].best_mv;
    f265_mv center_mv = {{ (orig_mv.x+2) & ~3, (orig_mv.y+2) & ~3 }};
    int bmx = center_mv.x, bmy = center_mv.y;

    #else
    // Keep the best candidate before the search.
    int orig_cost = bc;
    f265_mv orig_mv = me->ref[0].best_mv;

    // Center the motion vector on a halfpel position.
    int bmx = orig_mv.x & ~1, bmy = orig_mv.y & ~1;

    // The motion vector was not aligned on a halfpel position. Get the new
    // cost.
    if (orig_mv.p & 0x10001)
    {
        f265_pix *ref = venc_get_halfpel_plane_from_mv(ref_planes, bmx, bmy, stride);
        bc = sad(src, stride, ref, stride, packed_dims);
    }
    #endif

    #define SAD4() sad4(costs, src, stride, refs, stride, packed_dims)

    #define MVC(off, ox, oy)\
    do\
    {\
        f265_mv tmv = {{ (bmx+(ox)), (bmy+(oy)) }};\
        costs[off] += venc_me_mv_cost_test(me, tmv, 0);\
    } while (0)

    #define REF4(ref0, ref1, ref2, ref3, ox0, oy0, ox1, oy1, ox2, oy2, ox3, oy3)\
    do\
    {\
        refs[0] = ref0; refs[1] = ref1; refs[2] = ref2; refs[3] = ref3; SAD4();\
        MVC(0, ox0, oy0); MVC(1, ox1, oy1); MVC(2, ox2, oy2); MVC(3, ox3, oy3);\
    } while(0)

    // Xdia unweighted.
    if (1)
    {
        // X search.
        {
            bc <<= 3;
            f265_pix *nw = venc_get_halfpel_plane_from_mv(ref_planes, bmx - 2, bmy - 2, stride);
            REF4(nw+stride, nw+stride+1, nw, nw+1, -2,2, 2,2, -2,-2, 2,-2);
            bc = venc_update_best_cost4(bc, costs); idx = bc&7; bc &= ~7;
            bmx += x_pos[idx][0]; bmy += x_pos[idx][1];
        }

        // Diamond search.
        {
            f265_pix *tp = venc_get_halfpel_plane_from_mv(ref_planes, bmx, bmy - 2, stride);
            f265_pix *lp = venc_get_halfpel_plane_from_mv(ref_planes, bmx - 2, bmy, stride);
            REF4(tp, tp+stride, lp, lp+1, 0,-2, 0,2, -2,0, 2,0);
            bc = venc_update_best_cost4(bc, costs); idx = bc&7; bc >>= 3;
            bmx += dia_pos[idx][0]; bmy += dia_pos[idx][1];
        }
    }

    me->ref[0].best_mv.x = bmx;
    me->ref[0].best_mv.y = bmy;
    me->best_cost = bc;

    // Restore the non-centered position as needed.
    #if 0
    if (me->ref[0].best_mv.p == center_mv.p) me->ref[0].best_mv = orig_mv;
    #else
    if (orig_cost < bc)
    {
        me->ref[0].best_mv = orig_mv;
        me->best_cost = orig_cost;
    }
    #endif

    #undef SAD4
    #undef MVC
    #undef REF4
}

// Experimental code.
static void venc_fast_quarterpel_search(f265_enc_thread *t)
{
    // The first position is used for unification.
    const int8_t x_pos[5][2] = { {0, 0}, {-1, 1}, {1, 1}, {-1, -1}, {1, -1} };
    const int8_t dia_pos[5][2] = { {0, 0}, {0, -1}, {0, 1}, { -1, 0}, {1, 0} };

    f265_me_ctx *me = &t->me;
    f265_pix *luma_src = me->src_planes[0];
    int plane_stride = t->plane_stride;
    int costs[4];
    int idx;
    int bc = me->best_cost;
    int bmx = me->ref[0].best_mv.x, bmy = me->ref[0].best_mv.y;
    int packed_dims = me->packed_dims[0];
    int width = (packed_dims>>8)&0xff, height = packed_dims&0xff, awi = packed_dims>>24;
    venc_fsad_func luma_sad = venc_fsad[awi];
    int aligned_block_size = F265_ALIGN_VAL(width*height, 64);
    int alloc_size = aligned_block_size*sizeof(f265_pix);
    f265_pix *tmp_buf = (f265_pix*)t->store;
    t->store += alloc_size;

    #define PATTERN_SEARCH(pattern_pos)\
    do\
    {\
        for (int i = 0; i < 4; i++)\
        {\
            int mx = bmx + pattern_pos[i+1][0], my = bmy + pattern_pos[i+1][1], ref_stride;\
            f265_mv tmv = {{ mx, my }};\
            f265_pix *ref = venc_fake_luma_ref_p(t, &ref_stride, tmp_buf, mx, my);\
            costs[i] = luma_sad(luma_src, plane_stride, ref, ref_stride, packed_dims);\
            costs[i] += venc_me_mv_cost_test(me, tmv, 0);\
        }\
    } while (0)

    // Xdia unweighted.
    if (1)
    {
        // X search.
        {
            bc <<= 3;
            PATTERN_SEARCH(x_pos);
            bc = venc_update_best_cost4(bc, costs); idx = bc&7; bc &= ~7;
            bmx += x_pos[idx][0]; bmy += x_pos[idx][1];
        }

        // Diamond search.
        {
            PATTERN_SEARCH(dia_pos);
            bc = venc_update_best_cost4(bc, costs); idx = bc&7; bc >>= 3;
            bmx += dia_pos[idx][0]; bmy += dia_pos[idx][1];
        }
    }

    me->ref[0].best_mv.x = bmx;
    me->ref[0].best_mv.y = bmy;
    me->best_cost = bc;
    t->store -= alloc_size;

    #undef PATTERN_SEARCH
}

// Eliminate the duplicate and null MVs in place. Return the number of remaining MVs.
int venc_me_remove_duplicate(f265_mv mv_list[], int nb_mv)
{
    f265_mv pmv = mv_list[0];
    int32_t nb_out = 1;
    for (int32_t i = 1; i < nb_mv; i++)
        if (mv_list[i].p && mv_list[i].p != pmv.p)
            mv_list[nb_out++].p = mv_list[i].p;

    return nb_out;
}

// Search for the best MV using a three pass technique.
// Each pass can use different search algo, iteration and distortion function.
// The threshold argument is a pointer to the best hpel MV cost from all references.
// When set, the best hpel cost is compared against the threshold. If the hpel cost
// is above 8/7 of the threshold value, the reference is considered unlikely to
// yield a better MV, and qpel analysis is skipped. The threshold is updated if the
// current reference yields a better hpel cost. If threshold is not set, qpel is
// always tested.
// This was only tested with P-prediction, no support for B-prediction.
// TODO The code should support weighting, but that is untested.
// FIXME: pass the thread object here for now on.
void venc_early_me(f265_me_ctx *me,
                   f265_mv mv_candidate_list[],
                   int nb_candidate,
                   int *threshold)
{
    f265_enc_thread *t = me->t;

    // Set up the distortion function and the chroma_flag.
    int8_t *dist_id = me->early_me_params->dist_func_ids;
    me->dist_func_id = dist_id[0] >> 1;
    me->chroma_flag = dist_id[0] & 1;

    // Optimization: if only the PMV is present and it is equal to its clipped
    // value, keep its current cost. Do this better eventually.
    if (nb_candidate > 1 || f265_clip_mv(mv_candidate_list[0], me->me_bounds64).p != me->ref[0].best_mv.p)
    {
        // Reset the cost.
        me->best_cost = F265_MAX_SAD;

        // Get the best clipped candidate. FIXME: handle clipping better.
        nb_candidate = venc_me_remove_duplicate(mv_candidate_list, nb_candidate);
        for (int i = 0; i < nb_candidate; i++)
            venc_me_keep_best_mv(f265_clip_mv(mv_candidate_list[i], me->me_bounds64), me);
    }

    // Save the best candidate.
    f265_mv best_fpel_mv = me->ref[0].best_mv;
    int best_cand_cost = me->best_cost;

    // Round the best candidate to fpel.
    me->ref[0].best_mv.x = (me->ref[0].best_mv.x + 2) & ~3;
    me->ref[0].best_mv.y = (me->ref[0].best_mv.y + 2) & ~3;
    me->best_cost = venc_me_mv_total_cost(me->ref[0].best_mv, me);

    // Test the null candidate.
    f265_mv zero_cand; zero_cand.p = 0;
    venc_me_keep_best_mv(zero_cand, me);


    // If the first search function is not set, the job is done.
    if (!me->early_me_params->search_funcs[0]) return;

    // Search for the best fpel MV.
    if (t->enc->gd.algo&0x400) venc_fast_fullpel_search(t);
    else (*me->early_me_params->search_funcs[0])(me->early_me_params->nb_iters[0], 2, me);

    // Keep the best candidate MV if it has a better cost than the best fpel MV.
    if (best_cand_cost < me->best_cost)
    {
        me->best_cost = best_cand_cost;
        me->ref[0].best_mv = best_fpel_mv;
    }


    // If the second search function is not set, the job is done.
    if (!me->early_me_params->search_funcs[1]) return;

    // If a different cost metric is used, find the best MV cost in the new metric.
    if ((me->dist_func_id<<1) + me->chroma_flag != dist_id[1])
    {
        // Set up the new metric.
        me->dist_func_id = dist_id[1] >> 1;
        me->chroma_flag = dist_id[1] & 1;

        // Get the new cost.
        me->best_cost = venc_me_mv_total_cost(me->ref[0].best_mv, me);
    }

    // Search for the best hpel MV.
    if (t->enc->gd.algo&0x400) venc_fast_halfpel_search(t);
    else (*me->early_me_params->search_funcs[1])(me->early_me_params->nb_iters[1], 1, me);

    // If the third search function is not set, the job is done.
    if (!me->early_me_params->search_funcs[2]) return;

    // If a different cost metric is used, find the best MV cost in the new metric.
    if ((me->dist_func_id<<1) + me->chroma_flag != dist_id[2])
    {
        me->dist_func_id = dist_id[2] >> 1;
        me->chroma_flag = dist_id[2] & 1;
        me->best_cost = venc_me_mv_total_cost(me->ref[0].best_mv, me);
    }

    // Shortcut to speed up the search. If the best cost is not under 8/7 of a threshold
    // (the other reference best MV), don't take the time to search at qpel level.
    if (threshold)
    {
        int thres = (me->best_cost >> 3) * 7;
        if (thres > *threshold) return;

        // Update the threshold.
        if (me->best_cost < *threshold)
            *threshold = me->best_cost;
    }

    if (t->enc->gd.algo&0x400) venc_fast_quarterpel_search(t);
    else (*me->early_me_params->search_funcs[2])(me->early_me_params->nb_iters[2], 0, me);
}

// A wrapper around the SAD function.
// The SAD function takes one more argument than SATD.
int32_t venc_sad_wrap(f265_pix *src0, int32_t stride0, f265_pix *src1, int32_t stride1,
                      int32_t width, int32_t height, int32_t bitdepth)
{
    return venc_sad(src0, stride0, src1, stride1, width, height, 0, bitdepth);
}

// Set the ME parameter information.
void venc_me_set_params(f265_enc *enc, f265_me_ctx *me)
{
    f265_gen_data *gd = &enc->gd;
    me->enc = enc;
    me->bit_depth[0] = gd->bit_depth[0];
    me->bit_depth[1] = gd->bit_depth[1];
    me->csf[0] = gd->csf[0];
    me->csf[1] = gd->csf[1];
    me->early_me_params = gd->early_me_params;
    me->ref_stride = gd->stride;
    me->dist[0] = &venc_sad_wrap;
    me->dist[1] = &venc_satd;
}

// Set the ME data for the partition specified.
// TODO: A 4:2:0 color space is assumed.
void venc_me_set_partition(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    f265_me_ctx *me = &t->me;

    // Set the block location data.
    int part_size[2], part_off[2];
    venc_get_part_loc(part_size, part_off, cb->inter_part, part_idx, 1<<cb->lg_bs);

    for (int i = 0; i < 2; i++)
    {
        me->pos[i] = t->ctb_off[i] + cb->cb_off[i] + part_off[i];
        me->dim[i] = part_size[i];
    }

    me->plane_off[0] = me->pos[1]*me->ref_stride + me->pos[0];
    me->plane_off[1] = (me->pos[1]>>1)*me->ref_stride + (me->pos[0]>>1);

    int awi_luma = f265_width_to_awi[me->dim[0]];
    me->packed_dims[0] = awi_luma<<24 | me->dim[0]<<8 | me->dim[1];
    me->packed_dims[1] = (awi_luma-1)<<24 | (me->dim[0]>>1)<<8 | (me->dim[1]>>1);

    // Set the source plane pointers.
    me->src_planes[0] = t->src_frame->src_planes[0] + me->plane_off[0];
    me->src_planes[1] = t->src_frame->src_planes[1] + me->plane_off[1];
    me->src_planes[2] = t->src_frame->src_planes[2] + me->plane_off[1];

    // Cache the neighbour motion vectors.
    venc_get_neighbour_mvs(t, cb, part_idx, t->ib.neighbour_mv);
}

// Set the ME data for the reference index specified. The ref_id parameter is
// not the ref_idx. It is "0" for unidirectional prediction and the first
// reference in bi-prediction. It is "1" for the second reference in
// bi-prediction. TODO: A 4:2:0 color space is assumed.
void venc_me_set_ref(f265_enc_thread *t, f265_ref_ctx *ref_ctx, int ref_id)
{
    f265_me_ctx *me = &t->me;
    me->ref[ref_id].ref_planes[0] = ref_ctx->planes[0] + me->plane_off[0];
    me->ref[ref_id].ref_planes[1] = ref_ctx->planes[4] + me->plane_off[1];
    me->ref[ref_id].ref_planes[2] = ref_ctx->planes[5] + me->plane_off[1];
    me->ref[ref_id].weights = ref_ctx->weights;

    // New interface.
    for (int i = 0; i < 4; i++) me->ref[ref_id].planes[i] = ref_ctx->planes[i] + me->plane_off[0];
    me->ref[ref_id].planes[4] = ref_ctx->planes[4] + me->plane_off[1];
    me->ref[ref_id].planes[5] = ref_ctx->planes[5] + me->plane_off[1];
}

// Set the current PMV and center the mv_cost_table on the PMV.
void venc_me_set_pmv(f265_me_ctx *me, f265_mv pmv, int qp, int ref_id)
{
    me->ref[ref_id].pmv = pmv;
    me->mv_cost_table = f265_mv_costs[qp];
    int mv_offs = F265_NB_MV_COSTS >> 1;
    for (int i = 0; i < 2; i++)
        me->ref[ref_id].mv_costs[i] = me->mv_cost_table + mv_offs - pmv.v[i];
}

