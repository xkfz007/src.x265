// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include "f265/bdi.h"

// Include the auto-generated dispatch code.
#include "f265/asm.c"

// Dispatch table for handling bit depths. The index is the high bit depth flag.
static struct
{
    void (*analyze_params)(f265_enc_params *params);
    void (*init_enc)(f265_enc_params *params, uint8_t *buf, f265_enc **enc_handle, char **error_handle);
    void (*deinit_enc)(f265_enc *enc);
    int (*process_enc_req)(f265_enc *enc, f265_enc_req *req);
    void (*abort_enc)(f265_enc *enc);

} f265_bd_dispatch[2];

// Set the dispatch table for the bit depth.
static void f265_set_bit_depth_dispatch()
{
    // Link the dispatch table for each bit depth.
    #define DISPATCH(IDX, BD) \
        void f265_##BD##_analyze_params(f265_enc_params *params);\
        f265_bd_dispatch[IDX].analyze_params = f265_##BD##_analyze_params;\
        void f265_##BD##_init_enc(f265_enc_params *params, uint8_t *buf, f265_enc **enc_handle, char **error_handle);\
        f265_bd_dispatch[IDX].init_enc = f265_##BD##_init_enc;\
        void f265_##BD##_deinit_enc(f265_enc *enc);\
        f265_bd_dispatch[IDX].deinit_enc = f265_##BD##_deinit_enc;\
        int f265_##BD##_process_enc_req(f265_enc *enc, f265_enc_req *req);\
        f265_bd_dispatch[IDX].process_enc_req = f265_##BD##_process_enc_req;
    #ifdef F265_HAVE_LBD
    DISPATCH(0, lbd);
    #endif
    #ifdef F265_HAVE_HBD
    DISPATCH(1, hbd);
    #endif
    #undef DISPATCH
}

// Detect the assembly architecture (i.e. a set of features that the processor
// support). The return value is a bitfield, as follow:
//   Bit 1: "SSE41" (SSE 4.1).
//   Bit 2: "AVX2" (AVX2, popcnt, movbe, BMI2).
//
// The reliability of this code is not guaranteed since the API is crazy.
// We're mostly using the method suggested by Agner Fog.
static int f265_detect_asm_arch()
{
    // Use the CPUID instruction to identify the supported processor features.
    #define do_cpuid(func, eax, ebx, ecx, edx) \
        __asm__ __volatile__ ("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (func), "c" (0));

    // Use the XGETBV instruction to identify which features the OS supports.
    #ifdef __APPLE__
    // GCC on OS X complains that xgetbv is an unknown instruction.
    #define do_xgetbv(xcr0) __asm__ __volatile__ (".byte 0x0f, 0x01, 0xd0" : "=a" (xcr0) : "c" (0) : "%edx");
    #else
    #define do_xgetbv(xcr0) __asm__ __volatile__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
    #endif

    // Features detected.
    int sse41_flag = 0, movbe_flag = 0, popcnt_flag = 0, osxsave_flag = 0, avx_flag = 0, avx2_flag = 0, bmi2_flag = 0,
        ymm_on_flag = 0;

    // Register content.
    uint32_t eax, ebx, ecx, edx, xcr0;

    // CPUID function 1.
    do_cpuid(1, eax , ebx , ecx, edx);      // Assuming function 1 is available on 64-bits processors.
    sse41_flag =   !!(ecx & 0x00080000);    // Bit 19.
    movbe_flag =   !!(ecx & 0x00400000);    // Bit 22.
    popcnt_flag =  !!(ecx & 0x00800000);    // Bit 23.
    osxsave_flag = !!(ecx & 0x08000000);    // Bit 27.
    avx_flag =     !!(ecx & 0x10000000);    // Bit 28.

    // Safe to execute XGETBV to check if the YMM registers are supported by the OS.
    if (osxsave_flag && avx_flag)
    {
        do_xgetbv(xcr0);
        ymm_on_flag = (xcr0&6) == 6;
    }

    // Safe to execute CPUID function 7.
    if (ymm_on_flag)
    {
        do_cpuid(7, eax , ebx , ecx, edx);
        avx2_flag =    !!(ebx & 0x00000020);    // Bit 5.
        bmi2_flag =    !!(ebx & 0x00000100);    // Bit 8.
    }

    int ret = 0;
    ret |= (sse41_flag) << 0;
    ret |= (avx2_flag && popcnt_flag && movbe_flag && bmi2_flag) << 1;
    return ret;
}

// Set the dispatch tables for the assembly architecture.
static void f265_set_asm_dispatch(int no_asm_flag)
{
    int asm_arch_bf = no_asm_flag ? 0 : f265_detect_asm_arch();
    f265_link_asm(!!(asm_arch_bf&2));
}

// Populate f265_mv_costs table. The table yields the cost of a PMV in the SAD domain.
static void f265_generate_mv_cost_table()
{
    // Number of negative (or positive) motion vector cost offsets.
    int mv_offs = F265_NB_MV_COSTS >> 1;

    // Compute the logs.
    float cabac_magic = 0.718;
    float cabac_logs[mv_offs + 1];
    cabac_logs[0] = cabac_magic;
    for (int i = 1; i <= mv_offs; i++) cabac_logs[i] = 2.0 * log2f(i + 1) + 1.0 + cabac_magic;

    // Set up the motion vector costs.
    for (int qp = 0; qp < 52; qp++)
    {
        int lambda = f265_lambdas[qp];
        uint16_t *p = f265_mv_costs[qp] + mv_offs;
        for (int i = 0; i <= mv_offs; i++)
            p[-i] = p[i] = ((int)(lambda * cabac_logs[i] + 8.0)) >> 4;
    }
}

void f265_init_global(int no_asm_flag)
{
    f265_set_bit_depth_dispatch();
    f265_set_asm_dispatch(no_asm_flag);
    f265_generate_mv_cost_table();
}

void f265_set_default_params(f265_enc_params *p)
{
    memset(p, 0, sizeof(f265_enc_params));
    p->nb_workers[0] = -1;
    p->chroma_format = 1;
    for (int i = 0; i < 4; i++) p->bit_depth[i] = 8;
    p->cb_range[0] = 3;
    p->cb_range[1] = 6;
    p->pcm_range[0] = 3;
    p->pcm_range[1] = 5;
    p->tb_range[0] = 2;
    p->tb_range[1] = 5;
    p->tb_depth[0] = 4;
    p->tb_depth[1] = 4;
    p->qg_log = -1;
    p->nb_refs = 1;
    p->merge_cand = 5;
    p->parallel_merge_level = 2;
    p->deblock_flag = 1;
    p->smooth_intra_flag = 1;
    p->tmv_flag = 1;
    p->key_interval = 120;
    p->frame_rate_num = 30;
    p->frame_rate_den = 1;
    p->bitrate[0] = 500;
    p->bitrate[1] = 0;
    p->bitrate[2] = 0;
    p->rc_method = F265_RCM_CQP;
    p->qp_bounds[0] = 0;
    p->qp_bounds[1] = 51;
    p->lt_conv_min = 10;
    p->lt_conv_exp = 0.6;
    p->qp = 30;
    p->me_dist[0] = p->me_dist[1] = 0;
    p->me_dist[2] = 0;
    p->me_algo[0] = 2;
    p->me_algo[1] = p->me_algo[2] = 1;
    p->me_iter[0] = 16;
    p->me_iter[1] = p->me_iter[2] = 1;
}

void f265_normalize_params(f265_enc_params *p)
{
    #define CL(V, min, max) (V) = F265_CLAMP((V), (min), (max))
    #define CLF(V) (V) = !!(V)
    #define CLP(V) (V) = F265_MAX(0, (V))
    for (int i = 0; i < 2; i++) CL(p->clip_dim[i], 1, 16384);
    CL(p->mt_mode, 0, 2);
    CL(p->chroma_format, 0, 3);
    for (int i = 0; i < 4; i++) CL(p->bit_depth[i], 8, 14);
    CL(p->cb_range[0], 3, 6);
    CL(p->cb_range[1], p->cb_range[0], 6);
    CL(p->pcm_range[0], -1, 5);
    CL(p->pcm_range[1], p->pcm_range[0], 5);
    CL(p->tb_range[0], 2, 5);
    CL(p->tb_range[1], p->tb_range[0], 5);
    int max_tb_depth = p->cb_range[1] - p->tb_range[0];
    CL(p->tb_depth[0], 0, max_tb_depth);
    CL(p->tb_depth[1], 0, max_tb_depth);
    CL(p->qg_log, -1, 6);
    CL(p->nb_refs, 0, 16);
    CL(p->nb_b_frames, 0, 16);
    CL(p->chroma_qp_idx_off, -12, 12);
    for (int i = 0; i < 2; i++) CL(p->deblock_off[i], -6, 6);
    CL(p->merge_cand, 1, 5);
    CL(p->parallel_merge_level, 2, 6);
    for (int i = 0; i < 4; i++) CL(p->me_algo[i], 0, 3);
    if (p->me_algo[0] == 3) p->me_algo[0] = 2;
    for (int i = 0; i < 4; i++) CL(p->me_iter[i], 0, 16);
    for (int i = 0; i < 4; i++) CL(p->me_dist[i], 0, 1);
    CLF(p->wpp_flag);
    CLF(p->deblock_flag);
    CLF(p->sao_flag);
    CLF(p->scaling_flag);
    CLF(p->transquant_bypass_flag);
    CLF(p->sign_hiding_flag);
    CLF(p->transform_skip_flag);
    CLF(p->smooth_intra_flag);
    CLF(p->amp_flag);
    CLF(p->tmv_flag);
    CLF(p->weight_flag);
    CLF(p->chroma_me_flag);
    CL(p->rdo_level, 0, 1);
    CLF(p->hm_me_flag);
    CLF(p->all_intra_flag);
    CLF(p->nullify_inter_tb_flag);
    CLF(p->la_flag);
    CLF(p->reformat_flag);
    CL(p->la_decision_delay, 0, 1024);
    CL(p->la_process_delay, 0, 8);
    CLP(p->key_interval);
    CLP(p->frame_rate_num);
    p->frame_rate_den = F265_MAX(p->frame_rate_den, 1);
    CL(p->qp, 0, 51);
    #undef CL
    #undef CLF
    #undef CLP

    // Make sure there is room for B frames in the lookahead.
    p->la_decision_delay = F265_MAX(p->la_decision_delay, p->nb_b_frames);

    // Normalize the number of threads. FIXME.
    p->nb_workers[0] = p->nb_workers[1] = 0;
}

// Dispatch by the parameters/encoder bit depth. We rely on the bit depth being
// the first field of the encoder object.
#define DISPATCH f265_bd_dispatch[params->hbd_flag]
void f265_analyze_params(f265_enc_params *params) { DISPATCH.analyze_params(params); }
void f265_init_enc(f265_enc_params *params, uint8_t *buf, f265_enc **enc_handle, char **error_handle)
    { DISPATCH.init_enc(params, buf, enc_handle, error_handle); }
#undef DISPATCH
#define DISPATCH f265_bd_dispatch[(*(uint8_t*)enc)]
void f265_deinit_enc(f265_enc *enc) { DISPATCH.deinit_enc(enc); }
int f265_process_enc_req(f265_enc *enc, f265_enc_req *req) { return DISPATCH.process_enc_req(enc, req); }
void f265_abort_enc(f265_enc *enc) {}
#undef DISPATCH

