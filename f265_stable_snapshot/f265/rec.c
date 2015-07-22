// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Block encoding and reconstruction.

#include "f265/enc.h"

// Compile the encoding functions in RDO mode.
#define F265_RDO
#include "entropy.c"

extern int venc_trace_syntax_flag;

// Compute the residual of the block specified (src-pred). The residual is
// stored contiguously.
static void venc_get_block_residual(int16_t *dst, f265_pix *src, int src_stride, f265_pix *pred, int pred_stride,
                                    int bs)
{
    for (int y = 0; y < bs; y++, dst += bs, src += src_stride, pred += pred_stride)
        for (int x = 0; x < bs; x++)
            dst[x] = src[x] - pred[x];
}

// Compute the reconstruction of the block specified (pred+residual).
static void venc_get_block_reconstruction(f265_pix *dst, int dst_stride, f265_pix *pred, int pred_stride,
                                          int16_t *coeffs, int bs, int bd)
{
    int pix_clamp = (1<<bd)-1;
    for (int y = 0; y < bs; y++, dst += dst_stride, pred += pred_stride, coeffs += bs)
        for (int x = 0; x < bs; x++)
            dst[x] = F265_CLAMP(pred[x] + coeffs[x], 0, pix_clamp);
}

// DCT/IDCT helper functions.

// 4x4 DST, must be named this way for compatibility.
static void venc_dct_dst_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 4; i++, dst++, src += 4)
    {
        int a = src[0]+src[3];
        int b = src[0]+src[1];
        int c = src[0]-src[1];
        int d = src[1]+src[3];
        int e = 74*(b - src[3]);
        int f = 74*src[2];

        dst[0]  = (29*a + 55*d + f + add) >> shift;
        dst[4]  = (e + add) >> shift;
        dst[8]  = (29*c + 55*a - f + add) >> shift;
        dst[12] = (55*c - 29*d + f + add) >> shift;
    }
}

static void venc_idct_dst_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 4; i++, dst += 4, src++)
    {
        int a = src[0]+src[8];
        int b = src[0]-src[8];
        int c = src[0]-src[12];
        int d = src[8]+src[12];
        int e = 74*(b + src[12]);
        int f = 74*src[4];

        dst[0] = F265_CLAMP((29*a + 55*d + f + add) >> shift, -32768, 32767);
        dst[1] = F265_CLAMP((55*c - 29*d + f + add) >> shift, -32768, 32767);
        dst[2] = F265_CLAMP((e + add) >> shift, -32768, 32767);
        dst[3] = F265_CLAMP((55*a + 29*c - f + add) >> shift, -32768, 32767);
    }
}

static void venc_dct_4_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 4; i++, dst++, src += 4)
    {
        int add_0_3 = src[0]+src[3], add_1_2 = src[1]+src[2];
        int sub_0_3 = src[0]-src[3], sub_1_2 = src[1]-src[2];

        dst[0]  = (64*add_0_3 + 64*add_1_2 + add) >> shift;
        dst[4]  = (83*sub_0_3 + 36*sub_1_2 + add) >> shift;
        dst[8]  = (64*add_0_3 - 64*add_1_2 + add) >> shift;
        dst[12] = (36*sub_0_3 - 83*sub_1_2 + add) >> shift;
    }
}

static void venc_idct_4_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 4; i++, dst += 4, src++)
    {
        int add_0 = 64*src[0] + 64*src[8];
        int add_1 = 83*src[4] + 36*src[12];
        int sub_0 = 64*src[0] - 64*src[8];
        int sub_1 = 36*src[4] - 83*src[12];

        dst[0] = F265_CLAMP((add_0 + add_1 + add) >> shift, -32768, 32767);
        dst[1] = F265_CLAMP((sub_0 + sub_1 + add) >> shift, -32768, 32767);
        dst[2] = F265_CLAMP((sub_0 - sub_1 + add) >> shift, -32768, 32767);
        dst[3] = F265_CLAMP((add_0 - add_1 + add) >> shift, -32768, 32767);
    }
}

static void venc_dct_8_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 8; i++, dst++, src += 8)
    {
        int add_0_7 = src[0]+src[7], add_1_6 = src[1]+src[6], add_2_5 = src[2]+src[5], add_3_4 = src[3]+src[4];
        int sub_0_7 = src[0]-src[7], sub_1_6 = src[1]-src[6], sub_2_5 = src[2]-src[5], sub_3_4 = src[3]-src[4];

        dst[0]  = (64*add_0_7 + 64*add_1_6 + 64*add_2_5 + 64*add_3_4 + add) >> shift;
        dst[8]  = (89*sub_0_7 + 75*sub_1_6 + 50*sub_2_5 + 18*sub_3_4 + add) >> shift;
        dst[16] = (83*add_0_7 + 36*add_1_6 - 36*add_2_5 - 83*add_3_4 + add) >> shift;
        dst[24] = (75*sub_0_7 - 18*sub_1_6 - 89*sub_2_5 - 50*sub_3_4 + add) >> shift;
        dst[32] = (64*add_0_7 - 64*add_1_6 - 64*add_2_5 + 64*add_3_4 + add) >> shift;
        dst[40] = (50*sub_0_7 - 89*sub_1_6 + 18*sub_2_5 + 75*sub_3_4 + add) >> shift;
        dst[48] = (36*add_0_7 - 83*add_1_6 + 83*add_2_5 - 36*add_3_4 + add) >> shift;
        dst[56] = (18*sub_0_7 - 50*sub_1_6 + 75*sub_2_5 - 89*sub_3_4 + add) >> shift;
    }
}

static void venc_idct_8_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 8; i++, dst += 8, src++)
    {
        int add_0 = 64*src[0] + 64*src[32] + 83*src[16] + 36*src[48];
        int add_1 = 64*src[0] - 64*src[32] + 36*src[16] - 83*src[48];
        int add_2 = 64*src[0] - 64*src[32] - 36*src[16] + 83*src[48];
        int add_3 = 64*src[0] + 64*src[32] - 83*src[16] - 36*src[48];

        int sub_0 = 89*src[8] + 75*src[24] + 50*src[40] + 18*src[56];
        int sub_1 = 75*src[8] - 18*src[24] - 89*src[40] - 50*src[56];
        int sub_2 = 50*src[8] - 89*src[24] + 18*src[40] + 75*src[56];
        int sub_3 = 18*src[8] - 50*src[24] + 75*src[40] - 89*src[56];

        dst[0] = F265_CLAMP((add_0 + sub_0 + add) >> shift, -32768, 32767);
        dst[1] = F265_CLAMP((add_1 + sub_1 + add) >> shift, -32768, 32767);
        dst[2] = F265_CLAMP((add_2 + sub_2 + add) >> shift, -32768, 32767);
        dst[3] = F265_CLAMP((add_3 + sub_3 + add) >> shift, -32768, 32767);

        dst[4] = F265_CLAMP((add_3 - sub_3 + add) >> shift, -32768, 32767);
        dst[5] = F265_CLAMP((add_2 - sub_2 + add) >> shift, -32768, 32767);
        dst[6] = F265_CLAMP((add_1 - sub_1 + add) >> shift, -32768, 32767);
        dst[7] = F265_CLAMP((add_0 - sub_0 + add) >> shift, -32768, 32767);
    }
}

static void venc_dct_16_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 16; dst++, src += 16, i++)
    {
        int a_0_15 = src[0] + src[15], a_1_14 = src[1] + src[14], a_2_13 = src[2] + src[13];
        int a_3_12 = src[3] + src[12], a_4_11 = src[4] + src[11], a_5_10 = src[5] + src[10];
        int a_6_9 = src[6] + src[9], a_7_8 = src[7] + src[8];

        int s_0_15 = src[0] - src[15], s_1_14 = src[1] - src[14], s_2_13 = src[2] - src[13];
        int s_3_12 = src[3] - src[12], s_4_11 = src[4] - src[11], s_5_10 = src[5] - src[10];
        int s_6_9 = src[6] - src[9], s_7_8 = src[7] - src[8];

        dst[0] = (64*a_0_15 + 64*a_1_14 + 64*a_2_13 + 64*a_3_12 +
                  64*a_4_11 + 64*a_5_10 + 64*a_6_9 + 64*a_7_8 + add) >> shift;

        dst[16] = (90*s_0_15 + 87*s_1_14 + 80*s_2_13 + 70*s_3_12 +
                   57*s_4_11 + 43*s_5_10 + 25*s_6_9 + 9*s_7_8 + add) >> shift;

        dst[32] = (89*(a_0_15 - a_7_8) + 75*(a_1_14 - a_6_9) +
                   50*(a_2_13 - a_5_10) + 18*(a_3_12 - a_4_11) + add) >> shift;

        dst[48] = (87*s_0_15 + 57*s_1_14 + 9*s_2_13 - 43*s_3_12 -
                   80*s_4_11 - 90*s_5_10 - 70*s_6_9 - 25*s_7_8 + add) >> shift;

        dst[64] = (83*(a_0_15 + a_7_8 - a_3_12 - a_4_11) +
                   36*(a_1_14 + a_6_9 - a_2_13 - a_5_10) + add) >> shift;

        dst[80] = (80*s_0_15 + 9*s_1_14 - 70*s_2_13 - 87*s_3_12 -
                   25*s_4_11 + 57*s_5_10 + 90*s_6_9 + 43*s_7_8 + add) >> shift;

        dst[96] = (75*(a_0_15 - a_7_8) - 18*(a_1_14 - a_6_9) -
                   89*(a_2_13 - a_5_10) - 50*(a_3_12 - a_4_11) + add) >> shift;

        dst[112] = (70*s_0_15 - 43*s_1_14 - 87*s_2_13 + 9*s_3_12 +
                    90*s_4_11 + 25*s_5_10 - 80*s_6_9 - 57*s_7_8 + add) >> shift;

        dst[128] = (64*(a_0_15 + a_7_8 + a_3_12 + a_4_11) -
                    64*(a_1_14 + a_6_9 + a_2_13 + a_5_10) + add) >> shift;

        dst[144] = (57*s_0_15 - 80*s_1_14 - 25*s_2_13 + 90*s_3_12 -
                    9*s_4_11 - 87*s_5_10 + 43*s_6_9 + 70*s_7_8 + add) >> shift;

        dst[160] = (50*(a_0_15 - a_7_8) - 89*(a_1_14 - a_6_9) +
                    18*(a_2_13 - a_5_10) + 75*(a_3_12 - a_4_11) + add) >> shift;

        dst[176] = (43*s_0_15 - 90*s_1_14 + 57*s_2_13 + 25*s_3_12 -
                    87*s_4_11 + 70*s_5_10 + 9*s_6_9 - 80*s_7_8 + add) >> shift;

        dst[192] = (36*(a_0_15 + a_7_8 - a_3_12 - a_4_11) -
                    83*(a_1_14 + a_6_9 - a_2_13 - a_5_10) + add) >> shift;

        dst[208] = (25*s_0_15 - 70*s_1_14 + 90*s_2_13 - 80*s_3_12 +
                    43*s_4_11 + 9*s_5_10 - 57*s_6_9 + 87*s_7_8 + add) >> shift;

        dst[224] = (18*(a_0_15 - a_7_8) - 50*(a_1_14 - a_6_9) +
                    75*(a_2_13 - a_5_10) - 89*(a_3_12 - a_4_11) + add) >> shift;

        dst[240] = (9*s_0_15 - 25*s_1_14 + 43*s_2_13 - 57*s_3_12 +
                    70*s_4_11 - 80*s_5_10 + 87*s_6_9 - 90*s_7_8 + add) >> shift;
    }
}

static void venc_idct_16_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 16; dst += 16, src++, i++)
    {
        int s0 = 90*src[16] + 87*src[48] + 80*src[80] + 70*src[112] +
                 57*src[144] + 43*src[176] + 25*src[208] + 9*src[240];
        int s1 = 87*src[16] + 57*src[48] + 9*src[80] - 43*src[112] -
                 80*src[144] - 90*src[176] - 70*src[208] - 25*src[240];
        int s2 = 80*src[16] + 9*src[48] - 70*src[80] - 87*src[112] -
                 25*src[144] + 57*src[176] + 90*src[208] + 43*src[240];
        int s3 = 70*src[16] - 43*src[48] - 87*src[80] + 9*src[112] +
                 90*src[144] + 25*src[176] - 80*src[208] - 57*src[240];
        int s4 = 57*src[16] - 80*src[48] - 25*src[80] + 90*src[112] -
                 9*src[144] - 87*src[176] + 43*src[208] + 70*src[240];
        int s5 = 43*src[16] - 90*src[48] + 57*src[80] + 25*src[112] -
                 87*src[144] + 70*src[176] + 9*src[208] - 80*src[240];
        int s6 = 25*src[16] - 70*src[48] + 90*src[80] - 80*src[112] +
                 43*src[144] + 9*src[176] - 57*src[208] + 87*src[240];
        int s7 = 9*src[16] - 25*src[48] + 43*src[80] - 57*src[112] +
                 70*src[144] - 80*src[176] + 87*src[208] - 90*src[240];

        int a0 = 64*src[0] + 89*src[32] + 83*src[64] + 75*src[96] +
                 64*src[128] + 50*src[160] + 36*src[192] + 18*src[224];
        int a1 = 64*src[0] + 75*src[32] + 36*src[64] - 18*src[96] -
                 64*src[128] - 89*src[160] - 83*src[192] - 50*src[224];
        int a2 = 64*src[0] + 50*src[32] - 36*src[64] - 89*src[96] -
                 64*src[128] + 18*src[160] + 83*src[192] + 75*src[224];
        int a3 = 64*src[0] + 18*src[32] - 83*src[64] - 50*src[96] +
                 64*src[128] + 75*src[160] - 36*src[192] - 89*src[224];
        int a4 = 64*src[0] - 18*src[32] - 83*src[64] + 50*src[96] +
                 64*src[128] - 75*src[160] - 36*src[192] + 89*src[224];
        int a5 = 64*src[0] - 50*src[32] - 36*src[64] + 89*src[96] -
                 64*src[128] - 18*src[160] + 83*src[192] - 75*src[224];
        int a6 = 64*src[0] - 75*src[32] + 36*src[64] + 18*src[96] -
                 64*src[128] + 89*src[160] - 83*src[192] + 50*src[224];
        int a7 = 64*src[0] - 89*src[32] + 83*src[64] - 75*src[96] +
                 64*src[128] - 50*src[160] + 36*src[192] - 18*src[224];

        dst[0] = F265_CLAMP((a0 + s0 + add) >> shift, -32768, 32767);
        dst[1] = F265_CLAMP((a1 + s1 + add) >> shift, -32768, 32767);
        dst[2] = F265_CLAMP((a2 + s2 + add) >> shift, -32768, 32767);
        dst[3] = F265_CLAMP((a3 + s3 + add) >> shift, -32768, 32767);
        dst[4] = F265_CLAMP((a4 + s4 + add) >> shift, -32768, 32767);
        dst[5] = F265_CLAMP((a5 + s5 + add) >> shift, -32768, 32767);
        dst[6] = F265_CLAMP((a6 + s6 + add) >> shift, -32768, 32767);
        dst[7] = F265_CLAMP((a7 + s7 + add) >> shift, -32768, 32767);
        dst[8] = F265_CLAMP((a7 - s7 + add) >> shift, -32768, 32767);
        dst[9] = F265_CLAMP((a6 - s6 + add) >> shift, -32768, 32767);
        dst[10] = F265_CLAMP((a5 - s5 + add) >> shift, -32768, 32767);
        dst[11] = F265_CLAMP((a4 - s4 + add) >> shift, -32768, 32767);
        dst[12] = F265_CLAMP((a3 - s3 + add) >> shift, -32768, 32767);
        dst[13] = F265_CLAMP((a2 - s2 + add) >> shift, -32768, 32767);
        dst[14] = F265_CLAMP((a1 - s1 + add) >> shift, -32768, 32767);
        dst[15] = F265_CLAMP((a0 - s0 + add) >> shift, -32768, 32767);
    }
}

static void venc_dct_32_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 32; dst++, src += 32, i++)
    {
        int a_0_31 = src[0] + src[31], a_1_30 = src[1] + src[30], a_2_29 = src[2] + src[29], a_3_28 = src[3] + src[28];
        int a_4_27 = src[4] + src[27], a_5_26 = src[5] + src[26], a_6_25 = src[6] + src[25], a_7_24 = src[7] + src[24];
        int a_8_23 = src[8] + src[23], a_9_22 = src[9] + src[22], a_10_21 = src[10] + src[21];
        int a_11_20 = src[11] + src[20], a_12_19 = src[12] + src[19], a_13_18 = src[13] + src[18];
        int a_14_17 = src[14] + src[17], a_15_16 = src[15] + src[16];

        int s_0_31 = src[0] - src[31], s_1_30 = src[1] - src[30], s_2_29 = src[2] - src[29], s_3_28 = src[3] - src[28];
        int s_4_27 = src[4] - src[27], s_5_26 = src[5] - src[26], s_6_25 = src[6] - src[25], s_7_24 = src[7] - src[24];
        int s_8_23 = src[8] - src[23], s_9_22 = src[9] - src[22], s_10_21 = src[10] - src[21];
        int s_11_20 = src[11] - src[20], s_12_19 = src[12] - src[19], s_13_18 = src[13] - src[18];
        int s_14_17 = src[14] - src[17], s_15_16 = src[15] - src[16];

        int a0 = a_0_31 + a_15_16, a1 = a_1_30 + a_14_17, a2 = a_2_29 + a_13_18, a3 = a_3_28 + a_12_19;
        int a4 = a_4_27 + a_11_20, a5 = a_5_26 + a_10_21, a6 = a_6_25 + a_9_22, a7 = a_7_24 + a_8_23;

        int s0 = a_0_31 - a_15_16, s1 = a_1_30 - a_14_17, s2 = a_2_29 - a_13_18, s3 = a_3_28 - a_12_19;
        int s4 = a_4_27 - a_11_20, s5 = a_5_26 - a_10_21, s6 = a_6_25 - a_9_22, s7 = a_7_24 - a_8_23;

        dst[0]   = (64*a0 + 64*a1 + 64*a2 + 64*a3 + 64*a4 + 64*a5 + 64*a6 + 64*a7 + add) >> shift;
        dst[128] = (89*a0 + 75*a1 + 50*a2 + 18*a3 - 18*a4 - 50*a5 - 75*a6 - 89*a7 + add) >> shift;
        dst[256] = (83*a0 + 36*a1 - 36*a2 - 83*a3 - 83*a4 - 36*a5 + 36*a6 + 83*a7 + add) >> shift;
        dst[384] = (75*a0 - 18*a1 - 89*a2 - 50*a3 + 50*a4 + 89*a5 + 18*a6 - 75*a7 + add) >> shift;
        dst[512] = (64*a0 - 64*a1 - 64*a2 + 64*a3 + 64*a4 - 64*a5 - 64*a6 + 64*a7 + add) >> shift;
        dst[640] = (50*a0 - 89*a1 + 18*a2 + 75*a3 - 75*a4 - 18*a5 + 89*a6 - 50*a7 + add) >> shift;
        dst[768] = (36*a0 - 83*a1 + 83*a2 - 36*a3 - 36*a4 + 83*a5 - 83*a6 + 36*a7 + add) >> shift;
        dst[896] = (18*a0 - 50*a1 + 75*a2 - 89*a3 + 89*a4 - 75*a5 + 50*a6 - 18*a7 + add) >> shift;

        dst[64]  = (90*s0 + 87*s1 + 80*s2 + 70*s3 + 57*s4 + 43*s5 + 25*s6 + 9*s7 + add) >> shift;
        dst[192] = (87*s0 + 57*s1 + 9*s2 - 43*s3 - 80*s4 - 90*s5 - 70*s6 - 25*s7 + add) >> shift;
        dst[320] = (80*s0 + 9*s1 - 70*s2 - 87*s3 - 25*s4 + 57*s5 + 90*s6 + 43*s7 + add) >> shift;
        dst[448] = (70*s0 - 43*s1 - 87*s2 + 9*s3 + 90*s4 + 25*s5 - 80*s6 - 57*s7 + add) >> shift;
        dst[576] = (57*s0 - 80*s1 - 25*s2 + 90*s3 - 9*s4 - 87*s5 + 43*s6 + 70*s7 + add) >> shift;
        dst[704] = (43*s0 - 90*s1 + 57*s2 + 25*s3 - 87*s4 + 70*s5 + 9*s6 - 80*s7 + add) >> shift;
        dst[832] = (25*s0 - 70*s1 + 90*s2 - 80*s3 + 43*s4 + 9*s5 - 57*s6 + 87*s7 + add) >> shift;
        dst[960] = (9*s0 - 25*s1 + 43*s2 - 57*s3 + 70*s4 - 80*s5 + 87*s6 - 90*s7 + add) >> shift;

        dst[32] = (90*s_0_31 + 90*s_1_30 + 88*s_2_29 + 85*s_3_28 + 82*s_4_27 + 78*s_5_26 + 73*s_6_25 + 67*s_7_24 +
                   61*s_8_23 + 54*s_9_22 + 46*s_10_21 + 38*s_11_20 + 31*s_12_19 + 22*s_13_18 + 13*s_14_17 +
                   4*s_15_16 + add) >> shift;

        dst[96] = (90*s_0_31 + 82*s_1_30 + 67*s_2_29 + 46*s_3_28 + 22*s_4_27 - 4*s_5_26 - 31*s_6_25 - 54*s_7_24 -
                   73*s_8_23 - 85*s_9_22 - 90*s_10_21 - 88*s_11_20 - 78*s_12_19 - 61*s_13_18 - 38*s_14_17 -
                   13*s_15_16 + add) >> shift;

        dst[160] = (88*s_0_31 + 67*s_1_30 + 31*s_2_29 - 13*s_3_28 - 54*s_4_27 - 82*s_5_26 - 90*s_6_25 - 78*s_7_24 -
                    46*s_8_23 - 4*s_9_22 + 38*s_10_21 + 73*s_11_20 + 90*s_12_19 + 85*s_13_18 + 61*s_14_17 +
                    22*s_15_16 + add) >> shift;

        dst[224] = (85*s_0_31 + 46*s_1_30 - 13*s_2_29 - 67*s_3_28 - 90*s_4_27 - 73*s_5_26 - 22*s_6_25 + 38*s_7_24 +
                    82*s_8_23 + 88*s_9_22 + 54*s_10_21 - 4*s_11_20 - 61*s_12_19 - 90*s_13_18 - 78*s_14_17 -
                    31*s_15_16 + add) >> shift;

        dst[288] = (82*s_0_31 + 22*s_1_30 - 54*s_2_29 - 90*s_3_28 - 61*s_4_27 + 13*s_5_26 + 78*s_6_25 + 85*s_7_24 +
                    31*s_8_23 - 46*s_9_22 - 90*s_10_21 - 67*s_11_20 + 4*s_12_19 + 73*s_13_18 + 88*s_14_17 +
                    38*s_15_16 + add) >> shift;

        dst[352] = (78*s_0_31 - 4*s_1_30 - 82*s_2_29 - 73*s_3_28 + 13*s_4_27 + 85*s_5_26 + 67*s_6_25 - 22*s_7_24 -
                    88*s_8_23 - 61*s_9_22 + 31*s_10_21 + 90*s_11_20 + 54*s_12_19 - 38*s_13_18 - 90*s_14_17 -
                    46*s_15_16 +add) >> shift;

        dst[416] = (73*s_0_31 - 31*s_1_30 - 90*s_2_29 - 22*s_3_28 + 78*s_4_27 + 67*s_5_26 - 38*s_6_25 - 90*s_7_24 -
                    13*s_8_23 + 82*s_9_22 + 61*s_10_21 - 46*s_11_20 - 88*s_12_19 - 4*s_13_18 + 85*s_14_17 +
                    54*s_15_16 + add) >> shift;

        dst[480] = (67*s_0_31 - 54*s_1_30 - 78*s_2_29 + 38*s_3_28 + 85*s_4_27 - 22*s_5_26 - 90*s_6_25 + 4*s_7_24 +
                    90*s_8_23 + 13*s_9_22 - 88*s_10_21 - 31*s_11_20 + 82*s_12_19 + 46*s_13_18 - 73*s_14_17 -
                    61*s_15_16 + add) >> shift;

        dst[544] = (61*s_0_31 - 73*s_1_30 - 46*s_2_29 + 82*s_3_28 + 31*s_4_27 - 88*s_5_26 - 13*s_6_25 + 90*s_7_24 -
                    4*s_8_23 - 90*s_9_22 + 22*s_10_21 + 85*s_11_20 - 38*s_12_19 - 78*s_13_18 + 54*s_14_17 +
                    67*s_15_16 + add) >> shift;

        dst[608] = (54*s_0_31 - 85*s_1_30 - 4*s_2_29 + 88*s_3_28 - 46*s_4_27 - 61*s_5_26 + 82*s_6_25 + 13*s_7_24 -
                    90*s_8_23 + 38*s_9_22 + 67*s_10_21 - 78*s_11_20 - 22*s_12_19 + 90*s_13_18 - 31*s_14_17 -
                    73*s_15_16 + add) >> shift;

        dst[672] = (46*s_0_31 - 90*s_1_30 + 38*s_2_29 + 54*s_3_28 - 90*s_4_27 + 31*s_5_26 + 61*s_6_25 - 88*s_7_24 +
                    22*s_8_23 + 67*s_9_22 - 85*s_10_21 + 13*s_11_20 + 73*s_12_19 - 82*s_13_18 + 4*s_14_17 +
                    78*s_15_16 + add) >> shift;

        dst[736] = (38*s_0_31 - 88*s_1_30 + 73*s_2_29 - 4*s_3_28 - 67*s_4_27 + 90*s_5_26 - 46*s_6_25 - 31*s_7_24 +
                    85*s_8_23 - 78*s_9_22 + 13*s_10_21 + 61*s_11_20 - 90*s_12_19 + 54*s_13_18 + 22*s_14_17 -
                    82*s_15_16 + add) >> shift;

        dst[800] = (31*s_0_31 - 78*s_1_30 + 90*s_2_29 - 61*s_3_28 + 4*s_4_27 + 54*s_5_26 - 88*s_6_25 + 82*s_7_24 -
                    38*s_8_23 - 22*s_9_22 + 73*s_10_21 - 90*s_11_20 + 67*s_12_19 - 13*s_13_18 - 46*s_14_17 +
                    85*s_15_16 + add) >> shift;

        dst[864] = (22*s_0_31 - 61*s_1_30 + 85*s_2_29 - 90*s_3_28 + 73*s_4_27 - 38*s_5_26 - 4*s_6_25 + 46*s_7_24 -
                    78*s_8_23 + 90*s_9_22 - 82*s_10_21 + 54*s_11_20 - 13*s_12_19 - 31*s_13_18 + 67*s_14_17 -
                    88*s_15_16 + add) >> shift;

        dst[928] = (13*s_0_31 - 38*s_1_30 + 61*s_2_29 - 78*s_3_28 + 88*s_4_27 - 90*s_5_26 + 85*s_6_25 - 73*s_7_24 +
                    54*s_8_23 - 31*s_9_22 + 4*s_10_21 + 22*s_11_20 - 46*s_12_19 + 67*s_13_18 - 82*s_14_17 +
                    90*s_15_16 + add) >> shift;

        dst[992] = (4*s_0_31 - 13*s_1_30 + 22*s_2_29 - 31*s_3_28 + 38*s_4_27 - 46*s_5_26 + 54*s_6_25 - 61*s_7_24 +
                    67*s_8_23 - 73*s_9_22 + 78*s_10_21 - 82*s_11_20 + 85*s_12_19 - 88*s_13_18 + 90*s_14_17 -
                    90*s_15_16 + add) >> shift;
    }
}

static void venc_idct_32_1d(int16_t *dst, int16_t *src, int shift)
{
    int add = 1 << (shift - 1);

    for (int i = 0; i < 32; dst += 32, src++, i++)
    {
    	int s0 = 90*src[32] + 90*src[96] + 88*src[160] + 85*src[224] + 82*src[288] + 78*src[352] + 73*src[416] +
                 67*src[480] + 61*src[544] + 54*src[608] + 46*src[672] + 38*src[736] + 31*src[800] + 22*src[864] +
                 13*src[928] + 4*src[992];

    	int s1 = 90*src[32] + 82*src[96] + 67*src[160] + 46*src[224] + 22*src[288] - 4*src[352] - 31*src[416] -
                 54*src[480] - 73*src[544] - 85*src[608] - 90*src[672] - 88*src[736] - 78*src[800] - 61*src[864] -
                 38*src[928] - 13*src[992];

    	int s2 = 88*src[32] + 67*src[96] + 31*src[160] - 13*src[224] - 54*src[288] - 82*src[352] - 90*src[416] -
                 78*src[480] - 46*src[544] - 4*src[608] + 38*src[672] + 73*src[736] + 90*src[800] + 85*src[864] +
                 61*src[928] + 22*src[992];

    	int s3 = 85*src[32] + 46*src[96] - 13*src[160] - 67*src[224] - 90*src[288] - 73*src[352] - 22*src[416] +
                 38*src[480] + 82*src[544] + 88*src[608] + 54*src[672] - 4*src[736] - 61*src[800] - 90*src[864] -
                 78*src[928] - 31*src[992];

    	int s4 = 82*src[32] + 22*src[96] - 54*src[160] - 90*src[224] - 61*src[288] + 13*src[352] + 78*src[416] +
                 85*src[480] + 31*src[544] - 46*src[608] - 90*src[672] - 67*src[736] + 4*src[800] + 73*src[864] +
                 88*src[928] + 38*src[992];

    	int s5 = 78*src[32] - 4*src[96] - 82*src[160] - 73*src[224] + 13*src[288] + 85*src[352] + 67*src[416] -
                 22*src[480] - 88*src[544] - 61*src[608] + 31*src[672] + 90*src[736] + 54*src[800] - 38*src[864] -
                 90*src[928] - 46*src[992];

    	int s6 = 73*src[32] - 31*src[96] - 90*src[160] - 22*src[224] + 78*src[288] + 67*src[352] - 38*src[416] -
                 90*src[480] - 13*src[544] + 82*src[608] + 61*src[672] - 46*src[736] - 88*src[800] - 4*src[864] +
                 85*src[928] + 54*src[992];

    	int s7 = 67*src[32] - 54*src[96] - 78*src[160] + 38*src[224] + 85*src[288] - 22*src[352] - 90*src[416] +
                 4*src[480] + 90*src[544] + 13*src[608] - 88*src[672] - 31*src[736] + 82*src[800] + 46*src[864] -
                 73*src[928] - 61*src[992];

    	int s8 = 61*src[32] - 73*src[96] - 46*src[160] + 82*src[224] + 31*src[288] - 88*src[352] - 13*src[416] +
                 90*src[480] - 4*src[544] - 90*src[608] + 22*src[672] + 85*src[736] - 38*src[800] - 78*src[864] +
                 54*src[928] + 67*src[992];

    	int s9 = 54*src[32] - 85*src[96] - 4*src[160] + 88*src[224] - 46*src[288] - 61*src[352] + 82*src[416] +
                 13*src[480] - 90*src[544] + 38*src[608] + 67*src[672] - 78*src[736] - 22*src[800] + 90*src[864] -
                 31*src[928] - 73*src[992];

    	int s10 = 46*src[32] - 90*src[96] + 38*src[160] + 54*src[224] - 90*src[288] + 31*src[352] + 61*src[416] -
                  88*src[480] + 22*src[544] + 67*src[608] - 85*src[672] + 13*src[736] + 73*src[800] - 82*src[864] +
                  4*src[928] + 78*src[992];

    	int s11 = 38*src[32] - 88*src[96] + 73*src[160] - 4*src[224] - 67*src[288] + 90*src[352] - 46*src[416] -
                  31*src[480] + 85*src[544] - 78*src[608] + 13*src[672] + 61*src[736] - 90*src[800] + 54*src[864] +
                  22*src[928] - 82*src[992];

    	int s12 = 31*src[32] - 78*src[96] + 90*src[160] - 61*src[224] + 4*src[288] + 54*src[352] - 88*src[416] +
                  82*src[480] - 38*src[544] - 22*src[608] + 73*src[672] - 90*src[736] + 67*src[800] - 13*src[864] -
                  46*src[928] + 85*src[992];

    	int s13 = 22*src[32] - 61*src[96] + 85*src[160] - 90*src[224] + 73*src[288] - 38*src[352] - 4*src[416] +
                  46*src[480] - 78*src[544] + 90*src[608] - 82*src[672] + 54*src[736] - 13*src[800] - 31*src[864] +
                  67*src[928] - 88*src[992];

    	int s14 = 13*src[32] - 38*src[96] + 61*src[160] - 78*src[224] + 88*src[288] - 90*src[352] + 85*src[416] -
                  73*src[480] + 54*src[544] - 31*src[608] + 4*src[672] + 22*src[736] - 46*src[800] + 67*src[864] -
                  82*src[928] + 90*src[992];

    	int s15 = 4*src[32] - 13*src[96] + 22*src[160] - 31*src[224] + 38*src[288] - 46*src[352] + 54*src[416] -
                  61*src[480] + 67*src[544] - 73*src[608] + 78*src[672] - 82*src[736] + 85*src[800] - 88*src[864] +
                  90*src[928] - 90*src[992];

    	int e_0 = 90*src[64] + 87*src[192] + 80*src[320] + 70*src[448] +
                  57*src[576] + 43*src[704] + 25*src[832] + 9*src[960];
    	int e_1 = 87*src[64] + 57*src[192] + 9*src[320] - 43*src[448] -
                  80*src[576] - 90*src[704] - 70*src[832] - 25*src[960];
    	int e_2 = 80*src[64] + 9*src[192] - 70*src[320] - 87*src[448] -
                  25*src[576] + 57*src[704] + 90*src[832] + 43*src[960];
    	int e_3 = 70*src[64] - 43*src[192] - 87*src[320] + 9*src[448] +
                  90*src[576] + 25*src[704] - 80*src[832] - 57*src[960];
    	int e_4 = 57*src[64] - 80*src[192] - 25*src[320] + 90*src[448] -
                  9*src[576] - 87*src[704] + 43*src[832] + 70*src[960];
    	int e_5 = 43*src[64] - 90*src[192] + 57*src[320] + 25*src[448] -
                  87*src[576] + 70*src[704] + 9*src[832] - 80*src[960];
    	int e_6 = 25*src[64] - 70*src[192] + 90*src[320] - 80*src[448] +
                  43*src[576] + 9*src[704] - 57*src[832] + 87*src[960];
    	int e_7 = 9*src[64] - 25*src[192] + 43*src[320] - 57*src[448] +
                  70*src[576] - 80*src[704] + 87*src[832] - 90*src[960];

    	int a_0 = 64*src[0] + 64*src[512] + 83*src[256] + 36*src[768];
    	int a_1 = 64*src[0] - 64*src[512] + 36*src[256] - 83*src[768];
    	int a_2 = 64*src[0] - 64*src[512] - 36*src[256] + 83*src[768];
    	int a_3 = 64*src[0] + 64*src[512] - 83*src[256] - 36*src[768];
    	int a_4 = 89*src[128] + 75*src[384] + 50*src[640] + 18*src[896];
    	int a_5 = 75*src[128] - 18*src[384] - 89*src[640] - 50*src[896];
    	int a_6 = 50*src[128] - 89*src[384] + 18*src[640] + 75*src[896];
    	int a_7 = 18*src[128] - 50*src[384] + 75*src[640] - 89*src[896];

    	int a0 = a_0 + a_4 + e_0;
    	int a1 = a_1 + a_5 + e_1;
    	int a2 = a_2 + a_6 + e_2;
    	int a3 = a_3 + a_7 + e_3;
    	int a4 = a_3 - a_7 + e_4;
    	int a5 = a_2 - a_6 + e_5;
    	int a6 = a_1 - a_5 + e_6;
    	int a7 = a_0 - a_4 + e_7;
    	int a8 = a_0 - a_4 - e_7;
    	int a9 = a_1 - a_5 - e_6;
    	int a10 = a_2 - a_6 - e_5;
    	int a11 = a_3 - a_7 - e_4;
    	int a12 = a_3 + a_7 - e_3;
    	int a13 = a_2 + a_6 - e_2;
    	int a14 = a_1 + a_5 - e_1;
    	int a15 = a_0 + a_4 - e_0;

    	dst[0] = F265_CLAMP((a0 + s0 + add) >> shift, -32768, 32767);
    	dst[1] = F265_CLAMP((a1 + s1 + add) >> shift, -32768, 32767);
    	dst[2] = F265_CLAMP((a2 + s2 + add) >> shift, -32768, 32767);
     	dst[3] = F265_CLAMP((a3 + s3 + add) >> shift, -32768, 32767);
    	dst[4] = F265_CLAMP((a4 + s4 + add) >> shift, -32768, 32767);
    	dst[5] = F265_CLAMP((a5 + s5 + add) >> shift, -32768, 32767);
    	dst[6] = F265_CLAMP((a6 + s6 + add) >> shift, -32768, 32767);
    	dst[7] = F265_CLAMP((a7 + s7 + add) >> shift, -32768, 32767);
    	dst[8] = F265_CLAMP((a8 + s8 + add) >> shift, -32768, 32767);
    	dst[9] = F265_CLAMP((a9 + s9 + add) >> shift, -32768, 32767);
    	dst[10] = F265_CLAMP((a10 + s10 + add) >> shift, -32768, 32767);
    	dst[11] = F265_CLAMP((a11 + s11 + add) >> shift, -32768, 32767);
    	dst[12] = F265_CLAMP((a12 + s12 + add) >> shift, -32768, 32767);
    	dst[13] = F265_CLAMP((a13 + s13 + add) >> shift, -32768, 32767);
    	dst[14] = F265_CLAMP((a14 + s14 + add) >> shift, -32768, 32767);
    	dst[15] = F265_CLAMP((a15 + s15 + add) >> shift, -32768, 32767);

    	dst[16] = F265_CLAMP((a15 - s15 + add) >> shift, -32768, 32767);
    	dst[17] = F265_CLAMP((a14 - s14 + add) >> shift, -32768, 32767);
    	dst[18] = F265_CLAMP((a13 - s13 + add) >> shift, -32768, 32767);
    	dst[19] = F265_CLAMP((a12 - s12 + add) >> shift, -32768, 32767);
    	dst[20] = F265_CLAMP((a11 - s11 + add) >> shift, -32768, 32767);
    	dst[21] = F265_CLAMP((a10 - s10 + add) >> shift, -32768, 32767);
    	dst[22] = F265_CLAMP((a9 - s9 + add) >> shift, -32768, 32767);
    	dst[23] = F265_CLAMP((a8 - s8 + add) >> shift, -32768, 32767);
    	dst[24] = F265_CLAMP((a7 - s7 + add) >> shift, -32768, 32767);
    	dst[25] = F265_CLAMP((a6 - s6 + add) >> shift, -32768, 32767);
    	dst[26] = F265_CLAMP((a5 - s5 + add) >> shift, -32768, 32767);
    	dst[27] = F265_CLAMP((a4 - s4 + add) >> shift, -32768, 32767);
    	dst[28] = F265_CLAMP((a3 - s3 + add) >> shift, -32768, 32767);
    	dst[29] = F265_CLAMP((a2 - s2 + add) >> shift, -32768, 32767);
    	dst[30] = F265_CLAMP((a1 - s1 + add) >> shift, -32768, 32767);
    	dst[31] = F265_CLAMP((a0 - s0 + add) >> shift, -32768, 32767);
    }
}

// Define the DCT and IDCT C stubs. The DCT computes the coefficients from the
// source and the prediction. The IDCT adds the residual to the prediction at
// the destination. The coefficients stride is the block size.
#define DECLARE_DCT_LBD(LG_BS, NAME)\
void venc_dct_##NAME##_c(int16_t *dst, f265_pix *src, int src_stride, f265_pix *pred, int pred_stride,\
                         uint8_t *spill)\
{\
    int lg_bs = LG_BS, bd = 8;\
    int bs = 1<<lg_bs, bs2 = 1<<(lg_bs<<1);\
    int shift1 = lg_bs + bd - 9, shift2 = lg_bs + 6;\
    int16_t diff[bs2], tmp[bs2];\
    venc_get_block_residual(diff, src, src_stride, pred, pred_stride, bs);\
    venc_dct_##NAME##_1d(tmp, diff, shift1);\
    venc_dct_##NAME##_1d(dst, tmp, shift2);\
}\
void venc_idct_##NAME##_c(f265_pix *dst, int dst_stride, f265_pix *pred, int pred_stride, int16_t *coeffs,\
                          uint8_t *spill)\
{\
    int lg_bs = LG_BS, bd = 8;\
    int bs = 1<<lg_bs, bs2 = 1<<(lg_bs<<1);\
    int shift1 = 7;\
    int shift2 = 20 - bd;\
    int16_t tmp[bs2], idct[bs2];\
    venc_idct_##NAME##_1d(tmp, coeffs, shift1);\
    venc_idct_##NAME##_1d(idct, tmp, shift2);\
    venc_get_block_reconstruction(dst, dst_stride, pred, pred_stride, idct, bs, bd);\
}
DECLARE_DCT_LBD(2, dst);
DECLARE_DCT_LBD(2, 4);
DECLARE_DCT_LBD(3, 8);
DECLARE_DCT_LBD(4, 16);
DECLARE_DCT_LBD(5, 32);
#undef DECLARE_DCT_LBD

// Perform the 1D DCT/IDCT (or DST) on the source coefficients. 'idx_map'
// specifies the (i,j,k) indices used to index the source and DCT factor arrays
// (idx[0]*bs + idx[1]). The DCT multiplication result is shifted by 'shift' and
// clipped as requested. The destination is stored as i*bs+j (raster scan).
// LEGACY CODE.
void venc_do_dct_1d(int16_t *dst, int16_t *src, int lg_bs, int dst_flag, int idx_map[2][2], int shift, int clip_flag)
{
    const int8_t *mats[4] = { f265_dct_mat_4[0], f265_dct_mat_8[0], f265_dct_mat_16[0], f265_dct_mat_32[0] };
    int8_t *mat = dst_flag ? (int8_t*)f265_dst_mat : (int8_t*)mats[lg_bs-2];
    int bs = 1<<lg_bs;

    // Naive matrix-like multiply.
    for (int i = 0; i < bs; i++)
    {
        for (int j = 0; j < bs; j++)
        {
            int res = 0;
            for (int k = 0; k < bs; k++)
            {
                int idx[3] = { i, j, k };
                int src_idx = idx[idx_map[0][0]]*bs + idx[idx_map[0][1]];
                int mat_idx = idx[idx_map[1][0]]*bs + idx[idx_map[1][1]];
                res += src[src_idx] * mat[mat_idx];
            }
            res = (res + (1<<(shift-1)))>>shift;
            if (clip_flag) res = F265_CLAMP(res, -32768, 32767);
            dst[i*bs + j] = res;
        }
    }
}

// Quantize the coefficients. Return true if a coefficient is non-null.
int venc_quant_c(int16_t *dst, int16_t *src, int bs, int mult, int add, int shift)
{
    int ret = 0;
    for (int i = 0; i < bs*bs; i++)
    {
        // The deadzone for quantization is uneven so we must cater for the
        // sign.
        int coeff = src[i];
        int sign = (coeff < 0) ? -1 : 1;
        coeff = ((F265_ABS(coeff)*mult + add)>>shift)*sign;
        dst[i] = F265_CLAMP(coeff, -32768, 32767);
        ret |= !!coeff;
    }
    return ret;
}

// Dequantize the coefficients.
void venc_dequant_c(int16_t *dst, int16_t *src, int bs, int mult, int add, int shift)
{
    for (int i = 0; i < bs*bs; i++)
        dst[i] = F265_CLAMP((src[i]*mult + add)>>shift, -32768, 32767);
}

// Initialize the transform tree.
void venc_init_transform_tree(f265_enc_thread *t)
{
    f265_tt_enc *tt = &t->tt;
    tt->tn = t->tmap;
    tt->tb = t->tb;
    tt->sb = t->sb;
    tt->levels = t->levels;
}

// Get the subblock non-zero flags of a transform block.
static void venc_get_subblock_flags(f265_tb_enc *tb, int16_t *qc)
{
    int lg_bs = tb->lg_bs;
    int bs = 1<<lg_bs;
    int sb_width = bs>>2;

    // Zero the subblock flags.
    for (int i = 0; i < 3; i++) tb->nz_flags[i] = 0;

    // Compute the subblock flags in raster scan.
    uint64_t rs_nz_flags = 0;
    for (int y = 0; y < sb_width; y++)
        for (int x = 0; x < sb_width; x++)
        {
            int16_t coeffs[16];
            uint64_t nz_flag = 0;
            venc_copy_block_s16(coeffs, 4, qc + 4*(y*bs + x), bs, 4, 4);
            for (int i = 0; i < 16; i++) if (coeffs[i]) { nz_flag = 1; break; }
            rs_nz_flags |= nz_flag<<(y*sb_width + x);
        }

    // Compute the subblock flags in encoding order.
    const uint8_t *sb_scan = f265_scan_map_data + f265_scan_map_idx[lg_bs-2][tb->order];
    for (int i = 0; i < sb_width*sb_width; i++)
        tb->nz_flags[0] |= ((rs_nz_flags>>sb_scan[i])&1)<<i;

    // Compute the subblock flags of the neighbours.
    for (int y = 0; y < sb_width; y++)
        for (int x = 0; x < sb_width; x++)
        {
            int pos = y*sb_width + x;
            if (x < sb_width - 1) tb->nz_flags[1] |= ((rs_nz_flags>>(pos + 1))&1)<<pos;
            if (y < sb_width - 1) tb->nz_flags[2] |= ((rs_nz_flags>>(pos + sb_width))&1)<<pos;
        }
}

// Preprocess the quantized coefficients in a non-empty transform block for
// encoding.
void venc_preprocess_tb(f265_tt_enc *tt, int16_t *qc)
{
    f265_tb_enc *tb = tt->tb;
    f265_sb_enc *sb = tt->sb;
    int16_t *levels = tt->levels;
    int lg_bs = tb->lg_bs;
    int lg_sb = lg_bs-2;
    int bs = 1<<lg_bs;
    int nb_sb = 1<<(lg_sb<<1);
    const uint8_t *coeff_scan = f265_scan_map_data + f265_scan_map_idx[2][tb->order];
    const uint8_t *coeff_off = f265_scan_map_data + 256 + f265_scan_map_idx[lg_sb][tb->order];

    // Loop over the subblocks in encoding order.
    for (int sb_idx = 0; sb_idx < nb_sb; sb_idx++)
    {
        // Skip empty subblocks.
        if (!((tb->nz_flags[0]>>sb_idx)&1)) continue;

        // Extract the subblock coefficients in raster scan.
        int16_t coeffs[16];
        venc_copy_block_s16(coeffs, 4, qc + (coeff_off[sb_idx]<<2), bs, 4, 4);

        // Process the coefficients in encoding order.
        int nb_nz = 0, nz_flags = 0, signs = 0, gt1_flags = 0, gt2_flag = 0, remain_flags = 0;
        for (int i = 0; i < 16; i++)
        {
            // Get the next non-zero coefficient in encoding order.
            int coeff = coeffs[coeff_scan[i]];
            if (!coeff) continue;

            // The coefficient value is inferred when the actual level is lower
            // than the base level.
            int level = F265_ABS(coeff);
            int base_level = 1;

            // Update the flags.
            nz_flags |= 1<<i;
            signs |= (coeff<0)<<nb_nz;
            if (nb_nz < 8)
            {
                base_level++;
                int gt1_flag = level>1;
                if (gt1_flag && !gt1_flags) { base_level++; gt2_flag = level>2; }
                gt1_flags |= gt1_flag<<nb_nz;
            }
            remain_flags |= (level>=base_level)<<i;
            nb_nz++;

            // Write the coefficient level tentatively.
            levels[i] = level;
        }

        // Reverse the signs for bypass encoding.
        {
            int tmp = signs;
            signs = 0;
            for (int i = 0; i < nb_nz; i++) signs |= ((tmp>>(nb_nz-i-1))&1)<<i;
        }

        // Commit the subblock and the levels.
        sb->nz_flags = nz_flags;
        sb->signs = signs;
        sb->remain_flags = remain_flags;
        sb->gt1_flags = gt1_flags;
        sb->packed_data = nb_nz | (gt2_flag<<5);
        sb++;
        if (remain_flags) levels += 16;
    }

    // Add a null subblock entry when the first subblock is empty because it is
    // encoded unconditionally.
    if (!((tb->nz_flags[0]>>(nb_sb-1))&1)) { f265_sb_enc z = {0}; *sb++ = z; }

    tt->tb++;
    tt->sb = sb;
    tt->levels = levels;
}

// Reconstruct a block from its prediction and return the transform block
// non-zero flag. At this time we pass a single parameter to control its
// behavior. The interface will change once we know what we actually need.
//
// Assuming low-bit depth at this time.
int venc_rec_block(f265_rec_params *rp)
{
    f265_enc_thread *t = rp->t;

    // Not storing the reconstruction is no longer supported.
    assert(rp->rec);

    int lg_bs = rp->lg_bs;
    int bs = 1<<lg_bs;
    int dst_flag = rp->dst_flag;
    int dct_idx = lg_bs-2 + 4*dst_flag;
    int order = rp->order;
    int zero_flag = rp->zero_flag;
    int qp = rp->qp;
    int bd = rp->bd;
    f265_tt_enc *tt = rp->tt;

    // DCT coefficients, quantized coefficients, dequantized coefficients.
    F265_ALIGN64 int16_t dct[32*32], quant[32*32], dequant[32*32];

    // No residual, the reconstruction is the prediction.
    if (zero_flag)
    {
        venc_copy_block(rp->rec, rp->rec_stride, rp->pred, rp->pred_stride, bs, bs);
        return 0;
    }

    // DCT.
    venc_dct[dct_idx](dct, rp->src, rp->src_stride, rp->pred, rp->pred_stride, t->store);

    // Quantization.
    int quant_mult = rp->quant_mult = f265_quant_mult[qp%6];
    int quant_shift = rp->quant_shift = 29 + qp/6 - bd - lg_bs;
    int quant_add = rp->quant_add = rp->rdoq_flag ? 1 << (quant_shift - 1) :
                                                    (rp->iframe_flag ? 171 : 85) << (quant_shift - 9);
    int nz_flag;
    if (rp->rdoq_flag) nz_flag = venc_do_rdoq(rp, quant, dct);
    else nz_flag = venc_quant[lg_bs-2](quant, dct, bs, quant_mult, quant_add, quant_shift);

    // No residual, the reconstruction in the prediction.
    if (!nz_flag)
    {
        venc_copy_block(rp->rec, rp->rec_stride, rp->pred, rp->pred_stride, bs, bs);
        return 0;
    }

    // Dequantization.
    int dequant_mult = f265_dequant_mult[qp%6] << (qp/6);
    int dequant_shift = bd + lg_bs - 9;
    int dequant_add = 1<<(dequant_shift-1);
    venc_dequant[lg_bs-2](dequant, quant, bs, dequant_mult, dequant_add, dequant_shift);

    // Inverse DCT.
    venc_idct[dct_idx](rp->rec, rp->rec_stride, rp->pred, rp->pred_stride, dequant, t->store);

    // Kludge to compute the distortion prior to clamping the samples.
    if (rp->compute_ssd_from_idct2)
    {
        int idct_idx_map[2][2] = { { 2, 0 }, { 2, 1 } };
        int idct1_shift = 7;
        int idct2_shift = 20 - bd;
        F265_ALIGN64 int16_t diff[32*32], idct1[32*32], idct2[32*32];
        venc_do_dct_1d(idct1, dequant, lg_bs, dst_flag, idct_idx_map, idct1_shift, 1);
        venc_do_dct_1d(idct2, idct1, lg_bs, dst_flag, idct_idx_map, idct2_shift, 0);
        venc_get_block_residual(diff, rp->src, rp->src_stride, rp->pred, rp->pred_stride, bs);
        rp->idct2_dist = venc_ssd16(diff, bs, idct2, bs, bs, bs, 0, rp->bd);
    }

    // Update the transform tree.
    if (tt && nz_flag)
    {
        // Set up the transform block.
        f265_tb_enc *tb = tt->tb;
        tb->lg_bs = lg_bs;
        tb->order = order;

        if (!rp->rdoq_flag)
            tb->sign_hiding_flags = 0;

        // Get the subblock flags. In assembly this should be done during the
        // quantization.
        venc_get_subblock_flags(tb, quant);

        // Preprocess the coefficients. Should be done with one assembly
        // function.
        venc_preprocess_tb(tt, quant);
    }

    return 1;
}

// Reconstruct the transform block. Return the non-zero flag.
int venc_rec_tb(f265_enc_thread *t, f265_pix *pred, int pred_stride, int comp, int lg_bs, int dst_flag, int order,
                int zero_flag, int ct_ox, int ct_oy, int depth, int intra_flag, int final_enc_flag)
{
    int chroma_flag = !!comp;
    int bd = t->enc->gd.bit_depth[chroma_flag];
    int qp = t->qp[chroma_flag];
    int stride = t->me.ref_stride;
    int plane_off = venc_get_ctb_block_plane_off(t, comp, ct_ox, ct_oy);
    f265_rec_params rp;
    rp.t = t;
    rp.src = t->src_frame->src_planes[comp] + plane_off;
    rp.src_stride = stride;
    rp.pred = pred;
    rp.pred_stride = pred_stride;
    rp.rec = t->src_frame->rec_planes[comp ? 3+comp : 0] + plane_off;
    rp.rec_stride = stride;
    rp.tt = &t->tt;
    rp.cabac_contexts = t->cbs.contexts;
    rp.comp = comp;
    rp.lg_bs = lg_bs;
    rp.tb_depth = depth;
    rp.dst_flag = dst_flag;
    rp.order = order;
    rp.zero_flag = zero_flag;
    rp.sign_hiding_flag = F265_GET_FLAG(t->enc->gd.eflags, F265_PF_SIGN_HIDING);
    rp.rdoq_flag = final_enc_flag & F265_GET_FLAG(t->enc->gd.eflags, F265_PF_RDOQ);
    rp.qp = qp;
    rp.bd = bd;
    rp.iframe_flag = t->src_frame->frame_type == F265_FRAME_I;
    rp.dump_flag = 0;
    rp.intra_flag = intra_flag;
    rp.lambda = t->hm_lambda[chroma_flag];

    // Kludge to compute the distortion prior to clamping the samples.
    rp.compute_ssd_from_idct2 = t->an.compute_ssd_from_idct2;
    rp.idct2_dist = 0;

    int nz_flag = venc_rec_block(&rp);

    // Kludge to compute the distortion prior to clamping the samples.
    t->an.idct2_dist = rp.idct2_dist;

    return nz_flag;
}

// Encode the CB and update the CABAC contexts.
static int venc_rdoq_enc_tt(f265_enc_thread *t, f265_pix *pred, int pred_stride, int comp, int lg_bs, int dst_flag,
                            int order, int zero_flag, int ct_ox, int ct_oy, int depth, int intra_flag)
{
    // Save the TT state. This is needed to update the CABAC context.
    f265_tt_enc pre_process_tt = t->tt;

    // Encode the CB.
    int yuv_flag = venc_rec_tb(t, pred, pred_stride, comp, lg_bs, dst_flag, order, zero_flag, ct_ox, ct_oy, depth,
                               intra_flag, 1);

    // In HM, intra contexts are updated only for the luma block.
    if (intra_flag && comp) return yuv_flag;

    #ifdef VAN_TRACE_SYNTAX
    // Disable the trace of the RDOQ process.
    venc_trace_syntax_flag = 0;
    #endif

    f265_tt_enc current_tt = t->tt;
    t->tt = pre_process_tt;

    // Update the RDOQ CABAC contexts.
    venc_encode_rdoq_sub_part(t, yuv_flag, comp, lg_bs, depth);

    // Restore the TT to the state it was after the recording process.
    t->tt = current_tt;

    #ifdef VAN_TRACE_SYNTAX
    venc_trace_syntax_flag = 1;
    #endif

    return yuv_flag;
}

// Reconstruct the intra transform block (as zero if requested). Return the non-zero flag.
int venc_rec_intra_tb(f265_enc_thread *t, int comp, int lg_bs, int mode, int zero_flag, int ct_ox, int ct_oy, int depth)
{
    f265_pix pred[32*32];
    int dst_flag, order;
    venc_predict_intra(t, pred, comp, lg_bs, mode, ct_ox, ct_oy);
    venc_get_intra_encode_flags(&dst_flag, &order, comp, lg_bs, mode);

    if (F265_GET_FLAG(t->enc->gd.eflags, F265_PF_RDOQ))
        return venc_rdoq_enc_tt(t, pred, 1<<lg_bs, comp, lg_bs, dst_flag, order, zero_flag, ct_ox, ct_oy, depth, 1);
    else
        return venc_rec_tb(t, pred, 1<<lg_bs, comp, lg_bs, dst_flag, order, zero_flag, ct_ox, ct_oy, depth, 1, 1);
}

// Reconstruct the intra transform tree. Return the YUV non-zero flags.
int venc_rec_intra_tt(f265_enc_thread *t, f265_cb *cb, int split_part_flag, int part_idx, int lg_bs,
                      int ct_ox, int ct_oy)
{
    int depth = cb->lg_bs - lg_bs;
    uint8_t *tn = t->tt.tn++;
    int split_tt_flag = !!(*tn&8);
    int yuv_mask = *tn&7;
    int yuv_flags = 0;

    // Split the transform tree.
    if (split_tt_flag)
        for (int i = 0, pi = part_idx, sbs = 1<<(lg_bs-1); i < 4; i++, pi += split_part_flag)
            yuv_flags |= venc_rec_intra_tt(t, cb, 0, pi, lg_bs-1, ct_ox + (i&1)*sbs, ct_oy + (i>>1)*sbs);

    // Reconstruct the current luma transform block.
    else
    {
        yuv_flags = venc_rec_intra_tb(t, 0, lg_bs, cb->intra_luma_mode[part_idx], !(yuv_mask&1), ct_ox, ct_oy, depth);
    }

    // Reconstruct the current chroma transform block.
    if (lg_bs == 3 || (lg_bs > 3 && !split_tt_flag))
        for (int i = 1; i < 3; i++)
            yuv_flags |= (venc_rec_intra_tb(t, i, lg_bs-1, cb->intra_chroma_mode, !(yuv_mask&(1<<i)),
                                            ct_ox>>1, ct_oy>>1, depth + (lg_bs == 3 && split_tt_flag))<<i);

    // Set the transform node.
    *tn = (split_tt_flag<<3)|yuv_flags;

    return yuv_flags;
}

// Reconstruct the inter transform block (as zero if requested). Return the non-zero flag.
int venc_rec_inter_tb(f265_enc_thread *t, f265_pix *pred, int pred_stride, int comp, int lg_bs, int dst_flag,
                      int order, int zero_flag, int ct_ox, int ct_oy, int depth)
{
    if (F265_GET_FLAG(t->enc->gd.eflags, F265_PF_RDOQ))
        return venc_rdoq_enc_tt(t, pred, pred_stride, comp, lg_bs, dst_flag, order, zero_flag, ct_ox, ct_oy, depth, 0);

    else
        return venc_rec_tb(t, pred, pred_stride, comp, lg_bs, dst_flag, order, zero_flag, ct_ox, ct_oy, depth, 0, 1);
}

// Reconstruct the inter transform tree. Return the YUV non-zero flags.
int venc_rec_inter_tt(f265_enc_thread *t, f265_cb *cb, f265_pix pred[3][64*64], int lg_bs, int cb_ox, int cb_oy)
{
    int depth = cb->lg_bs- lg_bs;
    uint8_t *tn = t->tt.tn++;
    int split_tt_flag = !!(*tn&8);
    int yuv_mask = *tn&7;
    int yuv_flags = 0;
    int cb_bs = 1<<cb->lg_bs;
    int ct_ox = cb->cb_off[0] + cb_ox, ct_oy = cb->cb_off[1] + cb_oy;

    // Split the transform tree.
    if (split_tt_flag)
    {
        const int nb_luma_ctx = F265_CO_CBF_CHROMA - F265_CO_CBF_LUMA;
        const int nb_chroma_ctx = F265_CO_MVD_GREATER0 - F265_CO_CBF_CHROMA;
        const int nb_total_ctx = nb_chroma_ctx + nb_luma_ctx;
        int rdoq_flag = F265_GET_FLAG(t->enc->gd.eflags, F265_PF_RDOQ);

        // Save CBF contexts.
        uint8_t ctx[nb_total_ctx];

        if (rdoq_flag)
            memcpy(ctx, t->cbs.contexts + F265_CO_CBF_LUMA, nb_total_ctx);

        for (int i = 0, sbs = 1<<(lg_bs-1); i < 4; i++)
            yuv_flags |= venc_rec_inter_tt(t, cb, pred, lg_bs-1, cb_ox + (i&1)*sbs, cb_oy + (i>>1)*sbs);

        // If the luma subtree is empty, undo the cbf_luma context update.
        if (rdoq_flag & !(yuv_flags & 1))
            memcpy(t->cbs.contexts + F265_CO_CBF_LUMA, ctx, nb_luma_ctx);

        // If the chroma subtrees is empty, undo the cbf_chroma context update.
        if (rdoq_flag & !(yuv_flags & 6))
            memcpy(t->cbs.contexts + F265_CO_CBF_CHROMA, ctx + nb_luma_ctx, nb_chroma_ctx);

    }

    // Reconstruct the current luma transform block.
    else
        yuv_flags = venc_rec_inter_tb(t, pred[0] + cb_bs*cb_oy + cb_ox, cb_bs, 0, lg_bs, 0, 0, !(yuv_mask&1), ct_ox,
                                      ct_oy, depth);

    // Reconstruct the current chroma transform block.
    if (lg_bs == 3 || (lg_bs > 3 && !split_tt_flag))
        for (int i = 1, sbs=cb_bs>>1; i < 3; i++)
            yuv_flags |= (venc_rec_inter_tb(t, pred[i] + ((sbs*cb_oy + cb_ox)>>1), sbs, i, lg_bs-1, 0, 0,
                                            !(yuv_mask&(1<<i)), ct_ox>>1, ct_oy>>1,
                                            depth + (lg_bs == 3 && split_tt_flag))<<i);

    // Set the transform node.
    *tn = (split_tt_flag<<3)|yuv_flags;

    return yuv_flags;
}

// Reconstruct an inter coding block.
void venc_rec_inter_cb(f265_enc_thread *t, f265_cb *cb)
{
    // Do the motion compensation.
    f265_pix pred[3][64*64];
    for (int i = 0; i < 3; i++) venc_mc_cb(t, cb, pred[i], 1<<(cb->lg_bs-!!i), i);

    // Skip mode. Copy the prediction.
    if (cb->flags&F265_CB_SKIP)
    {
        // FIXME: can result in 64x64 TB.
        *t->tt.tn++ = 0;
        for (int i = 0; i < 3; i++)
        {
            int chroma_flag = !!i;
            int bs = 1<<(cb->lg_bs-chroma_flag);
            int plane_off = venc_get_ctb_block_plane_off(t, i, cb->cb_off[0]>>chroma_flag, cb->cb_off[1]>>chroma_flag);
            f265_pix *rec = t->src_frame->rec_planes[i ? 3+i : 0] + plane_off;
            venc_copy_block(rec, t->me.ref_stride, pred[i], bs, bs, bs);
        }
    }

    // Reconstruct the transform tree.
    else
        venc_rec_inter_tt(t, cb, pred, cb->lg_bs, 0, 0);
}

#ifdef VAN_DUMP_INTRA
// Helper function.
static void venc_set_tt0_blueprint_helper(f265_enc_thread *t, int split_part_flag, int lg_bs)
{
    // For inter, we assume max_transform_hierarchy_depth_inter == 0. Thus, we
    // split the transform tree if we're splitting in partitions or if the
    // current transform block is larger than the maximum transform block.
    //
    // For intra, we split the transform tree if we're splitting in partitions
    // or if the current transform block is larger than the maximum transform
    // block.
    int split_tt_flag = split_part_flag | lg_bs>t->enc->gd.tb_range[1];

    // Do not zero blocks.
    int yuv_flags = 7;

    // Set the node.
    *t->tt.tn++ = (split_tt_flag<<3)|yuv_flags;

    // Split the transform tree.
    if (split_tt_flag)
        for (int i = 0; i < 4; i++)
            venc_set_tt0_blueprint_helper(t, 0, lg_bs-1);
}

// Create a transform tree blueprint for intra/inter transform depth 0.
void venc_set_tt0_blueprint(f265_enc_thread *t, f265_cb *cb)
{
    int split_part_flag = (cb->flags&F265_CB_INTRA) ? cb->intra_luma_mode[1] != -1 : cb->inter_part != F265_PART_UN;
    uint8_t *tn = t->tt.tn;
    venc_set_tt0_blueprint_helper(t, split_part_flag, cb->lg_bs);
    t->tt.tn = tn;
}
#endif

// Reconstruct the CB quadtree recursively. We update the YUV flags of the
// transform nodes but we do not change the layout of the transform tree.
void venc_rec_cb(f265_enc_thread *t, f265_cb *cb)
{
    // Skip absent CB.
    if (!(cb->flags&F265_CB_PRESENT)) return;

    // Handle split CB.
    if (cb->flags&F265_CB_SPLIT)
    {
        for (int i = 0; i < 4; i++) venc_rec_cb(t, t->cb + cb->child_idx + i);
        return;
    }

    // Set to emulate the tt0 layout. Bit-exact with HM in case of intra.
    #ifdef VAN_DUMP_INTRA
    venc_set_tt0_blueprint(t, cb);
    #endif

    // Inter.
    if (!(cb->flags&F265_CB_INTRA))
        venc_rec_inter_cb(t, cb);

    // No PCM support yet.
    else if (cb->intra_luma_mode[0] == 35)
    {
        printf("PCM support not implemented.\n");
        exit(1);
    }

    // Intra.
    else
        venc_rec_intra_tt(t, cb, cb->intra_luma_mode[1] != -1, 0, cb->lg_bs, cb->cb_off[0], cb->cb_off[1]);
}

// Reconstruct the CTB with the transform tree blueprint. The blocks identified
// as zero in the blueprint are encoded as zero. Otherwise, the blueprint is
// updated with the actual YUV flags. The split flags are unaffected.
void venc_rec_ctb(f265_enc_thread *t)
{
    venc_init_transform_tree(t);

    // Kludge: backup the raw CABAC object.
    f265_cabac_bs saved_cbs = t->cbs;

    venc_rec_cb(t, t->cb);

    // Restore the CABAC object.
    t->cbs = saved_cbs;

    #ifdef VAN_VALIDATE_MODE_PRED
    venc_hm_validate(t);
    #endif
}


///////////////////////////////////////////////////////////////////////////////
// HM-compatible RDOQ implementation.
//
// FIXME: this code needs to be factorized.

// Clear the block specified.
static void venc_clear_block_s16(int16_t *dst, int dst_stride, int width, int height)
{
    for (int y = 0; y < height; y++, dst += dst_stride)
        for (int x = 0; x < width; x++)
            dst[x] = 0;
}

// Retrieve from table the CABAC entropy value of a bin.
static uint32_t venc_retrieve_context_bin_rdoq(uint8_t *contexts, int ctx_idx, int bin)
{
    // Retrieve the context following the context index.
    int context = contexts[ctx_idx];

    // Retrieve the entropy value following the context and the bin value.
    return f265_cabac_entropy_table[context^bin];
}

// Initialize the last significant position rate table.
static void venc_init_last_context(f265_rec_params *rp, int *last_bits_x, int *last_bits_y)
{
    int lg_bs = rp->lg_bs;
    int bs = 1<<lg_bs;
    int chroma_flag = !!rp->comp;
    uint8_t *contexts = rp->cabac_contexts;
    int bits, ctx;

    // Calculate the context index offset and the context shift.
    int last_ctx_idx = F265_CO_LAST_SIG_COEFF + (chroma_flag ? 15 : 3*(lg_bs-2) + ((lg_bs-1)>>2));
    int last_ctx_shift = chroma_flag ? lg_bs - 2 : (lg_bs+1)>>2;

    // For every bin needed to encode a 'X coordinate', compute its bit rate.
    for (bits = 0, ctx = 0; ctx < (f265_last_coeff_table[bs-1]&15); ctx++)
    {
        int last_ctx = last_ctx_idx + (ctx >> last_ctx_shift);
        last_bits_x[ctx] = bits + venc_retrieve_context_bin_rdoq(contexts, last_ctx, 0);
        bits += venc_retrieve_context_bin_rdoq(contexts, last_ctx, 1);
    }
    last_bits_x[ctx] = bits;

    // Move the context offset to the 'Y coordinate' index.
    last_ctx_idx += 18;

    // Same for Y.
    // FIXME: copy-pasted code.
    for (bits = 0, ctx = 0; ctx < (f265_last_coeff_table[bs-1]&15); ctx++)
    {
        int last_ctx = last_ctx_idx + (ctx >> last_ctx_shift);
        last_bits_y[ctx] = bits + venc_retrieve_context_bin_rdoq(contexts, last_ctx, 0);
        bits += venc_retrieve_context_bin_rdoq(contexts, last_ctx, 1);
    }
    last_bits_y[ctx] = bits;
}

// This function is called during the computation of the best RDO level. It
// returns the RD cost of coding this coefficient.
static uint32_t venc_retrieve_rate_rdoq(uint8_t *contexts, int coeff, int gt1_ctx, int gt2_ctx, int rice,
                                        int nb_gt1_flags, int nb_gt2_flags, int sign_flag)
{
    uint32_t rate;

    // If needed, start cost evaluation with an equiprobable cost to represent
    // the sign bit flag CABAC rate.
    rate = (sign_flag)? 32768 : 0;

    // Retrieve the inferred baseLevel from greater1/2 flags.
    uint32_t baselevel = (nb_gt1_flags < 8) ? (2 + (nb_gt2_flags < 1)) : 1;

    // If the coeff level cannot be inferred, calculate the remaining value to
    // encode.
    if (coeff >= baselevel)
    {
        uint32_t delta = coeff - baselevel;
        uint32_t length;

        // The rice parameter is high enough to encode the delta level.
        if (delta < (3 << rice))
        {
            length = delta>>rice;
            rate += (length+1+rice) << 15;
        }

        // The rice parameter is too low.
        else
        {
            length = rice;
            delta  = delta - (3 << rice);
            while (delta >= (1<<length))
                delta -=  (1<<(length++));

            rate += (3+length+1-rice+length) << 15;
        }

        // If less than eight greater1 flags were used, a greater1 flag is encoded.
        if (nb_gt1_flags < 8)
        {
            rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 1);

            // If less than one greater2 flag is used, a greater2 flag is encoded.
            if (nb_gt2_flags < 1)
                rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER2 + gt2_ctx, 1);
        }
    }

    // If the level is zero, there's no rate cost.
    else if (coeff == 0)
    {
        rate = 0;
    }

    // If the level is one, a greater1 flag is signalled as 0.
    else if (coeff == 1)
    {
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 0);
    }

    // If the level is two, a greater1 flag is signalled as 1 and the greater2
    // flag as 0.
    else if (coeff == 2)
    {
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 1);
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER2 + gt2_ctx, 0);
    }

    // Impossible case.
    else
    {
        assert(0);
    }

    return rate;
}

// This function is called during the computation of the best
// RDO level. It returns the RD cost of coding this coefficient.
static uint32_t venc_retrieve_rate_sbh(uint8_t *contexts, int coeff, int gt1_ctx, int gt2_ctx, int rice,
                                       int nb_gt1_flags, int nb_gt2_flags)
{
    // FIXME: copy-pasted code from above.

    // Lookup tables retrieved from HM.
    const uint32_t gorice_range[5] = {7, 14, 26, 46, 78};
    const uint32_t gorice_prefix_length[5] = { 8, 7, 6, 5, 4};

    uint32_t rate = 0;

    // Retrieve the inferred baseLevel from greater1/2 flags.
    uint32_t baselevel = (nb_gt1_flags < 8) ? (2 + (nb_gt2_flags < 1)) : 1;

    // If the coeff level cannot be inferred, calculate the remaining value to
    // encode.
    if (coeff >= baselevel)
    {
        uint32_t delta = coeff - baselevel;
        uint32_t max_vlc = gorice_range[rice];
        uint32_t exp_golomb = delta > max_vlc;

        if (exp_golomb)
        {
            coeff = delta - max_vlc;
            uint32_t exp = 1;
            for (uint32_t mav_val = 2; coeff >= mav_val; mav_val <<= 1, exp += 2) {}
            rate += exp << 15;
            delta = F265_MIN(delta, (max_vlc + 1));
        }

        uint16_t prefix_length = (delta >> rice) + 1;
        uint16_t nbr_bits = F265_MIN(prefix_length, gorice_prefix_length[rice]) + rice;

        rate += nbr_bits << 15;

        // If less than eight greater1 flags were used, a greater1 flag is encoded.
        if (nb_gt1_flags < 8)
        {
            rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 1);

            // If less than one greater2 flag is used, a greater2 flag is encoded.
            if (nb_gt2_flags < 1)
                rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER2 + gt2_ctx, 1);
        }
    }

    // If the level is zero, there's no rate cost.
    else if (coeff == 0)
    {
        rate = 0;
    }

    // If the level is one, a greater1 flag is signalled as 0.
    else if (coeff == 1)
    {
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 0);
    }

    // If the level is two, a greater1 flag is signalled as 1 and the greater2
    // flag as 0.
    else if (coeff == 2)
    {
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1 + gt1_ctx, 1);
        rate += venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER2 + gt2_ctx, 0);
    }

    // Impossible case.
    else
    {
        assert(0);
    }

    return rate;
}

// This function is called during the quantization process of a coefficient.
// It returns the level with the lowest RD cost for the current coefficient.
static int venc_encode_level(f265_rec_params *rp, int co_pxy, int32_t mcoeff, int coeff, int sig_ctx, int gt1_ctx,
                             int gt2_ctx, int rice, int nb_gt1_flags, int nb_gt2_flags, int last_coeff,
                             double *dist_cost_qc0, double *rate_cost_qc0, double *dist_cost_qc, double *rate_cost_qc,
                             double *sig_rate_cost_qc)
{
    int shift = rp->quant_shift;
    double temp = rp->quant_temp;
    uint8_t *contexts = rp->cabac_contexts;
    double lambda = rp->lambda;

    // Set default values as if the level is zero-ified.
    double sig_flag_rate_cost_qc = 0;
    double best_rd_cost_qc = DBL_MAX;
    int best_level = 0;

    // Set its costs as zero-ified coefficient.
    rate_cost_qc0[co_pxy] = venc_retrieve_context_bin_rdoq(contexts, F265_CO_SIG_COEFF + sig_ctx, 0);
    rate_cost_qc[co_pxy]  = rate_cost_qc0[co_pxy];
    dist_cost_qc[co_pxy]  = dist_cost_qc0[co_pxy];

    // Not the last coefficient of the TB and its level is lower than three.
    if ((!last_coeff)&&(coeff<3))
    {
        // Use the zero-ified cost as the best RD cost.
        best_rd_cost_qc = dist_cost_qc0[co_pxy] + lambda*rate_cost_qc0[co_pxy];

        // If the quantized coefficient is zero, don't process it.
        if (coeff == 0) return 0;
    }

    // If coeff isn't LastSigCoeff, there is a significance flag.
    if (!last_coeff)
    {
        // Set default sig rate with not significant flag.
        sig_rate_cost_qc[co_pxy] = venc_retrieve_context_bin_rdoq(contexts, F265_CO_SIG_COEFF + sig_ctx, 0);

        // Retrieve significant flag for further processing.
        sig_flag_rate_cost_qc = venc_retrieve_context_bin_rdoq(contexts, F265_CO_SIG_COEFF + sig_ctx, 1);
    }

    // For each level from "coeff" to "min_level", calculate its distortion and rate costs.
    int min_level = ((coeff>1) ? (coeff-1) : 1);
    for (int curr_level = coeff; curr_level >= min_level; curr_level--)
    {
        // Calculate the distortion between the original level and the current
        // level.
        double delta = (mcoeff - (curr_level << shift));
        double curr_dist_cost_qc = delta * delta * temp;

        // Calculate the rate of the current level.
        double curr_rate_cost_qc =
            sig_flag_rate_cost_qc +
            venc_retrieve_rate_rdoq(contexts, curr_level, gt1_ctx, gt2_ctx, rice, nb_gt1_flags, nb_gt2_flags, 1);

        // Sum these costs as the current coefficient RD cost.
        double curr_rd_cost_qc = curr_dist_cost_qc + lambda*curr_rate_cost_qc;

        // If the current coefficient cost is better than the previous candidate cost,
        // remember its level and cost as the best coefficient level.
        if (curr_rd_cost_qc < best_rd_cost_qc)
        {
            best_level = curr_level;
            dist_cost_qc[co_pxy] = curr_dist_cost_qc;
            rate_cost_qc[co_pxy] = curr_rate_cost_qc;
            sig_rate_cost_qc[co_pxy] = sig_flag_rate_cost_qc;
            best_rd_cost_qc = curr_rd_cost_qc;
        }
    }

    return best_level;
}

// Perform the quantization/dequantization on the coefficients. Return true if a
// coefficient is non-null.
// FIXME: split this function in sub-functions.
int venc_do_rdoq(f265_rec_params *rp, int16_t *dst, int16_t *src)
{
    int lg_bs = rp->lg_bs;
    int bs = 1<<lg_bs;
    int lg_sb = lg_bs-2;
    int sb_width = bs>>2;
    int nb_sb = 1<<(lg_sb<<1);
    int chroma_flag = !!rp->comp;
    uint8_t *contexts = rp->cabac_contexts;
    int order = rp->order;
    int qp = rp->qp;
    int mult = rp->quant_mult;
    int add = rp->quant_add;
    int shift = rp->quant_shift;

    // Set the quant_temp (HM compatibility).
    {
        int bd = 8;
        int quant_scale = 1 << (15 - ((15-bd-lg_bs)<<1));
        rp->quant_temp = (double)quant_scale / (double)mult / (double)mult;
    }
    double lambda = rp->lambda;
    double temp = rp->quant_temp;

    // Sum of RD cost of the chosen coefficients for this transform block.
    double base_rd_cost = 0;

    // Sum of distortion of every zeroified coefficient for this transform block.
    double sum_dist_cost_qc0 = 0;

    // General rule:      rd_cost = dist_cost + lambda*rate_cost.
    // dist_cost_qc0[]    contains the distortion cost of zero-ified coefficient level.
    // rate_cost_qc0[]    contains the rate cost of zero-ified coefficient level.
    // dist_cost_qc[]     contains the distortion cost of quantified 'qc' coefficient level.
    // rate_cost_qc[]     contains the rate cost of quantified 'qc' coefficient level.
    // sig_rate_cost_qc[] contains the significant flag rate cost of quantified 'qc' coefficient level.
    F265_ALIGN64 double dist_cost_qc0[32*32];
    F265_ALIGN64 double rate_cost_qc0[32*32];
    F265_ALIGN64 double dist_cost_qc[32*32];
    F265_ALIGN64 double rate_cost_qc[32*32];
    F265_ALIGN64 double sig_rate_cost_qc[32*32];
    memset(dist_cost_qc0, 0, sizeof(dist_cost_qc0));
    memset(rate_cost_qc0, 0, sizeof(rate_cost_qc0));
    memset(dist_cost_qc, 0, sizeof(dist_cost_qc));
    memset(rate_cost_qc, 0, sizeof(rate_cost_qc));
    memset(sig_rate_cost_qc, 0, sizeof(sig_rate_cost_qc));

    // dist_cost_sb0[]    contains the sum of distortion cost of zero-ified sub-block.
    // rate_cost_sb0[]    contains the rate cost of zero-ified sub-block.
    // dist_cost_sb[]     contains the sum of distortion cost of quantified sub-block.
    // rate_cost_sb[]     contains the sun of rate cost of quantified sub-block.
    // sig_rate_cost_sb[] contains the significant flag rate cost of quantified sub-block.
    F265_ALIGN64 double dist_cost_sb0[8*8];
    F265_ALIGN64 double rate_cost_sb0[8*8];
    F265_ALIGN64 double dist_cost_sb[8*8];
    F265_ALIGN64 double rate_cost_sb[8*8];
    F265_ALIGN64 double sig_rate_cost_sb[8*8];
    memset(dist_cost_sb0, 0, sizeof(dist_cost_sb0));
    memset(rate_cost_sb0, 0, sizeof(rate_cost_sb0));
    memset(dist_cost_sb, 0, sizeof(dist_cost_sb));
    memset(rate_cost_sb, 0, sizeof(rate_cost_sb));
    memset(sig_rate_cost_sb, 0, sizeof(sig_rate_cost_sb));

    // SBH stuff.
    // delta_rate_up[]          contains the rate cost to increment a coefficient level.
    // delta_rate_down[]        contains the rate cost to decrement a coefficient level.
    // delta_sig_rate[]         contains the rate cost of a significant coefficient flag.
    // delta_dist[]             contains the distortion between the quantized coefficient and the source coefficient.
    int delta_rate_up[32*32];
    int delta_rate_down[32*32];
    int delta_sig_rate[32*32];
    int delta_dist[32*32];
    memset(delta_rate_up, 0, sizeof(delta_rate_up));
    memset(delta_rate_down, 0, sizeof(delta_rate_down));
    memset(delta_sig_rate, 0, sizeof(delta_sig_rate));
    memset(delta_dist, 0, sizeof(delta_dist));

    // coeff_offset  (co_oxy) specify, inside the sub-block, the position of current coefficient.
    // sub_pos       (sb_pxy) specify, inside the transform block, the position of the first coefficient of a sub-block.
    // sub_offset    (sb_oxy) specify, inside the transform block, the position of the current sub-block.
    const uint8_t *coeff_offset = f265_scan_map_data       + f265_scan_map_idx[2][order];
    const uint8_t *sub_pos      = f265_scan_map_data + 256 + f265_scan_map_idx[lg_sb][order];
    const uint8_t *sub_offset   = f265_scan_map_data       + f265_scan_map_idx[lg_sb][order];

    // Sub-block significance flag.
    uint64_t nz_sb_flags = 0;

    // Position the [best] last coefficient and its sub-block.
    int last_co_pxy = -1, last_sb_idx = -1;
    int best_lco_idx = 0, best_last_sb_idx = -1;

    // Multiplied coefficient level before downscaling.
    int32_t mc[32*32];

    // For every coefficient in the transform block:
    // - compute its multiplied level.
    // - compute the distortion cost of zero-ified coefficient level.
    // - compute the quantized coefficient level.
    for (int i = 0; i < bs*bs; i++)
    {
        int coeff = F265_ABS(src[i]);

        int m_coeff = coeff*mult;
        mc[i] = m_coeff;

        // Distortion evaluation of q -> 0.
        // (This 'q' is the multiplied coefficient, before losing precision with downscale).
        double curr_dist_cost_qc0 = (double)((double)m_coeff * (double)m_coeff) * temp;
        sum_dist_cost_qc0 += curr_dist_cost_qc0;
        dist_cost_qc0[i] = curr_dist_cost_qc0;

        int q_coeff = F265_CLAMP((m_coeff + add) >> shift, 0, 32767);
        dst[i] = q_coeff;
    }

    // Set the non-zero coefficient context information.
    const uint8_t *coeff_nz_ctx_idx = f265_coeff_nz_ctx_table;
    int coeff_nz_base_ctx = 0;
    if (lg_bs == 2) coeff_nz_ctx_idx += 4*16;
    else if (lg_bs == 3) coeff_nz_base_ctx += (!chroma_flag && order) ? 15 : 9;
    else coeff_nz_base_ctx += chroma_flag ? 12 : 21;

    // Initialize the last greater-than-1 context counter.
    int gt1_ctx_counter = 1;

    // Compute the best RDOQ level for every coefficient.
    // Loop over the sub-blocks in encoding order.
    for (int sb_idx = 0; sb_idx < nb_sb; sb_idx++)
    {
        // Extract the sub-block coefficients in raster scan.
        int sb_pxy = (sub_pos[sb_idx] << 2);
        int sb_oxy = sub_offset[sb_idx];
        int16_t coeffs[16];
        venc_copy_block_s16(coeffs, 4, dst + sb_pxy, bs, 4, 4);

        // Calculate the sub-block coordinates.
        uint8_t sb_ox = sb_oxy&((1<<lg_sb)-1);
        uint8_t sb_oy = sb_oxy>>lg_sb;

        // Update the context set for the greater-than-1 flags.
        int gt1_ctx_set = (((!chroma_flag && sb_pxy)<<1) + !gt1_ctx_counter)<<2;
        gt1_ctx_counter = 1;

        // Current rice parameter.
        int rice = 0;

        // Set the non-zero coefficient context information.
        int coeff_neighbour_idx = 0;
        {
            int sig_right_sb = 0;
            int sig_lower_sb = 0;

            if (sb_ox < sb_width - 1)
                sig_right_sb = (nz_sb_flags>>(sb_oxy + 1) & 1);

            if (sb_oy < sb_width - 1)
                sig_lower_sb = ((nz_sb_flags>>(sb_oxy + sb_width)) & 1);

            coeff_neighbour_idx = (sig_right_sb|sig_lower_sb<<1)<<4;
        }

        // Evaluate the coefficient RD costs (significant, greater1, greater2, remaining, sign).
        int nb_nz_flags = 0, nb_gt1_flags = 0, nb_gt2_flags = 0;
        int nz_sb_flag = 0, nz_qc_flags = 0;
        for (int coeff_idx = 0; coeff_idx < 16; coeff_idx++)
        {
            int co_oxy = coeff_offset[coeff_idx];
            int coeff = coeffs[co_oxy];
            int co_px = (sb_ox<<2) + (co_oxy&3);
            int co_py = (sb_oy<<2) + (co_oxy>>2);
            int co_pxy = (co_py<<lg_bs) + co_px;

            // Compute the distortion cost of a zero-ified sub-block
            // by adding every zero-ified coefficient distortion cost.
            dist_cost_sb0[sb_oxy] += dist_cost_qc0[co_pxy];

            // If this coefficient is the LastSigCoeff of the transform block, remember its position.
            int last_coeff = 0;
            if (coeff && (last_co_pxy < 0))
            {
                last_co_pxy = co_pxy;
                last_sb_idx = sb_idx;
                last_coeff = 1;
            }

            // The LastSigCoeff position was encountered.
            if (last_co_pxy >= 0)
            {
                // Calculate greater1/greater2 context indexes.
                int gt1_ctx = gt1_ctx_set + gt1_ctx_counter + 16*chroma_flag;
                int gt2_ctx = (gt1_ctx_set>>2) + !gt1_ctx_counter + 4*chroma_flag;

                // Calculate the significance context index.
                int sig_ctx = coeff_nz_base_ctx + coeff_nz_ctx_idx[coeff_neighbour_idx + co_oxy];
                if (sb_pxy) sig_ctx += chroma_flag ? 0 : 3;
                else if (!co_oxy) sig_ctx = 0;
                sig_ctx += 27*chroma_flag;
                if (last_coeff) sig_ctx = 0;

                // Find the best RDO-quantized level for this coefficient.
                int best_level = venc_encode_level(rp, co_pxy, mc[co_pxy], coeff, sig_ctx, gt1_ctx, gt2_ctx, rice,
                                                   nb_gt1_flags, nb_gt2_flags, last_coeff, dist_cost_qc0,
                                                   rate_cost_qc0, dist_cost_qc, rate_cost_qc, sig_rate_cost_qc);

                dst[co_pxy] = best_level;

                // Add the rate and distortion costs of this coefficient to the RD cost of this transform block.
                base_rd_cost += dist_cost_qc[co_pxy] + lambda*rate_cost_qc[co_pxy];

                // Add the costs of this coefficient to the costs of this sub-block.
                dist_cost_sb[sb_oxy] += dist_cost_qc[co_pxy];
                rate_cost_sb[sb_oxy] += rate_cost_qc[co_pxy];

                // Compute the coefficient distortion.
                delta_dist[co_pxy] = (mc[co_pxy] - (best_level<<shift)) >> (shift-8);

                // Compute the rate cost difference if we set this coefficient
                // to zero.
                if (!last_coeff)
                    delta_sig_rate[co_pxy] = venc_retrieve_context_bin_rdoq(contexts, F265_CO_SIG_COEFF+sig_ctx, 1) -
                                             venc_retrieve_context_bin_rdoq(contexts, F265_CO_SIG_COEFF+sig_ctx, 0);

                // Not a null coefficient.
                if (best_level > 0)
                {
                    // Get the current cost.
                    int curr_rate = venc_retrieve_rate_sbh(contexts, best_level, gt1_ctx, gt2_ctx, rice, nb_gt1_flags,
                                                           nb_gt2_flags);

                    // Get the cost difference if increasing or decreasing the
                    // coefficient level.
                    delta_rate_up[co_pxy] = venc_retrieve_rate_sbh(contexts, best_level+1, gt1_ctx, gt2_ctx, rice,
                                                                   nb_gt1_flags, nb_gt2_flags) - curr_rate;
                    delta_rate_down[co_pxy] = venc_retrieve_rate_sbh(contexts, best_level-1, gt1_ctx, gt2_ctx, rice,
                                                                     nb_gt1_flags, nb_gt2_flags) - curr_rate;
                }

                // Null coefficient. Get the cost of increasing the coefficient
                // level.
                else
                    delta_rate_up[co_pxy] = venc_retrieve_context_bin_rdoq(contexts, F265_CO_COEFF_GREATER1+gt1_ctx, 0);

                // Compute the value of the Sig/Greater1/Greater2 flags.
                int sig_flag = !!best_level;

                // Update the flags.
                //   nz_qc_flags  Bin array with the significance flag of the coefficients for this sub-block.
                //   nz_sb_flag   Non-zero sub-block flag.
                nz_qc_flags |= (sig_flag << co_oxy);
                nb_nz_flags += sig_flag;
                nz_sb_flag |= sig_flag;
                nz_sb_flags |= ((uint64_t)sig_flag)<<sb_oxy;

                // Update the counters.
                //   nb_gt1_flags   Number of greater1 flags for this sub-block.
                //   nb_gt2_flags   Number of greater2 flags for this sub-block.
                nb_gt1_flags += (nb_gt1_flags<8) ? (best_level>=1) : 0;
                nb_gt2_flags += (nb_gt2_flags<1) ? (best_level>=2) : 0;

                // If a greater1 flag is to be encoded, update the gt1 context.
                if (best_level>=1)
                {
                    // The current coefficient is the last of this sub-block.
                    if (!co_oxy)
                        gt1_ctx_counter = f265_gt1_ctx_counter_table[(nb_gt2_flags<<2)|gt1_ctx_counter];

                    else
                        gt1_ctx_counter = f265_gt1_ctx_counter_table[((best_level>1)<<2)|gt1_ctx_counter];
                }

                // Update the rice parameter if the level busts the threshold.
                rice += (best_level > (3<<rice)) && rice != 4;
            }

            else
            {
                // LastSigCoeff position wasn't encountered yet (trailing coefficient).
                // Add the zero-ified coefficient distortion cost to...
                // - the distortion cost of sub-block.
                // - the RD cost of the transform block.
                base_rd_cost += dist_cost_qc0[co_pxy];
                dist_cost_sb[sb_oxy] += dist_cost_qc0[co_pxy];
            }
        }

        // Evaluate the sub-block RD cost.
        // Check if zero-ified sub-block is better than keeping a quantized sub-block.
        if (last_sb_idx >= 0)
        {
            // Set the non-zero sub-block flag context index.
            int sb_ctx = 0;
            {
                int sig_right_sb = 0;
                int sig_lower_sb = 0;

                if (sb_ox < sb_width - 1)
                    sig_right_sb = (nz_sb_flags>>(sb_oxy + 1) & 1);

                if (sb_oy < sb_width - 1)
                    sig_lower_sb = (nz_sb_flags>>(sb_oxy + sb_width) & 1);

                sb_ctx = (sig_right_sb|sig_lower_sb) + 2*chroma_flag;
            }

            // Compute the rate cost of the sub-block flags.
            double curr_rate_cost_sb = venc_retrieve_context_bin_rdoq(contexts, F265_CO_CODED_SUB_BLOCK + sb_ctx, 1);
            double curr_rate_cost_sb0 = venc_retrieve_context_bin_rdoq(contexts, F265_CO_CODED_SUB_BLOCK + sb_ctx, 0);

            // Add these to sub-block rate costs.
            rate_cost_sb[sb_oxy] += curr_rate_cost_sb;
            rate_cost_sb0[sb_oxy] += curr_rate_cost_sb0;

            // Not the first sub-block.
            if (sb_pxy)
            {
                // This sub-block has a significant flag.
                if ((nz_sb_flags>>sb_oxy)&1)
                {
                    // Not the last sub-block.
                    if (sb_idx != last_sb_idx)
                    {
                        // No significant coefficient above position '0'.
                        if ((nz_qc_flags>>1) == 0)
                        {
                            // Remove the rate cost of the first coefficient of this SB.
                            base_rd_cost -= (lambda*sig_rate_cost_qc[sb_pxy]);
                            rate_cost_sb[sb_oxy] -= sig_rate_cost_qc[sb_pxy];
                        }

                        // Add the flag rate cost of this sub-block to the RD cost of the transform block.
                        // Note : RD cost of coefficient of this SB was previously set to base cost.
                        base_rd_cost += lambda*curr_rate_cost_sb;

                        // Compute the RD cost of this sub-block.
                        //   rd_cost_sb    RD cost with quantized coefficients.
                        //   rd_cost_sb0   RD cost with zero-ified coefficients.
                        double rd_cost_sb = dist_cost_sb[sb_oxy] + (lambda*rate_cost_sb[sb_oxy]);
                        double rd_cost_sb0 = dist_cost_sb0[sb_oxy] + (lambda*rate_cost_sb0[sb_oxy]);

                        // Remember the rate cost of the non-zero sub-block flag.
                        sig_rate_cost_sb[sb_oxy] = curr_rate_cost_sb;

                        // Zero-ified SB has better RD cost than quantized SB.
                        if (rd_cost_sb0 < rd_cost_sb)
                        {
                            // Clear the sub-block flag.
                            nz_sb_flags ^= ((uint64_t)1<<sb_oxy);

                            // Update the base cost by replacing quantized cost by zero-ified cost.
                            base_rd_cost -= rd_cost_sb;
                            base_rd_cost += rd_cost_sb0;

                            // Replace SB costs by zero-ified costs.
                            rate_cost_sb[sb_oxy] = rate_cost_sb0[sb_oxy];
                            dist_cost_sb[sb_oxy] = dist_cost_sb0[sb_oxy];
                            sig_rate_cost_sb[sb_oxy] = curr_rate_cost_sb0;

                            // Zero-ified the coefficients of this sub-block.
                            venc_clear_block_s16(dst + (sb_pxy), bs, 4, 4);
                        }
                    }
                }

                // No significant flag.
                else
                {
                    // Add the rate of non-significant SB flag.
                    base_rd_cost += lambda*curr_rate_cost_sb0;

                    // Remove the rate of coefficient flags of this sub-block
                    // (less the significant SB flag previously added).
                    base_rd_cost -= lambda*(rate_cost_sb[sb_oxy] - curr_rate_cost_sb);

                    // Remember the rate cost of the not significant sub-block flag.
                    sig_rate_cost_sb[sb_oxy] = curr_rate_cost_sb0;

                    // Clear the nz flag.
                    nz_sb_flags &= ~((uint64_t)1<<sb_oxy);
                }
            }

            // First sub-block.
            else
            {
                // The flag is inferred (no rate cost) for the last sub-block.
                nz_sb_flags |= 1;
            }
        }
    }

    // If LastSigCoeff position has not been encountered (i.e. the whole block is zero), return to caller function.
    if (last_co_pxy < 0) return 0;

    // Evaluate the RD cost of the CBF. Identify the context for the "coded
    // block flag" (CBF), calculate the RD cost of a non-significant TB (i.e.
    // cost of zero-ified transform block with a not significant coding block
    // flag), and calculate the RD cost of a significant transform block.
    double rd_cost_tb0 = 0;
    {
        uint8_t tb_depth = rp->tb_depth;
        int cbf_ctx = 0;

        // Inter luma root coding block.
        if (!(rp->intra_flag) && !(chroma_flag) && !tb_depth)
        {
            rd_cost_tb0 = sum_dist_cost_qc0 + lambda*venc_retrieve_context_bin_rdoq(contexts, F265_CO_ROOT_CBF, 0);
            base_rd_cost += lambda*venc_retrieve_context_bin_rdoq(contexts, F265_CO_ROOT_CBF, 1);
        }
        else
        {
            if (chroma_flag)
                cbf_ctx = F265_CO_CBF_CHROMA + tb_depth;
            else
                cbf_ctx = F265_CO_CBF_LUMA + !tb_depth;

            rd_cost_tb0 = sum_dist_cost_qc0 + lambda*venc_retrieve_context_bin_rdoq(contexts, cbf_ctx, 0);
            base_rd_cost += lambda*venc_retrieve_context_bin_rdoq(contexts, cbf_ctx, 1);
        }
    }

    // Preprocessed arrays of encoding rate for the last coefficient positions.
    int last_bits_x[10];
    int last_bits_y[10];
    memset(last_bits_x, 0, sizeof(last_bits_x));
    memset(last_bits_y, 0, sizeof(last_bits_y));

    // Compute the encoding rates for the last coefficient positions.
    venc_init_last_context(rp, last_bits_x, last_bits_y);

    double best_rd_cost = rd_cost_tb0;
    int best_coeff_pxy = -1;
    int best_last_found = 0;

    // Pass every sub-block.
    for (int sb_idx = last_sb_idx; sb_idx < nb_sb; sb_idx++)
    {
        // Extract the sub-block coefficients in raster scan.
        int sb_pxy = sub_pos[sb_idx];
        int sb_oxy = sub_offset[sb_idx];

        // This sub-block is now considered to be the last sub-block.
        // Remove the rate cost of its NZ sub-block flag (inferred flag).
        base_rd_cost -= lambda*sig_rate_cost_sb[sb_oxy];

        // Process the SB if it is significant.
        if ((nz_sb_flags>>sb_oxy)&1)
        {
            uint8_t sb_ox = sb_oxy&((1<<lg_sb)-1);
            uint8_t sb_oy = sb_oxy>>lg_sb;

            int16_t coeffs[16];
            venc_copy_block_s16(coeffs, 4, dst + (sb_pxy<<2), bs, 4, 4);

            // For every coefficient in this sub-block.
            for (int coeff_idx = 0; coeff_idx < 16; coeff_idx++)
            {
                int co_oxy = coeff_offset[coeff_idx];
                int coeff = coeffs[co_oxy];
                int co_px = (sb_ox<<2) + (co_oxy&3);
                int co_py = (sb_oy<<2) + (co_oxy>>2);
                int co_pxy = (co_py<<lg_bs) + co_px;

                // Significant coefficient.
                if (coeff)
                {
                    // Swap the order for vertical.
                    if (unlikely(order == 2)) F265_SWAP(int, co_px, co_py);

                    // Encode 'X' coordinate.
                    int tmp;
                    tmp = f265_last_coeff_table[co_px];
                    int x_prefix_ones = tmp&15;
                    int x_suffix_len = tmp>>4;

                    // Encode 'Y' coordinate.
                    tmp = f265_last_coeff_table[co_py];
                    int y_prefix_ones = tmp&15;
                    int y_suffix_len = tmp>>4;

                    // Compute the rate cost of this position as the last significant position.
                    double last_rate_cost;
                    last_rate_cost = last_bits_x[x_prefix_ones] + last_bits_y[y_prefix_ones];
                    last_rate_cost += 0x8000*(x_suffix_len + y_suffix_len);

                    // Compute the total RD cost of this TB.
                    //   Add the base cost minus the sig_coeff_rate
                    //   (this flag is now inferred by the last rate position).
                    //   Add the last position rate cost.
                    double total_rd_cost = base_rd_cost + lambda*(last_rate_cost - sig_rate_cost_qc[co_pxy]);

                    // If the TB with the current LastSigCoeff position has a better cost than the previous
                    // best RD cost, remember this coefficient position and the current TB cost for further processing.
                    if (total_rd_cost < best_rd_cost)
                    {
                        best_rd_cost = total_rd_cost;
                        best_coeff_pxy = co_pxy;
                        best_lco_idx = coeff_idx;
                        best_last_sb_idx = sb_idx;
                    }

                    // If the current coefficient level is higher than 1, stop the search for the best last position.
                    if (coeff > 1)
                    {
                        best_last_found = 1;
                        break;
                    }

                    // Set the quantized coefficient RD cost to the zero-ified
                    // coefficient RD cost.
                    base_rd_cost -= (dist_cost_qc[co_pxy] + lambda*rate_cost_qc[co_pxy]);
                    base_rd_cost += dist_cost_qc0[co_pxy];
                }
                else
                {
                    // This non-significant coeff becomes a trailing coeff.
                    // Remove its rate from the base cost.
                    base_rd_cost -= lambda*rate_cost_qc[co_pxy];
                }

            }

            if (best_last_found)
                break;
        }
    }

    // Number of NZ coeffs in this TB after the RDOQ process.
    int nb_nz_tb = 0;

    // If a best last position exists, do this process.
    if (best_coeff_pxy >= 0)
    {
        // Clear every sub-block after the best last SB.
        for (int sb_idx = last_sb_idx; sb_idx < best_last_sb_idx; sb_idx++)
        {
            int sb_pxy = sub_pos[sb_idx];
            venc_clear_block_s16(dst + (sb_pxy<<2), bs, 4, 4);
        }

        // Restore the sign to every coefficient in the remaining sub-blocks.
        int16_t src_coeffs[16], dst_coeffs[16];
        int sb_pxy = 0;
        for (int sb_idx = nb_sb-1; sb_idx >= best_last_sb_idx; sb_idx--)
        {
            sb_pxy = sub_pos[sb_idx];
            venc_copy_block_s16(src_coeffs, 4, src + (sb_pxy<<2), bs, 4, 4);
            venc_copy_block_s16(dst_coeffs, 4, dst + (sb_pxy<<2), bs, 4, 4);

            for (int coeff_idx = 0; coeff_idx < 16; coeff_idx++)
            {
                int level = dst_coeffs[coeff_idx];
                if (level)
                {
                    dst_coeffs[coeff_idx] = (src_coeffs[coeff_idx]<0)? -level : level;
                    nb_nz_tb += level;
                }
            }
            venc_copy_block_s16(dst + (sb_pxy<<2), bs, dst_coeffs, 4, 4, 4);
        }

        // Clear coefficients after the best last coeff in the best last SB.
        for (int coeff_idx = 0; coeff_idx < best_lco_idx; coeff_idx++)
        {
            int co_oxy = coeff_offset[coeff_idx];
            if (dst_coeffs[co_oxy])
            {
                nb_nz_tb -= F265_ABS(dst_coeffs[co_oxy]);
                dst_coeffs[co_oxy] = 0;
            }
        }
        venc_copy_block_s16(dst + (sb_pxy<<2), bs, dst_coeffs, 4, 4, 4);
    }

    // The best TB is a zero-ified TB.
    else
    {
        memset(dst, 0, sizeof(int16_t) * bs*bs);
    }

    // Clear the sign hiding flags. They will be set if required.
    uint64_t shb_flags = 0;

    // Sign bit hiding (SBH) is enabled by the encoder and the level sum is > 1.
    if (rp->sign_hiding_flag && nb_nz_tb > 1)
    {
        // Calculate the RD factor to convert a distortion as a rate cost.
        int64_t rd_factor = (int64_t)
                            (f265_dequant_mult[qp%6] * f265_dequant_mult[qp%6] * (1<<(2*(qp/6)))
                             / lambda / 16  + 0.5);

        int32_t last_sb = -1;
        int32_t sum_levels = 0;

        // For every sub-block.
        for (int sb_idx = 0; sb_idx < nb_sb; sb_idx++)
        {
            int sb_pxy = (sub_pos[sb_idx] << 2);
            int sb_oxy = sub_offset[sb_idx];
            uint8_t sb_ox = sb_oxy&((1<<lg_sb)-1);
            uint8_t sb_oy = sb_oxy>>lg_sb;

            int first_coeff_idx = 0, last_coeff_idx = 16;
            sum_levels = 0;

            int16_t coeffs[16];
            venc_copy_block_s16(coeffs, 4, dst + sb_pxy, bs, 4, 4);

            // Find the last significant coefficient index.
            for (int coeff_idx = 0; coeff_idx <= 15; coeff_idx++)
            {
                int co_oxy = coeff_offset[coeff_idx];
                if (coeffs[co_oxy])
                {
                    last_coeff_idx = coeff_idx;
                    break;
                }
            }

            // Find the first significant coefficient index.
            for (int coeff_idx = 15; coeff_idx >= 0; coeff_idx--)
            {
                int co_oxy = coeff_offset[coeff_idx];
                if (coeffs[co_oxy])
                {
                    first_coeff_idx = coeff_idx;
                    break;
                }
            }

            // Sum every level between the first and the last coefficients.
            for (int coeff_idx = first_coeff_idx; coeff_idx >= last_coeff_idx; coeff_idx--)
            {
                int co_oxy = coeff_offset[coeff_idx];
                sum_levels += coeffs[co_oxy];
            }

            // Set a flag is current sub-block is the last sub-block.
            if (last_coeff_idx <= 15 && last_sb == -1)
                last_sb = 1;

            // The distance between the first and the last coefficients is >= 4.
            if (first_coeff_idx-last_coeff_idx >= 4)
            {
                // Mark the transform block as using SBH.
                shb_flags |= ((uint64_t)1<<sb_oxy);

                // If the parity of the sum doesn't match with the sign of first coefficient, tune it.
                // Else, continue with the next sub-block.
                int signbit = (coeffs[coeff_offset[first_coeff_idx]]>0?0:1);
                if (signbit != (sum_levels&0x1))
                {
                    int64_t min_cost = LLONG_MAX, cur_cost = LLONG_MAX;
                    int min_co_pxy = -1, final_change = 0, cur_change = 0;
                    int sign_pos = -1;

                    for (int coeff_idx = (last_sb==1?last_coeff_idx:0); coeff_idx < 16; coeff_idx++)
                    {
                        int co_oxy = coeff_offset[coeff_idx];
                        int coeff = coeffs[co_oxy];
                        int co_px = (sb_ox<<2) + (co_oxy&3);
                        int co_py = (sb_oy<<2) + (co_oxy>>2);
                        int co_pxy = (co_py<<lg_bs) + co_px;

                        // Non-null coefficient.
                        if (coeff != 0)
                        {
                            // Evaluate the rate cost to adjust this coefficient
                            // by either incrementing or decrementing.
                            int64_t cost_up   = rd_factor * (-delta_dist[co_pxy]) + delta_rate_up[co_pxy];
                            int64_t cost_down = rd_factor * ( delta_dist[co_pxy]) + delta_rate_down[co_pxy] -
                                                (abs(coeff)==1?((1<<15)+delta_sig_rate[co_pxy]):0);

                            // If this coefficient is the LastSigCoeff of the TB and its absolute level is 1, reduce the
                            // decrementing rate cost of 0x20000 (advantage this candidate).
                            if (last_sb == 1 && last_coeff_idx == coeff_idx && abs(coeff) == 1)
                                cost_down -= (4<<15);

                            // Incrementing the cost is cheaper.
                            if (cost_up < cost_down)
                            {
                                cur_cost = cost_up;
                                cur_change =  1;
                            }

                            // Decrementing the cost is cheaper.
                            else
                            {
                                cur_cost = cost_down;
                                cur_change = -1;

                                // If this coefficient is the first non-zero of this sub-block and its
                                // absolute level is 1, eliminate this candidate.
                                if (coeff_idx == first_coeff_idx && abs(coeff) == 1)
                                    cur_cost = LLONG_MAX;
                            }
                        }

                        // Null coefficient.
                        else
                        {
                            // Evaluate the rate cost to adjust this coefficient
                            // by incrementing, and remember it as cur_cost.
                            cur_cost = rd_factor * (-(abs(delta_dist[co_pxy]))) + (1<<15) + delta_rate_up[co_pxy] +
                                       delta_sig_rate[co_pxy];
                            cur_change = 1;

                            // If this coefficient is a leading zero of this sub-block, and the sign of the original
                            // coefficient is different from the sign of the parity sum (i.e. the tuning would modify
                            // the coefficient sign, creating a greater distortion), eliminate this candidate.
                            if (coeff_idx > first_coeff_idx)
                            {
                                int this_sign_bit = (src[co_pxy]>=0?0:1);
                                if (this_sign_bit != signbit)
                                    cur_cost = LLONG_MAX;
                            }
                        }

                        // If the current cost is better than the previous cost, remember the current cost,
                        // the candidate position, and the change to do (increment or decrement).
                        if (cur_cost < min_cost)
                        {
                            min_cost = cur_cost;
                            final_change = cur_change;
                            min_co_pxy = co_pxy;
                            sign_pos = coeff_idx;
                        }
                    }

                    // If the coefficient level to tune is saturated, force a decrement tuning.
                    if (dst[min_co_pxy] == 32767 || dst[min_co_pxy] == -32768)
                        final_change = -1;

                    // Apply the adjustment in the absolute way (respecting its sign).
                    if (src[min_co_pxy] >= 0)
                        dst[min_co_pxy] += final_change;
                    else
                        dst[min_co_pxy] -= final_change;

                    // Special case: if the modified coefficient is at the range limit, and its new value is 0,
                    // the SBH range check might now be invalid. Re-test it.
                    if (sign_pos == first_coeff_idx || sign_pos == last_coeff_idx)
                    {
                        // Find the last significant coefficient index.
                        for (int coeff_idx = 0; coeff_idx < 16; coeff_idx++)
                        {
                            int co_oxy = coeff_offset[coeff_idx];
                            int co_px = (sb_ox<<2) + (co_oxy&3);
                            int co_py = (sb_oy<<2) + (co_oxy>>2);
                            int co_pxy = (co_py<<lg_bs) + co_px;
                            if (dst[co_pxy])
                            {
                                last_coeff_idx = coeff_idx;
                                break;
                            }
                        }

                        // Find the first significant coefficient index.
                        for (int coeff_idx = 15; coeff_idx >= 0; coeff_idx--)
                        {
                            int co_oxy = coeff_offset[coeff_idx];
                            int co_px = (sb_ox<<2) + (co_oxy&3);
                            int co_py = (sb_oy<<2) + (co_oxy>>2);
                            int co_pxy = (co_py<<lg_bs) + co_px;
                            if (dst[co_pxy])
                            {
                                first_coeff_idx = coeff_idx;
                                break;
                            }
                        }

                        if (first_coeff_idx - last_coeff_idx < 4)
                            shb_flags &= ~((uint64_t)1 << sb_oxy);
                    }
                }
            }

            if (last_sb == 1)
                last_sb = 0;
        }
    }

    // Save the SBH flags.
    rp->tt->tb->sign_hiding_flags = shb_flags;

    // Return the significance flag of this TB.
    return !!nb_nz_tb;
}
///////////////////////////////////////////////////////////////////////////////

