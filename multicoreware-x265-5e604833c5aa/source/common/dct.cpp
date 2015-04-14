/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Mandar Gurav <mandar@multicorewareinc.com>
 *          Deepthi Devaki Akkoorath <deepthidevaki@multicorewareinc.com>
 *          Mahesh Pittala <mahesh@multicorewareinc.com>
 *          Rajesh Paulraj <rajesh@multicorewareinc.com>
 *          Min Chen <min.chen@multicorewareinc.com>
 *          Praveen Kumar Tiwari <praveen@multicorewareinc.com>
 *          Nabajit Deka <nabajit@multicorewareinc.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com.
 *****************************************************************************/

#include "common.h"
#include "primitives.h"

using namespace x265;

#if _MSC_VER
#pragma warning(disable: 4127) // conditional expression is constant, typical for templated functions
#endif

namespace {
// anonymous file-static namespace

// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm
// give identical results
void fastForwardDst(int16_t *block, int16_t *coeff, int shift)  // input block, output coeff
{
    int c[4];
    int rnd_factor = 1 << (shift - 1);

    for (int i = 0; i < 4; i++)
    {
        // Intermediate Variables
        c[0] = block[4 * i + 0] + block[4 * i + 3];
        c[1] = block[4 * i + 1] + block[4 * i + 3];
        c[2] = block[4 * i + 0] - block[4 * i + 1];
        c[3] = 74 * block[4 * i + 2];

        coeff[i] =      (int16_t)((29 * c[0] + 55 * c[1]  + c[3] + rnd_factor) >> shift);
        coeff[4 + i] =  (int16_t)((74 * (block[4 * i + 0] + block[4 * i + 1] - block[4 * i + 3]) + rnd_factor) >> shift);
        coeff[8 + i] =  (int16_t)((29 * c[2] + 55 * c[0]  - c[3] + rnd_factor) >> shift);
        coeff[12 + i] = (int16_t)((55 * c[2] - 29 * c[1] + c[3] + rnd_factor) >> shift);
    }
}

void inversedst(int16_t *tmp, int16_t *block, int shift)  // input tmp, output block
{
    int i, c[4];
    int rnd_factor = 1 << (shift - 1);

    for (i = 0; i < 4; i++)
    {
        // Intermediate Variables
        c[0] = tmp[i] + tmp[8 + i];
        c[1] = tmp[8 + i] + tmp[12 + i];
        c[2] = tmp[i] - tmp[12 + i];
        c[3] = 74 * tmp[4 + i];

        block[4 * i + 0] = (int16_t)Clip3(-32768, 32767, (29 * c[0] + 55 * c[1]     + c[3]               + rnd_factor) >> shift);
        block[4 * i + 1] = (int16_t)Clip3(-32768, 32767, (55 * c[2] - 29 * c[1]     + c[3]               + rnd_factor) >> shift);
        block[4 * i + 2] = (int16_t)Clip3(-32768, 32767, (74 * (tmp[i] - tmp[8 + i]  + tmp[12 + i])      + rnd_factor) >> shift);
        block[4 * i + 3] = (int16_t)Clip3(-32768, 32767, (55 * c[0] + 29 * c[2]     - c[3]               + rnd_factor) >> shift);
    }
}

void partialButterfly16(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* E and O */
        for (k = 0; k < 8; k++)
        {
            E[k] = src[k] + src[15 - k];
            O[k] = src[k] - src[15 - k];
        }

        /* EE and EO */
        for (k = 0; k < 4; k++)
        {
            EE[k] = E[k] + E[7 - k];
            EO[k] = E[k] - E[7 - k];
        }

        /* EEE and EEO */
        EEE[0] = EE[0] + EE[3];
        EEO[0] = EE[0] - EE[3];
        EEE[1] = EE[1] + EE[2];
        EEO[1] = EE[1] - EE[2];

        dst[0] = (int16_t)((g_t16[0][0] * EEE[0] + g_t16[0][1] * EEE[1] + add) >> shift);
        dst[8 * line] = (int16_t)((g_t16[8][0] * EEE[0] + g_t16[8][1] * EEE[1] + add) >> shift);
        dst[4 * line] = (int16_t)((g_t16[4][0] * EEO[0] + g_t16[4][1] * EEO[1] + add) >> shift);
        dst[12 * line] = (int16_t)((g_t16[12][0] * EEO[0] + g_t16[12][1] * EEO[1] + add) >> shift);

        for (k = 2; k < 16; k += 4)
        {
            dst[k * line] = (int16_t)((g_t16[k][0] * EO[0] + g_t16[k][1] * EO[1] + g_t16[k][2] * EO[2] +
                                       g_t16[k][3] * EO[3] + add) >> shift);
        }

        for (k = 1; k < 16; k += 2)
        {
            dst[k * line] =  (int16_t)((g_t16[k][0] * O[0] + g_t16[k][1] * O[1] + g_t16[k][2] * O[2] + g_t16[k][3] * O[3] +
                                        g_t16[k][4] * O[4] + g_t16[k][5] * O[5] + g_t16[k][6] * O[6] + g_t16[k][7] * O[7] +
                                        add) >> shift);
        }

        src += 16;
        dst++;
    }
}

void partialButterfly32(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* E and O*/
        for (k = 0; k < 16; k++)
        {
            E[k] = src[k] + src[31 - k];
            O[k] = src[k] - src[31 - k];
        }

        /* EE and EO */
        for (k = 0; k < 8; k++)
        {
            EE[k] = E[k] + E[15 - k];
            EO[k] = E[k] - E[15 - k];
        }

        /* EEE and EEO */
        for (k = 0; k < 4; k++)
        {
            EEE[k] = EE[k] + EE[7 - k];
            EEO[k] = EE[k] - EE[7 - k];
        }

        /* EEEE and EEEO */
        EEEE[0] = EEE[0] + EEE[3];
        EEEO[0] = EEE[0] - EEE[3];
        EEEE[1] = EEE[1] + EEE[2];
        EEEO[1] = EEE[1] - EEE[2];

        dst[0] = (int16_t)((g_t32[0][0] * EEEE[0] + g_t32[0][1] * EEEE[1] + add) >> shift);
        dst[16 * line] = (int16_t)((g_t32[16][0] * EEEE[0] + g_t32[16][1] * EEEE[1] + add) >> shift);
        dst[8 * line] = (int16_t)((g_t32[8][0] * EEEO[0] + g_t32[8][1] * EEEO[1] + add) >> shift);
        dst[24 * line] = (int16_t)((g_t32[24][0] * EEEO[0] + g_t32[24][1] * EEEO[1] + add) >> shift);
        for (k = 4; k < 32; k += 8)
        {
            dst[k * line] = (int16_t)((g_t32[k][0] * EEO[0] + g_t32[k][1] * EEO[1] + g_t32[k][2] * EEO[2] +
                                       g_t32[k][3] * EEO[3] + add) >> shift);
        }

        for (k = 2; k < 32; k += 4)
        {
            dst[k * line] = (int16_t)((g_t32[k][0] * EO[0] + g_t32[k][1] * EO[1] + g_t32[k][2] * EO[2] +
                                       g_t32[k][3] * EO[3] + g_t32[k][4] * EO[4] + g_t32[k][5] * EO[5] +
                                       g_t32[k][6] * EO[6] + g_t32[k][7] * EO[7] + add) >> shift);
        }

        for (k = 1; k < 32; k += 2)
        {
            dst[k * line] = (int16_t)((g_t32[k][0] * O[0] + g_t32[k][1] * O[1] + g_t32[k][2] * O[2] + g_t32[k][3] * O[3] +
                                       g_t32[k][4] * O[4] + g_t32[k][5] * O[5] + g_t32[k][6] * O[6] + g_t32[k][7] * O[7] +
                                       g_t32[k][8] * O[8] + g_t32[k][9] * O[9] + g_t32[k][10] * O[10] + g_t32[k][11] *
                                       O[11] + g_t32[k][12] * O[12] + g_t32[k][13] * O[13] + g_t32[k][14] * O[14] +
                                       g_t32[k][15] * O[15] + add) >> shift);
        }

        src += 32;
        dst++;
    }
}

void partialButterfly8(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* E and O*/
        for (k = 0; k < 4; k++)
        {
            E[k] = src[k] + src[7 - k];
            O[k] = src[k] - src[7 - k];
        }

        /* EE and EO */
        EE[0] = E[0] + E[3];
        EO[0] = E[0] - E[3];
        EE[1] = E[1] + E[2];
        EO[1] = E[1] - E[2];

        dst[0] = (int16_t)((g_t8[0][0] * EE[0] + g_t8[0][1] * EE[1] + add) >> shift);
        dst[4 * line] = (int16_t)((g_t8[4][0] * EE[0] + g_t8[4][1] * EE[1] + add) >> shift);
        dst[2 * line] = (int16_t)((g_t8[2][0] * EO[0] + g_t8[2][1] * EO[1] + add) >> shift);
        dst[6 * line] = (int16_t)((g_t8[6][0] * EO[0] + g_t8[6][1] * EO[1] + add) >> shift);

        dst[line] = (int16_t)((g_t8[1][0] * O[0] + g_t8[1][1] * O[1] + g_t8[1][2] * O[2] + g_t8[1][3] * O[3] + add) >> shift);
        dst[3 * line] = (int16_t)((g_t8[3][0] * O[0] + g_t8[3][1] * O[1] + g_t8[3][2] * O[2] + g_t8[3][3] * O[3] + add) >> shift);
        dst[5 * line] = (int16_t)((g_t8[5][0] * O[0] + g_t8[5][1] * O[1] + g_t8[5][2] * O[2] + g_t8[5][3] * O[3] + add) >> shift);
        dst[7 * line] = (int16_t)((g_t8[7][0] * O[0] + g_t8[7][1] * O[1] + g_t8[7][2] * O[2] + g_t8[7][3] * O[3] + add) >> shift);

        src += 8;
        dst++;
    }
}

void partialButterflyInverse4(int16_t *src, int16_t *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        O[0] = g_t4[1][0] * src[line] + g_t4[3][0] * src[3 * line];
        O[1] = g_t4[1][1] * src[line] + g_t4[3][1] * src[3 * line];
        E[0] = g_t4[0][0] * src[0] + g_t4[2][0] * src[2 * line];
        E[1] = g_t4[0][1] * src[0] + g_t4[2][1] * src[2 * line];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        dst[0] = (int16_t)(Clip3(-32768, 32767, (E[0] + O[0] + add) >> shift));
        dst[1] = (int16_t)(Clip3(-32768, 32767, (E[1] + O[1] + add) >> shift));
        dst[2] = (int16_t)(Clip3(-32768, 32767, (E[1] - O[1] + add) >> shift));
        dst[3] = (int16_t)(Clip3(-32768, 32767, (E[0] - O[0] + add) >> shift));

        src++;
        dst += 4;
    }
}

void partialButterflyInverse8(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 4; k++)
        {
            O[k] = g_t8[1][k] * src[line] + g_t8[3][k] * src[3 * line] + g_t8[5][k] * src[5 * line] + g_t8[7][k] * src[7 * line];
        }

        EO[0] = g_t8[2][0] * src[2 * line] + g_t8[6][0] * src[6 * line];
        EO[1] = g_t8[2][1] * src[2 * line] + g_t8[6][1] * src[6 * line];
        EE[0] = g_t8[0][0] * src[0] + g_t8[4][0] * src[4 * line];
        EE[1] = g_t8[0][1] * src[0] + g_t8[4][1] * src[4 * line];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];
        for (k = 0; k < 4; k++)
        {
            dst[k] = (int16_t)Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 4] = (int16_t)Clip3(-32768, 32767, (E[3 - k] - O[3 - k] + add) >> shift);
        }

        src++;
        dst += 8;
    }
}

void partialButterflyInverse16(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 8; k++)
        {
            O[k] = g_t16[1][k] * src[line] + g_t16[3][k] * src[3 * line] + g_t16[5][k] * src[5 * line] + g_t16[7][k] * src[7 * line] +
                g_t16[9][k] * src[9 * line] + g_t16[11][k] * src[11 * line] + g_t16[13][k] * src[13 * line] + g_t16[15][k] * src[15 * line];
        }

        for (k = 0; k < 4; k++)
        {
            EO[k] = g_t16[2][k] * src[2 * line] + g_t16[6][k] * src[6 * line] + g_t16[10][k] * src[10 * line] + g_t16[14][k] * src[14 * line];
        }

        EEO[0] = g_t16[4][0] * src[4 * line] + g_t16[12][0] * src[12 * line];
        EEE[0] = g_t16[0][0] * src[0] + g_t16[8][0] * src[8 * line];
        EEO[1] = g_t16[4][1] * src[4 * line] + g_t16[12][1] * src[12 * line];
        EEE[1] = g_t16[0][1] * src[0] + g_t16[8][1] * src[8 * line];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        for (k = 0; k < 2; k++)
        {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 2] = EEE[1 - k] - EEO[1 - k];
        }

        for (k = 0; k < 4; k++)
        {
            E[k] = EE[k] + EO[k];
            E[k + 4] = EE[3 - k] - EO[3 - k];
        }

        for (k = 0; k < 8; k++)
        {
            dst[k]   = (int16_t)Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 8] = (int16_t)Clip3(-32768, 32767, (E[7 - k] - O[7 - k] + add) >> shift);
        }

        src++;
        dst += 16;
    }
}

void partialButterflyInverse32(int16_t *src, int16_t *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 16; k++)
        {
            O[k] = g_t32[1][k] * src[line] + g_t32[3][k] * src[3 * line] + g_t32[5][k] * src[5 * line] + g_t32[7][k] * src[7 * line] +
                g_t32[9][k] * src[9 * line] + g_t32[11][k] * src[11 * line] + g_t32[13][k] * src[13 * line] + g_t32[15][k] * src[15 * line] +
                g_t32[17][k] * src[17 * line] + g_t32[19][k] * src[19 * line] + g_t32[21][k] * src[21 * line] + g_t32[23][k] * src[23 * line] +
                g_t32[25][k] * src[25 * line] + g_t32[27][k] * src[27 * line] + g_t32[29][k] * src[29 * line] + g_t32[31][k] * src[31 * line];
        }

        for (k = 0; k < 8; k++)
        {
            EO[k] = g_t32[2][k] * src[2 * line] + g_t32[6][k] * src[6 * line] + g_t32[10][k] * src[10 * line] + g_t32[14][k] * src[14 * line] +
                g_t32[18][k] * src[18 * line] + g_t32[22][k] * src[22 * line] + g_t32[26][k] * src[26 * line] + g_t32[30][k] * src[30 * line];
        }

        for (k = 0; k < 4; k++)
        {
            EEO[k] = g_t32[4][k] * src[4 * line] + g_t32[12][k] * src[12 * line] + g_t32[20][k] * src[20 * line] + g_t32[28][k] * src[28 * line];
        }

        EEEO[0] = g_t32[8][0] * src[8 * line] + g_t32[24][0] * src[24 * line];
        EEEO[1] = g_t32[8][1] * src[8 * line] + g_t32[24][1] * src[24 * line];
        EEEE[0] = g_t32[0][0] * src[0] + g_t32[16][0] * src[16 * line];
        EEEE[1] = g_t32[0][1] * src[0] + g_t32[16][1] * src[16 * line];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];
        for (k = 0; k < 4; k++)
        {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 4] = EEE[3 - k] - EEO[3 - k];
        }

        for (k = 0; k < 8; k++)
        {
            E[k] = EE[k] + EO[k];
            E[k + 8] = EE[7 - k] - EO[7 - k];
        }

        for (k = 0; k < 16; k++)
        {
            dst[k] = (int16_t)Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 16] = (int16_t)Clip3(-32768, 32767, (E[15 - k] - O[15 - k] + add) >> shift);
        }

        src++;
        dst += 32;
    }
}

void partialButterfly4(int16_t *src, int16_t *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1 << (shift - 1);

    for (j = 0; j < line; j++)
    {
        /* E and O */
        E[0] = src[0] + src[3];
        O[0] = src[0] - src[3];
        E[1] = src[1] + src[2];
        O[1] = src[1] - src[2];

        dst[0] = (int16_t)((g_t4[0][0] * E[0] + g_t4[0][1] * E[1] + add) >> shift);
        dst[2 * line] = (int16_t)((g_t4[2][0] * E[0] + g_t4[2][1] * E[1] + add) >> shift);
        dst[line] = (int16_t)((g_t4[1][0] * O[0] + g_t4[1][1] * O[1] + add) >> shift);
        dst[3 * line] = (int16_t)((g_t4[3][0] * O[0] + g_t4[3][1] * O[1] + add) >> shift);

        src += 4;
        dst++;
    }
}

void dst4_c(int16_t *src, int32_t *dst, intptr_t stride)
{
    const int shift_1st = 1 + X265_DEPTH - 8;
    const int shift_2nd = 8;

    ALIGN_VAR_32(int16_t, coef[4 * 4]);
    ALIGN_VAR_32(int16_t, block[4 * 4]);

    for (int i = 0; i < 4; i++)
    {
        memcpy(&block[i * 4], &src[i * stride], 4 * sizeof(int16_t));
    }

    fastForwardDst(block, coef, shift_1st);
    fastForwardDst(coef, block, shift_2nd);

#define N (4)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dst[i * N + j] = block[i * N + j];
        }
    }

#undef N
}

void dct4_c(int16_t *src, int32_t *dst, intptr_t stride)
{
    const int shift_1st = 1 + X265_DEPTH - 8;
    const int shift_2nd = 8;

    ALIGN_VAR_32(int16_t, coef[4 * 4]);
    ALIGN_VAR_32(int16_t, block[4 * 4]);

    for (int i = 0; i < 4; i++)
    {
        memcpy(&block[i * 4], &src[i * stride], 4 * sizeof(int16_t));
    }

    partialButterfly4(block, coef, shift_1st, 4);
    partialButterfly4(coef, block, shift_2nd, 4);
#define N (4)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dst[i * N + j] = block[i * N + j];
        }
    }

#undef N
}

void dct8_c(int16_t *src, int32_t *dst, intptr_t stride)
{
    const int shift_1st = 2 + X265_DEPTH - 8;
    const int shift_2nd = 9;

    ALIGN_VAR_32(int16_t, coef[8 * 8]);
    ALIGN_VAR_32(int16_t, block[8 * 8]);

    for (int i = 0; i < 8; i++)
    {
        memcpy(&block[i * 8], &src[i * stride], 8 * sizeof(int16_t));
    }

    partialButterfly8(block, coef, shift_1st, 8);
    partialButterfly8(coef, block, shift_2nd, 8);

#define N (8)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dst[i * N + j] = block[i * N + j];
        }
    }

#undef N
}

void dct16_c(int16_t *src, int32_t *dst, intptr_t stride)
{
    const int shift_1st = 3 + X265_DEPTH - 8;
    const int shift_2nd = 10;

    ALIGN_VAR_32(int16_t, coef[16 * 16]);
    ALIGN_VAR_32(int16_t, block[16 * 16]);

    for (int i = 0; i < 16; i++)
    {
        memcpy(&block[i * 16], &src[i * stride], 16 * sizeof(int16_t));
    }

    partialButterfly16(block, coef, shift_1st, 16);
    partialButterfly16(coef, block, shift_2nd, 16);

#define N (16)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dst[i * N + j] = block[i * N + j];
        }
    }

#undef N
}

void dct32_c(int16_t *src, int32_t *dst, intptr_t stride)
{
    const int shift_1st = 4 + X265_DEPTH - 8;
    const int shift_2nd = 11;

    ALIGN_VAR_32(int16_t, coef[32 * 32]);
    ALIGN_VAR_32(int16_t, block[32 * 32]);

    for (int i = 0; i < 32; i++)
    {
        memcpy(&block[i * 32], &src[i * stride], 32 * sizeof(int16_t));
    }

    partialButterfly32(block, coef, shift_1st, 32);
    partialButterfly32(coef, block, shift_2nd, 32);

#define N (32)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dst[i * N + j] = block[i * N + j];
        }
    }

#undef N
}

void idst4_c(int32_t *src, int16_t *dst, intptr_t stride)
{
    const int shift_1st = 7;
    const int shift_2nd = 12 - (X265_DEPTH - 8);

    ALIGN_VAR_32(int16_t, coef[4 * 4]);
    ALIGN_VAR_32(int16_t, block[4 * 4]);

#define N (4)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            block[i * N + j] = (int16_t)src[i * N + j];
        }
    }

#undef N

    inversedst(block, coef, shift_1st); // Forward DST BY FAST ALGORITHM, block input, coef output
    inversedst(coef, block, shift_2nd); // Forward DST BY FAST ALGORITHM, coef input, coeff output

    for (int i = 0; i < 4; i++)
    {
        memcpy(&dst[i * stride], &block[i * 4], 4 * sizeof(int16_t));
    }
}

void idct4_c(int32_t *src, int16_t *dst, intptr_t stride)
{
    const int shift_1st = 7;
    const int shift_2nd = 12 - (X265_DEPTH - 8);

    ALIGN_VAR_32(int16_t, coef[4 * 4]);
    ALIGN_VAR_32(int16_t, block[4 * 4]);

#define N (4)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            block[i * N + j] = (int16_t)src[i * N + j];
        }
    }

#undef N

    partialButterflyInverse4(block, coef, shift_1st, 4); // Forward DST BY FAST ALGORITHM, block input, coef output
    partialButterflyInverse4(coef, block, shift_2nd, 4); // Forward DST BY FAST ALGORITHM, coef input, coeff output

    for (int i = 0; i < 4; i++)
    {
        memcpy(&dst[i * stride], &block[i * 4], 4 * sizeof(int16_t));
    }
}

void idct8_c(int32_t *src, int16_t *dst, intptr_t stride)
{
    const int shift_1st = 7;
    const int shift_2nd = 12 - (X265_DEPTH - 8);

    ALIGN_VAR_32(int16_t, coef[8 * 8]);
    ALIGN_VAR_32(int16_t, block[8 * 8]);

#define N (8)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            block[i * N + j] = (int16_t)src[i * N + j];
        }
    }

#undef N

    partialButterflyInverse8(block, coef, shift_1st, 8);
    partialButterflyInverse8(coef, block, shift_2nd, 8);
    for (int i = 0; i < 8; i++)
    {
        memcpy(&dst[i * stride], &block[i * 8], 8 * sizeof(int16_t));
    }
}

void idct16_c(int32_t *src, int16_t *dst, intptr_t stride)
{
    const int shift_1st = 7;
    const int shift_2nd = 12 - (X265_DEPTH - 8);

    ALIGN_VAR_32(int16_t, coef[16 * 16]);
    ALIGN_VAR_32(int16_t, block[16 * 16]);

#define N (16)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            block[i * N + j] = (int16_t)src[i * N + j];
        }
    }

#undef N

    partialButterflyInverse16(block, coef, shift_1st, 16);
    partialButterflyInverse16(coef, block, shift_2nd, 16);
    for (int i = 0; i < 16; i++)
    {
        memcpy(&dst[i * stride], &block[i * 16], 16 * sizeof(int16_t));
    }
}

void idct32_c(int32_t *src, int16_t *dst, intptr_t stride)
{
    const int shift_1st = 7;
    const int shift_2nd = 12 - (X265_DEPTH - 8);

    ALIGN_VAR_32(int16_t, coef[32 * 32]);
    ALIGN_VAR_32(int16_t, block[32 * 32]);

#define N (32)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            block[i * N + j] = (int16_t)src[i * N + j];
        }
    }

#undef N

    partialButterflyInverse32(block, coef, shift_1st, 32);
    partialButterflyInverse32(coef, block, shift_2nd, 32);

    for (int i = 0; i < 32; i++)
    {
        memcpy(&dst[i * stride], &block[i * 32], 32 * sizeof(int16_t));
    }
}

void dequant_normal_c(const int16_t* quantCoef, int32_t* coef, int num, int scale, int shift)
{
#if HIGH_BIT_DEPTH
    X265_CHECK(scale < 32768 || ((scale & 3) == 0 && shift > 2), "dequant invalid scale %d\n", scale);
#else
    // NOTE: maximum of scale is (72 * 256)
    X265_CHECK(scale < 32768, "dequant invalid scale %d\n", scale);
#endif
    X265_CHECK(num <= 32 * 32, "dequant num %d too large\n", num);
    X265_CHECK((num % 8) == 0, "dequant num %d not multiple of 8\n", num);
    X265_CHECK(shift <= 10, "shift too large %d\n", shift);
    X265_CHECK(((intptr_t)coef & 31) == 0, "dequant coef buffer not aligned\n");

    int add, coeffQ;

    add = 1 << (shift - 1);

    for (int n = 0; n < num; n++)
    {
        coeffQ = (quantCoef[n] * scale + add) >> shift;
        coef[n] = Clip3(-32768, 32767, coeffQ);
    }
}

void dequant_scaling_c(const int16_t* quantCoef, const int32_t *deQuantCoef, int32_t* coef, int num, int per, int shift)
{
    X265_CHECK(num <= 32 * 32, "dequant num %d too large\n", num);

    int add, coeffQ;

    shift += 4;

    if (shift > per)
    {
        add = 1 << (shift - per - 1);

        for (int n = 0; n < num; n++)
        {
            coeffQ = ((quantCoef[n] * deQuantCoef[n]) + add) >> (shift - per);
            coef[n] = Clip3(-32768, 32767, coeffQ);
        }
    }
    else
    {
        for (int n = 0; n < num; n++)
        {
            coeffQ   = Clip3(-32768, 32767, quantCoef[n] * deQuantCoef[n]);
            coef[n] = Clip3(-32768, 32767, coeffQ << (per - shift));
        }
    }
}

uint32_t quant_c(int32_t* coef, int32_t* quantCoeff, int32_t* deltaU, int16_t* qCoef, int qBits, int add, int numCoeff)
{
    X265_CHECK(qBits >= 8, "qBits less than 8\n");
    X265_CHECK((numCoeff % 16) == 0, "numCoeff must be multiple of 16\n");
    int qBits8 = qBits - 8;
    uint32_t numSig = 0;

    for (int blockpos = 0; blockpos < numCoeff; blockpos++)
    {
        int level = coef[blockpos];
        int sign  = (level < 0 ? -1 : 1);

        int tmplevel = abs(level) * quantCoeff[blockpos];
        level = ((tmplevel + add) >> qBits);
        deltaU[blockpos] = ((tmplevel - (level << qBits)) >> qBits8);
        if (level)
            ++numSig;
        level *= sign;
        qCoef[blockpos] = (int16_t)Clip3(-32768, 32767, level);
    }

    return numSig;
}

uint32_t nquant_c(int32_t* coef, int32_t* quantCoeff, int16_t* qCoef, int qBits, int add, int numCoeff)
{
    X265_CHECK((numCoeff % 16) == 0, "number of quant coeff is not multiple of 4x4\n");
    X265_CHECK((uint32_t)add < ((uint32_t)1 << qBits), "2 ^ qBits less than add\n");
    X265_CHECK(((intptr_t)quantCoeff & 31) == 0, "quantCoeff buffer not aligned\n");

    uint32_t numSig = 0;

    for (int blockpos = 0; blockpos < numCoeff; blockpos++)
    {
        int level = coef[blockpos];
        int sign  = (level < 0 ? -1 : 1);

        int tmplevel = abs(level) * quantCoeff[blockpos];
        level = ((tmplevel + add) >> qBits);
        if (level)
            ++numSig;
        level *= sign;
        qCoef[blockpos] = (int16_t)Clip3(-32768, 32767, level);
    }

    return numSig;
}

int  count_nonzero_c(const int16_t *quantCoeff, int numCoeff)
{
    X265_CHECK(((intptr_t)quantCoeff & 15) == 0, "quant buffer not aligned\n");
    X265_CHECK(numCoeff > 0 && (numCoeff & 15) == 0, "numCoeff invalid %d\n", numCoeff);

    int count = 0;

    for (int i = 0; i < numCoeff; i++)
    {
        count += quantCoeff[i] != 0;
    }

    return count;
}

template<int trSize>
uint32_t copy_count(int16_t* coeff, int16_t* residual, intptr_t stride)
{
    uint32_t numSig = 0;
    for (int k = 0; k < trSize; k++)
    {
        for (int j = 0; j < trSize; j++)
        {
            coeff[k * trSize + j] = residual[k * stride + j];
            numSig += (residual[k * stride + j] != 0);
        }
    }

    return numSig;
}

void denoiseDct_c(int32_t* dctCoef, uint32_t* resSum, uint16_t* offset, int numCoeff)
{
    for (int i = 0; i < numCoeff; i++)
    {
        int level = dctCoef[i];
        int sign = level >> 31;
        level = (level + sign) ^ sign;
        resSum[i] += level;
        level -= offset[i];
        dctCoef[i] = level < 0 ? 0 : (level ^ sign) - sign;
    }
}

}  // closing - anonymous file-static namespace

namespace x265 {
// x265 private namespace

void Setup_C_DCTPrimitives(EncoderPrimitives& p)
{
    p.dequant_scaling = dequant_scaling_c;
    p.dequant_normal = dequant_normal_c;
    p.quant = quant_c;
    p.nquant = nquant_c;
    p.dct[DST_4x4] = dst4_c;
    p.dct[DCT_4x4] = dct4_c;
    p.dct[DCT_8x8] = dct8_c;
    p.dct[DCT_16x16] = dct16_c;
    p.dct[DCT_32x32] = dct32_c;
    p.idct[IDST_4x4] = idst4_c;
    p.idct[IDCT_4x4] = idct4_c;
    p.idct[IDCT_8x8] = idct8_c;
    p.idct[IDCT_16x16] = idct16_c;
    p.idct[IDCT_32x32] = idct32_c;
    p.count_nonzero = count_nonzero_c;
    p.denoiseDct = denoiseDct_c;

    p.copy_cnt[BLOCK_4x4] = copy_count<4>;
    p.copy_cnt[BLOCK_8x8] = copy_count<8>;
    p.copy_cnt[BLOCK_16x16] = copy_count<16>;
    p.copy_cnt[BLOCK_32x32] = copy_count<32>;
}
}
