/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
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
#include "mbdstharness.h"

using namespace x265;

struct DctConf
{
    const char *name;
    int width;
};

const DctConf dctInfo[] =
{
    { "dst4x4\t",    4 },
    { "dct4x4\t",    4 },
    { "dct8x8\t",    8 },
    { "dct16x16",   16 },
    { "dct32x32",   32 },
};

const DctConf idctInfo[] =
{
    { "idst4x4\t",    4 },
    { "idct4x4\t",    4 },
    { "idct8x8\t",    8 },
    { "idct16x16",   16 },
    { "idct32x32",   32 },
};

MBDstHarness::MBDstHarness()
{
    const int idct_max = (1 << (BIT_DEPTH + 4)) - 1;

    /* [0] --- Random values
     * [1] --- Minimum
     * [2] --- Maximum */
    for (int i = 0; i < TEST_BUF_SIZE; i++)
    {
        short_test_buff[0][i]    = (rand() & PIXEL_MAX) - (rand() & PIXEL_MAX);
        int_test_buff[0][i]      = rand() % PIXEL_MAX;
        int_idct_test_buff[0][i] = (rand() % (SHORT_MAX - SHORT_MIN)) - SHORT_MAX;
        int_denoise_test_buff1[0][i] = int_denoise_test_buff2[0][i] = (rand() & UNSIGNED_SHORT_MAX) - (rand() & UNSIGNED_SHORT_MAX);

        short_test_buff[1][i]    = -PIXEL_MAX;
        int_test_buff[1][i]      = -PIXEL_MAX;
        int_idct_test_buff[1][i] = SHORT_MIN;
        int_denoise_test_buff1[1][i] = int_denoise_test_buff2[1][i] = -UNSIGNED_SHORT_MAX;

        short_test_buff[2][i]    = PIXEL_MAX;
        int_test_buff[2][i]      = PIXEL_MAX;
        int_idct_test_buff[2][i] = SHORT_MAX;
        int_denoise_test_buff1[2][i] = int_denoise_test_buff2[2][i] = UNSIGNED_SHORT_MAX;

        mbuf1[i] = rand() & PIXEL_MAX;
        mbufdct[i] = (rand() & PIXEL_MAX) - (rand() & PIXEL_MAX);
        mbufidct[i] = (rand() & idct_max);
    }

#if _DEBUG
    memset(mshortbuf2, 0, MAX_TU_SIZE * sizeof(int16_t));
    memset(mshortbuf3, 0, MAX_TU_SIZE * sizeof(int16_t));

    memset(mintbuf1, 0, MAX_TU_SIZE * sizeof(int));
    memset(mintbuf2, 0, MAX_TU_SIZE * sizeof(int));
    memset(mintbuf3, 0, MAX_TU_SIZE * sizeof(int));
    memset(mintbuf4, 0, MAX_TU_SIZE * sizeof(int));
#endif // if _DEBUG
}

bool MBDstHarness::check_dct_primitive(dct_t ref, dct_t opt, intptr_t width)
{
    int j = 0;
    intptr_t cmp_size = sizeof(int) * width * width;

    for (int i = 0; i < ITERS; i++)
    {
        int index = rand() % TEST_CASES;

        ref(short_test_buff[index] + j, mintbuf3, width);
        checked(opt, short_test_buff[index] + j, mintbuf4, width);

        if (memcmp(mintbuf3, mintbuf4, cmp_size))
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_idct_primitive(idct_t ref, idct_t opt, intptr_t width)
{
    int j = 0;
    intptr_t cmp_size = sizeof(int16_t) * width * width;

    for (int i = 0; i < ITERS; i++)
    {
        int index = rand() % TEST_CASES;

        ref(int_idct_test_buff[index] + j, mshortbuf2, width);
        checked(opt, int_idct_test_buff[index] + j, mshortbuf3, width);

        if (memcmp(mshortbuf2, mshortbuf3, cmp_size))
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_dequant_primitive(dequant_normal_t ref, dequant_normal_t opt)
{
    int j = 0;

    for (int i = 0; i < ITERS; i++)
    {
        int index = rand() % TEST_CASES;
        int log2TrSize = (rand() % 4) + 2;

        int width = (1 << log2TrSize);
        int height = width;
        int qp = rand() % (QP_MAX_SPEC + QP_BD_OFFSET + 1);
        int per = qp / 6;
        int rem = qp % 6;
        static const int invQuantScales[6] = { 40, 45, 51, 57, 64, 72 };
        int scale = invQuantScales[rem] << per;
        int transformShift = MAX_TR_DYNAMIC_RANGE - X265_DEPTH - log2TrSize;
        int shift = QUANT_IQUANT_SHIFT - QUANT_SHIFT - transformShift;

        ref(short_test_buff[index] + j, mintbuf3, width * height, scale, shift);
        checked(opt, short_test_buff[index] + j, mintbuf4, width * height, scale, shift);

        if (memcmp(mintbuf3, mintbuf4, sizeof(int) * height * width))
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_dequant_primitive(dequant_scaling_t ref, dequant_scaling_t opt)
{
    int j = 0;

    for (int i = 0; i < ITERS; i++)
    {
        int log2TrSize = (rand() % 4) + 2;

        int width = (1 << log2TrSize);
        int height = width;

        int qp = rand() % (QP_MAX_SPEC + QP_BD_OFFSET + 1);
        int per = qp / 6;
        int transformShift = MAX_TR_DYNAMIC_RANGE - X265_DEPTH - log2TrSize;
        int shift = QUANT_IQUANT_SHIFT - QUANT_SHIFT - transformShift;

        int cmp_size = sizeof(int) * height * width;
        int index1 = rand() % TEST_CASES;

        ref(short_test_buff[index1] + j, mintbuf3, mintbuf1, width * height, per, shift);
        checked(opt, short_test_buff[index1] + j, mintbuf4, mintbuf2, width * height, per, shift);

        if (memcmp(mintbuf1, mintbuf2, cmp_size))
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_quant_primitive(quant_t ref, quant_t opt)
{
    int j = 0;

    for (int i = 0; i < ITERS; i++)
    {
        int width = (rand() % 4 + 1) * 4;
        int height = width;

        uint32_t optReturnValue = 0;
        uint32_t refReturnValue = 0;

        int bits = (rand() % 24) + 8;
        int valueToAdd = rand() % (1 << bits);
        int cmp_size = sizeof(int) * height * width;
        int cmp_size1 = sizeof(short) * height * width;
        int numCoeff = height * width;

        int index1 = rand() % TEST_CASES;
        int index2 = rand() % TEST_CASES;

        refReturnValue = ref(int_test_buff[index1] + j, int_test_buff[index2] + j, mintbuf1, mshortbuf2, bits, valueToAdd, numCoeff);
        optReturnValue = (uint32_t)checked(opt, int_test_buff[index1] + j, int_test_buff[index2] + j, mintbuf3, mshortbuf3, bits, valueToAdd, numCoeff);

        if (memcmp(mintbuf1, mintbuf3, cmp_size))
            return false;

        if (memcmp(mshortbuf2, mshortbuf3, cmp_size1))
            return false;

        if (optReturnValue != refReturnValue)
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_nquant_primitive(nquant_t ref, nquant_t opt)
{
    int j = 0;

    for (int i = 0; i < ITERS; i++)
    {
        int width = (rand() % 4 + 1) * 4;
        int height = width;

        uint32_t optReturnValue = 0;
        uint32_t refReturnValue = 0;

        int bits = rand() % 32;
        int valueToAdd = rand() % (1 << bits);
        int cmp_size = sizeof(short) * height * width;
        int numCoeff = height * width;

        int index1 = rand() % TEST_CASES;
        int index2 = rand() % TEST_CASES;

        refReturnValue = ref(int_test_buff[index1] + j, int_test_buff[index2] + j, mshortbuf2, bits, valueToAdd, numCoeff);
        optReturnValue = (uint32_t)checked(opt, int_test_buff[index1] + j, int_test_buff[index2] + j, mshortbuf3, bits, valueToAdd, numCoeff);

        if (memcmp(mshortbuf2, mshortbuf3, cmp_size))
            return false;

        if (optReturnValue != refReturnValue)
            return false;

        reportfail();
        j += INCR;
    }

    return true;
}

bool MBDstHarness::check_count_nonzero_primitive(count_nonzero_t ref, count_nonzero_t opt)
{
    ALIGN_VAR_32(int16_t, qcoeff[32 * 32]);

    for (int i = 0; i < 4; i++)
    {
        int log2TrSize = i + 2;
        int num = 1 << (log2TrSize * 2);
        int mask = num - 1;

        for (int n = 0; n <= num; n++)
        {
            memset(qcoeff, 0, num * sizeof(int16_t));

            for (int j = 0; j < n; j++)
            {
                int k = rand() & mask;
                while (qcoeff[k])
                {
                    k = (k + 11) & mask;
                }

                qcoeff[k] = (int16_t)rand() - RAND_MAX / 2;
            }

            int refval = ref(qcoeff, num);
            int optval = (int)checked(opt, qcoeff, num);

            if (refval != optval)
                return false;

            reportfail();
        }
    }

    return true;
}

bool MBDstHarness::check_denoise_dct_primitive(denoiseDct_t ref, denoiseDct_t opt)
{
    int j = 0;

    for (int s = 0; s < 4; s++)
    {
        int log2TrSize = s + 2;
        int num = 1 << (log2TrSize * 2);
        int cmp_size = sizeof(int) * num;

        for (int i = 0; i < ITERS; i++)
        {
            memset(mubuf1, 0, num * sizeof(uint32_t));
            memset(mubuf2, 0, num * sizeof(uint32_t));
            memset(mushortbuf1, 0,  num * sizeof(uint16_t));

            for (int k = 0; k < num; k++)
                mushortbuf1[k] = rand() % UNSIGNED_SHORT_MAX;

            int index = rand() % TEST_CASES;

            ref(int_denoise_test_buff1[index] + j, mubuf1, mushortbuf1, num);
            checked(opt, int_denoise_test_buff2[index] + j, mubuf2, mushortbuf1, num);

            if (memcmp(int_denoise_test_buff1[index] + j, int_denoise_test_buff2[index] + j, cmp_size))
                return false;

            if (memcmp(mubuf1, mubuf2, cmp_size))
                return false;

            reportfail();
            j += INCR;
        }
        j = 0;
    }

    return true;
}


bool MBDstHarness::testCorrectness(const EncoderPrimitives& ref, const EncoderPrimitives& opt)
{
    for (int i = 0; i < NUM_DCTS; i++)
    {
        if (opt.dct[i])
        {
            if (!check_dct_primitive(ref.dct[i], opt.dct[i], dctInfo[i].width))
            {
                printf("\n%s failed\n", dctInfo[i].name);
                return false;
            }
        }
    }

    for (int i = 0; i < NUM_IDCTS; i++)
    {
        if (opt.idct[i])
        {
            if (!check_idct_primitive(ref.idct[i], opt.idct[i], idctInfo[i].width))
            {
                printf("%s failed\n", idctInfo[i].name);
                return false;
            }
        }
    }

    if (opt.dequant_normal)
    {
        if (!check_dequant_primitive(ref.dequant_normal, opt.dequant_normal))
        {
            printf("dequant: Failed!\n");
            return false;
        }
    }

    if (opt.dequant_scaling)
    {
        if (!check_dequant_primitive(ref.dequant_scaling, opt.dequant_scaling))
        {
            printf("dequant_scaling: Failed!\n");
            return false;
        }
    }

    if (opt.quant)
    {
        if (!check_quant_primitive(ref.quant, opt.quant))
        {
            printf("quant: Failed!\n");
            return false;
        }
    }

    if (opt.nquant)
    {
        if (!check_nquant_primitive(ref.nquant, opt.nquant))
        {
            printf("nquant: Failed!\n");
            return false;
        }
    }

    if (opt.count_nonzero)
    {
        if (!check_count_nonzero_primitive(ref.count_nonzero, opt.count_nonzero))
        {
            printf("count_nonzero: Failed!\n");
            return false;
        }
    }

    if (opt.dequant_scaling)
    {
        if (!check_dequant_primitive(ref.dequant_scaling, opt.dequant_scaling))
        {
            printf("dequant_scaling: Failed!\n");
            return false;
        }
    }

    if (opt.denoiseDct)
    {
        if (!check_denoise_dct_primitive(ref.denoiseDct, opt.denoiseDct))
        {
            printf("denoiseDct: Failed!\n");
            return false;
        }
    }

    return true;
}

void MBDstHarness::measureSpeed(const EncoderPrimitives& ref, const EncoderPrimitives& opt)
{
    for (int value = 0; value < NUM_DCTS; value++)
    {
        if (opt.dct[value])
        {
            printf("%s\t", dctInfo[value].name);
            REPORT_SPEEDUP(opt.dct[value], ref.dct[value], mbuf1, mintbuf3, dctInfo[value].width);
        }
    }

    for (int value = 0; value < NUM_IDCTS; value++)
    {
        if (opt.idct[value])
        {
            printf("%s\t", idctInfo[value].name);
            REPORT_SPEEDUP(opt.idct[value], ref.idct[value], mbufidct, mshortbuf2, idctInfo[value].width);
        }
    }

    if (opt.dequant_normal)
    {
        printf("dequant_normal\t");
        REPORT_SPEEDUP(opt.dequant_normal, ref.dequant_normal, short_test_buff[0], mintbuf3, 32 * 32, 70, 1);
    }

    if (opt.dequant_scaling)
    {
        printf("dequant_scaling\t");
        REPORT_SPEEDUP(opt.dequant_scaling, ref.dequant_scaling, short_test_buff[0], mintbuf3, mintbuf4, 32 * 32, 5, 1);
    }

    if (opt.quant)
    {
        printf("quant\t\t");
        REPORT_SPEEDUP(opt.quant, ref.quant, int_test_buff[0], int_test_buff[1], mintbuf3, mshortbuf2, 23, 23785, 32 * 32);
    }

    if (opt.nquant)
    {
        printf("nquant\t\t");
        REPORT_SPEEDUP(opt.nquant, ref.nquant, int_test_buff[0], int_test_buff[1], mshortbuf2, 23, 23785, 32 * 32);
    }

    if (opt.count_nonzero)
    {
        for (int i = 4; i <= 32; i <<= 1)
        {
            printf("count_nonzero[%dx%d]", i, i);
            REPORT_SPEEDUP(opt.count_nonzero, ref.count_nonzero, mbuf1, i * i)
        }
    }

    if (opt.denoiseDct)
    {
        printf("denoiseDct\t");
        REPORT_SPEEDUP(opt.denoiseDct, ref.denoiseDct, int_denoise_test_buff1[0], mubuf1, mushortbuf1, 32 * 32);
    }

}
