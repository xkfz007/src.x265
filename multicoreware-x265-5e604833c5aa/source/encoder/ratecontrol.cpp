/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Sumalatha Polureddy <sumalatha@multicorewareinc.com>
 *          Aarthi Priya Thirumalai <aarthi@multicorewareinc.com>
 *          Xun Xu, PPLive Corporation <xunxu@pptv.com>
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
#include "param.h"
#include "frame.h"
#include "framedata.h"
#include "picyuv.h"

#include "encoder.h"
#include "slicetype.h"
#include "ratecontrol.h"
#include "sei.h"

#define BR_SHIFT  6
#define CPB_SHIFT 4

using namespace x265;

/* Amortize the partial cost of I frames over the next N frames */
const double RateControl::s_amortizeFraction = 0.85;
const int RateControl::s_amortizeFrames = 75;
const int RateControl::s_slidingWindowFrames = 20;
const char *RateControl::s_defaultStatFileName = "x265_2pass.log";

namespace {
#define CMP_OPT_FIRST_PASS(opt, param_val)\
{\
    bErr = 0;\
    p = strstr(opts, opt "=");\
    char* q = strstr(opts, "no-"opt);\
    if (p && sscanf(p, opt "=%d" , &i) && param_val != i)\
        bErr = 1;\
    else if (!param_val && !q && !p)\
        bErr = 1;\
    else if (param_val && (q || !strstr(opts, opt)))\
        bErr = 1;\
    if (bErr)\
    {\
        x265_log(m_param, X265_LOG_ERROR, "different " opt " setting than first pass (%d vs %d)\n", param_val, i);\
        return false;\
    }\
}

inline int calcScale(uint32_t x)
{
    static uint8_t lut[16] = {4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};
    int y, z = (((x & 0xffff) - 1) >> 27) & 16;
    x >>= z;
    z += y = (((x & 0xff) - 1) >> 28) & 8;
    x >>= y;
    z += y = (((x & 0xf) - 1) >> 29) & 4;
    x >>= y;
    return z + lut[x&0xf];
}

inline int calcLength(uint32_t x)
{
    static uint8_t lut[16] = {4, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    int y, z = (((x >> 16) - 1) >> 27) & 16;
    x >>= z ^ 16;
    z += y = ((x - 0x100) >> 28) & 8;
    x >>= y ^ 8;
    z += y = ((x - 0x10) >> 29) & 4;
    x >>= y ^ 4;
    return z + lut[x];
}

inline void reduceFraction(int* n, int* d)
{
    int a = *n;
    int b = *d;
    int c;
    if (!a || !b)
        return;
    c = a % b;
    while (c)
    {
        a = b;
        b = c;
        c = a % b;
    }
    *n /= b;
    *d /= b;
}

inline char *strcatFilename(const char *input, const char *suffix)
{
    char *output = X265_MALLOC(char, strlen(input) + strlen(suffix) + 1);
    if (!output)
    {
        x265_log(NULL, X265_LOG_ERROR, "unable to allocate memory for filename\n");
        return NULL;
    }
    strcpy(output, input);
    strcat(output, suffix);
    return output;
}

inline double qScale2bits(RateControlEntry *rce, double qScale)
{
    if (qScale < 0.1)
        qScale = 0.1;
    return (rce->coeffBits + .1) * pow(rce->qScale / qScale, 1.1)
           + rce->mvBits * pow(X265_MAX(rce->qScale, 1) / X265_MAX(qScale, 1), 0.5)
           + rce->miscBits;
}

inline void copyRceData(RateControlEntry* rce, RateControlEntry* rce2Pass)
{
    rce->coeffBits = rce2Pass->coeffBits;
    rce->mvBits = rce2Pass->mvBits;
    rce->miscBits = rce2Pass->miscBits;
    rce->iCuCount = rce2Pass->iCuCount;
    rce->pCuCount = rce2Pass->pCuCount;
    rce->skipCuCount = rce2Pass->skipCuCount;
    rce->keptAsRef = rce2Pass->keptAsRef;
    rce->qScale = rce2Pass->qScale;
    rce->newQScale = rce2Pass->newQScale;
    rce->expectedBits = rce2Pass->expectedBits;
    rce->expectedVbv = rce2Pass->expectedVbv;
    rce->blurredComplexity = rce2Pass->blurredComplexity;
    rce->sliceType = rce2Pass->sliceType;
}

}  // end anonymous namespace
/* Compute variance to derive AC energy of each block */
static inline uint32_t acEnergyVar(Frame *curFrame, uint64_t sum_ssd, int shift, int i)
{
    uint32_t sum = (uint32_t)sum_ssd;
    uint32_t ssd = (uint32_t)(sum_ssd >> 32);

    curFrame->m_lowres.wp_sum[i] += sum;
    curFrame->m_lowres.wp_ssd[i] += ssd;
    return ssd - ((uint64_t)sum * sum >> shift);
}

/* Find the energy of each block in Y/Cb/Cr plane */
static inline uint32_t acEnergyPlane(Frame *curFrame, pixel* src, intptr_t srcStride, int bChroma, int colorFormat)
{
    if ((colorFormat != X265_CSP_I444) && bChroma)
    {
        ALIGN_VAR_8(pixel, pix[8 * 8]);
        primitives.luma_copy_pp[LUMA_8x8](pix, 8, src, srcStride);
        return acEnergyVar(curFrame, primitives.var[BLOCK_8x8](pix, 8), 6, bChroma);
    }
    else
        return acEnergyVar(curFrame, primitives.var[BLOCK_16x16](src, srcStride), 8, bChroma);
}

/* Find the total AC energy of each block in all planes */
uint32_t RateControl::acEnergyCu(Frame* curFrame, uint32_t block_x, uint32_t block_y)
{
    intptr_t stride = curFrame->m_origPicYuv->m_stride;
    intptr_t cStride = curFrame->m_origPicYuv->m_strideC;
    intptr_t blockOffsetLuma = block_x + (block_y * stride);
    int colorFormat = m_param->internalCsp;
    int hShift = CHROMA_H_SHIFT(colorFormat);
    int vShift = CHROMA_V_SHIFT(colorFormat);
    intptr_t blockOffsetChroma = (block_x >> hShift) + ((block_y >> vShift) * cStride);

    uint32_t var;

    var  = acEnergyPlane(curFrame, curFrame->m_origPicYuv->m_picOrg[0] + blockOffsetLuma, stride, 0, colorFormat);
    var += acEnergyPlane(curFrame, curFrame->m_origPicYuv->m_picOrg[1] + blockOffsetChroma, cStride, 1, colorFormat);
    var += acEnergyPlane(curFrame, curFrame->m_origPicYuv->m_picOrg[2] + blockOffsetChroma, cStride, 2, colorFormat);
    x265_emms();
    return var;
}

void RateControl::calcAdaptiveQuantFrame(Frame *curFrame)
{
    /* Actual adaptive quantization */
    int maxCol = curFrame->m_origPicYuv->m_picWidth;
    int maxRow = curFrame->m_origPicYuv->m_picHeight;

    for (int y = 0; y < 3; y++)
    {
        curFrame->m_lowres.wp_ssd[y] = 0;
        curFrame->m_lowres.wp_sum[y] = 0;
    }

    /* Calculate Qp offset for each 16x16 block in the frame */
    int block_xy = 0;
    int block_x = 0, block_y = 0;
    double strength = 0.f;
    if (m_param->rc.aqMode == X265_AQ_NONE || m_param->rc.aqStrength == 0)
    {
        /* Need to init it anyways for CU tree */
        int cuWidth = ((maxCol / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
        int cuHeight = ((maxRow / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
        int cuCount = cuWidth * cuHeight;

        if (m_param->rc.aqMode && m_param->rc.aqStrength == 0)
        {
            memset(curFrame->m_lowres.qpCuTreeOffset, 0, cuCount * sizeof(double));
            memset(curFrame->m_lowres.qpAqOffset, 0, cuCount * sizeof(double));
            for (int cuxy = 0; cuxy < cuCount; cuxy++)
                curFrame->m_lowres.invQscaleFactor[cuxy] = 256;
        }

        /* Need variance data for weighted prediction */
        if (m_param->bEnableWeightedPred || m_param->bEnableWeightedBiPred)
        {
            for (block_y = 0; block_y < maxRow; block_y += 16)
                for (block_x = 0; block_x < maxCol; block_x += 16)
                    acEnergyCu(curFrame, block_x, block_y);
        }
    }
    else
    {
        block_xy = 0;
        double avg_adj_pow2 = 0, avg_adj = 0, qp_adj = 0;
        if (m_param->rc.aqMode == X265_AQ_AUTO_VARIANCE)
        {
            double bit_depth_correction = pow(1 << (X265_DEPTH - 8), 0.5);
            for (block_y = 0; block_y < maxRow; block_y += 16)
            {
                for (block_x = 0; block_x < maxCol; block_x += 16)
                {
                    uint32_t energy = acEnergyCu(curFrame, block_x, block_y);
                    qp_adj = pow(energy + 1, 0.1);
                    curFrame->m_lowres.qpCuTreeOffset[block_xy] = qp_adj;
                    avg_adj += qp_adj;
                    avg_adj_pow2 += qp_adj * qp_adj;
                    block_xy++;
                }
            }

            avg_adj /= m_ncu;
            avg_adj_pow2 /= m_ncu;
            strength = m_param->rc.aqStrength * avg_adj / bit_depth_correction;
            avg_adj = avg_adj - 0.5f * (avg_adj_pow2 - (11.f * bit_depth_correction)) / avg_adj;
        }
        else
            strength = m_param->rc.aqStrength * 1.0397f;

        block_xy = 0;
        for (block_y = 0; block_y < maxRow; block_y += 16)
        {
            for (block_x = 0; block_x < maxCol; block_x += 16)
            {
                if (m_param->rc.aqMode == X265_AQ_AUTO_VARIANCE)
                {
                    qp_adj = curFrame->m_lowres.qpCuTreeOffset[block_xy];
                    qp_adj = strength * (qp_adj - avg_adj);
                }
                else
                {
                    uint32_t energy = acEnergyCu(curFrame, block_x, block_y);
                    qp_adj = strength * (X265_LOG2(X265_MAX(energy, 1)) - (14.427f + 2 * (X265_DEPTH - 8)));
                }
                curFrame->m_lowres.qpAqOffset[block_xy] = qp_adj;
                curFrame->m_lowres.qpCuTreeOffset[block_xy] = qp_adj;
                curFrame->m_lowres.invQscaleFactor[block_xy] = x265_exp2fix8(qp_adj);
                block_xy++;
            }
        }
    }

    if (m_param->bEnableWeightedPred || m_param->bEnableWeightedBiPred)
    {
        int hShift = CHROMA_H_SHIFT(m_param->internalCsp);
        int vShift = CHROMA_V_SHIFT(m_param->internalCsp);
        maxCol = ((maxCol + 8) >> 4) << 4;
        maxRow = ((maxRow + 8) >> 4) << 4;
        int width[3]  = { maxCol, maxCol >> hShift, maxCol >> hShift };
        int height[3] = { maxRow, maxRow >> vShift, maxRow >> vShift };

        for (int i = 0; i < 3; i++)
        {
            uint64_t sum, ssd;
            sum = curFrame->m_lowres.wp_sum[i];
            ssd = curFrame->m_lowres.wp_ssd[i];
            curFrame->m_lowres.wp_ssd[i] = ssd - (sum * sum + (width[i] * height[i]) / 2) / (width[i] * height[i]);
        }
    }
}

RateControl::RateControl(x265_param *p)
{
    m_param = p;
    int lowresCuWidth = ((m_param->sourceWidth / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
    int lowresCuHeight = ((m_param->sourceHeight / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
    m_ncu = lowresCuWidth * lowresCuHeight;

    if (m_param->rc.cuTree)
        m_qCompress = 1;
    else
        m_qCompress = m_param->rc.qCompress;

    // validate for param->rc, maybe it is need to add a function like x265_parameters_valiate()
    m_residualFrames = 0;
    m_partialResidualFrames = 0;
    m_residualCost = 0;
    m_partialResidualCost = 0;
    m_rateFactorMaxIncrement = 0;
    m_rateFactorMaxDecrement = 0;
    m_fps = m_param->fpsNum / m_param->fpsDenom;
    m_startEndOrder.set(0);
    m_bTerminated = false;
    m_finalFrameCount = 0;
    m_numEntries = 0;
    if (m_param->rc.rateControlMode == X265_RC_CRF)
    {
        m_param->rc.qp = (int)m_param->rc.rfConstant;
        m_param->rc.bitrate = 0;

        double baseCplx = m_ncu * (m_param->bframes ? 120 : 80);
        double mbtree_offset = m_param->rc.cuTree ? (1.0 - m_param->rc.qCompress) * 13.5 : 0;
        m_rateFactorConstant = pow(baseCplx, 1 - m_qCompress) /
            x265_qp2qScale(m_param->rc.rfConstant + mbtree_offset);
        if (m_param->rc.rfConstantMax)
        {
            m_rateFactorMaxIncrement = m_param->rc.rfConstantMax - m_param->rc.rfConstant;
            if (m_rateFactorMaxIncrement <= 0)
            {
                x265_log(m_param, X265_LOG_WARNING, "CRF max must be greater than CRF\n");
                m_rateFactorMaxIncrement = 0;
            }
        }
        if (m_param->rc.rfConstantMin)
            m_rateFactorMaxDecrement = m_param->rc.rfConstant - m_param->rc.rfConstantMin;
    }
    m_isAbr = m_param->rc.rateControlMode != X265_RC_CQP && !m_param->rc.bStatRead;
    m_2pass = m_param->rc.rateControlMode == X265_RC_ABR && m_param->rc.bStatRead;
    m_bitrate = m_param->rc.bitrate * 1000;
    m_frameDuration = (double)m_param->fpsDenom / m_param->fpsNum;
    m_qp = m_param->rc.qp;
    m_lastRceq = 1; /* handles the cmplxrsum when the previous frame cost is zero */
    m_shortTermCplxSum = 0;
    m_shortTermCplxCount = 0;
    m_lastNonBPictType = I_SLICE;
    m_isAbrReset = false;
    m_lastAbrResetPoc = -1;
    m_statFileOut = NULL;
    m_cutreeStatFileOut = m_cutreeStatFileIn = NULL;
    m_rce2Pass = NULL;

    // vbv initialization
    m_param->rc.vbvBufferSize = Clip3(0, 2000000, m_param->rc.vbvBufferSize);
    m_param->rc.vbvMaxBitrate = Clip3(0, 2000000, m_param->rc.vbvMaxBitrate);
    m_param->rc.vbvBufferInit = Clip3(0.0, 2000000.0, m_param->rc.vbvBufferInit);
    m_singleFrameVbv = 0;
    if (m_param->rc.vbvBufferSize)
    {
        if (m_param->rc.rateControlMode == X265_RC_CQP)
        {
            x265_log(m_param, X265_LOG_WARNING, "VBV is incompatible with constant QP, ignored.\n");
            m_param->rc.vbvBufferSize = 0;
            m_param->rc.vbvMaxBitrate = 0;
        }
        else if (m_param->rc.vbvMaxBitrate == 0)
        {
            if (m_param->rc.rateControlMode == X265_RC_ABR)
            {
                x265_log(m_param, X265_LOG_WARNING, "VBV maxrate unspecified, assuming CBR\n");
                m_param->rc.vbvMaxBitrate = m_param->rc.bitrate;
            }
            else
            {
                x265_log(m_param, X265_LOG_WARNING, "VBV bufsize set but maxrate unspecified, ignored\n");
                m_param->rc.vbvBufferSize = 0;
            }
        }
        else if (m_param->rc.vbvMaxBitrate < m_param->rc.bitrate &&
                 m_param->rc.rateControlMode == X265_RC_ABR)
        {
            x265_log(m_param, X265_LOG_WARNING, "max bitrate less than average bitrate, assuming CBR\n");
            m_param->rc.bitrate = m_param->rc.vbvMaxBitrate;
        }
    }
    else if (m_param->rc.vbvMaxBitrate)
    {
        x265_log(m_param, X265_LOG_WARNING, "VBV maxrate specified, but no bufsize, ignored\n");
        m_param->rc.vbvMaxBitrate = 0;
    }
    m_isVbv = m_param->rc.vbvMaxBitrate > 0 && m_param->rc.vbvBufferSize > 0;
    if (m_param->bEmitHRDSEI && !m_isVbv)
    {
        x265_log(m_param, X265_LOG_WARNING, "NAL HRD parameters require VBV parameters, ignored\n");
        m_param->bEmitHRDSEI = 0;
    }

    m_isCbr = m_param->rc.rateControlMode == X265_RC_ABR && m_isVbv && !m_2pass && m_param->rc.vbvMaxBitrate <= m_param->rc.bitrate;
    m_leadingBframes = m_param->bframes;
    m_bframeBits = 0;
    m_leadingNoBSatd = 0;
    m_ipOffset = 6.0 * X265_LOG2(m_param->rc.ipFactor);
    m_pbOffset = 6.0 * X265_LOG2(m_param->rc.pbFactor);

    /* Adjust the first frame in order to stabilize the quality level compared to the rest */
#define ABR_INIT_QP_MIN (24)
#define ABR_INIT_QP_MAX (40)
#define CRF_INIT_QP (int)m_param->rc.rfConstant
    for (int i = 0; i < 3; i++)
        m_lastQScaleFor[i] = x265_qp2qScale(m_param->rc.rateControlMode == X265_RC_CRF ? CRF_INIT_QP : ABR_INIT_QP_MIN);

    if (m_param->rc.rateControlMode == X265_RC_CQP)
    {
        if (m_qp && !m_param->bLossless)
        {
            m_qpConstant[P_SLICE] = m_qp;
            m_qpConstant[I_SLICE] = Clip3(0, QP_MAX_MAX, (int)(m_qp - m_ipOffset + 0.5));
            m_qpConstant[B_SLICE] = Clip3(0, QP_MAX_MAX, (int)(m_qp + m_pbOffset + 0.5));
        }
        else
        {
            m_qpConstant[P_SLICE] = m_qpConstant[I_SLICE] = m_qpConstant[B_SLICE] = m_qp;
        }
    }

    /* qstep - value set as encoder specific */
    m_lstep = pow(2, m_param->rc.qpStep / 6.0);

    for (int i = 0; i < 2; i++)
        m_cuTreeStats.qpBuffer[i] = NULL;
}

bool RateControl::init(const SPS *sps)
{
    if (m_isVbv)
    {
        /* We don't support changing the ABR bitrate right now,
         * so if the stream starts as CBR, keep it CBR. */
        if (m_param->rc.vbvBufferSize < (int)(m_param->rc.vbvMaxBitrate / m_fps))
        {
            m_param->rc.vbvBufferSize = (int)(m_param->rc.vbvMaxBitrate / m_fps);
            x265_log(m_param, X265_LOG_WARNING, "VBV buffer size cannot be smaller than one frame, using %d kbit\n",
                     m_param->rc.vbvBufferSize);
        }
        int vbvBufferSize = m_param->rc.vbvBufferSize * 1000;
        int vbvMaxBitrate = m_param->rc.vbvMaxBitrate * 1000;

        if (m_param->bEmitHRDSEI)
        {
            const HRDInfo* hrd = &sps->vuiParameters.hrdParameters;
            vbvBufferSize = hrd->cpbSizeValue << (hrd->cpbSizeScale + CPB_SHIFT);
            vbvMaxBitrate = hrd->bitRateValue << (hrd->bitRateScale + BR_SHIFT);
        }
        m_bufferRate = vbvMaxBitrate / m_fps;
        m_vbvMaxRate = vbvMaxBitrate;
        m_bufferSize = vbvBufferSize;
        m_singleFrameVbv = m_bufferRate * 1.1 > m_bufferSize;

        if (m_param->rc.vbvBufferInit > 1.)
            m_param->rc.vbvBufferInit = Clip3(0.0, 1.0, m_param->rc.vbvBufferInit / m_param->rc.vbvBufferSize);
        m_param->rc.vbvBufferInit = Clip3(0.0, 1.0, X265_MAX(m_param->rc.vbvBufferInit, m_bufferRate / m_bufferSize));
        m_bufferFillFinal = m_bufferSize * m_param->rc.vbvBufferInit;
    }

    m_totalBits = 0;
    m_framesDone = 0;
    m_residualCost = 0;
    m_partialResidualCost = 0;
    for (int i = 0; i < s_slidingWindowFrames; i++)
    {
        m_satdCostWindow[i] = 0;
        m_encodedBitsWindow[i] = 0;
    }
    m_sliderPos = 0;

    /* 720p videos seem to be a good cutoff for cplxrSum */
    double tuneCplxFactor = (m_param->rc.cuTree && m_ncu > 3600) ? 2.5 : 1;

    /* estimated ratio that produces a reasonable QP for the first I-frame */
    m_cplxrSum = .01 * pow(7.0e5, m_qCompress) * pow(m_ncu, 0.5) * tuneCplxFactor;
    m_wantedBitsWindow = m_bitrate * m_frameDuration;
    m_accumPNorm = .01;
    m_accumPQp = (m_param->rc.rateControlMode == X265_RC_CRF ? CRF_INIT_QP : ABR_INIT_QP_MIN) * m_accumPNorm;

    /* Frame Predictors and Row predictors used in vbv */
    for (int i = 0; i < 5; i++)
    {
        m_pred[i].coeff = 2.0;
        m_pred[i].count = 1.0;
        m_pred[i].decay = 0.5;
        m_pred[i].offset = 0.0;
    }
    m_predBfromP = m_pred[0];
    if (!m_statFileOut && (m_param->rc.bStatWrite || m_param->rc.bStatRead))
    {
        /* If the user hasn't defined the stat filename, use the default value */
        const char *fileName = m_param->rc.statFileName;
        if (!fileName)
            fileName = s_defaultStatFileName;
        /* Load stat file and init 2pass algo */
        if (m_param->rc.bStatRead)
        {
            m_expectedBitsSum = 0;
            char *p, *statsIn, *statsBuf;
            /* read 1st pass stats */
            statsIn = statsBuf = x265_slurp_file(fileName);
            if (!statsBuf)
                return false;
            if (m_param->rc.cuTree)
            {
                char *tmpFile = strcatFilename(fileName, ".cutree");
                if (!tmpFile)
                    return false;
                m_cutreeStatFileIn = fopen(tmpFile, "rb");
                X265_FREE(tmpFile);
                if (!m_cutreeStatFileIn)
                {
                    x265_log(m_param, X265_LOG_ERROR, "can't open stats file %s\n", tmpFile);
                    return false;
                }
            }

            /* check whether 1st pass options were compatible with current options */
            if (strncmp(statsBuf, "#options:", 9))
            {
                x265_log(m_param, X265_LOG_ERROR,"options list in stats file not valid\n");
                return false;
            }
            {
                int i, j;
                uint32_t k , l;
                bool bErr = false;
                char *opts = statsBuf;
                statsIn = strchr(statsBuf, '\n');
                if (!statsIn)
                {
                    x265_log(m_param, X265_LOG_ERROR, "Malformed stats file\n");
                    return false;
                }
                *statsIn = '\0';
                statsIn++;
                if (sscanf(opts, "#options: %dx%d", &i, &j) != 2)
                {
                    x265_log(m_param, X265_LOG_ERROR, "Resolution specified in stats file not valid\n");
                    return false;
                }
                if ((p = strstr(opts, " fps=")) == 0 || sscanf(p, " fps=%u/%u", &k, &l) != 2)
                {
                    x265_log(m_param, X265_LOG_ERROR, "fps specified in stats file not valid\n");
                    return false;
                }
                if (k != m_param->fpsNum || l != m_param->fpsDenom)
                {
                    x265_log(m_param, X265_LOG_ERROR, "fps mismatch with 1st pass (%u/%u vs %u/%u)\n",
                              m_param->fpsNum, m_param->fpsDenom, k, l);
                    return false;
                }
                CMP_OPT_FIRST_PASS("bitdepth", m_param->internalBitDepth);
                CMP_OPT_FIRST_PASS("weightp", m_param->bEnableWeightedPred);
                CMP_OPT_FIRST_PASS("bframes", m_param->bframes);
                CMP_OPT_FIRST_PASS("b-pyramid", m_param->bBPyramid);
                CMP_OPT_FIRST_PASS("open-gop", m_param->bOpenGOP);
                CMP_OPT_FIRST_PASS("keyint", m_param->keyframeMax);
                CMP_OPT_FIRST_PASS("scenecut", m_param->scenecutThreshold);

                if ((p = strstr(opts, "b-adapt=")) != 0 && sscanf(p, "b-adapt=%d", &i) && i >= X265_B_ADAPT_NONE && i <= X265_B_ADAPT_TRELLIS)
                {
                    m_param->bFrameAdaptive = i;
                }
                else if (m_param->bframes)
                {
                    x265_log(m_param, X265_LOG_ERROR, "b-adapt method specified in stats file not valid\n");
                    return false;
                }

                if ((p = strstr(opts, "rc-lookahead=")) != 0 && sscanf(p, "rc-lookahead=%d", &i))
                    m_param->lookaheadDepth = i;
            }
            /* find number of pics */
            p = statsIn;
            int numEntries;
            for (numEntries = -1; p; numEntries++)
                p = strchr(p + 1, ';');
            if (!numEntries)
            {
                x265_log(m_param, X265_LOG_ERROR, "empty stats file\n");
                return false;
            }
            m_numEntries = numEntries;

            if (m_param->totalFrames < m_numEntries && m_param->totalFrames > 0)
            {
                x265_log(m_param, X265_LOG_WARNING, "2nd pass has fewer frames than 1st pass (%d vs %d)\n",
                         m_param->totalFrames, m_numEntries);
            }
            if (m_param->totalFrames > m_numEntries)
            {
                x265_log(m_param, X265_LOG_ERROR, "2nd pass has more frames than 1st pass (%d vs %d)\n",
                         m_param->totalFrames, m_numEntries);
                return false;
            }

            m_rce2Pass = X265_MALLOC(RateControlEntry, m_numEntries);
            if (!m_rce2Pass)
            {
                 x265_log(m_param, X265_LOG_ERROR, "Rce Entries for 2 pass cannot be allocated\n");
                 return false;
            }
            /* init all to skipped p frames */
            for (int i = 0; i < m_numEntries; i++)
            {
                RateControlEntry *rce = &m_rce2Pass[i];
                rce->sliceType = P_SLICE;
                rce->qScale = rce->newQScale = x265_qp2qScale(20);
                rce->miscBits = m_ncu + 10;
                rce->newQp = 0;
            }
            /* read stats */
            p = statsIn;
            double totalQpAq = 0;
            for (int i = 0; i < m_numEntries; i++)
            {
                RateControlEntry *rce;
                int frameNumber;
                char picType;
                int e;
                char *next;
                double qpRc, qpAq;
                next = strstr(p, ";");
                if (next)
                    *next++ = 0;
                e = sscanf(p, " in:%d ", &frameNumber);
                if (frameNumber < 0 || frameNumber >= m_numEntries)
                {
                    x265_log(m_param, X265_LOG_ERROR, "bad frame number (%d) at stats line %d\n", frameNumber, i);
                    return false;
                }
                rce = &m_rce2Pass[frameNumber];
                e += sscanf(p, " in:%*d out:%*d type:%c q:%lf q-aq:%lf tex:%d mv:%d misc:%d icu:%lf pcu:%lf scu:%lf",
                       &picType, &qpRc, &qpAq, &rce->coeffBits,
                       &rce->mvBits, &rce->miscBits, &rce->iCuCount, &rce->pCuCount,
                       &rce->skipCuCount);
                rce->keptAsRef = true;
                if (picType == 'b' || picType == 'p')
                    rce->keptAsRef = false;
                if (picType == 'I' || picType == 'i')
                    rce->sliceType = I_SLICE;
                else if (picType == 'P' || picType == 'p')
                    rce->sliceType = P_SLICE;
                else if (picType == 'B' || picType == 'b')
                    rce->sliceType = B_SLICE;
                else
                    e = -1;
                if (e < 10)
                {
                    x265_log(m_param, X265_LOG_ERROR, "statistics are damaged at line %d, parser out=%d\n", i, e);
                    return false;
                }
                rce->qScale = x265_qp2qScale(qpRc);
                totalQpAq += qpAq;
                p = next;
            }
            X265_FREE(statsBuf);

            if (m_param->rc.rateControlMode == X265_RC_ABR)
            {
                if (!initPass2())
                    return false;
            } /* else we're using constant quant, so no need to run the bitrate allocation */
        }
        /* Open output file */
        /* If input and output files are the same, output to a temp file
         * and move it to the real name only when it's complete */
        if (m_param->rc.bStatWrite)
        {
            char *p, *statFileTmpname;
            statFileTmpname = strcatFilename(fileName, ".temp");
            if (!statFileTmpname)
                return false;
            m_statFileOut = fopen(statFileTmpname, "wb");
            X265_FREE(statFileTmpname);
            if (!m_statFileOut)
            {
                x265_log(m_param, X265_LOG_ERROR, "can't open stats file %s\n", statFileTmpname);
                return false;
            }
            p = x265_param2string(m_param);
            if (p)
                fprintf(m_statFileOut, "#options: %s\n", p);
            X265_FREE(p);
            if (m_param->rc.cuTree && !m_param->rc.bStatRead)
            {
                statFileTmpname = strcatFilename(fileName, ".cutree.temp");
                if (!statFileTmpname)
                    return false;
                m_cutreeStatFileOut = fopen(statFileTmpname, "wb");
                X265_FREE(statFileTmpname);
                if (!m_cutreeStatFileOut)
                {
                    x265_log(m_param, X265_LOG_ERROR, "can't open mbtree stats file %s\n", statFileTmpname);
                    return false;
                }
            }
        }
        if (m_param->rc.cuTree)
        {
            m_cuTreeStats.qpBuffer[0] = X265_MALLOC(uint16_t, m_ncu * sizeof(uint16_t));
            if (m_param->bBPyramid && m_param->rc.bStatRead)
                m_cuTreeStats.qpBuffer[1] = X265_MALLOC(uint16_t, m_ncu * sizeof(uint16_t));
            m_cuTreeStats.qpBufPos = -1;
        }
    }
    return true;
}

void RateControl::initHRD(SPS *sps)
{
    int vbvBufferSize = m_param->rc.vbvBufferSize * 1000;
    int vbvMaxBitrate = m_param->rc.vbvMaxBitrate * 1000;

    // Init HRD
    HRDInfo* hrd = &sps->vuiParameters.hrdParameters;
    hrd->cbrFlag = m_isCbr;

    // normalize HRD size and rate to the value / scale notation
    hrd->bitRateScale = Clip3(0, 15, calcScale(vbvMaxBitrate) - BR_SHIFT);
    hrd->bitRateValue = (vbvMaxBitrate >> (hrd->bitRateScale + BR_SHIFT));

    hrd->cpbSizeScale = Clip3(0, 15, calcScale(vbvBufferSize) - CPB_SHIFT);
    hrd->cpbSizeValue = (vbvBufferSize >> (hrd->cpbSizeScale + CPB_SHIFT));
    int bitRateUnscale = hrd->bitRateValue << (hrd->bitRateScale + BR_SHIFT);
    int cpbSizeUnscale = hrd->cpbSizeValue << (hrd->cpbSizeScale + CPB_SHIFT);

    // arbitrary
    #define MAX_DURATION 0.5

    TimingInfo *time = &sps->vuiParameters.timingInfo;
    int maxCpbOutputDelay = (int)(X265_MIN(m_param->keyframeMax * MAX_DURATION * time->timeScale / time->numUnitsInTick, INT_MAX));
    int maxDpbOutputDelay = (int)(sps->maxDecPicBuffering * MAX_DURATION * time->timeScale / time->numUnitsInTick);
    int maxDelay = (int)(90000.0 * cpbSizeUnscale / bitRateUnscale + 0.5);

    hrd->initialCpbRemovalDelayLength = 2 + Clip3(4, 22, 32 - calcLength(maxDelay));
    hrd->cpbRemovalDelayLength = Clip3(4, 31, 32 - calcLength(maxCpbOutputDelay));
    hrd->dpbOutputDelayLength = Clip3(4, 31, 32 - calcLength(maxDpbOutputDelay));

    #undef MAX_DURATION
}

bool RateControl::initPass2()
{
    uint64_t allConstBits = 0;
    uint64_t allAvailableBits = uint64_t(m_param->rc.bitrate * 1000. * m_numEntries * m_frameDuration);
    double rateFactor, stepMult;
    double qBlur = m_param->rc.qblur;
    double cplxBlur = m_param->rc.complexityBlur;
    const int filterSize = (int)(qBlur * 4) | 1;
    double expectedBits;
    double *qScale, *blurredQscale;
    double baseCplx = m_ncu * (m_param->bframes ? 120 : 80);
    double clippedDuration = CLIP_DURATION(m_frameDuration) / BASE_FRAME_DURATION;

    /* find total/average complexity & const_bits */
    for (int i = 0; i < m_numEntries; i++)
        allConstBits += m_rce2Pass[i].miscBits;

    if (allAvailableBits < allConstBits)
    {
        x265_log(m_param, X265_LOG_ERROR, "requested bitrate is too low. estimated minimum is %d kbps\n",
                 (int)(allConstBits * m_fps / m_numEntries * 1000.));
        return false;
    }

    /* Blur complexities, to reduce local fluctuation of QP.
     * We don't blur the QPs directly, because then one very simple frame
     * could drag down the QP of a nearby complex frame and give it more
     * bits than intended. */
    for (int i = 0; i < m_numEntries; i++)
    {
        double weightSum = 0;
        double cplxSum = 0;
        double weight = 1.0;
        double gaussianWeight;
        /* weighted average of cplx of future frames */
        for (int j = 1; j < cplxBlur * 2 && j < m_numEntries - i; j++)
        {
            RateControlEntry *rcj = &m_rce2Pass[i + j];
            weight *= 1 - pow(rcj->iCuCount / m_ncu, 2);
            if (weight < 0.0001)
                break;
            gaussianWeight = weight * exp(-j * j / 200.0);
            weightSum += gaussianWeight;
            cplxSum += gaussianWeight * (qScale2bits(rcj, 1) - rcj->miscBits) / clippedDuration;
        }
        /* weighted average of cplx of past frames */
        weight = 1.0;
        for (int j = 0; j <= cplxBlur * 2 && j <= i; j++)
        {
            RateControlEntry *rcj = &m_rce2Pass[i - j];
            gaussianWeight = weight * exp(-j * j / 200.0);
            weightSum += gaussianWeight;
            cplxSum += gaussianWeight * (qScale2bits(rcj, 1) - rcj->miscBits) / clippedDuration;
            weight *= 1 - pow(rcj->iCuCount / m_ncu, 2);
            if (weight < .0001)
                break;
        }
        m_rce2Pass[i].blurredComplexity = cplxSum / weightSum;
    }

    CHECKED_MALLOC(qScale, double, m_numEntries);
    if (filterSize > 1)
    {
        CHECKED_MALLOC(blurredQscale, double, m_numEntries);
    }
    else
        blurredQscale = qScale;

    /* Search for a factor which, when multiplied by the RCEQ values from
     * each frame, adds up to the desired total size.
     * There is no exact closed-form solution because of VBV constraints and
     * because qscale2bits is not invertible, but we can start with the simple
     * approximation of scaling the 1st pass by the ratio of bitrates.
     * The search range is probably overkill, but speed doesn't matter here. */

    expectedBits = 1;
    for (int i = 0; i < m_numEntries; i++)
    {
        RateControlEntry* rce = &m_rce2Pass[i];
        double q = getQScale(rce, 1.0);
        expectedBits += qScale2bits(rce, q);
        m_lastQScaleFor[rce->sliceType] = q;
    }
    stepMult = allAvailableBits / expectedBits;

    rateFactor = 0;
    for (double step = 1E4 * stepMult; step > 1E-7 * stepMult; step *= 0.5)
    {
        expectedBits = 0;
        rateFactor += step;

        m_lastNonBPictType = -1;
        m_lastAccumPNorm = 1;
        m_accumPNorm = 0;

        m_lastQScaleFor[0] = m_lastQScaleFor[1] =
        m_lastQScaleFor[2] = pow(baseCplx, 1 - m_qCompress) / rateFactor;

        /* find qscale */
        for (int i = 0; i < m_numEntries; i++)
        {
            RateControlEntry *rce = &m_rce2Pass[i];
            qScale[i] = getQScale(rce, rateFactor);
            m_lastQScaleFor[rce->sliceType] = qScale[i];
        }

        /* fixed I/B qscale relative to P */
        for (int i = m_numEntries - 1; i >= 0; i--)
        {
            qScale[i] = getDiffLimitedQScale(&m_rce2Pass[i], qScale[i]);
            X265_CHECK(qScale[i] >= 0, "qScale became negative\n");
        }

        /* smooth curve */
        if (filterSize > 1)
        {
            X265_CHECK(filterSize % 2 == 1, "filterSize not an odd number\n");
            for (int i = 0; i < m_numEntries; i++)
            {
                double q = 0.0, sum = 0.0;

                for (int j = 0; j < filterSize; j++)
                {
                    int idx = i + j - filterSize / 2;
                    double d = idx - i;
                    double coeff = qBlur == 0 ? 1.0 : exp(-d * d / (qBlur * qBlur));
                    if (idx < 0 || idx >= m_numEntries)
                        continue;
                    if (m_rce2Pass[i].sliceType != m_rce2Pass[idx].sliceType)
                        continue;
                    q += qScale[idx] * coeff;
                    sum += coeff;
                }
                blurredQscale[i] = q / sum;
            }
        }

        /* find expected bits */
        for (int i = 0; i < m_numEntries; i++)
        {
            RateControlEntry *rce = &m_rce2Pass[i];
            rce->newQScale = clipQscale(NULL, rce, blurredQscale[i]); // check if needed
            X265_CHECK(rce->newQScale >= 0, "new Qscale is negative\n");
            expectedBits += qScale2bits(rce, rce->newQScale);
        }

        if (expectedBits > allAvailableBits)
            rateFactor -= step;
    }

    X265_FREE(qScale);
    if (filterSize > 1)
        X265_FREE(blurredQscale);

    if (m_isVbv)
        if (!vbv2Pass(allAvailableBits))
            return false;
    expectedBits = countExpectedBits();

    if (fabs(expectedBits / allAvailableBits - 1.0) > 0.01)
    {
        double avgq = 0;
        for (int i = 0; i < m_numEntries; i++)
            avgq += m_rce2Pass[i].newQScale;
        avgq = x265_qScale2qp(avgq / m_numEntries);

        if (expectedBits > allAvailableBits || !m_isVbv)
            x265_log(m_param, X265_LOG_WARNING, "Error: 2pass curve failed to converge\n");
        x265_log(m_param, X265_LOG_WARNING, "target: %.2f kbit/s, expected: %.2f kbit/s, avg QP: %.4f\n",
                 (double)m_param->rc.bitrate,
                 expectedBits * m_fps / (m_numEntries * 1000.),
                 avgq);
        if (expectedBits < allAvailableBits && avgq < QP_MIN + 2)
        {
            x265_log(m_param, X265_LOG_WARNING, "try reducing target bitrate\n");
        }
        else if (expectedBits > allAvailableBits && avgq > QP_MAX_SPEC - 2)
        {
            x265_log(m_param, X265_LOG_WARNING, "try increasing target bitrate\n");
        }
        else if (!(m_2pass && m_isVbv))
            x265_log(m_param, X265_LOG_WARNING, "internal error\n");
    }

    return true;

fail:
    x265_log(m_param, X265_LOG_WARNING, "two-pass ABR initialization failed\n");
    return false;
}

bool RateControl::vbv2Pass(uint64_t allAvailableBits)
{
    /* for each interval of bufferFull .. underflow, uniformly increase the qp of all
     * frames in the interval until either buffer is full at some intermediate frame or the
     * last frame in the interval no longer underflows.  Recompute intervals and repeat.
     * Then do the converse to put bits back into overflow areas until target size is met */

    double *fills;
    double expectedBits = 0;
    double adjustment;
    double prevBits = 0;
    int t0, t1;
    int iterations = 0 , adjMin, adjMax;
    CHECKED_MALLOC(fills, double, m_numEntries + 1);
    fills++;

    /* adjust overall stream size */
    do
    {
        iterations++;
        prevBits = expectedBits;

        if (expectedBits)
        {   /* not first iteration */
            adjustment = X265_MAX(X265_MIN(expectedBits / allAvailableBits, 0.999), 0.9);
            fills[-1] = m_bufferSize * m_param->rc.vbvBufferInit;
            t0 = 0;
            /* fix overflows */
            adjMin = 1;
            while (adjMin && findUnderflow(fills, &t0, &t1, 1))
            {
                adjMin = fixUnderflow(t0, t1, adjustment, MIN_QPSCALE, MAX_MAX_QPSCALE);
                t0 = t1;
            }
        }

        fills[-1] = m_bufferSize * (1. - m_param->rc.vbvBufferInit);
        t0 = 0;
        /* fix underflows -- should be done after overflow, as we'd better undersize target than underflowing VBV */
        adjMax = 1;
        while (adjMax && findUnderflow(fills, &t0, &t1, 0))
            adjMax = fixUnderflow(t0, t1, 1.001, MIN_QPSCALE, MAX_MAX_QPSCALE );

        expectedBits = countExpectedBits();
    }
    while ((expectedBits < .995 * allAvailableBits) && ((int64_t)(expectedBits+.5) > (int64_t)(prevBits+.5)));

    if (!adjMax)
        x265_log(m_param, X265_LOG_WARNING, "vbv-maxrate issue, qpmax or vbv-maxrate too low\n");

    /* store expected vbv filling values for tracking when encoding */
    for (int i = 0; i < m_numEntries; i++)
        m_rce2Pass[i].expectedVbv = m_bufferSize - fills[i];

    X265_FREE(fills - 1);
    return true;

fail:
    x265_log(m_param, X265_LOG_ERROR, "malloc failure in two-pass VBV init\n");
    return false;
}

/* In 2pass, force the same frame types as in the 1st pass */
int RateControl::rateControlSliceType(int frameNum)
{
    if (m_param->rc.bStatRead)
    {
        if (frameNum >= m_numEntries)
        {
            /* We could try to initialize everything required for ABR and
             * adaptive B-frames, but that would be complicated.
             * So just calculate the average QP used so far. */
            m_param->rc.qp = (m_accumPQp < 1) ? ABR_INIT_QP_MAX : (int)(m_accumPQp + 0.5);
            m_qpConstant[P_SLICE] = Clip3(0, QP_MAX_MAX, m_param->rc.qp);
            m_qpConstant[I_SLICE] = Clip3(0, QP_MAX_MAX, (int)(m_param->rc.qp - m_ipOffset + 0.5));
            m_qpConstant[B_SLICE] = Clip3(0, QP_MAX_MAX, (int)(m_param->rc.qp + m_pbOffset + 0.5));

            x265_log(m_param, X265_LOG_ERROR, "2nd pass has more frames than 1st pass (%d)\n", m_numEntries);
            x265_log(m_param, X265_LOG_ERROR, "continuing anyway, at constant QP=%d\n", m_param->rc.qp);
            if (m_param->bFrameAdaptive)
                x265_log(m_param, X265_LOG_ERROR, "disabling adaptive B-frames\n");

            m_isAbr = 0;
            m_2pass = 0;
            m_param->rc.rateControlMode = X265_RC_CQP;
            m_param->rc.bStatRead = 0;
            m_param->bFrameAdaptive = 0;
            m_param->scenecutThreshold = 0;
            m_param->rc.cuTree = 0;
            if (m_param->bframes > 1)
                m_param->bframes = 1;
            return X265_TYPE_AUTO;
        }
        int frameType = m_rce2Pass[frameNum].sliceType == I_SLICE ? (frameNum > 0 && m_param->bOpenGOP ? X265_TYPE_I : X265_TYPE_IDR)
                            : m_rce2Pass[frameNum].sliceType == P_SLICE ? X265_TYPE_P
                            : (m_rce2Pass[frameNum].sliceType == B_SLICE && m_rce2Pass[frameNum].keptAsRef? X265_TYPE_BREF : X265_TYPE_B);
        return frameType;
    }
    else
        return X265_TYPE_AUTO;
}

int RateControl::rateControlStart(Frame* curFrame, RateControlEntry* rce, Encoder* enc)
{
    int orderValue = m_startEndOrder.get();
    int startOrdinal = rce->encodeOrder * 2;

    while (orderValue < startOrdinal && !m_bTerminated)
       orderValue = m_startEndOrder.waitForChange(orderValue);

    if (!curFrame)
    {
        // faked rateControlStart calls when the encoder is flushing
        m_startEndOrder.incr();
        return 0;
    }

    FrameData& curEncData = *curFrame->m_encData;
    m_curSlice = curEncData.m_slice;
    m_sliceType = m_curSlice->m_sliceType;
    rce->sliceType = m_sliceType;
    rce->poc = m_curSlice->m_poc;
    if (m_param->rc.bStatRead)
    {
        X265_CHECK(rce->poc >= 0 && rce->poc < m_numEntries, "bad encode ordinal\n");
        copyRceData(rce, &m_rce2Pass[rce->poc]);
    }
    rce->isActive = true;
    if (m_sliceType == B_SLICE)
        rce->bframes = m_leadingBframes;
    else
        m_leadingBframes = curFrame->m_lowres.leadingBframes;

    rce->bLastMiniGopBFrame = curFrame->m_lowres.bLastMiniGopBFrame;
    rce->bufferRate = m_bufferRate;
    rce->rowCplxrSum = 0.0;
    rce->rowTotalBits = 0;
    if (m_isVbv)
    {
        if (rce->rowPreds[0][0].count == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    rce->rowPreds[i][j].coeff = 0.25;
                    rce->rowPreds[i][j].count = 1.0;
                    rce->rowPreds[i][j].decay = 0.5;
                    rce->rowPreds[i][j].offset = 0.0;
                }
            }
        }
        rce->rowPred[0] = &rce->rowPreds[m_sliceType][0];
        rce->rowPred[1] = &rce->rowPreds[m_sliceType][1];
        m_predictedBits = m_totalBits;
        updateVbvPlan(enc);
        rce->bufferFill = m_bufferFill;

        int mincr = enc->m_vps.ptl.minCrForLevel;
        /* Profiles above Main10 don't require maxAU size check, so just set the maximum to a large value. */
        if (enc->m_vps.ptl.profileIdc > Profile::MAIN10 || enc->m_vps.ptl.levelIdc == Level::NONE)
            rce->frameSizeMaximum = 1e9;
        else
        {
            /* The spec has a special case for the first frame. */
            if (rce->encodeOrder == 0)
            {
                /* 1.5 * (Max( PicSizeInSamplesY, fR * MaxLumaSr) + MaxLumaSr * (AuCpbRemovalTime[ 0 ] -AuNominalRemovalTime[ 0 ])) ? MinCr */
                double fr = 1. / 300;
                int picSizeInSamplesY = m_param->sourceWidth * m_param->sourceHeight;
                rce->frameSizeMaximum = 8 * 1.5 * X265_MAX(picSizeInSamplesY, fr * enc->m_vps.ptl.maxLumaSrForLevel) / mincr;
            }
            else
            {
                /* 1.5 * MaxLumaSr * (AuCpbRemovalTime[ n ] - AyCpbRemovalTime[ n - 1 ]) ? MinCr */
                rce->frameSizeMaximum = 8 * 1.5 * enc->m_vps.ptl.maxLumaSrForLevel * m_frameDuration / mincr;
            }
        }
    }
    if (m_isAbr || m_2pass) // ABR,CRF
    {
        if (m_isAbr || m_isVbv)
        {
            m_currentSatd = curFrame->m_lowres.satdCost >> (X265_DEPTH - 8);
            /* Update rce for use in rate control VBV later */
            rce->lastSatd = m_currentSatd;
        }
        double q = x265_qScale2qp(rateEstimateQscale(curFrame, rce));
        q = Clip3((double)QP_MIN, (double)QP_MAX_MAX, q);
        m_qp = int(q + 0.5);
        rce->qpaRc = curEncData.m_avgQpRc = curEncData.m_avgQpAq = q;
        /* copy value of lastRceq into thread local rce struct *to be used in RateControlEnd() */
        rce->qRceq = m_lastRceq;
        accumPQpUpdate();
    }
    else // CQP
    {
        if (m_sliceType == B_SLICE && IS_REFERENCED(curFrame))
            m_qp = (m_qpConstant[B_SLICE] + m_qpConstant[P_SLICE]) / 2;
        else
            m_qp = m_qpConstant[m_sliceType];
        curEncData.m_avgQpAq = curEncData.m_avgQpRc = m_qp;
    }
    if (m_sliceType != B_SLICE)
    {
        m_lastNonBPictType = m_sliceType;
        m_leadingNoBSatd = m_currentSatd;
    }
    rce->leadingNoBSatd = m_leadingNoBSatd;
    if (curFrame->m_forceqp)
    {
        m_qp = int32_t(curFrame->m_forceqp + 0.5) - 1;
        m_qp = Clip3(QP_MIN, QP_MAX_MAX, m_qp);
        rce->qpaRc = curEncData.m_avgQpRc = curEncData.m_avgQpAq = m_qp;
    }
    // Do not increment m_startEndOrder here. Make rateControlEnd of previous thread
    // to wait until rateControlUpdateStats of this frame is called
    m_framesDone++;
    return m_qp;
}

void RateControl::accumPQpUpdate()
{
    m_accumPQp   *= .95;
    m_accumPNorm *= .95;
    m_accumPNorm += 1;
    if (m_sliceType == I_SLICE)
        m_accumPQp += m_qp + m_ipOffset;
    else
        m_accumPQp += m_qp;
}

double RateControl::getDiffLimitedQScale(RateControlEntry *rce, double q)
{
    // force I/B quants as a function of P quants
    const double lastPqScale    = m_lastQScaleFor[P_SLICE];
    const double lastNonBqScale = m_lastQScaleFor[m_lastNonBPictType];
    if (rce->sliceType == I_SLICE)
    {
        double iq = q;
        double pq = x265_qp2qScale(m_accumPQp / m_accumPNorm);
        double ipFactor = fabs(m_param->rc.ipFactor);
        /* don't apply ipFactor if the following frame is also I */
        if (m_accumPNorm <= 0)
            q = iq;
        else if (m_param->rc.ipFactor < 0)
            q = iq / ipFactor;
        else if (m_accumPNorm >= 1)
            q = pq / ipFactor;
        else
            q = m_accumPNorm * pq / ipFactor + (1 - m_accumPNorm) * iq;
    }
    else if (rce->sliceType == B_SLICE)
    {
        if (m_param->rc.pbFactor > 0)
            q = lastNonBqScale;
        if (!rce->keptAsRef)
            q *= fabs(m_param->rc.pbFactor);
    }
    else if (rce->sliceType == P_SLICE
             && m_lastNonBPictType == P_SLICE
             && rce->coeffBits == 0)
    {
        q = lastPqScale;
    }

    /* last qscale / qdiff stuff */
    if (m_lastNonBPictType == rce->sliceType &&
        (rce->sliceType != I_SLICE || m_lastAccumPNorm < 1))
    {
        double maxQscale = m_lastQScaleFor[rce->sliceType] * m_lstep;
        double minQscale = m_lastQScaleFor[rce->sliceType] / m_lstep;
        q = Clip3(minQscale, maxQscale, q);
    }

    m_lastQScaleFor[rce->sliceType] = q;
    if (rce->sliceType != B_SLICE)
        m_lastNonBPictType = rce->sliceType;
    if (rce->sliceType == I_SLICE)
    {
        m_lastAccumPNorm = m_accumPNorm;
        m_accumPNorm = 0;
        m_accumPQp = 0;
    }
    if (rce->sliceType == P_SLICE)
    {
        double mask = 1 - pow(rce->iCuCount / m_ncu, 2);
        m_accumPQp   = mask * (x265_qScale2qp(q) + m_accumPQp);
        m_accumPNorm = mask * (1 + m_accumPNorm);
    }

    return q;
}

double RateControl::countExpectedBits()
{
    double expectedBits = 0;
    for( int i = 0; i < m_numEntries; i++ )
    {
        RateControlEntry *rce = &m_rce2Pass[i];
        rce->expectedBits = (uint64_t)expectedBits;
        expectedBits += qScale2bits(rce, rce->newQScale);
    }
    return expectedBits;
}

bool RateControl::findUnderflow(double *fills, int *t0, int *t1, int over)
{
    /* find an interval ending on an overflow or underflow (depending on whether
     * we're adding or removing bits), and starting on the earliest frame that
     * can influence the buffer fill of that end frame. */
    const double bufferMin = .1 * m_bufferSize;
    const double bufferMax = .9 * m_bufferSize;
    double fill = fills[*t0 - 1];
    double parity = over ? 1. : -1.;
    int start = -1, end = -1;
    for (int i = *t0; i < m_numEntries; i++)
    {
        fill += (m_frameDuration * m_vbvMaxRate -
                 qScale2bits(&m_rce2Pass[i], m_rce2Pass[i].newQScale)) * parity;
        fill = Clip3(0.0, m_bufferSize, fill);
        fills[i] = fill;
        if (fill <= bufferMin || i == 0)
        {
            if (end >= 0)
                break;
            start = i;
        }
        else if (fill >= bufferMax && start >= 0)
            end = i;
    }
    *t0 = start;
    *t1 = end;
    return start >= 0 && end >= 0;
}

bool RateControl::fixUnderflow(int t0, int t1, double adjustment, double qscaleMin, double qscaleMax)
{
    double qscaleOrig, qscaleNew;
    bool adjusted = false;
    if (t0 > 0)
        t0++;
    for (int i = t0; i <= t1; i++)
    {
        qscaleOrig = m_rce2Pass[i].newQScale;
        qscaleOrig = Clip3(qscaleMin, qscaleMax, qscaleOrig);
        qscaleNew  = qscaleOrig * adjustment;
        qscaleNew  = Clip3(qscaleMin, qscaleMax, qscaleNew);
        m_rce2Pass[i].newQScale = qscaleNew;
        adjusted = adjusted || (qscaleNew != qscaleOrig);
    }
    return adjusted;
}

bool RateControl::cuTreeReadFor2Pass(Frame* frame)
{
    uint8_t sliceTypeActual = (uint8_t)m_rce2Pass[frame->m_poc].sliceType;

    if (m_rce2Pass[frame->m_poc].keptAsRef)
    {
        uint8_t type;
        if (m_cuTreeStats.qpBufPos < 0)
        {
            do
            {
                m_cuTreeStats.qpBufPos++;

                if (!fread(&type, 1, 1, m_cutreeStatFileIn))
                    goto fail;
                if (fread(m_cuTreeStats.qpBuffer[m_cuTreeStats.qpBufPos], sizeof(uint16_t), m_ncu, m_cutreeStatFileIn) != (size_t)m_ncu)
                    goto fail;

                if (type != sliceTypeActual && m_cuTreeStats.qpBufPos == 1)
                {
                    x265_log(m_param, X265_LOG_ERROR, "CU-tree frametype %d doesn't match actual frametype %d.\n", type, sliceTypeActual);
                    return false;
                }
            }
            while(type != sliceTypeActual);
        }
        for (int i = 0; i < m_ncu; i++)
        {
            int16_t qpFix8 = m_cuTreeStats.qpBuffer[m_cuTreeStats.qpBufPos][i];
            frame->m_lowres.qpCuTreeOffset[i] = (double)(qpFix8) / 256.0;
            frame->m_lowres.invQscaleFactor[i] = x265_exp2fix8(frame->m_lowres.qpCuTreeOffset[i]);
        }
        m_cuTreeStats.qpBufPos--;
    }
    else
        calcAdaptiveQuantFrame(frame);
    return true;

fail:
    x265_log(m_param, X265_LOG_ERROR, "Incomplete CU-tree stats file.\n");
    return false;
}

double RateControl::rateEstimateQscale(Frame* curFrame, RateControlEntry *rce)
{
    double q;

    if (m_2pass)
    {
        if (m_sliceType != rce->sliceType)
        {
            x265_log(m_param, X265_LOG_ERROR, "slice=%c but 2pass stats say %c\n",
                     g_sliceTypeToChar[m_sliceType], g_sliceTypeToChar[rce->sliceType]);
        }
    }
    else
    {
        if (m_isAbr)
        {
            double slidingWindowCplxSum = 0;
            int start = m_sliderPos > s_slidingWindowFrames ?  m_sliderPos : 0;
            for (int cnt = 0; cnt < s_slidingWindowFrames; cnt++, start++)
            {
                int pos = start % s_slidingWindowFrames;
                slidingWindowCplxSum *= 0.5;
                if (!m_satdCostWindow[pos])
                    break;
                slidingWindowCplxSum += m_satdCostWindow[pos] / (CLIP_DURATION(m_frameDuration) / BASE_FRAME_DURATION);
            }
            rce->movingAvgSum = slidingWindowCplxSum;
            m_satdCostWindow[m_sliderPos % s_slidingWindowFrames] = rce->lastSatd;
            m_sliderPos++;
        }
    }

    if (m_sliceType == B_SLICE)
    {
        /* B-frames don't have independent rate control, but rather get the
         * average QP of the two adjacent P-frames + an offset */
        Slice* prevRefSlice = m_curSlice->m_refPicList[0][0]->m_encData->m_slice;
        Slice* nextRefSlice = m_curSlice->m_refPicList[1][0]->m_encData->m_slice;
        double q0 = m_curSlice->m_refPicList[0][0]->m_encData->m_avgQpRc;
        double q1 = m_curSlice->m_refPicList[1][0]->m_encData->m_avgQpRc;
        bool i0 = prevRefSlice->m_sliceType == I_SLICE;
        bool i1 = nextRefSlice->m_sliceType == I_SLICE;
        int dt0 = abs(m_curSlice->m_poc - prevRefSlice->m_poc);
        int dt1 = abs(m_curSlice->m_poc - nextRefSlice->m_poc);

        // Skip taking a reference frame before the Scenecut if ABR has been reset.
        if (m_lastAbrResetPoc >= 0)
        {
            if (prevRefSlice->m_sliceType == P_SLICE && prevRefSlice->m_poc < m_lastAbrResetPoc)
            {
                i0 = i1;
                dt0 = dt1;
                q0 = q1;
            }
        }
        if (prevRefSlice->m_sliceType == B_SLICE && IS_REFERENCED(m_curSlice->m_refPicList[0][0]))
            q0 -= m_pbOffset / 2;
        if (nextRefSlice->m_sliceType == B_SLICE && IS_REFERENCED(m_curSlice->m_refPicList[1][0]))
            q1 -= m_pbOffset / 2;
        if (i0 && i1)
            q = (q0 + q1) / 2 + m_ipOffset;
        else if (i0)
            q = q1;
        else if (i1)
            q = q0;
        else
            q = (q0 * dt1 + q1 * dt0) / (dt0 + dt1);

        if (IS_REFERENCED(curFrame))
            q += m_pbOffset / 2;
        else
            q += m_pbOffset;
        rce->qpNoVbv = q;
        double qScale = x265_qp2qScale(q);

        if (!m_2pass && m_isVbv)
        {
            if (m_leadingBframes > 5)
            {
                qScale = clipQscale(curFrame, rce, qScale);
                m_lastQScaleFor[m_sliceType] = qScale;
            }
            rce->frameSizePlanned = predictSize(&m_predBfromP, qScale, (double)m_leadingNoBSatd);
        }
        else if (m_2pass && m_isVbv)
        {
            rce->frameSizePlanned = qScale2bits(rce, qScale);
        }
        /* Limit planned size by MinCR */
        if (m_isVbv)
            rce->frameSizePlanned = X265_MIN(rce->frameSizePlanned, rce->frameSizeMaximum);
        rce->frameSizeEstimated = rce->frameSizePlanned;
        rce->newQScale = qScale;
        return qScale;
    }
    else
    {
        double abrBuffer = 2 * m_param->rc.rateTolerance * m_bitrate;
        if (m_2pass)
        {
            int64_t diff;
            if (!m_isVbv)
            {
                m_predictedBits = m_totalBits;
                if (rce->encodeOrder < m_param->frameNumThreads)
                    m_predictedBits += (int64_t)(rce->encodeOrder * m_bitrate / m_fps);
                else
                    m_predictedBits += (int64_t)(m_param->frameNumThreads * m_bitrate / m_fps);
            }
            /* Adjust ABR buffer based on distance to the end of the video. */
            if (m_numEntries > rce->encodeOrder)
            {
                uint64_t finalBits = m_rce2Pass[m_numEntries - 1].expectedBits;
                double videoPos = (double)rce->expectedBits / finalBits;
                double scaleFactor = sqrt((1 - videoPos) * m_numEntries);
                abrBuffer *= 0.5 * X265_MAX(scaleFactor, 0.5);
            }
            diff = m_predictedBits - (int64_t)rce->expectedBits;
            q = rce->newQScale;
            q /= Clip3(0.5, 2.0, (double)(abrBuffer - diff) / abrBuffer);
            if (m_expectedBitsSum > 0)
            {
                /* Adjust quant based on the difference between
                 * achieved and expected bitrate so far */
                double curTime = (double)rce->encodeOrder / m_numEntries;
                double w = Clip3(0.0, 1.0, curTime * 100);
                q *= pow((double)m_totalBits / m_expectedBitsSum, w);
            }
            rce->qpNoVbv = x265_qScale2qp(q);
            if (m_isVbv)
            {
                /* Do not overflow vbv */
                double expectedSize = qScale2bits(rce, q);
                double expectedVbv = m_bufferFill + m_bufferRate - expectedSize;
                double expectedFullness = rce->expectedVbv / m_bufferSize;
                double qmax = q * (2 - expectedFullness);
                double sizeConstraint = 1 + expectedFullness;
                qmax = X265_MAX(qmax, rce->newQScale);
                if (expectedFullness < .05)
                    qmax = MAX_MAX_QPSCALE;
                qmax = X265_MIN(qmax, MAX_MAX_QPSCALE);
                while (((expectedVbv < rce->expectedVbv/sizeConstraint) && (q < qmax)) ||
                        ((expectedVbv < 0) && (q < MAX_MAX_QPSCALE)))
                {
                    q *= 1.05;
                    expectedSize = qScale2bits(rce, q);
                    expectedVbv = m_bufferFill + m_bufferRate - expectedSize;
                }
            }
            q = Clip3(MIN_QPSCALE, MAX_MAX_QPSCALE, q);
        }
        else
        {
            /* 1pass ABR */

            /* Calculate the quantizer which would have produced the desired
             * average bitrate if it had been applied to all frames so far.
             * Then modulate that quant based on the current frame's complexity
             * relative to the average complexity so far (using the 2pass RCEQ).
             * Then bias the quant up or down if total size so far was far from
             * the target.
             * Result: Depending on the value of rate_tolerance, there is a
             * tradeoff between quality and bitrate precision. But at large
             * tolerances, the bit distribution approaches that of 2pass. */

            double wantedBits, overflow = 1;

            m_shortTermCplxSum *= 0.5;
            m_shortTermCplxCount *= 0.5;
            m_shortTermCplxSum += m_currentSatd / (CLIP_DURATION(m_frameDuration) / BASE_FRAME_DURATION);
            m_shortTermCplxCount++;
            /* coeffBits to be used in 2-pass */
            rce->coeffBits = (int)m_currentSatd;
            rce->blurredComplexity = m_shortTermCplxSum / m_shortTermCplxCount;
            rce->mvBits = 0;
            rce->sliceType = m_sliceType;

            if (m_param->rc.rateControlMode == X265_RC_CRF)
            {
                q = getQScale(rce, m_rateFactorConstant);
            }
            else
            {
                if (!m_param->rc.bStatRead)
                    checkAndResetABR(rce, false);
                q = getQScale(rce, m_wantedBitsWindow / m_cplxrSum);

                /* ABR code can potentially be counterproductive in CBR, so just
                 * don't bother.  Don't run it if the frame complexity is zero
                 * either. */
                if (!m_isCbr && m_currentSatd)
                {
                    /* use framesDone instead of POC as poc count is not serial with bframes enabled */
                    double timeDone = (double)(m_framesDone - m_param->frameNumThreads + 1) * m_frameDuration;
                    wantedBits = timeDone * m_bitrate;
                    if (wantedBits > 0 && m_totalBits > 0 && !m_partialResidualFrames)
                    {
                        abrBuffer *= X265_MAX(1, sqrt(timeDone));
                        overflow = Clip3(.5, 2.0, 1.0 + (m_totalBits - wantedBits) / abrBuffer);
                        q *= overflow;
                    }
                }
            }

            if (m_sliceType == I_SLICE && m_param->keyframeMax > 1
                && m_lastNonBPictType != I_SLICE && !m_isAbrReset)
            {
                q = x265_qp2qScale(m_accumPQp / m_accumPNorm);
                q /= fabs(m_param->rc.ipFactor);
            }
            else if (m_framesDone > 0)
            {
                if (m_param->rc.rateControlMode != X265_RC_CRF)
                {
                    double lqmin = 0, lqmax = 0;
                    lqmin = m_lastQScaleFor[m_sliceType] / m_lstep;
                    lqmax = m_lastQScaleFor[m_sliceType] * m_lstep;
                    if (!m_partialResidualFrames)
                    {
                        if (overflow > 1.1 && m_framesDone > 3)
                            lqmax *= m_lstep;
                        else if (overflow < 0.9)
                            lqmin /= m_lstep;
                    }
                    q = Clip3(lqmin, lqmax, q);
                }
            }
            else if (m_qCompress != 1 && m_param->rc.rateControlMode == X265_RC_CRF)
            {
                q = x265_qp2qScale(CRF_INIT_QP) / fabs(m_param->rc.ipFactor);
            }
            else if (m_framesDone == 0 && !m_isVbv)
            {
                /* for ABR alone, clip the first I frame qp */
                double lqmax = x265_qp2qScale(ABR_INIT_QP_MAX) * m_lstep;
                q = X265_MIN(lqmax, q);
            }
            q = Clip3(MIN_QPSCALE, MAX_MAX_QPSCALE, q);
            rce->qpNoVbv = x265_qScale2qp(q);
            q = clipQscale(curFrame, rce, q);
        }
        m_lastQScaleFor[m_sliceType] = q;
        if ((m_curSlice->m_poc == 0 || m_lastQScaleFor[P_SLICE] < q) && !(m_2pass && !m_isVbv))
            m_lastQScaleFor[P_SLICE] = q * fabs(m_param->rc.ipFactor);

        if (m_2pass && m_isVbv)
            rce->frameSizePlanned = qScale2bits(rce, q);
        else
            rce->frameSizePlanned = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);

        /* Always use up the whole VBV in this case. */
        if (m_singleFrameVbv)
            rce->frameSizePlanned = m_bufferRate;
        /* Limit planned size by MinCR */
        if (m_isVbv)
            rce->frameSizePlanned = X265_MIN(rce->frameSizePlanned, rce->frameSizeMaximum);
        rce->frameSizeEstimated = rce->frameSizePlanned;
        rce->newQScale = q;
        return q;
    }
}

void RateControl::rateControlUpdateStats(RateControlEntry* rce)
{
    if (!m_param->rc.bStatWrite && !m_param->rc.bStatRead)
    {
        if (rce->sliceType == I_SLICE)
        {
            /* previous I still had a residual; roll it into the new loan */
            if (m_partialResidualFrames)
                rce->rowTotalBits += m_partialResidualCost * m_partialResidualFrames;

            m_partialResidualFrames = X265_MIN(s_amortizeFrames, m_param->keyframeMax);
            m_partialResidualCost = (int)((rce->rowTotalBits * s_amortizeFraction) /m_partialResidualFrames);
            rce->rowTotalBits -= m_partialResidualCost * m_partialResidualFrames;
        }
        else if (m_partialResidualFrames)
        {
             rce->rowTotalBits += m_partialResidualCost;
             m_partialResidualFrames--;
        }
    }
    if (rce->sliceType != B_SLICE)
        rce->rowCplxrSum = rce->rowTotalBits * x265_qp2qScale(rce->qpaRc) / rce->qRceq;
    else
        rce->rowCplxrSum = rce->rowTotalBits * x265_qp2qScale(rce->qpaRc) / (rce->qRceq * fabs(m_param->rc.pbFactor));

    m_cplxrSum += rce->rowCplxrSum;
    m_totalBits += rce->rowTotalBits;

    /* do not allow the next frame to enter rateControlStart() until this
     * frame has updated its mid-frame statistics */
    m_startEndOrder.incr();

    if (rce->encodeOrder < m_param->frameNumThreads - 1)
        m_startEndOrder.incr(); // faked rateControlEnd calls for negative frames
}

void RateControl::checkAndResetABR(RateControlEntry* rce, bool isFrameDone)
{
    double abrBuffer = 2 * m_param->rc.rateTolerance * m_bitrate;

    // Check if current Slice is a scene cut that follows low detailed/blank frames
    if (rce->lastSatd > 4 * rce->movingAvgSum)
    {
        if (!m_isAbrReset && rce->movingAvgSum > 0)
        {
            int64_t shrtTermWantedBits = (int64_t) (X265_MIN(m_sliderPos, s_slidingWindowFrames) * m_bitrate * m_frameDuration);
            int64_t shrtTermTotalBitsSum = 0;
            // Reset ABR if prev frames are blank to prevent further sudden overflows/ high bit rate spikes.
            for (int i = 0; i < s_slidingWindowFrames ; i++)
                shrtTermTotalBitsSum += m_encodedBitsWindow[i];
            double underflow = (shrtTermTotalBitsSum - shrtTermWantedBits) / abrBuffer;
            const double epsilon = 0.0001f;
            if (underflow < epsilon && !isFrameDone)
            {
                init(m_curSlice->m_sps);
                m_shortTermCplxSum = rce->lastSatd / (CLIP_DURATION(m_frameDuration) / BASE_FRAME_DURATION);
                m_shortTermCplxCount = 1;
                m_isAbrReset = true;
                m_lastAbrResetPoc = rce->poc;
            }
        }
        else
        {
            // Clear flag to reset ABR and continue as usual.
            m_isAbrReset = false;
        }
    }
}

void RateControl::hrdFullness(SEIBufferingPeriod *seiBP)
{
    const VUI* vui = &m_curSlice->m_sps->vuiParameters;
    const HRDInfo* hrd = &vui->hrdParameters;
    int num = 90000;
    int denom = hrd->bitRateValue << (hrd->bitRateScale + BR_SHIFT);
    reduceFraction(&num, &denom);
    int64_t cpbState = (int64_t)m_bufferFillFinal;
    int64_t cpbSize = (int64_t)hrd->cpbSizeValue << (hrd->cpbSizeScale + CPB_SHIFT);

    if (cpbState < 0 || cpbState > cpbSize)
    {
        x265_log(m_param, X265_LOG_WARNING, "CPB %s: %.0lf bits in a %.0lf-bit buffer\n",
                 cpbState < 0 ? "underflow" : "overflow", (float)cpbState/denom, (float)cpbSize/denom);
    }

    seiBP->m_initialCpbRemovalDelay = (uint32_t)(num * cpbState + denom) / denom;
    seiBP->m_initialCpbRemovalDelayOffset = (uint32_t)(num * cpbSize + denom) / denom - seiBP->m_initialCpbRemovalDelay;
}

void RateControl::updateVbvPlan(Encoder* enc)
{
    m_bufferFill = m_bufferFillFinal;
    enc->updateVbvPlan(this);
}

double RateControl::predictSize(Predictor *p, double q, double var)
{
    return (p->coeff * var + p->offset) / (q * p->count);
}

double RateControl::clipQscale(Frame* curFrame, RateControlEntry* rce, double q)
{
    // B-frames are not directly subject to VBV,
    // since they are controlled by referenced P-frames' QPs.
    double q0 = q;
    if (m_isVbv && m_currentSatd > 0 && curFrame)
    {
        if (m_param->lookaheadDepth || m_param->rc.cuTree ||
            m_param->scenecutThreshold ||
            (m_param->bFrameAdaptive && m_param->bframes))
        {
           /* Lookahead VBV: If lookahead is done, raise the quantizer as necessary
            * such that no frames in the lookahead overflow and such that the buffer
            * is in a reasonable state by the end of the lookahead. */
            int loopTerminate = 0;
            /* Avoid an infinite loop. */
            for (int iterations = 0; iterations < 1000 && loopTerminate != 3; iterations++)
            {
                double frameQ[3];
                double curBits;
                if (m_sliceType == B_SLICE)
                    curBits = predictSize(&m_predBfromP, q, (double)m_currentSatd);
                else
                    curBits = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);
                double bufferFillCur = m_bufferFill - curBits;
                double targetFill;
                double totalDuration = 0;
                frameQ[P_SLICE] = m_sliceType == I_SLICE ? q * m_param->rc.ipFactor : (m_sliceType == B_SLICE ? q / m_param->rc.pbFactor : q);
                frameQ[B_SLICE] = frameQ[P_SLICE] * m_param->rc.pbFactor;
                frameQ[I_SLICE] = frameQ[P_SLICE] / m_param->rc.ipFactor;
                /* Loop over the planned future frames. */
                for (int j = 0; bufferFillCur >= 0; j++)
                {
                    int type = curFrame->m_lowres.plannedType[j];
                    if (type == X265_TYPE_AUTO)
                        break;
                    totalDuration += m_frameDuration;
                    double wantedFrameSize = m_vbvMaxRate * m_frameDuration;
                    if (bufferFillCur + wantedFrameSize <= m_bufferSize)
                        bufferFillCur += wantedFrameSize;
                    int64_t satd = curFrame->m_lowres.plannedSatd[j] >> (X265_DEPTH - 8);
                    type = IS_X265_TYPE_I(type) ? I_SLICE : IS_X265_TYPE_B(type) ? B_SLICE : P_SLICE;
                    curBits = predictSize(&m_pred[type], frameQ[type], (double)satd);
                    bufferFillCur -= curBits;
                }

                /* Try to get the buffer at least 50% filled, but don't set an impossible goal. */
                targetFill = X265_MIN(m_bufferFill + totalDuration * m_vbvMaxRate * 0.5, m_bufferSize * 0.5);
                if (bufferFillCur < targetFill)
                {
                    q *= 1.01;
                    loopTerminate |= 1;
                    continue;
                }
                /* Try to get the buffer no more than 80% filled, but don't set an impossible goal. */
                targetFill = Clip3(m_bufferSize * 0.8, m_bufferSize, m_bufferFill - totalDuration * m_vbvMaxRate * 0.5);
                if (m_isCbr && bufferFillCur > targetFill)
                {
                    q /= 1.01;
                    loopTerminate |= 2;
                    continue;
                }
                break;
            }
            q = X265_MAX(q0 / 2, q);
        }
        else
        {
            /* Fallback to old purely-reactive algorithm: no lookahead. */
            if ((m_sliceType == P_SLICE || m_sliceType == B_SLICE ||
                    (m_sliceType == I_SLICE && m_lastNonBPictType == I_SLICE)) &&
                m_bufferFill / m_bufferSize < 0.5)
            {
                q /= Clip3(0.5, 1.0, 2.0 * m_bufferFill / m_bufferSize);
            }
            // Now a hard threshold to make sure the frame fits in VBV.
            // This one is mostly for I-frames.
            double bits = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);

            // For small VBVs, allow the frame to use up the entire VBV.
            double maxFillFactor;
            maxFillFactor = m_bufferSize >= 5 * m_bufferRate ? 2 : 1;
            // For single-frame VBVs, request that the frame use up the entire VBV.
            double minFillFactor = m_singleFrameVbv ? 1 : 2;

            for (int iterations = 0; iterations < 10; iterations++)
            {
                double qf = 1.0;
                if (bits > m_bufferFill / maxFillFactor)
                    qf = Clip3(0.2, 1.0, m_bufferFill / (maxFillFactor * bits));
                q /= qf;
                bits *= qf;
                if (bits < m_bufferRate / minFillFactor)
                    q *= bits * minFillFactor / m_bufferRate;
                bits = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);
            }

            q = X265_MAX(q0, q);
        }

        /* Apply MinCR restrictions */
        double pbits = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);
        if (pbits > rce->frameSizeMaximum)
            q *= pbits / rce->frameSizeMaximum;

        // Check B-frame complexity, and use up any bits that would
        // overflow before the next P-frame.
        if (m_leadingBframes <= 5 && m_sliceType == P_SLICE && !m_singleFrameVbv)
        {
            int nb = m_leadingBframes;
            double bits = predictSize(&m_pred[m_sliceType], q, (double)m_currentSatd);
            double bbits = predictSize(&m_predBfromP, q * m_param->rc.pbFactor, (double)m_currentSatd);
            double space;
            if (bbits > m_bufferRate)
                nb = 0;
            double pbbits = nb * bbits;

            space = m_bufferFill + (1 + nb) * m_bufferRate - m_bufferSize;
            if (pbbits < space)
                q *= X265_MAX(pbbits / space, bits / (0.5 * m_bufferSize));

            q = X265_MAX(q0 / 2, q);
        }

        if (!m_isCbr || (m_isAbr && m_currentSatd >= rce->movingAvgSum && q <= q0 / 2))
            q = X265_MAX(q0, q);

        if (m_rateFactorMaxIncrement)
        {
            double qpNoVbv = x265_qScale2qp(q0);
            double qmax = X265_MIN(MAX_MAX_QPSCALE,x265_qp2qScale(qpNoVbv + m_rateFactorMaxIncrement));
            return Clip3(MIN_QPSCALE, qmax, q);
        }
    }
    if (m_2pass)
    {
        double min = log(MIN_QPSCALE);
        double max = log(MAX_MAX_QPSCALE);
        q = (log(q) - min) / (max - min) - 0.5;
        q = 1.0 / (1.0 + exp(-4 * q));
        q = q*(max - min) + min;
        return exp(q);
    }
    return Clip3(MIN_QPSCALE, MAX_MAX_QPSCALE, q);
}

double RateControl::predictRowsSizeSum(Frame* curFrame, RateControlEntry* rce, double qpVbv, int32_t& encodedBitsSoFar)
{
    uint32_t rowSatdCostSoFar = 0, totalSatdBits = 0;
    encodedBitsSoFar = 0;

    double qScale = x265_qp2qScale(qpVbv);
    FrameData& curEncData = *curFrame->m_encData;
    int picType = curEncData.m_slice->m_sliceType;
    Frame* refFrame = curEncData.m_slice->m_refPicList[0][0];

    uint32_t maxRows = curEncData.m_slice->m_sps->numCuInHeight;
    uint32_t maxCols = curEncData.m_slice->m_sps->numCuInWidth;

    for (uint32_t row = 0; row < maxRows; row++)
    {
        encodedBitsSoFar += curEncData.m_rowStat[row].encodedBits;
        rowSatdCostSoFar = curEncData.m_rowStat[row].diagSatd;
        uint32_t satdCostForPendingCus = curEncData.m_rowStat[row].satdForVbv - rowSatdCostSoFar;
        satdCostForPendingCus >>= X265_DEPTH - 8;
        if (satdCostForPendingCus  > 0)
        {
            double pred_s = predictSize(rce->rowPred[0], qScale, satdCostForPendingCus);
            uint32_t refRowSatdCost = 0, refRowBits = 0, intraCost = 0;
            double refQScale = 0;

            if (picType != I_SLICE)
            {
                FrameData& refEncData = *refFrame->m_encData;
                uint32_t endCuAddr = maxCols * (row + 1);
                for (uint32_t cuAddr = curEncData.m_rowStat[row].numEncodedCUs + 1; cuAddr < endCuAddr; cuAddr++)
                {
                    refRowSatdCost += refEncData.m_cuStat[cuAddr].vbvCost;
                    refRowBits += refEncData.m_cuStat[cuAddr].totalBits;
                    intraCost += curEncData.m_cuStat[cuAddr].intraVbvCost;
                }

                refRowSatdCost >>= X265_DEPTH - 8;
                refQScale = refEncData.m_rowStat[row].diagQpScale;
            }

            if (picType == I_SLICE || qScale >= refQScale)
            {
                if (picType == P_SLICE 
                    && !refFrame 
                    && refFrame->m_encData->m_slice->m_sliceType == picType
                    && refQScale > 0
                    && refRowSatdCost > 0)
                {
                    if (abs(int32_t(refRowSatdCost - satdCostForPendingCus)) < (int32_t)satdCostForPendingCus / 2)
                    {
                        double predTotal = refRowBits * satdCostForPendingCus / refRowSatdCost * refQScale / qScale;
                        totalSatdBits += int32_t((pred_s + predTotal) * 0.5);
                        continue;
                    }
                }
                totalSatdBits += int32_t(pred_s);
            }
            else
            {
                /* Our QP is lower than the reference! */
                double pred_intra = predictSize(rce->rowPred[1], qScale, intraCost);
                /* Sum: better to overestimate than underestimate by using only one of the two predictors. */
                totalSatdBits += int32_t(pred_intra + pred_s);
            }
        }
    }

    return totalSatdBits + encodedBitsSoFar;
}

int RateControl::rowDiagonalVbvRateControl(Frame* curFrame, uint32_t row, RateControlEntry* rce, double& qpVbv)
{
    FrameData& curEncData = *curFrame->m_encData;
    double qScaleVbv = x265_qp2qScale(qpVbv);
    uint64_t rowSatdCost = curEncData.m_rowStat[row].diagSatd;
    double encodedBits = curEncData.m_rowStat[row].encodedBits;

    if (row == 1)
    {
        rowSatdCost += curEncData.m_rowStat[0].diagSatd;
        encodedBits += curEncData.m_rowStat[0].encodedBits;
    }
    rowSatdCost >>= X265_DEPTH - 8;
    updatePredictor(rce->rowPred[0], qScaleVbv, (double)rowSatdCost, encodedBits);
    if (curEncData.m_slice->m_sliceType == P_SLICE)
    {
        Frame* refFrame = curEncData.m_slice->m_refPicList[0][0];
        if (qpVbv < refFrame->m_encData->m_rowStat[row].diagQp)
        {
            uint64_t intraRowSatdCost = curEncData.m_rowStat[row].diagIntraSatd;
            if (row == 1)
                intraRowSatdCost += curEncData.m_rowStat[0].diagIntraSatd;

            updatePredictor(rce->rowPred[1], qScaleVbv, (double)intraRowSatdCost, encodedBits);
        }
    }

    int canReencodeRow = 1;
    /* tweak quality based on difference from predicted size */
    double prevRowQp = qpVbv;
    double qpAbsoluteMax = QP_MAX_MAX;
    double qpAbsoluteMin = QP_MIN;
    if (m_rateFactorMaxIncrement)
        qpAbsoluteMax = X265_MIN(qpAbsoluteMax, rce->qpNoVbv + m_rateFactorMaxIncrement);

    if (m_rateFactorMaxDecrement)
        qpAbsoluteMin = X265_MAX(qpAbsoluteMin, rce->qpNoVbv - m_rateFactorMaxDecrement);

    double qpMax = X265_MIN(prevRowQp + m_param->rc.qpStep, qpAbsoluteMax);
    double qpMin = X265_MAX(prevRowQp - m_param->rc.qpStep, qpAbsoluteMin);
    double stepSize = 0.5;
    double bufferLeftPlanned = rce->bufferFill - rce->frameSizePlanned;

    const SPS& sps = *curEncData.m_slice->m_sps;
    double maxFrameError = X265_MAX(0.05, 1.0 / sps.numCuInHeight);

    if (row < sps.numCuInHeight - 1)
    {
        /* B-frames shouldn't use lower QP than their reference frames. */
        if (rce->sliceType == B_SLICE)
        {
            Frame* refSlice1 = curEncData.m_slice->m_refPicList[0][0];
            Frame* refSlice2 = curEncData.m_slice->m_refPicList[1][0];
            qpMin = X265_MAX(qpMin, X265_MAX(refSlice1->m_encData->m_rowStat[row].diagQp, refSlice2->m_encData->m_rowStat[row].diagQp));
            qpVbv = X265_MAX(qpVbv, qpMin);
        }
        /* More threads means we have to be more cautious in letting ratecontrol use up extra bits. */
        double rcTol = bufferLeftPlanned / m_param->frameNumThreads * m_param->rc.rateTolerance;
        int32_t encodedBitsSoFar = 0;
        double accFrameBits = predictRowsSizeSum(curFrame, rce, qpVbv, encodedBitsSoFar);

        /* * Don't increase the row QPs until a sufficent amount of the bits of
         * the frame have been processed, in case a flat area at the top of the
         * frame was measured inaccurately. */
        if (encodedBitsSoFar < 0.05f * rce->frameSizePlanned)
            qpMax = qpAbsoluteMax = prevRowQp;

        if (rce->sliceType != I_SLICE)
            rcTol *= 0.5;

        if (!m_isCbr)
            qpMin = X265_MAX(qpMin, rce->qpNoVbv);

        while (qpVbv < qpMax
               && ((accFrameBits > rce->frameSizePlanned + rcTol) ||
                   (rce->bufferFill - accFrameBits < bufferLeftPlanned * 0.5) ||
                   (accFrameBits > rce->frameSizePlanned && qpVbv < rce->qpNoVbv)))
        {
            qpVbv += stepSize;
            accFrameBits = predictRowsSizeSum(curFrame, rce, qpVbv, encodedBitsSoFar);
        }

        while (qpVbv > qpMin
               && (qpVbv > curEncData.m_rowStat[0].diagQp || m_singleFrameVbv)
               && ((accFrameBits < rce->frameSizePlanned * 0.8f && qpVbv <= prevRowQp)
                   || accFrameBits < (rce->bufferFill - m_bufferSize + m_bufferRate) * 1.1))
        {
            qpVbv -= stepSize;
            accFrameBits = predictRowsSizeSum(curFrame, rce, qpVbv, encodedBitsSoFar);
        }

        /* avoid VBV underflow or MinCr violation */
        while ((qpVbv < qpAbsoluteMax)
               && ((rce->bufferFill - accFrameBits < m_bufferRate * maxFrameError) ||
                   (rce->frameSizeMaximum - accFrameBits < rce->frameSizeMaximum * maxFrameError)))
        {
            qpVbv += stepSize;
            accFrameBits = predictRowsSizeSum(curFrame, rce, qpVbv, encodedBitsSoFar);
        }

        rce->frameSizeEstimated = accFrameBits;

        /* If the current row was large enough to cause a large QP jump, try re-encoding it. */
        if (qpVbv > qpMax && prevRowQp < qpMax && canReencodeRow)
        {
            /* Bump QP to halfway in between... close enough. */
            qpVbv = Clip3(prevRowQp + 1.0f, qpMax, (prevRowQp + qpVbv) * 0.5);
            return -1;
        }

        if (m_param->rc.rfConstantMin)
        {
            if (qpVbv < qpMin && prevRowQp > qpMin && canReencodeRow)
            {
                qpVbv = Clip3(qpMin, prevRowQp, (prevRowQp + qpVbv) * 0.5);
                return -1;
            }
        }
    }
    else
    {
        int32_t encodedBitsSoFar = 0;
        rce->frameSizeEstimated = predictRowsSizeSum(curFrame, rce, qpVbv, encodedBitsSoFar);

        /* Last-ditch attempt: if the last row of the frame underflowed the VBV,
         * try again. */
        if ((rce->frameSizeEstimated > (rce->bufferFill - m_bufferRate * maxFrameError) &&
             qpVbv < qpMax && canReencodeRow))
        {
            qpVbv = qpMax;
            return -1;
        }
    }
    return 0;
}

/* modify the bitrate curve from pass1 for one frame */
double RateControl::getQScale(RateControlEntry *rce, double rateFactor)
{
    double q;

    if (m_param->rc.cuTree)
    {
        // Scale and units are obtained from rateNum and rateDenom for videos with fixed frame rates.
        double timescale = (double)m_param->fpsDenom / (2 * m_param->fpsNum);
        q = pow(BASE_FRAME_DURATION / CLIP_DURATION(2 * timescale), 1 - m_param->rc.qCompress);
    }
    else
        q = pow(rce->blurredComplexity, 1 - m_param->rc.qCompress);
    // avoid NaN's in the Rceq
    if (rce->coeffBits + rce->mvBits == 0)
        q = m_lastQScaleFor[rce->sliceType];
    else
    {
        m_lastRceq = q;
        q /= rateFactor;
    }
    return q;
}

void RateControl::updatePredictor(Predictor *p, double q, double var, double bits)
{
    if (var < 10)
        return;
    const double range = 1.5;
    double old_coeff = p->coeff / p->count;
    double new_coeff = bits * q / var;
    double new_coeff_clipped = Clip3(old_coeff / range, old_coeff * range, new_coeff);
    double new_offset = bits * q - new_coeff_clipped * var;
    if (new_offset >= 0)
        new_coeff = new_coeff_clipped;
    else
        new_offset = 0;
    p->count  *= p->decay;
    p->coeff  *= p->decay;
    p->offset *= p->decay;
    p->count++;
    p->coeff  += new_coeff;
    p->offset += new_offset;
}

void RateControl::updateVbv(int64_t bits, RateControlEntry* rce)
{
    if (rce->lastSatd >= m_ncu)
        updatePredictor(&m_pred[rce->sliceType], x265_qp2qScale(rce->qpaRc), (double)rce->lastSatd, (double)bits);
    if (!m_isVbv)
        return;

    m_bufferFillFinal -= bits;

    if (m_bufferFillFinal < 0)
        x265_log(m_param, X265_LOG_WARNING, "poc:%d, VBV underflow (%.0f bits)\n", rce->poc, m_bufferFillFinal);

    m_bufferFillFinal = X265_MAX(m_bufferFillFinal, 0);
    m_bufferFillFinal += m_bufferRate;
    m_bufferFillFinal = X265_MIN(m_bufferFillFinal, m_bufferSize);
}

/* After encoding one frame, update rate control state */
int RateControl::rateControlEnd(Frame* curFrame, int64_t bits, RateControlEntry* rce, FrameStats* stats)
{
    int orderValue = m_startEndOrder.get();
    int endOrdinal = (rce->encodeOrder + m_param->frameNumThreads) * 2 - 1;
    while (orderValue < endOrdinal && !m_bTerminated)
    {
        /* no more frames are being encoded, so fake the start event if we would
         * have blocked on it. Note that this does not enforce rateControlEnd()
         * ordering during flush, but this has no impact on the outputs */
        if (m_finalFrameCount && orderValue >= 2 * m_finalFrameCount)
            break;
        orderValue = m_startEndOrder.waitForChange(orderValue);
    }

    FrameData& curEncData = *curFrame->m_encData;
    int64_t actualBits = bits;
    Slice *slice = curEncData.m_slice;
    if (m_isAbr)
    {
        if (m_param->rc.rateControlMode == X265_RC_ABR && !m_param->rc.bStatRead)
            checkAndResetABR(rce, true);

        if (m_param->rc.rateControlMode == X265_RC_CRF)
        {
            if (int(curEncData.m_avgQpRc + 0.5) == slice->m_sliceQp)
                curEncData.m_rateFactor = m_rateFactorConstant;
            else
            {
                /* If vbv changed the frame QP recalculate the rate-factor */
                double baseCplx = m_ncu * (m_param->bframes ? 120 : 80);
                double mbtree_offset = m_param->rc.cuTree ? (1.0 - m_param->rc.qCompress) * 13.5 : 0;
                curEncData.m_rateFactor = pow(baseCplx, 1 - m_qCompress) /
                    x265_qp2qScale(int(curEncData.m_avgQpRc + 0.5) + mbtree_offset);
            }
        }
    }

    if (m_param->rc.aqMode || m_isVbv)
    {
        if (m_isVbv)
        {
            for (uint32_t i = 0; i < slice->m_sps->numCuInHeight; i++)
                curEncData.m_avgQpRc += curEncData.m_rowStat[i].sumQpRc;

            curEncData.m_avgQpRc /= slice->m_sps->numCUsInFrame;
            rce->qpaRc = curEncData.m_avgQpRc;

            // copy avg RC qp to m_avgQpAq. To print out the correct qp when aq/cutree is disabled.
            curEncData.m_avgQpAq = curEncData.m_avgQpRc;
        }

        if (m_param->rc.aqMode)
        {
            for (uint32_t i = 0; i < slice->m_sps->numCuInHeight; i++)
                curEncData.m_avgQpAq += curEncData.m_rowStat[i].sumQpAq;

            curEncData.m_avgQpAq /= slice->m_sps->numCUsInFrame;
        }
    }

    // Write frame stats into the stats file if 2 pass is enabled.
    if (m_param->rc.bStatWrite)
    {
        char cType = rce->sliceType == I_SLICE ? (rce->poc > 0 && m_param->bOpenGOP ? 'i' : 'I')
            : rce->sliceType == P_SLICE ? 'P'
            : IS_REFERENCED(curFrame) ? 'B' : 'b';
        if (fprintf(m_statFileOut,
                    "in:%d out:%d type:%c q:%.2f q-aq:%.2f tex:%d mv:%d misc:%d icu:%.2f pcu:%.2f scu:%.2f ;\n",
                    rce->poc, rce->encodeOrder,
                    cType, curEncData.m_avgQpRc, curEncData.m_avgQpAq,
                    stats->coeffBits,
                    stats->mvBits,
                    stats->miscBits,
                    stats->percentIntra * m_ncu,
                    stats->percentInter * m_ncu,
                    stats->percentSkip  * m_ncu) < 0)
            goto writeFailure;
        /* Don't re-write the data in multi-pass mode. */
        if (m_param->rc.cuTree && IS_REFERENCED(curFrame) && !m_param->rc.bStatRead)
        {
            uint8_t sliceType = (uint8_t)rce->sliceType;
            for (int i = 0; i < m_ncu; i++)
                    m_cuTreeStats.qpBuffer[0][i] = (uint16_t)(curFrame->m_lowres.qpCuTreeOffset[i] * 256.0);
            if (fwrite(&sliceType, 1, 1, m_cutreeStatFileOut) < 1)
                goto writeFailure;
            if (fwrite(m_cuTreeStats.qpBuffer[0], sizeof(uint16_t), m_ncu, m_cutreeStatFileOut) < (size_t)m_ncu)
                goto writeFailure;
        }
    }
    if (m_isAbr && !m_isAbrReset)
    {
        /* amortize part of each I slice over the next several frames, up to
         * keyint-max, to avoid over-compensating for the large I slice cost */
        if (!m_param->rc.bStatWrite && !m_param->rc.bStatRead)
        {
            if (rce->sliceType == I_SLICE)
            {
                /* previous I still had a residual; roll it into the new loan */
                if (m_residualFrames)
                    bits += m_residualCost * m_residualFrames;
                m_residualFrames = X265_MIN(s_amortizeFrames, m_param->keyframeMax);
                m_residualCost = (int)((bits * s_amortizeFraction) / m_residualFrames);
                bits -= m_residualCost * m_residualFrames;
            }
            else if (m_residualFrames)
            {
                bits += m_residualCost;
                m_residualFrames--;
            }
        }
        if (rce->sliceType != B_SLICE)
        {
            /* The factor 1.5 is to tune up the actual bits, otherwise the cplxrSum is scaled too low
                * to improve short term compensation for next frame. */
            m_cplxrSum += (bits * x265_qp2qScale(rce->qpaRc) / rce->qRceq) - (rce->rowCplxrSum);
        }
        else
        {
            /* Depends on the fact that B-frame's QP is an offset from the following P-frame's.
                * Not perfectly accurate with B-refs, but good enough. */
            m_cplxrSum += (bits * x265_qp2qScale(rce->qpaRc) / (rce->qRceq * fabs(m_param->rc.pbFactor))) - (rce->rowCplxrSum);
        }
        m_wantedBitsWindow += m_frameDuration * m_bitrate;
        m_totalBits += bits - rce->rowTotalBits;
        int pos = m_sliderPos - m_param->frameNumThreads;
        if (pos >= 0)
            m_encodedBitsWindow[pos % s_slidingWindowFrames] = actualBits;
    }

    if (m_2pass)
    {
        m_expectedBitsSum += qScale2bits(rce, x265_qp2qScale(rce->newQp));
        m_totalBits += bits - rce->rowTotalBits;
    }

    if (m_isVbv)
    {
        if (rce->sliceType == B_SLICE)
        {
            m_bframeBits += actualBits;
            if (rce->bLastMiniGopBFrame)
            {
                if (rce->bframes != 0)
                    updatePredictor(&m_predBfromP, x265_qp2qScale(rce->qpaRc), (double)rce->leadingNoBSatd, (double)m_bframeBits / rce->bframes);
                m_bframeBits = 0;
            }
        }
        updateVbv(actualBits, rce);

        if (m_param->bEmitHRDSEI)
        {
            const VUI *vui = &curEncData.m_slice->m_sps->vuiParameters;
            const HRDInfo *hrd = &vui->hrdParameters;
            const TimingInfo *time = &vui->timingInfo;
            if (!curFrame->m_poc)
            {
                // first access unit initializes the HRD
                rce->hrdTiming->cpbInitialAT = 0;
                rce->hrdTiming->cpbRemovalTime = m_nominalRemovalTime = (double)m_bufPeriodSEI.m_initialCpbRemovalDelay / 90000;
            }
            else
            {
                rce->hrdTiming->cpbRemovalTime = m_nominalRemovalTime + (double)rce->picTimingSEI->m_auCpbRemovalDelay * time->numUnitsInTick / time->timeScale;
                double cpbEarliestAT = rce->hrdTiming->cpbRemovalTime - (double)m_bufPeriodSEI.m_initialCpbRemovalDelay / 90000;
                if (!curFrame->m_lowres.bKeyframe)
                    cpbEarliestAT -= (double)m_bufPeriodSEI.m_initialCpbRemovalDelayOffset / 90000;

                rce->hrdTiming->cpbInitialAT = hrd->cbrFlag ? m_prevCpbFinalAT : X265_MAX(m_prevCpbFinalAT, cpbEarliestAT);
            }

            uint32_t cpbsizeUnscale = hrd->cpbSizeValue << (hrd->cpbSizeScale + CPB_SHIFT);
            rce->hrdTiming->cpbFinalAT = m_prevCpbFinalAT = rce->hrdTiming->cpbInitialAT + actualBits / cpbsizeUnscale;
            rce->hrdTiming->dpbOutputTime = (double)rce->picTimingSEI->m_picDpbOutputDelay * time->numUnitsInTick / time->timeScale + rce->hrdTiming->cpbRemovalTime;
        }
    }
    // Allow rateControlStart of next frame only when rateControlEnd of previous frame is over
    m_startEndOrder.incr();
    rce->isActive = false;
    return 0;

writeFailure:
    x265_log(m_param, X265_LOG_ERROR, "RatecontrolEnd: stats file write failure\n");
    return 1;
}

#if defined(_MSC_VER)
#pragma warning(disable: 4996) // POSIX function names are just fine, thank you
#endif

/* called when the encoder is flushing, and thus the final frame count is
 * unambiguously known */
void RateControl::setFinalFrameCount(int count)
{
    m_finalFrameCount = count;
    /* unblock waiting threads */
    m_startEndOrder.set(m_startEndOrder.get());
}

/* called when the encoder is closing, and no more frames will be output.
 * all blocked functions must finish so the frame encoder threads can be
 * closed */
void RateControl::terminate()
{
    m_bTerminated = true;
    /* unblock waiting threads */
    m_startEndOrder.set(m_startEndOrder.get());
}

void RateControl::destroy()
{
    const char *fileName = m_param->rc.statFileName;
    if (!fileName)
        fileName = s_defaultStatFileName;

    if (m_statFileOut)
    {
        fclose(m_statFileOut);
        char *tmpFileName = strcatFilename(fileName, ".temp");
        int bError = 1;
        if (tmpFileName)
        {
           unlink(fileName);
           bError = rename(tmpFileName, fileName);
        }
        if (bError)
        {
            x265_log(m_param, X265_LOG_ERROR, "failed to rename output stats file to \"%s\"\n",
                     fileName);
        }
        X265_FREE(tmpFileName);
    }

    if (m_cutreeStatFileOut)
    {
        fclose(m_cutreeStatFileOut);
        char *tmpFileName = strcatFilename(fileName, ".cutree.temp");
        char *newFileName = strcatFilename(fileName, ".cutree");
        int bError = 1;
        if (tmpFileName && newFileName)
        {
           unlink(newFileName);
           bError = rename(tmpFileName, newFileName);
        }
        if (bError)
        {
            x265_log(m_param, X265_LOG_ERROR, "failed to rename cutree output stats file to \"%s\"\n",
                     newFileName);
        }
        X265_FREE(tmpFileName);
        X265_FREE(newFileName);
    }

    if (m_cutreeStatFileIn)
        fclose(m_cutreeStatFileIn);

    X265_FREE(m_rce2Pass);
    for (int i = 0; i < 2; i++)
        X265_FREE(m_cuTreeStats.qpBuffer[i]);
}

