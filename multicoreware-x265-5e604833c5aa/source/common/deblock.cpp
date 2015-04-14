/*****************************************************************************
* Copyright (C) 2013 x265 project
*
* Author: Gopu Govindaswamy <gopu@multicorewareinc.com>
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
#include "deblock.h"
#include "framedata.h"
#include "picyuv.h"
#include "slice.h"
#include "mv.h"

using namespace x265;

#define DEBLOCK_SMALLEST_BLOCK  8
#define DEFAULT_INTRA_TC_OFFSET 2

void Deblock::deblockCTU(CUData* cu, int32_t dir)
{
    uint8_t blockingStrength[MAX_NUM_PARTITIONS];

    memset(blockingStrength, 0, sizeof(uint8_t) * m_numPartitions);

    deblockCU(cu, 0, 0, dir, blockingStrength);
}

/* Deblocking filter process in CU-based (the same function as conventional's)
 * param Edge the direction of the edge in block boundary (horizonta/vertical), which is added newly */
void Deblock::deblockCU(CUData* cu, uint32_t absPartIdx, uint32_t depth, const int32_t dir, uint8_t blockingStrength[])
{
    if (cu->m_partSize[absPartIdx] == SIZE_NONE)
        return;

    uint32_t curNumParts = NUM_CU_PARTITIONS >> (depth << 1);

    const SPS& sps = *cu->m_slice->m_sps;

    if (cu->m_cuDepth[absPartIdx] > depth)
    {
        uint32_t qNumParts   = curNumParts >> 2;
        uint32_t xmax = sps.picWidthInLumaSamples  - cu->m_cuPelX;
        uint32_t ymax = sps.picHeightInLumaSamples - cu->m_cuPelY;
        for (uint32_t partIdx = 0; partIdx < 4; partIdx++, absPartIdx += qNumParts)
            if (g_zscanToPelX[absPartIdx] < xmax && g_zscanToPelY[absPartIdx] < ymax)
                deblockCU(cu, absPartIdx, depth + 1, dir, blockingStrength);
        return;
    }

    const uint32_t widthInBaseUnits = sps.numPartInCUSize >> depth;
    Param params;
    setLoopfilterParam(cu, absPartIdx, &params);
    setEdgefilterPU(cu, absPartIdx, dir, blockingStrength, widthInBaseUnits);
    setEdgefilterTU(cu, absPartIdx, depth, dir, blockingStrength);
    setEdgefilterMultiple(cu, absPartIdx, dir, 0, (dir == EDGE_VER ? params.leftEdge : params.topEdge), blockingStrength, widthInBaseUnits);

    for (uint32_t partIdx = absPartIdx; partIdx < absPartIdx + curNumParts; partIdx++)
    {
        uint32_t bsCheck = !(partIdx & (1 << dir));

        if (bsCheck && blockingStrength[partIdx])
            getBoundaryStrengthSingle(cu, dir, partIdx, blockingStrength);
    }

    const uint32_t partIdxIncr = DEBLOCK_SMALLEST_BLOCK >> LOG2_UNIT_SIZE;
    uint32_t sizeInPU = sps.numPartInCUSize >> depth;
    uint32_t shiftFactor = (dir == EDGE_VER) ? cu->m_hChromaShift : cu->m_vChromaShift;
    uint32_t chromaMask = ((DEBLOCK_SMALLEST_BLOCK << shiftFactor) >> LOG2_UNIT_SIZE) - 1;
    uint32_t e0 = (dir == EDGE_VER ? g_zscanToPelX[absPartIdx] : g_zscanToPelY[absPartIdx]) >> LOG2_UNIT_SIZE;
        
    for (uint32_t e = 0; e < sizeInPU; e += partIdxIncr)
    {
        edgeFilterLuma(cu, absPartIdx, depth, dir, e, blockingStrength);
        if (!((e0 + e) & chromaMask))
            edgeFilterChroma(cu, absPartIdx, depth, dir, e, blockingStrength);
    }
}

static inline uint32_t calcBsIdx(CUData* cu, uint32_t absPartIdx, int32_t dir, int32_t edgeIdx, int32_t baseUnitIdx)
{
    uint32_t ctuWidthInBaseUnits = cu->m_slice->m_sps->numPartInCUSize;

    if (dir)
        return g_rasterToZscan[g_zscanToRaster[absPartIdx] + edgeIdx * ctuWidthInBaseUnits + baseUnitIdx];
    else
        return g_rasterToZscan[g_zscanToRaster[absPartIdx] + baseUnitIdx * ctuWidthInBaseUnits + edgeIdx];
}

void Deblock::setEdgefilterMultiple(CUData* cu, uint32_t scanIdx, int32_t dir, int32_t edgeIdx, uint8_t value, uint8_t blockingStrength[], uint32_t widthInBaseUnits)
{
    const uint32_t numElem = widthInBaseUnits;
    X265_CHECK(numElem > 0, "numElem edge filter check\n");
    for (uint32_t i = 0; i < numElem; i++)
    {
        const uint32_t bsidx = calcBsIdx(cu, scanIdx, dir, edgeIdx, i);
        blockingStrength[bsidx] = value;
    }
}

void Deblock::setEdgefilterTU(CUData* cu, uint32_t absPartIdx, uint32_t depth, int32_t dir, uint8_t blockingStrength[])
{
    if ((uint32_t)cu->m_tuDepth[absPartIdx] + cu->m_cuDepth[absPartIdx] > depth)
    {
        const uint32_t curNumParts = NUM_CU_PARTITIONS >> (depth << 1);
        const uint32_t qNumParts   = curNumParts >> 2;

        for (uint32_t partIdx = 0; partIdx < 4; partIdx++, absPartIdx += qNumParts)
            setEdgefilterTU(cu, absPartIdx, depth + 1, dir, blockingStrength);
        return;
    }

    uint32_t widthInBaseUnits  = 1 << (cu->m_log2CUSize[absPartIdx] - cu->m_tuDepth[absPartIdx] - LOG2_UNIT_SIZE);
    setEdgefilterMultiple(cu, absPartIdx, dir, 0, 2, blockingStrength, widthInBaseUnits);
}

void Deblock::setEdgefilterPU(CUData* cu, uint32_t absPartIdx, int32_t dir, uint8_t blockingStrength[], uint32_t widthInBaseUnits)
{
    const uint32_t hWidthInBaseUnits = widthInBaseUnits >> 1;
    const uint32_t qWidthInBaseUnits = widthInBaseUnits >> 2;

    switch (cu->m_partSize[absPartIdx])
    {
    case SIZE_2NxN:
        if (EDGE_HOR == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, hWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_Nx2N:
        if (EDGE_VER == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, hWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_NxN:
        setEdgefilterMultiple(cu, absPartIdx, dir, hWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_2NxnU:
        if (EDGE_HOR == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, qWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_nLx2N:
        if (EDGE_VER == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, qWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_2NxnD:
        if (EDGE_HOR == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, widthInBaseUnits - qWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;
    case SIZE_nRx2N:
        if (EDGE_VER == dir)
            setEdgefilterMultiple(cu, absPartIdx, dir, widthInBaseUnits - qWidthInBaseUnits, 1, blockingStrength, widthInBaseUnits);
        break;

    case SIZE_2Nx2N:
    default:
        break;
    }
}

void Deblock::setLoopfilterParam(CUData* cu, uint32_t absPartIdx, Param *params)
{
    uint32_t x = cu->m_cuPelX + g_zscanToPelX[absPartIdx];
    uint32_t y = cu->m_cuPelY + g_zscanToPelY[absPartIdx];

    const CUData* tempCU;
    uint32_t    tempPartIdx;

    if (!x)
        params->leftEdge = 0;
    else
    {
        tempCU = cu->getPULeft(tempPartIdx, absPartIdx);
        if (tempCU)
            params->leftEdge = 2;
        else
            params->leftEdge = 0;
    }

    if (!y)
        params->topEdge = 0;
    else
    {
        tempCU = cu->getPUAbove(tempPartIdx, absPartIdx);
        if (tempCU)
            params->topEdge = 2;
        else
            params->topEdge = 0;
    }
}

void Deblock::getBoundaryStrengthSingle(CUData* cu, int32_t dir, uint32_t absPartIdx, uint8_t blockingStrength[])
{
    const Slice* const slice = cu->m_slice;
    const uint32_t partQ = absPartIdx;
    CUData* const cuQ = cu;

    uint32_t partP;
    const CUData* cuP;
    uint8_t bs = 0;

    // Calculate block index
    if (dir == EDGE_VER)
        cuP = cuQ->getPULeft(partP, partQ);
    else // (dir == EDGE_HOR)
        cuP = cuQ->getPUAbove(partP, partQ);

    // Set BS for Intra MB : BS = 4 or 3
    if (cuP->isIntra(partP) || cuQ->isIntra(partQ))
        bs = 2;

    // Set BS for not Intra MB : BS = 2 or 1 or 0
    if (!cuP->isIntra(partP) && !cuQ->isIntra(partQ))
    {
        uint32_t nsPartQ = partQ;
        uint32_t nsPartP = partP;

        if (blockingStrength[absPartIdx] > 1 &&
            (cuQ->getCbf(nsPartQ, TEXT_LUMA, cuQ->m_tuDepth[nsPartQ]) ||
             cuP->getCbf(nsPartP, TEXT_LUMA, cuP->m_tuDepth[nsPartP])))
            bs = 1;
        else
        {
            if (dir == EDGE_HOR)
                cuP = cuQ->getPUAbove(partP, partQ);

            if (slice->isInterB() || cuP->m_slice->isInterB())
            {
                int32_t refIdx;
                Frame *refP0, *refP1, *refQ0, *refQ1;
                refIdx = cuP->m_refIdx[0][partP];
                refP0 = (refIdx < 0) ? NULL : cuP->m_slice->m_refPicList[0][refIdx];
                refIdx = cuP->m_refIdx[1][partP];
                refP1 = (refIdx < 0) ? NULL : cuP->m_slice->m_refPicList[1][refIdx];
                refIdx = cuQ->m_refIdx[0][partQ];
                refQ0 = (refIdx < 0) ? NULL : slice->m_refPicList[0][refIdx];
                refIdx = cuQ->m_refIdx[1][partQ];
                refQ1 = (refIdx < 0) ? NULL : slice->m_refPicList[1][refIdx];

                MV mvp0 = cuP->m_mv[0][partP];
                MV mvp1 = cuP->m_mv[1][partP];
                MV mvq0 = cuQ->m_mv[0][partQ];
                MV mvq1 = cuQ->m_mv[1][partQ];

                if (!refP0) mvp0 = 0;
                if (!refP1) mvp1 = 0;
                if (!refQ0) mvq0 = 0;
                if (!refQ1) mvq1 = 0;

                if (((refP0 == refQ0) && (refP1 == refQ1)) || ((refP0 == refQ1) && (refP1 == refQ0)))
                {
                    if (refP0 != refP1) // Different L0 & L1
                    {
                        if (refP0 == refQ0)
                        {
                            bs  = ((abs(mvq0.x - mvp0.x) >= 4) ||
                                   (abs(mvq0.y - mvp0.y) >= 4) ||
                                   (abs(mvq1.x - mvp1.x) >= 4) ||
                                   (abs(mvq1.y - mvp1.y) >= 4)) ? 1 : 0;
                        }
                        else
                        {
                            bs  = ((abs(mvq1.x - mvp0.x) >= 4) ||
                                   (abs(mvq1.y - mvp0.y) >= 4) ||
                                   (abs(mvq0.x - mvp1.x) >= 4) ||
                                   (abs(mvq0.y - mvp1.y) >= 4)) ? 1 : 0;
                        }
                    }
                    else // Same L0 & L1
                    {
                        bs  = ((abs(mvq0.x - mvp0.x) >= 4) ||
                               (abs(mvq0.y - mvp0.y) >= 4) ||
                               (abs(mvq1.x - mvp1.x) >= 4) ||
                               (abs(mvq1.y - mvp1.y) >= 4)) &&
                              ((abs(mvq1.x - mvp0.x) >= 4) ||
                               (abs(mvq1.y - mvp0.y) >= 4) ||
                               (abs(mvq0.x - mvp1.x) >= 4) ||
                               (abs(mvq0.y - mvp1.y) >= 4)) ? 1 : 0;
                    }
                }
                else // for all different Ref_Idx
                    bs = 1;
            }
            else // slice->isInterP()
            {
                int32_t refIdx;
                Frame *refp0, *refq0;
                refIdx = cuP->m_refIdx[0][partP];
                refp0 = (refIdx < 0) ? NULL : cuP->m_slice->m_refPicList[0][refIdx];
                refIdx = cuQ->m_refIdx[0][partQ];
                refq0 = (refIdx < 0) ? NULL : slice->m_refPicList[0][refIdx];
                MV mvp0 = cuP->m_mv[0][partP];
                MV mvq0 = cuQ->m_mv[0][partQ];

                if (!refp0) mvp0 = 0;
                if (!refq0) mvq0 = 0;

                bs = ((refp0 != refq0) ||
                      (abs(mvq0.x - mvp0.x) >= 4) ||
                      (abs(mvq0.y - mvp0.y) >= 4)) ? 1 : 0;
            }
        }
    }

    blockingStrength[absPartIdx] = bs;
}

static inline int32_t calcDP(pixel* src, intptr_t offset)
{
    return abs(static_cast<int32_t>(src[-offset * 3]) - 2 * src[-offset * 2] + src[-offset]);
}

static inline int32_t calcDQ(pixel* src, intptr_t offset)
{
    return abs(static_cast<int32_t>(src[0]) - 2 * src[offset] + src[offset * 2]);
}

static inline bool useStrongFiltering(intptr_t offset, int32_t beta, int32_t tc, pixel* src)
{
    int16_t m4     = (int16_t)src[0];
    int16_t m3     = (int16_t)src[-offset];
    int16_t m7     = (int16_t)src[offset * 3];
    int16_t m0     = (int16_t)src[-offset * 4];
    int32_t strong = abs(m0 - m3) + abs(m7 - m4);

    return (strong < (beta >> 3)) && (abs(m3 - m4) < ((tc * 5 + 1) >> 1));
}

/* Deblocking for the luminance component with strong or weak filter
 * \param src            pointer to picture data
 * \param offset         offset value for picture data
 * \param tc             tc value
 * \param partPNoFilter  indicator to disable filtering on partP
 * \param partQNoFilter  indicator to disable filtering on partQ
 * \param filterSecondP  decision weak filter/no filter for partP
 * \param filterSecondQ  decision weak filter/no filter for partQ */
static inline void pelFilterLumaStrong(pixel* src, intptr_t srcStep, intptr_t offset, int32_t tc, bool partPNoFilter, bool partQNoFilter)
{
    for (int32_t i = 0; i < UNIT_SIZE; i++, src += srcStep)
    {
        int16_t m4  = (int16_t)src[0];
        int16_t m3  = (int16_t)src[-offset];
        int16_t m5  = (int16_t)src[offset];
        int16_t m2  = (int16_t)src[-offset * 2];
        int32_t tc2 = 2 * tc;
        if (!partPNoFilter)
        {
            int16_t m1  = (int16_t)src[-offset * 3];
            int16_t m0  = (int16_t)src[-offset * 4];
            src[-offset * 3] = (pixel)(Clip3(-tc2, tc2, ((2 * m0 + 3 * m1 + m2 + m3 + m4 + 4) >> 3) - m1) + m1);
            src[-offset * 2] = (pixel)(Clip3(-tc2, tc2, ((m1 + m2 + m3 + m4 + 2) >> 2) - m2) + m2);
            src[-offset]     = (pixel)(Clip3(-tc2, tc2, ((m1 + 2 * m2 + 2 * m3 + 2 * m4 + m5 + 4) >> 3) - m3) + m3);
        }
        if (!partQNoFilter)
        {
            int16_t m6  = (int16_t)src[offset * 2];
            int16_t m7  = (int16_t)src[offset * 3];
            src[0]           = (pixel)(Clip3(-tc2, tc2, ((m2 + 2 * m3 + 2 * m4 + 2 * m5 + m6 + 4) >> 3) - m4) + m4);
            src[offset]      = (pixel)(Clip3(-tc2, tc2, ((m3 + m4 + m5 + m6 + 2) >> 2) - m5) + m5);
            src[offset * 2]  = (pixel)(Clip3(-tc2, tc2, ((m3 + m4 + m5 + 3 * m6 + 2 * m7 + 4) >> 3) - m6) + m6);
        }
    }
}

/* Weak filter */
static inline void pelFilterLuma(pixel* src, intptr_t srcStep, intptr_t offset, int32_t tc, bool partPNoFilter, bool partQNoFilter,
                                 bool filterSecondP, bool filterSecondQ)
{
    int32_t thrCut = tc * 10;

    for (int32_t i = 0; i < UNIT_SIZE; i++, src += srcStep)
    {
        int16_t m4  = (int16_t)src[0];
        int16_t m3  = (int16_t)src[-offset];
        int16_t m5  = (int16_t)src[offset];
        int16_t m2  = (int16_t)src[-offset * 2];

        int32_t delta = (9 * (m4 - m3) - 3 * (m5 - m2) + 8) >> 4;

        if (abs(delta) < thrCut)
        {
            delta = Clip3(-tc, tc, delta);

            int32_t tc2 = tc >> 1;
            if (!partPNoFilter)
            {
                src[-offset] = Clip(m3 + delta);
                if (filterSecondP)
                {
                    int16_t m1  = (int16_t)src[-offset * 3];
                    int32_t delta1 = Clip3(-tc2, tc2, ((((m1 + m3 + 1) >> 1) - m2 + delta) >> 1));
                    src[-offset * 2] = Clip(m2 + delta1);
                }
            }
            if (!partQNoFilter)
            {
                src[0] = Clip(m4 - delta);
                if (filterSecondQ)
                {
                    int16_t m6  = (int16_t)src[offset * 2];
                    int32_t delta2 = Clip3(-tc2, tc2, ((((m6 + m4 + 1) >> 1) - m5 - delta) >> 1));
                    src[offset] = Clip(m5 + delta2);
                }
            }
        }
    }
}

/* Deblocking of one line/column for the chrominance component
 * \param src            pointer to picture data
 * \param offset         offset value for picture data
 * \param tc             tc value
 * \param partPNoFilter  indicator to disable filtering on partP
 * \param partQNoFilter  indicator to disable filtering on partQ */
static inline void pelFilterChroma(pixel* src, intptr_t srcStep, intptr_t offset, int32_t tc, bool partPNoFilter, bool partQNoFilter)
{
    for (int32_t i = 0; i < UNIT_SIZE; i++, src += srcStep)
    {
        int16_t m4  = (int16_t)src[0];
        int16_t m3  = (int16_t)src[-offset];
        int16_t m5  = (int16_t)src[offset];
        int16_t m2  = (int16_t)src[-offset * 2];

        int32_t delta = Clip3(-tc, tc, ((((m4 - m3) << 2) + m2 - m5 + 4) >> 3));
        if (!partPNoFilter)
            src[-offset] = Clip(m3 + delta);
        if (!partQNoFilter)
            src[0] = Clip(m4 - delta);
    }
}

void Deblock::edgeFilterLuma(CUData* cu, uint32_t absPartIdx, uint32_t depth, int32_t dir, int32_t edge, const uint8_t blockingStrength[])
{
    PicYuv* reconYuv = cu->m_encData->m_reconPicYuv;
    pixel* src = reconYuv->getLumaAddr(cu->m_cuAddr, absPartIdx);

    intptr_t stride = reconYuv->m_stride;
    uint32_t numParts = cu->m_slice->m_sps->numPartInCUSize >> depth;

    intptr_t offset, srcStep;

    bool  partPNoFilter = false;
    bool  partQNoFilter = false;
    uint32_t  partP = 0;
    uint32_t  partQ = 0;
    const CUData* cuP = cu;
    const CUData* cuQ = cu;
    int32_t betaOffset = cuQ->m_slice->m_pps->deblockingFilterBetaOffsetDiv2 << 1;
    int32_t tcOffset = cuQ->m_slice->m_pps->deblockingFilterTcOffsetDiv2 << 1;

    if (dir == EDGE_VER)
    {
        offset = 1;
        srcStep = stride;
        src += (edge << LOG2_UNIT_SIZE);
    }
    else // (dir == EDGE_HOR)
    {
        offset = stride;
        srcStep = 1;
        src += (edge << LOG2_UNIT_SIZE) * stride;
    }

    for (uint32_t idx = 0; idx < numParts; idx++)
    {
        uint32_t unitOffset = idx << LOG2_UNIT_SIZE;
        uint32_t bsAbsIdx = calcBsIdx(cu, absPartIdx, dir, edge, idx);
        uint32_t bs = blockingStrength[bsAbsIdx];
        if (bs)
        {
            int32_t qpQ = cu->m_qp[bsAbsIdx];
            partQ = bsAbsIdx;

            // Derive neighboring PU index
            if (dir == EDGE_VER)
                cuP = cuQ->getPULeft(partP, partQ);
            else // (dir == EDGE_HOR)
                cuP = cuQ->getPUAbove(partP, partQ);

            int32_t qpP = cuP->m_qp[partP];
            int32_t qp = (qpP + qpQ + 1) >> 1;

            int32_t indexB = Clip3(0, QP_MAX_SPEC, qp + betaOffset);

            const int32_t bitdepthShift = X265_DEPTH - 8;
            int32_t beta = s_betaTable[indexB] << bitdepthShift;

            int32_t dp0 = calcDP(src + srcStep * (unitOffset + 0), offset);
            int32_t dq0 = calcDQ(src + srcStep * (unitOffset + 0), offset);
            int32_t dp3 = calcDP(src + srcStep * (unitOffset + 3), offset);
            int32_t dq3 = calcDQ(src + srcStep * (unitOffset + 3), offset);
            int32_t d0 = dp0 + dq0;
            int32_t d3 = dp3 + dq3;

            int32_t d =  d0 + d3;

            if (d < beta)
            {
                if (cu->m_slice->m_pps->bTransquantBypassEnabled)
                {
                    // check if each of PUs is lossless coded
                    partPNoFilter = !!cuP->m_tqBypass[partP];
                    partQNoFilter = !!cuQ->m_tqBypass[partQ];
                }

                int32_t indexTC = Clip3(0, QP_MAX_SPEC + DEFAULT_INTRA_TC_OFFSET, int32_t(qp + DEFAULT_INTRA_TC_OFFSET * (bs - 1) + tcOffset));
                int32_t tc = s_tcTable[indexTC] << bitdepthShift;

                bool sw = (2 * d0 < (beta >> 2) &&
                           2 * d3 < (beta >> 2) &&
                           useStrongFiltering(offset, beta, tc, src + srcStep * (unitOffset + 0)) &&
                           useStrongFiltering(offset, beta, tc, src + srcStep * (unitOffset + 3)));

                if (sw)
                    pelFilterLumaStrong(src + srcStep * unitOffset, srcStep, offset, tc, partPNoFilter, partQNoFilter);
                else
                {
                    int32_t sideThreshold = (beta + (beta >> 1)) >> 3;
                    int32_t dp = dp0 + dp3;
                    int32_t dq = dq0 + dq3;
                    bool filterP = (dp < sideThreshold);
                    bool filterQ = (dq < sideThreshold);

                    pelFilterLuma(src + srcStep * unitOffset, srcStep, offset, tc, partPNoFilter, partQNoFilter, filterP, filterQ);
                }
            }
        }
    }
}

void Deblock::edgeFilterChroma(CUData* cu, uint32_t absPartIdx, uint32_t depth, int32_t dir, int32_t edge, const uint8_t blockingStrength[])
{
    int32_t chFmt = cu->m_chromaFormat, chromaShift;
    intptr_t offset, srcStep;

    bool partPNoFilter = false;
    bool partQNoFilter = false;
    uint32_t partP;
    uint32_t partQ;
    const CUData* cuP;
    const CUData* cuQ = cu;
    int32_t tcOffset = cu->m_slice->m_pps->deblockingFilterTcOffsetDiv2 << 1;

    X265_CHECK(((dir == EDGE_VER)
                ? ((g_zscanToPelX[absPartIdx] + edge * UNIT_SIZE) >> cu->m_hChromaShift)
                : ((g_zscanToPelY[absPartIdx] + edge * UNIT_SIZE) >> cu->m_vChromaShift)) % DEBLOCK_SMALLEST_BLOCK == 0,
               "invalid edge\n");

    PicYuv* reconPic = cu->m_encData->m_reconPicYuv;
    intptr_t stride = reconPic->m_strideC;
    intptr_t srcOffset = reconPic->getChromaAddrOffset(cu->m_cuAddr, absPartIdx);

    if (dir == EDGE_VER)
    {
        chromaShift = cu->m_vChromaShift;
        srcOffset += (edge << (LOG2_UNIT_SIZE - cu->m_hChromaShift));
        offset     = 1;
        srcStep    = stride;
    }
    else // (dir == EDGE_HOR)
    {
        chromaShift = cu->m_hChromaShift;
        srcOffset += edge * stride << (LOG2_UNIT_SIZE - cu->m_vChromaShift);
        offset     = stride;
        srcStep    = 1;
    }

    pixel* srcChroma[2];
    srcChroma[0] = reconPic->m_picOrg[1] + srcOffset;
    srcChroma[1] = reconPic->m_picOrg[2] + srcOffset;

    uint32_t numUnits = cu->m_slice->m_sps->numPartInCUSize >> (depth + chromaShift);

    for (uint32_t idx = 0; idx < numUnits; idx++)
    {
        uint32_t unitOffset = idx << LOG2_UNIT_SIZE;
        uint32_t bsAbsIdx = calcBsIdx(cu, absPartIdx, dir, edge, idx << chromaShift);
        uint32_t bs = blockingStrength[bsAbsIdx];

        if (bs > 1)
        {
            int32_t qpQ = cu->m_qp[bsAbsIdx];
            partQ = bsAbsIdx;

            // Derive neighboring PU index
            if (dir == EDGE_VER)
                cuP = cuQ->getPULeft(partP, partQ);
            else // (dir == EDGE_HOR)
                cuP = cuQ->getPUAbove(partP, partQ);

            int32_t qpP = cuP->m_qp[partP];

            if (cu->m_slice->m_pps->bTransquantBypassEnabled)
            {
                // check if each of PUs is lossless coded
                partPNoFilter = !!cuP->m_tqBypass[partP];
                partQNoFilter = !!cuQ->m_tqBypass[partQ];
            }

            for (uint32_t chromaIdx = 0; chromaIdx < 2; chromaIdx++)
            {
                int32_t chromaQPOffset  = !chromaIdx ? cu->m_slice->m_pps->chromaCbQpOffset : cu->m_slice->m_pps->chromaCrQpOffset;
                int32_t qp = ((qpP + qpQ + 1) >> 1) + chromaQPOffset;
                if (qp >= 30)
                {
                    if (chFmt == X265_CSP_I420)
                        qp = g_chromaScale[qp];
                    else
                        qp = X265_MIN(qp, 51);
                }

                int32_t indexTC = Clip3(0, QP_MAX_SPEC + DEFAULT_INTRA_TC_OFFSET, int32_t(qp + DEFAULT_INTRA_TC_OFFSET + tcOffset));
                const int32_t bitdepthShift = X265_DEPTH - 8;
                int32_t tc = s_tcTable[indexTC] << bitdepthShift;
                pixel* srcC = srcChroma[chromaIdx];

                pelFilterChroma(srcC + srcStep * unitOffset, srcStep, offset, tc, partPNoFilter, partQNoFilter);
            }
        }
    }
}

const uint8_t Deblock::s_tcTable[54] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
    2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 11, 13, 14, 16, 18, 20, 22, 24
};

const uint8_t Deblock::s_betaTable[52] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64
};

