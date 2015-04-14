/*****************************************************************************
 * Copyright (C) 2014 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
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
#include "yuv.h"
#include "shortyuv.h"
#include "picyuv.h"
#include "primitives.h"

using namespace x265;

Yuv::Yuv()
{
    m_buf[0] = NULL;
    m_buf[1] = NULL;
    m_buf[2] = NULL;
}

bool Yuv::create(uint32_t size, int csp)
{
    m_csp = csp;
    m_hChromaShift = CHROMA_H_SHIFT(csp);
    m_vChromaShift = CHROMA_V_SHIFT(csp);

    // set width and height
    m_size  = size;
    m_csize = size >> m_hChromaShift;
    m_part = partitionFromSizes(size, size);

    size_t sizeL = size * size;
    size_t sizeC = sizeL >> (m_vChromaShift + m_hChromaShift);

    X265_CHECK((sizeC & 15) == 0, "invalid size");

    // memory allocation (padded for SIMD reads)
    CHECKED_MALLOC(m_buf[0], pixel, sizeL + sizeC * 2 + 8);
    m_buf[1] = m_buf[0] + sizeL;
    m_buf[2] = m_buf[0] + sizeL + sizeC;
    return true;

fail:
    return false;
}

void Yuv::destroy()
{
    X265_FREE(m_buf[0]);
}

void Yuv::copyToPicYuv(PicYuv& dstPic, uint32_t cuAddr, uint32_t absPartIdx) const
{
    pixel* dstY = dstPic.getLumaAddr(cuAddr, absPartIdx);

    primitives.luma_copy_pp[m_part](dstY, dstPic.m_stride, m_buf[0], m_size);

    pixel* dstU = dstPic.getCbAddr(cuAddr, absPartIdx);
    pixel* dstV = dstPic.getCrAddr(cuAddr, absPartIdx);
    primitives.chroma[m_csp].copy_pp[m_part](dstU, dstPic.m_strideC, m_buf[1], m_csize);
    primitives.chroma[m_csp].copy_pp[m_part](dstV, dstPic.m_strideC, m_buf[2], m_csize);
}

void Yuv::copyFromPicYuv(const PicYuv& srcPic, uint32_t cuAddr, uint32_t absPartIdx)
{
    /* We cheat with const_cast internally because the get methods are not capable of
     * returning const buffers and the primitives are not const aware, but we know
     * this function does not modify srcPic */
    PicYuv& srcPicSafe = const_cast<PicYuv&>(srcPic);
    pixel* srcY = srcPicSafe.getLumaAddr(cuAddr, absPartIdx);

    primitives.luma_copy_pp[m_part](m_buf[0], m_size, srcY, srcPic.m_stride);

    pixel* srcU = srcPicSafe.getCbAddr(cuAddr, absPartIdx);
    pixel* srcV = srcPicSafe.getCrAddr(cuAddr, absPartIdx);
    primitives.chroma[m_csp].copy_pp[m_part](m_buf[1], m_csize, srcU, srcPicSafe.m_strideC);
    primitives.chroma[m_csp].copy_pp[m_part](m_buf[2], m_csize, srcV, srcPicSafe.m_strideC);
}

void Yuv::copyFromYuv(const Yuv& srcYuv)
{
    X265_CHECK(m_size <= srcYuv.m_size, "invalid size\n");

    primitives.luma_copy_pp[m_part](m_buf[0], m_size, srcYuv.m_buf[0], srcYuv.m_size);
    primitives.chroma[m_csp].copy_pp[m_part](m_buf[1], m_csize, srcYuv.m_buf[1], srcYuv.m_csize);
    primitives.chroma[m_csp].copy_pp[m_part](m_buf[2], m_csize, srcYuv.m_buf[2], srcYuv.m_csize);
}

void Yuv::copyToPartYuv(Yuv& dstYuv, uint32_t absPartIdx) const
{
    pixel* dstY = dstYuv.getLumaAddr(absPartIdx);
    primitives.luma_copy_pp[m_part](dstY, dstYuv.m_size, m_buf[0], m_size);

    pixel* dstU = dstYuv.getCbAddr(absPartIdx);
    pixel* dstV = dstYuv.getCrAddr(absPartIdx);
    primitives.chroma[m_csp].copy_pp[m_part](dstU, dstYuv.m_csize, m_buf[1], m_csize);
    primitives.chroma[m_csp].copy_pp[m_part](dstV, dstYuv.m_csize, m_buf[2], m_csize);
}

void Yuv::copyPartToYuv(Yuv& dstYuv, uint32_t absPartIdx) const
{
    pixel* srcY = m_buf[0] + getAddrOffset(absPartIdx, m_size);
    pixel* dstY = dstYuv.m_buf[0];

    primitives.luma_copy_pp[dstYuv.m_part](dstY, dstYuv.m_size, srcY, m_size);

    pixel* srcU = m_buf[1] + getChromaAddrOffset(absPartIdx);
    pixel* srcV = m_buf[2] + getChromaAddrOffset(absPartIdx);
    pixel* dstU = dstYuv.m_buf[1];
    pixel* dstV = dstYuv.m_buf[2];
    primitives.chroma[m_csp].copy_pp[dstYuv.m_part](dstU, dstYuv.m_csize, srcU, m_csize);
    primitives.chroma[m_csp].copy_pp[dstYuv.m_part](dstV, dstYuv.m_csize, srcV, m_csize);
}

void Yuv::addClip(const Yuv& srcYuv0, const ShortYuv& srcYuv1, uint32_t log2SizeL)
{
    primitives.luma_add_ps[log2SizeL - 2](m_buf[0], m_size, srcYuv0.m_buf[0], srcYuv1.m_buf[0], srcYuv0.m_size, srcYuv1.m_size);
    primitives.chroma[m_csp].add_ps[log2SizeL - 2](m_buf[1], m_csize, srcYuv0.m_buf[1], srcYuv1.m_buf[1], srcYuv0.m_csize, srcYuv1.m_csize);
    primitives.chroma[m_csp].add_ps[log2SizeL - 2](m_buf[2], m_csize, srcYuv0.m_buf[2], srcYuv1.m_buf[2], srcYuv0.m_csize, srcYuv1.m_csize);
}

void Yuv::addAvg(const ShortYuv& srcYuv0, const ShortYuv& srcYuv1, uint32_t absPartIdx, uint32_t width, uint32_t height, bool bLuma, bool bChroma)
{
    int part = partitionFromSizes(width, height);

    if (bLuma)
    {
        int16_t* srcY0 = const_cast<ShortYuv&>(srcYuv0).getLumaAddr(absPartIdx);
        int16_t* srcY1 = const_cast<ShortYuv&>(srcYuv1).getLumaAddr(absPartIdx);
        pixel* dstY = getLumaAddr(absPartIdx);

        primitives.luma_addAvg[part](srcY0, srcY1, dstY, srcYuv0.m_size, srcYuv1.m_size, m_size);
    }
    if (bChroma)
    {
        int16_t* srcU0 = const_cast<ShortYuv&>(srcYuv0).getCbAddr(absPartIdx);
        int16_t* srcV0 = const_cast<ShortYuv&>(srcYuv0).getCrAddr(absPartIdx);
        int16_t* srcU1 = const_cast<ShortYuv&>(srcYuv1).getCbAddr(absPartIdx);
        int16_t* srcV1 = const_cast<ShortYuv&>(srcYuv1).getCrAddr(absPartIdx);
        pixel* dstU = getCbAddr(absPartIdx);
        pixel* dstV = getCrAddr(absPartIdx);

        primitives.chroma[m_csp].addAvg[part](srcU0, srcU1, dstU, srcYuv0.m_csize, srcYuv1.m_csize, m_csize);
        primitives.chroma[m_csp].addAvg[part](srcV0, srcV1, dstV, srcYuv0.m_csize, srcYuv1.m_csize, m_csize);
    }
}

void Yuv::copyPartToPartLuma(Yuv& dstYuv, uint32_t absPartIdx, uint32_t log2Size) const
{
    const pixel* src = getLumaAddr(absPartIdx);
    pixel* dst = dstYuv.getLumaAddr(absPartIdx);
    primitives.square_copy_pp[log2Size - 2](dst, dstYuv.m_size, const_cast<pixel*>(src), m_size);
}

void Yuv::copyPartToPartChroma(Yuv& dstYuv, uint32_t absPartIdx, uint32_t log2SizeL) const
{
    int part = partitionFromLog2Size(log2SizeL);
    const pixel* srcU = getCbAddr(absPartIdx);
    const pixel* srcV = getCrAddr(absPartIdx);
    pixel* dstU = dstYuv.getCbAddr(absPartIdx);
    pixel* dstV = dstYuv.getCrAddr(absPartIdx);

    primitives.chroma[m_csp].copy_pp[part](dstU, dstYuv.m_csize, const_cast<pixel*>(srcU), m_csize);
    primitives.chroma[m_csp].copy_pp[part](dstV, dstYuv.m_csize, const_cast<pixel*>(srcV), m_csize);
}
