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

#ifndef X265_DEBLOCK_H
#define X265_DEBLOCK_H

#include "common.h"

namespace x265 {
// private namespace

class CUData;

class Deblock
{
public:
    enum { EDGE_VER, EDGE_HOR };

    uint32_t m_numPartitions;

    Deblock() : m_numPartitions(0) {}

    void init() { m_numPartitions = 1 << (g_maxFullDepth * 2); }

    void deblockCTU(CUData* cu, int32_t dir);

protected:

    // CU-level deblocking function
    void deblockCU(CUData* cu, uint32_t absZOrderIdx, uint32_t depth, const int32_t Edge, uint8_t blockingStrength[]);

    struct Param
    {
        uint8_t leftEdge;
        uint8_t topEdge;
    };

    // set filtering functions
    void setLoopfilterParam(CUData* cu, uint32_t absZOrderIdx, Param *params);
    void setEdgefilterTU(CUData* cu, uint32_t absZOrderIdx, uint32_t depth, int32_t dir, uint8_t blockingStrength[]);
    void setEdgefilterPU(CUData* cu, uint32_t absZOrderIdx, int32_t dir, uint8_t blockingStrength[], uint32_t widthInBaseUnits);
    void setEdgefilterMultiple(CUData* cu, uint32_t absZOrderIdx, int32_t dir, int32_t edgeIdx, uint8_t value, uint8_t blockingStrength[], uint32_t widthInBaseUnits);

    // get filtering functions
    void getBoundaryStrengthSingle(CUData* cu, int32_t dir, uint32_t partIdx, uint8_t blockingStrength[]);

    // filter luma/chroma functions
    void edgeFilterLuma(CUData* cu, uint32_t absZOrderIdx, uint32_t depth, int32_t dir, int32_t edge, const uint8_t blockingStrength[]);
    void edgeFilterChroma(CUData* cu, uint32_t absZOrderIdx, uint32_t depth, int32_t dir, int32_t edge, const uint8_t blockingStrength[]);

    static const uint8_t s_tcTable[54];
    static const uint8_t s_betaTable[52];
};
}
#endif // ifndef X265_DEBLOCK_H
