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

#ifndef X265_CUDATA_H
#define X265_CUDATA_H

#include "common.h"
#include "slice.h"
#include "mv.h"

namespace x265 {
// private namespace

class FrameData;
class Slice;
struct TUEntropyCodingParameters;
struct CUDataMemPool;

enum PartSize
{
    SIZE_2Nx2N, // symmetric motion partition,  2Nx2N
    SIZE_2NxN,  // symmetric motion partition,  2Nx N
    SIZE_Nx2N,  // symmetric motion partition,   Nx2N
    SIZE_NxN,   // symmetric motion partition,   Nx N
    SIZE_2NxnU, // asymmetric motion partition, 2Nx( N/2) + 2Nx(3N/2)
    SIZE_2NxnD, // asymmetric motion partition, 2Nx(3N/2) + 2Nx( N/2)
    SIZE_nLx2N, // asymmetric motion partition, ( N/2)x2N + (3N/2)x2N
    SIZE_nRx2N, // asymmetric motion partition, (3N/2)x2N + ( N/2)x2N
    SIZE_NONE = 15
};

enum PredMode
{
    MODE_INTER,
    MODE_INTRA,
    MODE_NONE = 15
};

// motion vector predictor direction used in AMVP
enum MVP_DIR
{
    MD_LEFT = 0,    // MVP of left block
    MD_ABOVE,       // MVP of above block
    MD_ABOVE_RIGHT, // MVP of above right block
    MD_BELOW_LEFT,  // MVP of below left block
    MD_ABOVE_LEFT   // MVP of above left block
};

struct CUGeom
{
    enum {
        INTRA           = 1<<0, // CU is intra predicted
        PRESENT         = 1<<1, // CU is not completely outside the frame
        SPLIT_MANDATORY = 1<<2, // CU split is mandatory if CU is inside frame and can be split
        LEAF            = 1<<3, // CU is a leaf node of the CTU
        SPLIT           = 1<<4, // CU is currently split in four child CUs.
    };
    
    // (1 + 4 + 16 + 64) = 85.
    enum { MAX_GEOMS = 85 };

    uint32_t log2CUSize;    // Log of the CU size.
    uint32_t childOffset;   // offset of the first child CU from current CU
    uint32_t encodeIdx;     // Encoding index of this CU in terms of 4x4 blocks.
    uint32_t numPartitions; // Number of 4x4 blocks in the CU
    uint32_t depth;         // depth of this CU relative from CTU
    uint32_t flags;         // CU flags.
};

struct MVField
{
    MV  mv;
    int refIdx;
};

typedef void(*cucopy_t)(uint8_t* dst, uint8_t* src); // dst and src are aligned to MIN(size, 32)
typedef void(*cubcast_t)(uint8_t* dst, uint8_t val); // dst is aligned to MIN(size, 32)

// Partition count table, index represents partitioning mode.
const uint32_t nbPartsTable[8] = { 1, 2, 2, 4, 2, 2, 2, 2 };

// Holds part data for a CU of a given size, from an 8x8 CU to a CTU
class CUData
{
public:

    static cubcast_t s_partSet[NUM_FULL_DEPTH]; // pointer to broadcast set functions per absolute depth
    static uint32_t  s_numPartInCUSize;

    FrameData*    m_encData;
    const Slice*  m_slice;

    cucopy_t      m_partCopy;         // pointer to function that copies m_numPartitions elements
    cubcast_t     m_partSet;          // pointer to function that sets m_numPartitions elements
    cucopy_t      m_subPartCopy;      // pointer to function that copies m_numPartitions/4 elements, may be NULL
    cubcast_t     m_subPartSet;       // pointer to function that sets m_numPartitions/4 elements, may be NULL

    uint32_t      m_cuAddr;           // address of CTU within the picture in raster order
    uint32_t      m_absIdxInCTU;      // address of CU within its CTU in Z scan order
    uint32_t      m_cuPelX;           // CU position within the picture, in pixels (X)
    uint32_t      m_cuPelY;           // CU position within the picture, in pixels (Y)
    uint32_t      m_numPartitions;    // maximum number of 4x4 partitions within this CU

    int           m_chromaFormat;
    int           m_hChromaShift;
    int           m_vChromaShift;

    /* Per-part data, stored contiguously */
    char*         m_qp;               // array of QP values
    uint8_t*      m_log2CUSize;       // array of cu log2Size TODO: seems redundant to depth
    uint8_t*      m_partSize;         // array of partition sizes
    uint8_t*      m_predMode;         // array of prediction modes
    uint8_t*      m_lumaIntraDir;     // array of intra directions (luma)
    uint8_t*      m_tqBypass;         // array of CU lossless flags
    char*         m_refIdx[2];        // array of motion reference indices per list
    uint8_t*      m_cuDepth;          // array of depths
    uint8_t*      m_skipFlag;         // array of skip flags
    uint8_t*      m_mergeFlag;        // array of merge flags
    uint8_t*      m_interDir;         // array of inter directions
    uint8_t*      m_mvpIdx[2];        // array of motion vector predictor candidates or merge candidate indices [0]
    uint8_t*      m_tuDepth;          // array of transform indices
    uint8_t*      m_transformSkip[3]; // array of transform skipping flags per plane
    uint8_t*      m_cbf[3];           // array of coded block flags (CBF) per plane
    uint8_t*      m_chromaIntraDir;   // array of intra directions (chroma)
    enum { BytesPerPartition = 22 };  // combined sizeof() of all per-part data

    coeff_t*      m_trCoeff[3];       // transformed coefficient buffer per plane

    MV*           m_mv[2];            // array of motion vectors per list
    MV*           m_mvd[2];           // array of coded motion vector deltas per list
    enum { TMVP_UNIT_MASK = 0xF0 };  // mask for mapping index to into a compressed (reference) MV field

    const CUData* m_cuAboveLeft;      // pointer to above-left neighbor CTU
    const CUData* m_cuAboveRight;     // pointer to above-right neighbor CTU
    const CUData* m_cuAbove;          // pointer to above neighbor CTU
    const CUData* m_cuLeft;           // pointer to left neighbor CTU

    CUData();

    void     initialize(const CUDataMemPool& dataPool, uint32_t depth, int csp, int instance);
    void     calcCTUGeoms(uint32_t picWidth, uint32_t picHeight, uint32_t maxCUSize, CUGeom cuDataArray[CUGeom::MAX_GEOMS]) const;

    void     initCTU(const Frame& frame, uint32_t cuAddr, int qp);
    void     initSubCU(const CUData& ctu, const CUGeom& cuGeom);
    void     initLosslessCU(const CUData& cu, const CUGeom& cuGeom);

    void     copyPartFrom(const CUData& cu, const CUGeom& childGeom, uint32_t subPartIdx);
    void     setEmptyPart(const CUGeom& childGeom, uint32_t subPartIdx);
    void     copyToPic(uint32_t depth) const;

    /* RD-0 methods called only from encodeResidue */
    void     copyFromPic(const CUData& ctu, const CUGeom& cuGeom);
    void     updatePic(uint32_t depth) const;

    void     setPartSizeSubParts(PartSize size)    { m_partSet(m_partSize, (uint8_t)size); }
    void     setSkipFlagSubParts(uint8_t skipFlag) { m_partSet(m_skipFlag, skipFlag); }
    void     setPredModeSubParts(PredMode mode)    { m_partSet(m_predMode, (uint8_t)mode); }
    void     clearCbf()                            { m_partSet(m_cbf[0], 0); m_partSet(m_cbf[1], 0); m_partSet(m_cbf[2], 0); }

    /* these functions all take depth as an absolute depth from CTU, it is used to calculate the number of parts to copy */
    void     setQPSubParts(char qp, uint32_t absPartIdx, uint32_t depth)                      { s_partSet[depth]((uint8_t*)m_qp + absPartIdx, (uint8_t)qp); }
    void     setTUDepthSubParts(uint8_t tuDepth, uint32_t absPartIdx, uint32_t depth)         { s_partSet[depth](m_tuDepth + absPartIdx, tuDepth); }
    void     setLumaIntraDirSubParts(uint8_t dir, uint32_t absPartIdx, uint32_t depth)        { s_partSet[depth](m_lumaIntraDir + absPartIdx, dir); }
    void     setChromIntraDirSubParts(uint8_t dir, uint32_t absPartIdx, uint32_t depth)       { s_partSet[depth](m_chromaIntraDir + absPartIdx, dir); }
    void     setCbfSubParts(uint8_t cbf, TextType ttype, uint32_t absPartIdx, uint32_t depth) { s_partSet[depth](m_cbf[ttype] + absPartIdx, cbf); }
    void     setCbfPartRange(uint8_t cbf, TextType ttype, uint32_t absPartIdx, uint32_t coveredPartIdxes) { memset(m_cbf[ttype] + absPartIdx, cbf, coveredPartIdxes); }
    void     setTransformSkipSubParts(uint8_t tskip, TextType ttype, uint32_t absPartIdx, uint32_t depth) { s_partSet[depth](m_transformSkip[ttype] + absPartIdx, tskip); }
    void     setTransformSkipPartRange(uint8_t tskip, TextType ttype, uint32_t absPartIdx, uint32_t coveredPartIdxes) { memset(m_transformSkip[ttype] + absPartIdx, tskip, coveredPartIdxes); }

    bool     setQPSubCUs(char qp, uint32_t absPartIdx, uint32_t depth);

    void     setPUInterDir(uint8_t dir, uint32_t absPartIdx, uint32_t puIdx);
    void     setPUMv(int list, const MV& mv, int absPartIdx, int puIdx);
    void     setPURefIdx(int list, char refIdx, int absPartIdx, int puIdx);

    uint8_t  getCbf(uint32_t absPartIdx, TextType ttype, uint32_t trDepth) const { return (m_cbf[ttype][absPartIdx] >> trDepth) & 0x1; }
    uint8_t  getQtRootCbf(uint32_t absPartIdx) const                             { return m_cbf[0][absPartIdx] || m_cbf[1][absPartIdx] || m_cbf[2][absPartIdx]; }
    char     getRefQP(uint32_t currAbsIdxInCTU) const;
    uint32_t getInterMergeCandidates(uint32_t absPartIdx, uint32_t puIdx, MVField (*mvFieldNeighbours)[2], uint8_t* interDirNeighbours) const;
    void     clipMv(MV& outMV) const;
    int      fillMvpCand(uint32_t puIdx, uint32_t absPartIdx, int picList, int refIdx, MV* amvpCand, MV* mvc) const;
    void     getIntraTUQtDepthRange(uint32_t tuDepthRange[2], uint32_t absPartIdx) const;
    void     getInterTUQtDepthRange(uint32_t tuDepthRange[2], uint32_t absPartIdx) const;

    uint32_t getNumPartInter() const              { return nbPartsTable[(int)m_partSize[0]]; }
    bool     isIntra(uint32_t absPartIdx) const   { return m_predMode[absPartIdx] == MODE_INTRA; }
    bool     isSkipped(uint32_t absPartIdx) const { return !!m_skipFlag[absPartIdx]; }
    bool     isBipredRestriction() const          { return m_log2CUSize[0] == 3 && m_partSize[0] != SIZE_2Nx2N; }

    void     getPartIndexAndSize(uint32_t puIdx, uint32_t& absPartIdx, int& puWidth, int& puHeight) const;
    void     getMvField(const CUData* cu, uint32_t absPartIdx, int picList, MVField& mvField) const;

    void     getAllowedChromaDir(uint32_t absPartIdx, uint32_t* modeList) const;
    int      getIntraDirLumaPredictor(uint32_t absPartIdx, uint32_t* intraDirPred) const;
    void     deriveLeftRightTopIdxAdi(uint32_t& partIdxLT, uint32_t& partIdxRT, uint32_t partOffset, uint32_t partDepth) const;

    uint32_t getSCUAddr() const                  { return (m_cuAddr << g_maxFullDepth * 2) + m_absIdxInCTU; }
    uint32_t getCtxSplitFlag(uint32_t absPartIdx, uint32_t depth) const;
    uint32_t getCtxSkipFlag(uint32_t absPartIdx) const;
    ScanType getCoefScanIdx(uint32_t absPartIdx, uint32_t log2TrSize, bool bIsLuma, bool bIsIntra) const;
    void     getTUEntropyCodingParameters(TUEntropyCodingParameters &result, uint32_t absPartIdx, uint32_t log2TrSize, bool bIsLuma) const;

    const CUData* getPULeft(uint32_t& lPartUnitIdx, uint32_t curPartUnitIdx) const;
    const CUData* getPUAbove(uint32_t& aPartUnitIdx, uint32_t curPartUnitIdx, bool planarAtCTUBoundary = false) const;
    const CUData* getPUAboveLeft(uint32_t& alPartUnitIdx, uint32_t curPartUnitIdx) const;
    const CUData* getPUAboveRight(uint32_t& arPartUnitIdx, uint32_t curPartUnitIdx) const;
    const CUData* getPUBelowLeft(uint32_t& blPartUnitIdx, uint32_t curPartUnitIdx) const;

    const CUData* getQpMinCuLeft(uint32_t& lPartUnitIdx, uint32_t currAbsIdxInCTU) const;
    const CUData* getQpMinCuAbove(uint32_t& aPartUnitIdx, uint32_t currAbsIdxInCTU) const;

    const CUData* getPUAboveRightAdi(uint32_t& arPartUnitIdx, uint32_t curPartUnitIdx, uint32_t partUnitOffset) const;
    const CUData* getPUBelowLeftAdi(uint32_t& blPartUnitIdx, uint32_t curPartUnitIdx, uint32_t partUnitOffset) const;

protected:

    template<typename T>
    void setAllPU(T *p, const T& val, int absPartIdx, int puIdx);

    char getLastCodedQP(uint32_t absPartIdx) const;
    int  getLastValidPartIdx(int absPartIdx) const;

    bool hasEqualMotion(uint32_t absPartIdx, const CUData& candCU, uint32_t candAbsPartIdx) const;

    bool isDiffMER(int xN, int yN, int xP, int yP) const;

    // add possible motion vector predictor candidates
    bool addMVPCand(MV& mvp, int picList, int refIdx, uint32_t absPartIdx, MVP_DIR dir) const;
    bool addMVPCandOrder(MV& mvp, int picList, int refIdx, uint32_t absPartIdx, MVP_DIR dir) const;

    bool getColMVP(MV& outMV, int& outRefIdx, int picList, int cuAddr, int absPartIdx) const;

    void scaleMvByPOCDist(MV& outMV, const MV& inMV, int curPOC, int curRefPOC, int colPOC, int colRefPOC) const;

    void     deriveLeftRightTopIdx(uint32_t puIdx, uint32_t& partIdxLT, uint32_t& partIdxRT) const;

    uint32_t deriveCenterIdx(uint32_t puIdx) const;
    uint32_t deriveRightBottomIdx(uint32_t puIdx) const;
    uint32_t deriveLeftBottomIdx(uint32_t puIdx) const;
};

// TU settings for entropy encoding
struct TUEntropyCodingParameters
{
    const uint16_t *scan;
    const uint16_t *scanCG;
    ScanType        scanType;
    uint32_t        log2TrSizeCG;
    uint32_t        firstSignificanceMapContext;
};

struct CUDataMemPool
{
    uint8_t* charMemBlock;
    coeff_t* trCoeffMemBlock;
    MV*      mvMemBlock;

    CUDataMemPool() { charMemBlock = NULL; trCoeffMemBlock = NULL; mvMemBlock = NULL; }

    bool create(uint32_t depth, uint32_t csp, uint32_t numInstances)
    {
        uint32_t numPartition = NUM_CU_PARTITIONS >> (depth * 2);
        uint32_t cuSize = g_maxCUSize >> depth;
        uint32_t sizeL = cuSize * cuSize;
        uint32_t sizeC = sizeL >> (CHROMA_H_SHIFT(csp) + CHROMA_V_SHIFT(csp));
        CHECKED_MALLOC(trCoeffMemBlock, coeff_t, (sizeL + sizeC * 2) * numInstances);
        CHECKED_MALLOC(charMemBlock, uint8_t, numPartition * numInstances * CUData::BytesPerPartition);
        CHECKED_MALLOC(mvMemBlock, MV, numPartition * 4 * numInstances);
        return true;

    fail:
        return false;
    }

    void destroy()
    {
        X265_FREE(trCoeffMemBlock);
        X265_FREE(mvMemBlock);
        X265_FREE(charMemBlock);
    }
};
}

#endif // ifndef X265_CUDATA_H
