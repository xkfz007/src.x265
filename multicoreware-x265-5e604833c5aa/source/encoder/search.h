/*****************************************************************************
* Copyright (C) 2013 x265 project
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

#ifndef X265_SEARCH_H
#define X265_SEARCH_H

#include "common.h"
#include "predict.h"
#include "quant.h"
#include "bitcost.h"
#include "yuv.h"
#include "threadpool.h"

#include "rdcost.h"
#include "entropy.h"
#include "motion.h"

#define MVP_IDX_BITS 1
#define NUM_LAYERS 4

namespace x265 {
// private namespace

class Entropy;
struct ThreadLocalData;

/* All the CABAC contexts that Analysis needs to keep track of at each depth
 * and temp buffers for residual, coeff, and recon for use during residual
 * quad-tree depth recursion */
struct RQTData
{
    Entropy  cur;     /* starting context for current CU */

    /* these are indexed by qtLayer (log2size - 2) so nominally 0=4x4, 1=8x8, 2=16x16, 3=32x32
     * the coeffRQT and reconQtYuv are allocated to the max CU size at every depth. The parts
     * which are reconstructed at each depth are valid. At the end, the transform depth table
     * is walked and the coeff and recon at the final split depths are collected */
    Entropy  rqtRoot;      /* residual quad-tree start context */
    Entropy  rqtTemp;      /* residual quad-tree temp context */
    Entropy  rqtTest;      /* residual quad-tree test context */
    coeff_t* coeffRQT[3];  /* coeff storage for entire CTU for each RQT layer */
    Yuv      reconQtYuv;   /* recon storage for entire CTU for each RQT layer (intra) */
    ShortYuv resiQtYuv;    /* residual storage for entire CTU for each RQT layer (inter) */
    
    /* per-depth temp buffers for inter prediction */
    ShortYuv tmpResiYuv;
    Yuv      tmpPredYuv;
    Yuv      bidirPredYuv[2];
};

inline int getTUBits(int idx, int numIdx)
{
    return idx + (idx < numIdx - 1);
}

class Search : public JobProvider, public Predict
{
public:

    static const pixel   zeroPixel[MAX_CU_SIZE];
    static const int16_t zeroShort[MAX_CU_SIZE];

    MotionEstimate  m_me;
    Quant           m_quant;
    RDCost          m_rdCost;
    const x265_param* m_param;
    Frame*          m_frame;
    const Slice*    m_slice;

    Entropy         m_entropyCoder;
    RQTData         m_rqt[NUM_FULL_DEPTH];

    uint8_t*        m_qtTempCbf[3];
    uint8_t*        m_qtTempTransformSkipFlag[3];

    bool            m_bFrameParallel;
    bool            m_bEnableRDOQ;
    uint32_t        m_numLayers;
    uint32_t        m_refLagPixels;

    struct Mode
    {
        CUData     cu;
        const Yuv* fencYuv;
        Yuv        predYuv;
        Yuv        reconYuv;
        Entropy    contexts;

        uint64_t   rdCost;     // sum of partition (psy) RD costs          (sse(fenc, recon) + lambda2 * bits)
        uint64_t   sa8dCost;   // sum of partition sa8d distortion costs   (sa8d(fenc, pred) + lambda * bits)
        uint32_t   sa8dBits;   // signal bits used in sa8dCost calculation
        uint32_t   psyEnergy;  // sum of partition psycho-visual energy difference
        uint32_t   distortion; // sum of partition SSE distortion
        uint32_t   totalBits;  // sum of partition bits (mv + coeff)
        uint32_t   mvBits;     // Mv bits + Ref + block type (or intra mode)
        uint32_t   coeffBits;  // Texture bits (DCT Coeffs)

        void initCosts()
        {
            rdCost = 0;
            sa8dCost = 0;
            sa8dBits = 0;
            psyEnergy = 0;
            distortion = 0;
            totalBits = 0;
            mvBits = 0;
            coeffBits = 0;
        }

        void addSubCosts(const Mode& subMode)
        {
            rdCost += subMode.rdCost;
            sa8dCost += subMode.sa8dCost;
            sa8dBits += subMode.sa8dBits;
            psyEnergy += subMode.psyEnergy;
            distortion += subMode.distortion;
            totalBits += subMode.totalBits;
            mvBits += subMode.mvBits;
            coeffBits += subMode.coeffBits;
        }
    };

    struct MotionData
    {
        MV  mv;
        MV  mvp;
        int mvpIdx;
        int ref;
        uint32_t cost;
        int bits;
    };

    Search();
    ~Search();

    bool     initSearch(const x265_param& param, ScalingList& scalingList);
    void     setQP(const Slice& slice, int qp);

    // mark temp RD entropy contexts as uninitialized; useful for finding loads without stores
    void     invalidateContexts(int fromDepth);

    // full RD search of intra modes. if sharedModes is not NULL, it directly uses them
    void     checkIntra(Mode& intraMode, const CUGeom& cuGeom, PartSize partSize, uint8_t* sharedModes);

    // estimation inter prediction (non-skip)
    bool     predInterSearch(Mode& interMode, const CUGeom& cuGeom, bool bMergeOnly, bool bChroma);

    // encode residual and compute rd-cost for inter mode
    void     encodeResAndCalcRdInterCU(Mode& interMode, const CUGeom& cuGeom);
    void     encodeResAndCalcRdSkipCU(Mode& interMode);

    void     generateCoeffRecon(Mode& mode, const CUGeom& cuGeom);
    void     residualTransformQuantInter(Mode& mode, const CUGeom& cuGeom, uint32_t absPartIdx, uint32_t depth, uint32_t depthRange[2]);

    uint32_t getIntraRemModeBits(CUData & cu, uint32_t absPartIdx, uint32_t preds[3], uint64_t& mpms) const;

protected:

    /* motion estimation distribution */
    ThreadLocalData* m_tld;
    CUData*       m_curMECu;
    const CUGeom* m_curGeom;
    int           m_curPart;
    MotionData    m_bestME[2];
    uint32_t      m_listSelBits[3];
    int           m_totalNumME;
    volatile int  m_numAcquiredME;
    volatile int  m_numCompletedME;
    Event         m_meCompletionEvent;
    Lock          m_outputLock;
    bool          m_bJobsQueued;
    void     singleMotionEstimation(Search& master, const CUData& cu, const CUGeom& cuGeom, int part, int list, int ref);

    void     saveResidualQTData(CUData& cu, ShortYuv& resiYuv, uint32_t absPartIdx, uint32_t depth);

    // RDO search of luma intra modes; result is fully encoded luma. luma distortion is returned
    uint32_t estIntraPredQT(Mode &intraMode, const CUGeom& cuGeom, uint32_t depthRange[2], uint8_t* sharedModes);

    // RDO select best chroma mode from luma; result is fully encode chroma. chroma distortion is returned
    uint32_t estIntraPredChromaQT(Mode &intraMode, const CUGeom& cuGeom);

    void     codeSubdivCbfQTChroma(const CUData& cu, uint32_t trDepth, uint32_t absPartIdx,  uint32_t absPartIdxStep, uint32_t width, uint32_t height);
    void     codeCoeffQTChroma(const CUData& cu, uint32_t trDepth, uint32_t absPartIdx, TextType ttype);

    struct Cost
    {
        uint64_t rdcost;
        uint32_t bits;
        uint32_t distortion;
        uint32_t energy;
        Cost() { rdcost = 0; bits = 0; distortion = 0; energy = 0; }
    };

    void     estimateResidualQT(Mode& mode, const CUGeom& cuGeom, uint32_t absPartIdx, uint32_t depth, ShortYuv& resiYuv, Cost& costs, uint32_t depthRange[2]);

    void     encodeResidualQT(CUData& cu, uint32_t absPartIdx, uint32_t depth, bool bSubdivAndCbf, TextType ttype, uint32_t depthRange[2]);

    // generate prediction, generate residual and recon. if bAllowSplit, find optimal RQT splits
    void     codeIntraLumaQT(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t absPartIdx, bool bAllowSplit, Cost& costs, uint32_t depthRange[2]);
    void     codeIntraLumaTSkip(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t absPartIdx, Cost& costs);
    void     extractIntraResultQT(CUData& cu, Yuv& reconYuv, uint32_t trDepth, uint32_t absPartIdx);

    // generate chroma prediction, generate residual and recon
    uint32_t codeIntraChromaQt(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t absPartIdx, uint32_t& psyEnergy);
    uint32_t codeIntraChromaTSkip(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t trDepthC, uint32_t absPartIdx, uint32_t& psyEnergy);
    void     extractIntraResultChromaQT(CUData& cu, Yuv& reconYuv, uint32_t absPartIdx, uint32_t trDepth, bool tuQuad);

    void     residualTransformQuantIntra(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t absPartIdx, uint32_t depthRange[2]);
    void     residualQTIntraChroma(Mode& mode, const CUGeom& cuGeom, uint32_t trDepth, uint32_t absPartIdx);

    void     offsetSubTUCBFs(CUData& cu, TextType ttype, uint32_t trDepth, uint32_t absPartIdx);

    struct MergeData
    {
        /* merge candidate data, cached between calls to mergeEstimation */
        MVField  mvFieldNeighbours[MRG_MAX_NUM_CANDS][2];
        uint8_t  interDirNeighbours[MRG_MAX_NUM_CANDS];
        uint32_t maxNumMergeCand;

        /* data updated for each partition */
        uint32_t absPartIdx;
        int      width;
        int      height;

        /* outputs */
        MVField  mvField[2];
        uint32_t interDir;
        uint32_t index;
        uint32_t bits;
    };

    /* inter/ME helper functions */
    void     checkBestMVP(MV* amvpCand, MV cMv, MV& mvPred, int& mvpIdx, uint32_t& outBits, uint32_t& outCost) const;
    void     setSearchRange(const CUData& cu, MV mvp, int merange, MV& mvmin, MV& mvmax) const;
    uint32_t mergeEstimation(CUData& cu, const CUGeom& cuGeom, int partIdx, MergeData& m);
    static void getBlkBits(PartSize cuMode, bool bPSlice, int partIdx, uint32_t lastMode, uint32_t blockBit[3]);

    /* intra helper functions */
    enum { MAX_RD_INTRA_MODES = 16 };
    static void updateCandList(uint32_t mode, uint64_t cost, int maxCandCount, uint32_t* candModeList, uint64_t* candCostList);
    void     getBestIntraModeChroma(Mode& intraMode, const CUGeom& cuGeom);

    void updateModeCost(Mode& m) const { m.rdCost = m_rdCost.m_psyRd ? m_rdCost.calcPsyRdCost(m.distortion, m.totalBits, m.psyEnergy) : m_rdCost.calcRdCost(m.distortion, m.totalBits); }
};
}

#endif // ifndef X265_SEARCH_H
