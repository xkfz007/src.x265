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

#ifndef X265_SLICETYPE_H
#define X265_SLICETYPE_H

#include "common.h"
#include "slice.h"
#include "motion.h"
#include "piclist.h"
#include "threadpool.h"

namespace x265 {
// private namespace

struct Lowres;
class Frame;
class Lookahead;

#define LOWRES_COST_MASK  ((1 << 14) - 1)
#define LOWRES_COST_SHIFT 14

/* Thread local data for lookahead tasks */
struct LookaheadTLD
{
    MotionEstimate  me;
    ReferencePlanes weightedRef;
    pixel*          wbuffer[4];
    int             widthInCU;
    int             heightInCU;
    int             ncu;
    int             paddedLines;

#if DETAILED_CU_STATS
    int64_t         batchElapsedTime;
    int64_t         coopSliceElapsedTime;
    uint64_t        countBatches;
    uint64_t        countCoopSlices;
#endif

    LookaheadTLD()
    {
        me.setQP(X265_LOOKAHEAD_QP);
        me.init(X265_HEX_SEARCH, 1, X265_CSP_I400);
        for (int i = 0; i < 4; i++)
            wbuffer[i] = NULL;
        widthInCU = heightInCU = ncu = paddedLines = 0;

#if DETAILED_CU_STATS
        batchElapsedTime = 0;
        coopSliceElapsedTime = 0;
        countBatches = 0;
        countCoopSlices = 0;
#endif
    }

    void init(int w, int h, int n)
    {
        widthInCU = w;
        heightInCU = h;
        ncu = n;
    }

    ~LookaheadTLD() { X265_FREE(wbuffer[0]); }

    void calcAdaptiveQuantFrame(Frame *curFrame, x265_param* param);
    void lowresIntraEstimate(Lowres& fenc);

    void weightsAnalyse(Lowres& fenc, Lowres& ref);

protected:

    uint32_t acEnergyCu(Frame* curFrame, uint32_t blockX, uint32_t blockY, int csp);
    uint32_t weightCostLuma(Lowres& fenc, Lowres& ref, WeightParam& wp);
    bool     allocWeightedRef(Lowres& fenc);
};

class Lookahead : public JobProvider
{
public:

    PicList       m_inputQueue;      // input pictures in order received
    PicList       m_outputQueue;     // pictures to be encoded, in encode order
    Lock          m_inputLock;
    Lock          m_outputLock;

    /* pre-lookahead */
    Frame*        m_preframes[X265_LOOKAHEAD_MAX];
    int           m_preTotal, m_preAcquired, m_preCompleted;
    int           m_fullQueueSize;
    bool          m_isActive;
    bool          m_sliceTypeBusy;
    bool          m_bAdaptiveQuant;
    bool          m_outputSignalRequired;
    bool          m_bBatchMotionSearch;
    bool          m_bBatchFrameCosts;
    Lock          m_preLookaheadLock;
    Event         m_outputSignal;

    LookaheadTLD* m_tld;
    x265_param*   m_param;
    Lowres*       m_lastNonB;
    int*          m_scratch;         // temp buffer for cutree propagate
    
    int           m_histogram[X265_BFRAME_MAX + 1];
    int           m_lastKeyframe;
    int           m_8x8Width;
    int           m_8x8Height;
    int           m_8x8Blocks;
    int           m_numCoopSlices;
    int           m_numRowsPerSlice;
    bool          m_filled;

    Lookahead(x265_param *param, ThreadPool *pool);

#if DETAILED_CU_STATS
    int64_t       m_slicetypeDecideElapsedTime;
    int64_t       m_preLookaheadElapsedTime;
    uint64_t      m_countSlicetypeDecide;
    uint64_t      m_countPreLookahead;
    void          getWorkerStats(int64_t& batchElapsedTime, uint64_t& batchCount, int64_t& coopSliceElapsedTime, uint64_t& coopSliceCount);
#endif

    bool    create();
    void    destroy();
    void    stop();

    void    addPicture(Frame&, int sliceType);
    void    flush();
    Frame*  getDecidedPicture();

    void    getEstimatedPictureCost(Frame *pic);


protected:

    void    findJob(int workerThreadID);
    void    slicetypeDecide();
    void    slicetypeAnalyse(Lowres **frames, bool bKeyframe);

    /* called by slicetypeAnalyse() to make slice decisions */
    bool    scenecut(Lowres **frames, int p0, int p1, bool bRealScenecut, int numFrames, int maxSearch);
    bool    scenecutInternal(Lowres **frames, int p0, int p1, bool bRealScenecut);
    void    slicetypePath(Lowres **frames, int length, char(*best_paths)[X265_LOOKAHEAD_MAX + 1]);
    int64_t slicetypePathCost(Lowres **frames, char *path, int64_t threshold);
    int64_t vbvFrameCost(Lowres **frames, int p0, int p1, int b);
    void    vbvLookahead(Lowres **frames, int numFrames, int keyframes);

    /* called by slicetypeAnalyse() to effect cuTree adjustments to adaptive
     * quant offsets */
    void    cuTree(Lowres **frames, int numframes, bool bintra);
    void    estimateCUPropagate(Lowres **frames, double average_duration, int p0, int p1, int b, int referenced);
    void    cuTreeFinish(Lowres *frame, double averageDuration, int ref0Distance);

    /* called by getEstimatedPictureCost() to finalize cuTree costs */
    int64_t frameCostRecalculate(Lowres **frames, int p0, int p1, int b);
};

class CostEstimateGroup : public BondedTaskGroup
{
public:

    Lookahead& m_lookahead;
    Lowres**   m_frames;
    bool       m_batchMode;

    CostEstimateGroup(Lookahead& l, Lowres** f) : m_lookahead(l), m_frames(f), m_batchMode(false) {}

    /* Cooperative cost estimate using multiple slices of downscaled frame */
    struct Coop
    {
        int  p0, b, p1;
        bool bDoSearch[2];
    } m_coop;

    enum { MAX_COOP_SLICES = 32 };
    struct Slice
    {
        int  costEst;
        int  costEstAq;
        int  intraMbs;
    } m_slice[MAX_COOP_SLICES];

    int64_t singleCost(int p0, int p1, int b, bool intraPenalty = false);

    /* Batch cost estimates, using one worker thread per estimateFrameCost() call */
    enum { MAX_BATCH_SIZE = 512 };
    struct Estimate
    {
        int  p0, b, p1;
    } m_estimates[MAX_BATCH_SIZE];

    void add(int p0, int p1, int b);
    void finishBatch();

protected:

    static const int s_merange = 16;

    void    processTasks(int workerThreadID);

    int64_t estimateFrameCost(LookaheadTLD& tld, int p0, int p1, int b, bool intraPenalty);
    void    estimateCUCost(LookaheadTLD& tld, int cux, int cuy, int p0, int p1, int b, bool bDoSearch[2], bool lastRow, int slice);

    CostEstimateGroup& operator=(const CostEstimateGroup&);
};

}

#endif // ifndef X265_SLICETYPE_H