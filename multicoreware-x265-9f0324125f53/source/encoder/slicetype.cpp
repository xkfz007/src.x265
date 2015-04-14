/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Gopu Govindaswamy <gopu@multicorewareinc.com>
 *          Steve Borho <steve@borho.org>
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
#include "frame.h"
#include "framedata.h"
#include "picyuv.h"
#include "primitives.h"
#include "lowres.h"
#include "mv.h"

#include "slicetype.h"
#include "motion.h"
#include "ratecontrol.h"

#define NUM_CUS (m_widthInCU > 2 && m_heightInCU > 2 ? (m_widthInCU - 2) * (m_heightInCU - 2) : m_widthInCU * m_heightInCU)

using namespace x265;

static inline int16_t median(int16_t a, int16_t b, int16_t c)
{
    int16_t t = (a - b) & ((a - b) >> 31);

    a -= t;
    b += t;
    b -= (b - c) & ((b - c) >> 31);
    b += (a - b) & ((a - b) >> 31);
    return b;
}

static inline void median_mv(MV &dst, MV a, MV b, MV c)
{
    dst.x = median(a.x, b.x, c.x);
    dst.y = median(a.y, b.y, c.y);
}

Lookahead::Lookahead(x265_param *param, ThreadPool* pool)
    : JobProvider(pool)
    , m_est(pool)
{
    m_bReady = false;
    m_bBusy = false;
    m_param = param;
    m_lastKeyframe = -m_param->keyframeMax;
    m_lastNonB = NULL;
    m_bFilled = false;
    m_bFlushed = false;
    m_bFlush = false;
    m_widthInCU = ((m_param->sourceWidth / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
    m_heightInCU = ((m_param->sourceHeight / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
    m_scratch = (int*)x265_malloc(m_widthInCU * sizeof(int));
    memset(m_histogram, 0, sizeof(m_histogram));
}

Lookahead::~Lookahead() { }

void Lookahead::init()
{
    if (m_pool && m_pool->getThreadCount() >= 4 &&
        ((m_param->bFrameAdaptive && m_param->bframes) ||
         m_param->rc.cuTree || m_param->scenecutThreshold ||
         (m_param->lookaheadDepth && m_param->rc.vbvBufferSize)))
    {
        JobProvider::enqueue();
    }
    else
        m_pool = NULL; /* disable use of worker thread */
}

void Lookahead::stop()
{
    /* do not allow slicetypeDecide() to get started again */
    m_bReady = false;
    m_bFlushed = false;
    m_bFlush = false;
    m_bBusy = false;

    if (m_pool)
        JobProvider::flush(); // flush will dequeue, if it is necessary
}

void Lookahead::destroy()
{
    // these two queues will be empty unless the encode was aborted
    while (!m_inputQueue.empty())
    {
        Frame* curFrame = m_inputQueue.popFront();
        curFrame->destroy();
        delete curFrame;
    }

    while (!m_outputQueue.empty())
    {
        Frame* curFrame = m_outputQueue.popFront();
        curFrame->destroy();
        delete curFrame;
    }

    x265_free(m_scratch);
}

/* Called by API thread */
void Lookahead::addPicture(Frame *curFrame, int sliceType)
{
    {
        ProfileScopeEvent(prelookahead);
        PicYuv *orig = curFrame->m_fencPic;
        curFrame->m_lowres.init(orig, curFrame->m_poc, sliceType);
    }

    m_inputQueueLock.acquire();
    m_inputQueue.pushBack(*curFrame);

    if (m_inputQueue.size() >= m_param->lookaheadDepth)
    {
        if (m_pool)
        {
            m_bReady = !m_bBusy;
            m_inputQueueLock.release();
            m_pool->pokeIdleThread();
        }
        else
            slicetypeDecide();
    }
    else
        m_inputQueueLock.release();

    /* determine if the lookahead is (over) filled enough for frames to begin to
     * be consumed by frame encoders */
    if (!m_bFilled)
    {
        if (!m_param->bframes & !m_param->lookaheadDepth)
            m_bFilled = true; /* zero-latency */
        else if (curFrame->m_poc >= m_param->lookaheadDepth + 2 + m_param->bframes)
            m_bFilled = true; /* full capacity plus mini-gop lag */
    }
}

/* Called by API thread */
void Lookahead::flush()
{
    m_bFlush = true;
    m_bFilled = true;

    /* just in case the input queue is never allowed to fill */
    m_inputQueueLock.acquire();
    if (m_inputQueue.empty())
    {
        m_bFlushed = true;
        m_inputQueueLock.release();
    }
    else
    {
        if (m_pool)
        {
            m_bReady = !m_bBusy;
            m_inputQueueLock.release();
            m_pool->pokeIdleThread();
        }
        else
            slicetypeDecide();
    }
}

/* Called by API thread. If the lookahead queue has not yet been filled the
 * first time, it immediately returns NULL.  Else the function blocks until
 * outputs are available and then pops the first frame from the output queue. If
 * flush() has been called and the output queue is empty, NULL is returned. */
Frame* Lookahead::getDecidedPicture()
{
    if (!m_bFilled)
        return NULL;

    m_outputQueueLock.acquire();
    Frame *fenc = m_outputQueue.popFront();
    m_outputQueueLock.release();

    if (fenc || m_bFlushed)
        return fenc;

    do
    {
        m_outputAvailable.wait();

        m_outputQueueLock.acquire();
        fenc = m_outputQueue.popFront();
        m_outputQueueLock.release();
    }
    while (!fenc);

    return fenc;
}

/* Called by pool worker threads */
bool Lookahead::findJob(int)
{
    if (!m_bReady)
        return false;

    m_inputQueueLock.acquire();
    if (!m_bReady)
    {
        m_inputQueueLock.release();
        return false;
    }

    m_bReady = false;
    m_bBusy = true;

    do
    {
        slicetypeDecide(); // releases input queue lock

        m_inputQueueLock.acquire();

        if (!m_bBusy)
            break;
    }
    while (m_inputQueue.size() >= m_param->lookaheadDepth ||
           (m_bFlush && m_inputQueue.size()));

    m_bBusy = false;
    m_inputQueueLock.release();
    return true;
}

/* Called by rate-control to calculate the estimated SATD cost for a given
 * picture.  It assumes dpb->prepareEncode() has already been called for the
 * picture and all the references are established */
void Lookahead::getEstimatedPictureCost(Frame *curFrame)
{
    Lowres *frames[X265_LOOKAHEAD_MAX];

    // POC distances to each reference
    Slice *slice = curFrame->m_encData->m_slice;
    int p0 = 0, p1, b;
    int poc = slice->m_poc;
    int l0poc = slice->m_refPOCList[0][0];
    int l1poc = slice->m_refPOCList[1][0];

    switch (slice->m_sliceType)
    {
    case I_SLICE:
        frames[p0] = &curFrame->m_lowres;
        b = p1 = 0;
        break;

    case P_SLICE:
        b = p1 = poc - l0poc;
        frames[p0] = &slice->m_refPicList[0][0]->m_lowres;
        frames[b] = &curFrame->m_lowres;
        break;

    case B_SLICE:
        b = poc - l0poc;
        p1 = b + l1poc - poc;
        frames[p0] = &slice->m_refPicList[0][0]->m_lowres;
        frames[b] = &curFrame->m_lowres;
        frames[p1] = &slice->m_refPicList[1][0]->m_lowres;
        break;

    default:
        return;
    }

#if DEBUG_FRMCOST_RECALC&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
    fprintf(fp, "ENCODING:poc=%d,type=%d: p0=%d b=%d p1=%d\n", poc, slice->m_sliceType, p0, b, p1);
    fclose(fp);
  }
#endif
    if (m_param->rc.cuTree && !m_param->rc.bStatRead)
        /* update row satds based on cutree offsets */
        curFrame->m_lowres.satdCost = frameCostRecalculate(frames, p0, p1, b);
    else if (m_param->rc.aqMode)
        curFrame->m_lowres.satdCost = curFrame->m_lowres.costEstAq[b - p0][p1 - b];
    else
        curFrame->m_lowres.satdCost = curFrame->m_lowres.costEst[b - p0][p1 - b];

#if DEBUG_FRMCOST_RECALC&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
    fprintf(fp, "poc=%lld p0=%d b=%d p1=%d satd=%d\n", poc, p0, b, p1, curFrame->m_lowres.satdCost);
    fclose(fp);
  }
#endif
    if (m_param->rc.vbvBufferSize && m_param->rc.vbvMaxBitrate)
    {
        /* aggregate lowres row satds to CTU resolution */
        curFrame->m_lowres.lowresCostForRc = curFrame->m_lowres.lowresCosts[b - p0][p1 - b];
        uint32_t lowresRow = 0, lowresCol = 0, lowresCuIdx = 0, sum = 0;
        uint32_t scale = m_param->maxCUSize / (2 * X265_LOWRES_CU_SIZE);
        uint32_t numCuInHeight = (m_param->sourceHeight + g_maxCUSize - 1) / g_maxCUSize;
        uint32_t widthInLowresCu = (uint32_t)m_widthInCU, heightInLowresCu = (uint32_t)m_heightInCU;
        double *qp_offset = 0;
        /* Factor in qpoffsets based on Aq/Cutree in CU costs */
        if (m_param->rc.aqMode)
            qp_offset = (frames[b]->sliceType == X265_TYPE_B || !m_param->rc.cuTree) ? frames[b]->qpAqOffset : frames[b]->qpCuTreeOffset;

        for (uint32_t row = 0; row < numCuInHeight; row++)
        {
            lowresRow = row * scale;
            for (uint32_t cnt = 0; cnt < scale && lowresRow < heightInLowresCu; lowresRow++, cnt++)
            {
                sum = 0;
                lowresCuIdx = lowresRow * widthInLowresCu;
                for (lowresCol = 0; lowresCol < widthInLowresCu; lowresCol++, lowresCuIdx++)
                {
                    uint16_t lowresCuCost = curFrame->m_lowres.lowresCostForRc[lowresCuIdx] & LOWRES_COST_MASK;
                    if (qp_offset)
                    {
                        lowresCuCost = (uint16_t)((lowresCuCost * x265_exp2fix8(qp_offset[lowresCuIdx]) + 128) >> 8);
                        int32_t intraCuCost = curFrame->m_lowres.intraCost[lowresCuIdx]; 
                        curFrame->m_lowres.intraCost[lowresCuIdx] = (intraCuCost * x265_exp2fix8(qp_offset[lowresCuIdx]) + 128) >> 8;
                    }
                    curFrame->m_lowres.lowresCostForRc[lowresCuIdx] = lowresCuCost;
                    sum += lowresCuCost;
                }
                curFrame->m_encData->m_rowStat[row].satdForVbv += sum;
            }
#if DEBUG_FRMCOST_RECALC&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
        fprintf(fp, "row=%d row_satd=%d\n", row, curFrame->m_encData->m_rowStat[row].satdForVbv);
        fclose(fp);
      }
#endif
        }

#if DEBUG_FRMCOST_RECALC&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* pf = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
      int i_mb_xy = 0;
      //int m_widthInCU = ((m_param->sourceWidth / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
      //int m_heightInCU = ((m_param->sourceHeight / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
      for(int j = 0; j < m_heightInCU; j++) {
        for(int i = 0; i < m_widthInCU; i++) {
          int i_mb_cost = curFrame->m_lowres.lowresCostForRc[i_mb_xy] & LOWRES_COST_MASK;
          fprintf(pf, "%6d ", i_mb_cost);
          i_mb_xy++;
        }
        fprintf(pf, "\n");
        fflush(pf);
      }
      fclose(pf);

    }
#endif
    }
}

/* called by API thread or worker thread with inputQueueLock acquired */
void Lookahead::slicetypeDecide()
{
    ProfileScopeEvent(slicetypeDecideEV);

    Lowres *frames[X265_LOOKAHEAD_MAX];
    Frame *list[X265_LOOKAHEAD_MAX];
    int maxSearch = X265_MIN(m_param->lookaheadDepth, X265_LOOKAHEAD_MAX);

    memset(frames, 0, sizeof(frames));
    memset(list, 0, sizeof(list));
    {
        Frame *curFrame = m_inputQueue.first();
        int j;
        for (j = 0; j < m_param->bframes + 2; j++)
        {
            if (!curFrame) break;
            list[j] = curFrame;
            curFrame = curFrame->m_next;
        }

        curFrame = m_inputQueue.first();
        frames[0] = m_lastNonB;
        for (j = 0; j < maxSearch; j++)
        {
            if (!curFrame) break;
            frames[j + 1] = &curFrame->m_lowres;
            curFrame = curFrame->m_next;
        }

        maxSearch = j;
    }

    m_inputQueueLock.release();

    if (!m_est.m_rows && list[0])
        m_est.init(m_param, list[0]);

    if (m_lastNonB && !m_param->rc.bStatRead &&
        ((m_param->bFrameAdaptive && m_param->bframes) ||
         m_param->rc.cuTree || m_param->scenecutThreshold ||
         (m_param->lookaheadDepth && m_param->rc.vbvBufferSize)))
    {
        slicetypeAnalyse(frames, false);
    }

    int bframes, brefs;
    for (bframes = 0, brefs = 0;; bframes++)
    {
        Lowres& frm = list[bframes]->m_lowres;

        if (frm.sliceType == X265_TYPE_BREF && !m_param->bBPyramid && brefs == m_param->bBPyramid)
        {
            frm.sliceType = X265_TYPE_B;
            x265_log(m_param, X265_LOG_WARNING, "B-ref at frame %d incompatible with B-pyramid\n",
                     frm.frameNum);
        }

        /* pyramid with multiple B-refs needs a big enough dpb that the preceding P-frame stays available.
           smaller dpb could be supported by smart enough use of mmco, but it's easier just to forbid it.*/
        else if (frm.sliceType == X265_TYPE_BREF && m_param->bBPyramid && brefs &&
                 m_param->maxNumReferences <= (brefs + 3))
        {
            frm.sliceType = X265_TYPE_B;
            x265_log(m_param, X265_LOG_WARNING, "B-ref at frame %d incompatible with B-pyramid and %d reference frames\n",
                     frm.sliceType, m_param->maxNumReferences);
        }

        if ( /*(!param->intraRefresh || frm.frameNum == 0) && */ frm.frameNum - m_lastKeyframe >= m_param->keyframeMax)
        {
            if (frm.sliceType == X265_TYPE_AUTO || frm.sliceType == X265_TYPE_I)
                frm.sliceType = m_param->bOpenGOP && m_lastKeyframe >= 0 ? X265_TYPE_I : X265_TYPE_IDR;
            bool warn = frm.sliceType != X265_TYPE_IDR;
            if (warn && m_param->bOpenGOP)
                warn &= frm.sliceType != X265_TYPE_I;
            if (warn)
            {
                x265_log(m_param, X265_LOG_WARNING, "specified frame type (%d) at %d is not compatible with keyframe interval\n",
                         frm.sliceType, frm.frameNum);
                frm.sliceType = m_param->bOpenGOP && m_lastKeyframe >= 0 ? X265_TYPE_I : X265_TYPE_IDR;
            }
        }
        if (frm.sliceType == X265_TYPE_I && frm.frameNum - m_lastKeyframe >= m_param->keyframeMin)
        {
            if (m_param->bOpenGOP)
            {
                m_lastKeyframe = frm.frameNum;
                frm.bKeyframe = true;
            }
            else
                frm.sliceType = X265_TYPE_IDR;
        }
        if (frm.sliceType == X265_TYPE_IDR)
        {
            /* Closed GOP */
            m_lastKeyframe = frm.frameNum;
            frm.bKeyframe = true;
            if (bframes > 0)
            {
                list[bframes - 1]->m_lowres.sliceType = X265_TYPE_P;
                bframes--;
            }
        }
        if (bframes == m_param->bframes || !list[bframes + 1])
        {
            if (IS_X265_TYPE_B(frm.sliceType))
                x265_log(m_param, X265_LOG_WARNING, "specified frame type is not compatible with max B-frames\n");
            if (frm.sliceType == X265_TYPE_AUTO || IS_X265_TYPE_B(frm.sliceType))
                frm.sliceType = X265_TYPE_P;
        }
        if (frm.sliceType == X265_TYPE_BREF)
            brefs++;
        if (frm.sliceType == X265_TYPE_AUTO)
            frm.sliceType = X265_TYPE_B;
        else if (!IS_X265_TYPE_B(frm.sliceType))
            break;
    }

    if (bframes)
        list[bframes - 1]->m_lowres.bLastMiniGopBFrame = true;
    list[bframes]->m_lowres.leadingBframes = bframes;
    m_lastNonB = &list[bframes]->m_lowres;
    m_histogram[bframes]++;

    /* insert a bref into the sequence */
    if (m_param->bBPyramid && bframes > 1 && !brefs)
    {
        list[bframes / 2]->m_lowres.sliceType = X265_TYPE_BREF;
        brefs++;
    }
    /* calculate the frame costs ahead of time for estimateFrameCost while we still have lowres */
    if (m_param->rc.rateControlMode != X265_RC_CQP)
    {
        int p0, p1, b;
        /* For zero latency tuning, calculate frame cost to be used later in RC */
        if (!maxSearch)
        {
            for (int i = 0; i <= bframes; i++)
               frames[i + 1] = &list[i]->m_lowres;
        }

        /* estimate new non-B cost */
        p1 = b = bframes + 1;
        p0 = (IS_X265_TYPE_I(frames[bframes + 1]->sliceType)) ? b : 0;
        m_est.estimateFrameCost(frames, p0, p1, b, 0);

        if (bframes)
        {
            p0 = 0; // last nonb
            for (b = 1; b <= bframes; b++)
            {
                if (frames[b]->sliceType == X265_TYPE_B)
                    for (p1 = b; frames[p1]->sliceType == X265_TYPE_B; p1++)
                        ; // find new nonb or bref
                else
                    p1 = bframes + 1;

                m_est.estimateFrameCost(frames, p0, p1, b, 0);

                if (frames[b]->sliceType == X265_TYPE_BREF)
                    p0 = b;
            }
        }
    }

    m_inputQueueLock.acquire();

    /* dequeue all frames from inputQueue that are about to be enqueued
     * in the output queue. The order is important because Frame can
     * only be in one list at a time */
    int64_t pts[X265_BFRAME_MAX + 1];
    for (int i = 0; i <= bframes; i++)
    {
        Frame *curFrame;
        curFrame = m_inputQueue.popFront();
        pts[i] = curFrame->m_pts;
        maxSearch--;
    }

    m_inputQueueLock.release();

    m_outputQueueLock.acquire();
    /* add non-B to output queue */
    int idx = 0;
    list[bframes]->m_reorderedPts = pts[idx++];
    m_outputQueue.pushBack(*list[bframes]);

    /* Add B-ref frame next to P frame in output queue, the B-ref encode before non B-ref frame */
    if (bframes > 1 && m_param->bBPyramid)
    {
        for (int i = 0; i < bframes; i++)
        {
            if (list[i]->m_lowres.sliceType == X265_TYPE_BREF)
            {
                list[i]->m_reorderedPts = pts[idx++];
                m_outputQueue.pushBack(*list[i]);
            }
        }
    }

    /* add B frames to output queue */
    for (int i = 0; i < bframes; i++)
    {
        /* push all the B frames into output queue except B-ref, which already pushed into output queue*/
        if (list[i]->m_lowres.sliceType != X265_TYPE_BREF)
        {
            list[i]->m_reorderedPts = pts[idx++];
            m_outputQueue.pushBack(*list[i]);
        }
    }

    bool isKeyFrameAnalyse = (m_param->rc.cuTree || (m_param->rc.vbvBufferSize && m_param->lookaheadDepth)) && !m_param->rc.bStatRead;
    if (isKeyFrameAnalyse && IS_X265_TYPE_I(m_lastNonB->sliceType))
    {
        m_inputQueueLock.acquire();
        Frame *curFrame = m_inputQueue.first();
        frames[0] = m_lastNonB;
        int j;
        for (j = 0; j < maxSearch; j++)
        {
            frames[j + 1] = &curFrame->m_lowres;
            curFrame = curFrame->m_next;
        }

        frames[j + 1] = NULL;
        m_inputQueueLock.release();
        slicetypeAnalyse(frames, true);
    }

    m_outputQueueLock.release();
    m_outputAvailable.trigger();
}
#if !KEEP_AS265_SAME_WITH_X265
void Lookahead::vbvLookahead(Lowres **frames, int numFrames, int keyframe)
{
    int prevNonB = 0, curNonB = 1, idx = 0;
    while (curNonB < numFrames && frames[curNonB]->sliceType == X265_TYPE_B)
        curNonB++;
    int nextNonB = keyframe ? prevNonB : curNonB;
    int nextB = prevNonB + 1;
    int nextBRef = 0;
    int miniGopEnd = keyframe ? prevNonB : curNonB;
    while (curNonB < numFrames + !keyframe)
    {
        /* P/I cost: This shouldn't include the cost of nextNonB */
        if (nextNonB != curNonB)
        {
            int p0 = IS_X265_TYPE_I(frames[curNonB]->sliceType) ? curNonB : prevNonB;
            frames[nextNonB]->plannedSatd[idx] = vbvFrameCost(frames, p0, curNonB, curNonB);
            frames[nextNonB]->plannedType[idx] = frames[curNonB]->sliceType;
            /* Save the nextNonB Cost in each B frame of the current miniGop */
            if (curNonB > miniGopEnd)
            {
                for (int j = nextB; j < miniGopEnd; j++)
                {
                    frames[j]->plannedSatd[frames[j]->indB] = frames[nextNonB]->plannedSatd[idx];
                    frames[j]->plannedType[frames[j]->indB++] = frames[nextNonB]->plannedType[idx];
                
                }
            }
            idx++;
        }
        /* Handle the B-frames: coded order */
        if (m_param->bBPyramid && curNonB - prevNonB > 1)
            nextBRef = (prevNonB + curNonB + 1) / 2;

        for (int i = prevNonB + 1; i < curNonB; i++, idx++)
        {
            int64_t satdCost = 0; int type = X265_TYPE_B;
            if (nextBRef)
            {
                if (i == nextBRef)
                {
                    satdCost = vbvFrameCost(frames, prevNonB, curNonB, nextBRef);
                    type = X265_TYPE_BREF;
                }
                else if (i < nextBRef)
                    satdCost = vbvFrameCost(frames, prevNonB, nextBRef, i);
                else
                    satdCost = vbvFrameCost(frames, nextBRef, curNonB, i);
            }
            else
                satdCost = vbvFrameCost(frames, prevNonB, nextNonB, i);
            frames[nextNonB]->plannedSatd[idx] = satdCost;
            frames[nextNonB]->plannedType[idx] = type;
            /* Save the nextB Cost in each B frame of the current miniGop */

            for (int j = nextB; j < miniGopEnd; j++)
            {
                if (nextBRef && i == nextBRef)
                    break;
                if (j >= i && j !=nextBRef)
                    continue;
                frames[j]->plannedSatd[frames[j]->indB] = satdCost;
                frames[j]->plannedType[frames[j]->indB++] = X265_TYPE_B;
            }
        }
        prevNonB = curNonB;
        curNonB++;
        while (curNonB <= numFrames && frames[curNonB]->sliceType == X265_TYPE_B)
            curNonB++;
    }

    frames[nextNonB]->plannedType[idx] = X265_TYPE_AUTO;
}

#else

void Lookahead::vbvLookahead(Lowres** frames, int numFrames, int keyframe) {
  int prevNonB = 0, curNonB = 1, idx = 0;
  bool isNextNonB = false;

  while(curNonB < numFrames && frames[curNonB]->sliceType == X265_TYPE_B)
    curNonB++;

  int nextNonB = keyframe ? prevNonB : curNonB;
  int nextB = keyframe ? prevNonB + 1 : curNonB + 1;

#if DEBUG_VBV_LOOKAHEAD
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
    fprintf(fp, "last_nonb=%d cur_nonb=%d next_nonb=%d\n", prevNonB, curNonB, nextNonB);
    fclose(fp);
  }
#endif
  while(curNonB < numFrames + !keyframe) {
#if DEBUG_VBV_LOOKAHEAD
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
      fprintf(fp, "while:last_nonb=%d cur_nonb=%d next_nonb=%d\n", prevNonB, curNonB, nextNonB);
      fclose(fp);
    }
#endif
    /* P/I cost: This shouldn't include the cost of nextNonB */
    if(nextNonB != curNonB) {
#if DEBUG_VBV_LOOKAHEAD
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
        fprintf(fp, "if:cur_nonb poc=%lld frames[%d]->i_type=%d\n",
        frames[curNonB]->frameNum, curNonB, frames[curNonB]->sliceType);
        fclose(fp);
      }
#endif
      int p0 = IS_X265_TYPE_I(frames[curNonB]->sliceType) ? curNonB : prevNonB;
      frames[nextNonB]->plannedSatd[idx] = vbvFrameCost(frames, p0, curNonB, curNonB);
      frames[nextNonB]->plannedType[idx] = frames[curNonB]->sliceType;
#if DEBUG_VBV_LOOKAHEAD
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
        fprintf(fp, "if:next_nonb poc=%d idx=%d frames[%d]->i_planned_satd[%d]=%d ->i_planned_type[%d]=%d\n",
                frames[nextNonB]->frameNum, idx,
                nextNonB, idx, frames[nextNonB]->plannedSatd[idx],
                idx, frames[nextNonB]->plannedType[idx]);
        fclose(fp);
      }
#endif
      idx++;
    }
    /* Handle the B-frames: coded order */
    for(int i = prevNonB + 1; i < curNonB; i++, idx++) {
      frames[nextNonB]->plannedSatd[idx] = vbvFrameCost(frames, prevNonB, curNonB, i);
      frames[nextNonB]->plannedType[idx] = X265_TYPE_B;
#if DEBUG_VBV_LOOKAHEAD
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
        fprintf(fp, "for:i=%d poc=%d idx=%d frames[%d]->i_planned_satd[%d]=%d ->i_planned_type[%d]=%d\n",
                i, frames[nextNonB]->frameNum,
                idx, nextNonB, idx, frames[nextNonB]->plannedSatd[idx],
                idx, frames[nextNonB]->plannedType[idx]);
        fclose(fp);
      }
#endif
    }
#if 0
    for(int i = nextB; i <= curNonB; i++) {
      for(int j = frames[i]->indB + i + 1; j <= curNonB; j++, frames[i]->indB++) {
        if(j == curNonB) {
          if(isNextNonB) {
            int p0 = IS_X265_TYPE_I(frames[curNonB]->sliceType) ? curNonB : prevNonB;
            frames[i]->plannedSatd[frames[i]->indB] = vbvFrameCost(frames, p0, curNonB, curNonB);
            frames[i]->plannedType[frames[i]->indB] = frames[curNonB]->sliceType;
          }
        } else {
          frames[i]->plannedSatd[frames[i]->indB] = vbvFrameCost(frames, prevNonB, curNonB, j);
          frames[i]->plannedType[frames[i]->indB] = X265_TYPE_B;
        }
      }
      if(i == curNonB && !isNextNonB)
        isNextNonB = true;
    }
#endif
    prevNonB = curNonB;
    curNonB++;
    while(curNonB <= numFrames && frames[curNonB]->sliceType == X265_TYPE_B)
      curNonB++;
#if DEBUG_VBV_LOOKAHEAD
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
      fprintf(fp, "last_nonb=%d cur_nonb=%d next_nonb=%d\n", prevNonB, curNonB, nextNonB);
      fclose(fp);
    }
#endif
  }

  frames[nextNonB]->plannedType[idx] = X265_TYPE_AUTO;
#if DEBUG_VBV_LOOKAHEAD
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
    fprintf(fp, "last_nonb=%d cur_nonb=%d next_nonb=%d idx=%d\n\n", prevNonB, curNonB, nextNonB, idx);
    fclose(fp);
  }
#endif
}
#endif
int64_t Lookahead::vbvFrameCost(Lowres **frames, int p0, int p1, int b)
{
    int64_t cost = m_est.estimateFrameCost(frames, p0, p1, b, 0);

    if (m_param->rc.aqMode)
    {
        if (m_param->rc.cuTree)
            return frameCostRecalculate(frames, p0, p1, b);
        else
            return frames[b]->costEstAq[b - p0][p1 - b];
    }
    return cost;
}

void Lookahead::slicetypeAnalyse(Lowres **frames, bool bKeyframe)
{
    int numFrames, origNumFrames, keyintLimit, framecnt;
    int maxSearch = X265_MIN(m_param->lookaheadDepth, X265_LOOKAHEAD_MAX);
    int cuCount = NUM_CUS;
    int resetStart;
    bool bIsVbvLookahead = m_param->rc.vbvBufferSize && m_param->lookaheadDepth;

#if DEBUG_VBV_LOOKAHEAD&&KEEP_AS265_SAME_WITH_X265
  {
    static int kk = 0;
    FILE* fp;
    if(kk == 0) {
      fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "w");
      kk = 1;
      fclose(fp);
    }
  }
#endif
    /* count undecided frames */
    for (framecnt = 0; framecnt < maxSearch; framecnt++)
    {
        Lowres *fenc = frames[framecnt + 1];
        if (!fenc || fenc->sliceType != X265_TYPE_AUTO)
            break;
    }

    if (!framecnt)
    {
        if (m_param->rc.cuTree)
            cuTree(frames, 0, bKeyframe);
        return;
    }

    frames[framecnt + 1] = NULL;

    keyintLimit = m_param->keyframeMax - frames[0]->frameNum + m_lastKeyframe - 1;
    origNumFrames = numFrames = X265_MIN(framecnt, keyintLimit);

    if (bIsVbvLookahead)
        numFrames = framecnt;
    else if (m_param->bOpenGOP && numFrames < framecnt)
        numFrames++;
    else if (numFrames == 0)
    {
        frames[1]->sliceType = X265_TYPE_I;
        return;
    }

    int numBFrames = 0;
    int numAnalyzed = numFrames;
    if (m_param->scenecutThreshold && scenecut(frames, 0, 1, true, origNumFrames, maxSearch))
    {
        frames[1]->sliceType = X265_TYPE_I;
        return;
    }

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
  {
    static int kk = 0;
    FILE* fp;
    if(kk == 0) {
      fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "w");
      kk = 1;
      fclose(fp);
    }
  }
#endif
    if (m_param->bframes)
    {
        if (m_param->bFrameAdaptive == X265_B_ADAPT_TRELLIS)
        {
            if (numFrames > 1)
            {
                char best_paths[X265_BFRAME_MAX + 1][X265_LOOKAHEAD_MAX + 1] = { "", "P" };
                int best_path_index = numFrames % (X265_BFRAME_MAX + 1);

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
          fprintf(fp, "num_frames=%d best_path_index=%d\n", numFrames, best_path_index);
          fclose(fp);
        }
#endif
                /* Perform the frametype analysis. */
                for (int j = 2; j <= numFrames; j++)
                {
                    slicetypePath(frames, j, best_paths);
                }

                numBFrames = (int)strspn(best_paths[best_path_index], "B");

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
          fprintf(fp, "num_bframes=%d\n", numBFrames);
          for(int i = 1; i < numFrames; i++) {
            fprintf(fp, "Before:i=%d poc=%d i_type=%d\n", i, frames[i]->frameNum, frames[i]->sliceType);
          }
          fclose(fp);
        }
#endif
                /* Load the results of the analysis into the frame types. */
                for (int j = 1; j < numFrames; j++)
                {
                    frames[j]->sliceType = best_paths[best_path_index][j - 1] == 'B' ? X265_TYPE_B : X265_TYPE_P;
                }
#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
          for(int i = 1; i < numFrames; i++) {
            fprintf(fp, "After:i=%d poc=%d i_type=%d\n", i, frames[i]->frameNum, frames[i]->sliceType);
          }
          fclose(fp);
        }
#endif
            }
            frames[numFrames]->sliceType = X265_TYPE_P;
        }
        else if (m_param->bFrameAdaptive == X265_B_ADAPT_FAST)
        {
            int64_t cost1p0, cost2p0, cost1b1, cost2p1;

            for (int i = 0; i <= numFrames - 2; )
            {
                cost2p1 = m_est.estimateFrameCost(frames, i + 0, i + 2, i + 2, 1);
#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
          fprintf(fp, "i=%d cost2p1=%d intra_mbs=%d mb_count=%d\n", i, cost2p1, frames[i + 2]->intraMbs[2], cuCount);
          fclose(fp);
        }
#endif
                if (frames[i + 2]->intraMbs[2] > cuCount / 2)
                {
                    frames[i + 1]->sliceType = X265_TYPE_P;
                    frames[i + 2]->sliceType = X265_TYPE_P;
                    i += 2;
                    continue;
                }

                cost1b1 = m_est.estimateFrameCost(frames, i + 0, i + 2, i + 1, 0);
                cost1p0 = m_est.estimateFrameCost(frames, i + 0, i + 1, i + 1, 0);
                cost2p0 = m_est.estimateFrameCost(frames, i + 1, i + 2, i + 2, 0);

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
          fprintf(fp, "cost1b1=%d cost1p0=%d cost2p0=%d cost2p1=%d\n",
                  cost1b1, cost1p0, cost2p0, cost2p1);
          fclose(fp);
        }
#endif
                if (cost1p0 + cost2p0 < cost1b1 + cost2p1)
                {
                    frames[i + 1]->sliceType = X265_TYPE_P;
                    i += 1;
                    continue;
                }

// arbitrary and untuned
#define INTER_THRESH 300
#define P_SENS_BIAS (50 - m_param->bFrameBias)
                frames[i + 1]->sliceType = X265_TYPE_B;

                int j;
                for (j = i + 2; j <= X265_MIN(i + m_param->bframes, numFrames - 1); j++)
                {
                    int64_t pthresh = X265_MAX(INTER_THRESH - P_SENS_BIAS * (j - i - 1), INTER_THRESH / 10);
                    int64_t pcost = m_est.estimateFrameCost(frames, i + 0, j + 1, j + 1, 1);
#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
          {
            FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
            fprintf(fp, "j=%d pthresh=%lld pcost=%lld intra_mbs=%d\n", j, pthresh, pcost, frames[j + 1]->intraMbs[j - i + 1]);
            fclose(fp);
          }
#endif
                    if (pcost > pthresh * cuCount || frames[j + 1]->intraMbs[j - i + 1] > cuCount / 3)
                        break;
                    frames[j]->sliceType = X265_TYPE_B;
                }

                frames[j]->sliceType = X265_TYPE_P;
                i = j;
            }
            frames[numFrames]->sliceType = X265_TYPE_P;
            numBFrames = 0;
            while (numBFrames < numFrames && frames[numBFrames + 1]->sliceType == X265_TYPE_B)
            {
                numBFrames++;
            }
#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
        fprintf(fp, "num_bframes=%d\n", numBFrames);
        fclose(fp);
      }
#endif
        }
        else
        {
            numBFrames = X265_MIN(numFrames - 1, m_param->bframes);
            for (int j = 1; j < numFrames; j++)
            {
                frames[j]->sliceType = (j % (numBFrames + 1)) ? X265_TYPE_B : X265_TYPE_P;
            }

            frames[numFrames]->sliceType = X265_TYPE_P;
        }
        /* Check scenecut on the first minigop. */
        for (int j = 1; j < numBFrames + 1; j++)
        {
            if (m_param->scenecutThreshold && scenecut(frames, j, j + 1, false, origNumFrames, maxSearch))
            {
                frames[j]->sliceType = X265_TYPE_P;
                numAnalyzed = j;
                break;
            }
        }

        resetStart = bKeyframe ? 1 : X265_MIN(numBFrames + 2, numAnalyzed + 1);
    }
    else
    {
        for (int j = 1; j <= numFrames; j++)
        {
            frames[j]->sliceType = X265_TYPE_P;
        }

        resetStart = bKeyframe ? 1 : 2;
    }

#if DEBUG_VBV_LOOKAHEAD&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
    fprintf(fp, "num_frames=%d num_bframes=%d reset_start=%d\n", numFrames, numBFrames, resetStart);
    for(int i = 0; i <= numFrames; i++) {
      fprintf(fp, "i=%d poc=%d i_type=%d\n", i, frames[i]->frameNum, frames[i]->sliceType);
    }
    fclose(fp);
  }
#endif
    if (m_param->rc.cuTree)
        cuTree(frames, X265_MIN(numFrames, m_param->keyframeMax), bKeyframe);

    // if (!param->bIntraRefresh)
    for (int j = keyintLimit + 1; j <= numFrames; j += m_param->keyframeMax)
    {
        frames[j]->sliceType = X265_TYPE_I;
        resetStart = X265_MIN(resetStart, j + 1);
    }

    if (bIsVbvLookahead)
        vbvLookahead(frames, numFrames, bKeyframe);

    /* Restore frametypes for all frames that haven't actually been decided yet. */
    for (int j = resetStart; j <= numFrames; j++)
    {
        frames[j]->sliceType = X265_TYPE_AUTO;
    }
#if DEBUG_VBV_LOOKAHEAD&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_VBV_LOOKAHEAD), "a");
    fprintf(fp, "final:num_frames=%d num_bframes=%d reset_start=%d\n", numFrames, numBFrames, resetStart);
    for(int i = 0; i <= numFrames; i++) {
      fprintf(fp, "i=%d poc=%d i_type=%d\n", i, frames[i]->frameNum, frames[i]->sliceType);
    }
    fclose(fp);
  }
#endif
}

bool Lookahead::scenecut(Lowres **frames, int p0, int p1, bool bRealScenecut, int numFrames, int maxSearch)
{
    /* Only do analysis during a normal scenecut check. */
    if (bRealScenecut && m_param->bframes)
    {
        int origmaxp1 = p0 + 1;
        /* Look ahead to avoid coding short flashes as scenecuts. */
        if (m_param->bFrameAdaptive == X265_B_ADAPT_TRELLIS)
            /* Don't analyse any more frames than the trellis would have covered. */
            origmaxp1 += m_param->bframes;
        else
            origmaxp1++;
        int maxp1 = X265_MIN(origmaxp1, numFrames);

        /* Where A and B are scenes: AAAAAABBBAAAAAA
         * If BBB is shorter than (maxp1-p0), it is detected as a flash
         * and not considered a scenecut. */
        for (int cp1 = p1; cp1 <= maxp1; cp1++)
        {
            if (!scenecutInternal(frames, p0, cp1, false))
                /* Any frame in between p0 and cur_p1 cannot be a real scenecut. */
                for (int i = cp1; i > p0; i--)
                {
                    frames[i]->bScenecut = false;
                }
        }

        /* Where A-F are scenes: AAAAABBCCDDEEFFFFFF
         * If each of BB ... EE are shorter than (maxp1-p0), they are
         * detected as flashes and not considered scenecuts.
         * Instead, the first F frame becomes a scenecut.
         * If the video ends before F, no frame becomes a scenecut. */
        for (int cp0 = p0; cp0 <= maxp1; cp0++)
        {
            if (origmaxp1 > maxSearch || (cp0 < maxp1 && scenecutInternal(frames, cp0, maxp1, false)))
                /* If cur_p0 is the p0 of a scenecut, it cannot be the p1 of a scenecut. */
                frames[cp0]->bScenecut = false;
        }
    }

    /* Ignore frames that are part of a flash, i.e. cannot be real scenecuts. */
    if (!frames[p1]->bScenecut)
        return false;
    return scenecutInternal(frames, p0, p1, bRealScenecut);
}

bool Lookahead::scenecutInternal(Lowres **frames, int p0, int p1, bool bRealScenecut)
{
    Lowres *frame = frames[p1];

    m_est.estimateFrameCost(frames, p0, p1, p1, 0);

    int64_t icost = frame->costEst[0][0];
    int64_t pcost = frame->costEst[p1 - p0][0];
    int gopSize = frame->frameNum - m_lastKeyframe;
    float threshMax = (float)(m_param->scenecutThreshold / 100.0);

    /* magic numbers pulled out of thin air */
    float threshMin = (float)(threshMax * 0.25);
    float bias;

    if (m_param->keyframeMin == m_param->keyframeMax)
        threshMin = threshMax;
    if (gopSize <= m_param->keyframeMin / 4)
        bias = threshMin / 4;
    else if (gopSize <= m_param->keyframeMin)
        bias = threshMin * gopSize / m_param->keyframeMin;
    else
    {
        bias = threshMin
            + (threshMax - threshMin)
            * (gopSize - m_param->keyframeMin)
            / (m_param->keyframeMax - m_param->keyframeMin);
    }

    bool res = pcost >= (1.0 - bias) * icost;
    if (res && bRealScenecut)
    {
        int imb = frame->intraMbs[p1 - p0];
        int pmb = NUM_CUS - imb;
        x265_log(m_param, X265_LOG_DEBUG, "scene cut at %d Icost:%d Pcost:%d ratio:%.4f bias:%.4f gop:%d (imb:%d pmb:%d)\n",
                 frame->frameNum, icost, pcost, 1. - (double)pcost / icost, bias, gopSize, imb, pmb);
    }
    return res;
}

void Lookahead::slicetypePath(Lowres **frames, int length, char(*best_paths)[X265_LOOKAHEAD_MAX + 1])
{
    char paths[2][X265_LOOKAHEAD_MAX + 1];
    int num_paths = X265_MIN(m_param->bframes + 1, length);
    int64_t best_cost = 1LL << 62;
#if KEEP_AS265_SAME_WITH_X265
  best_cost = 1LL << 28;
#endif
    int idx = 0;

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
    fprintf(fp, "length=%d num_paths=%d best_cost=%lld\n", length, num_paths, best_cost);
    fclose(fp);
  }
#endif
    /* Iterate over all currently possible paths */
    for (int path = 0; path < num_paths; path++)
    {
        /* Add suffixes to the current path */
        int len = length - (path + 1);
        memcpy(paths[idx], best_paths[len % (X265_BFRAME_MAX + 1)], len);
        memset(paths[idx] + len, 'B', path);
        strcpy(paths[idx] + len + path, "P");

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
      fprintf(fp, "path=%d:idx=%d path=%s\n", path, idx, paths[idx]);
      fclose(fp);
    }
#endif
        /* Calculate the actual cost of the current path */
        int64_t cost = slicetypePathCost(frames, paths[idx], best_cost);
        if (cost < best_cost)
        {
            best_cost = cost;
            idx ^= 1;
        }

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
      fprintf(fp, "cost=%lld best_cost=%lld idx=%d\n", cost, best_cost, idx);
      fclose(fp);
    }
#endif
    }

    /* Store the best path. */
    memcpy(best_paths[length % (X265_BFRAME_MAX + 1)], paths[idx ^ 1], length);

#if DEBUG_B_FRAME_ADAPTIVE&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_B_FRAME_ADAPTIVE), "a");
    fprintf(fp, "best_paths=%s\n", best_paths[length % (X265_BFRAME_MAX + 1)]);
    fclose(fp);
  }
#endif
}

int64_t Lookahead::slicetypePathCost(Lowres **frames, char *path, int64_t threshold)
{
    int64_t cost = 0;
    int loc = 1;
    int cur_p = 0;

    path--; /* Since the 1st path element is really the second frame */
    while (path[loc])
    {
        int next_p = loc;
        /* Find the location of the next P-frame. */
        while (path[next_p] != 'P')
        {
            next_p++;
        }

        /* Add the cost of the P-frame found above */
        cost += m_est.estimateFrameCost(frames, cur_p, next_p, next_p, 0);
        /* Early terminate if the cost we have found is larger than the best path cost so far */
        if (cost > threshold)
            break;

        if (m_param->bBPyramid && next_p - cur_p > 2)
        {
            int middle = cur_p + (next_p - cur_p) / 2;
            cost += m_est.estimateFrameCost(frames, cur_p, next_p, middle, 0);
            for (int next_b = loc; next_b < middle && cost < threshold; next_b++)
            {
                cost += m_est.estimateFrameCost(frames, cur_p, middle, next_b, 0);
            }

            for (int next_b = middle + 1; next_b < next_p && cost < threshold; next_b++)
            {
                cost += m_est.estimateFrameCost(frames, middle, next_p, next_b, 0);
            }
        }
        else
        {
            for (int next_b = loc; next_b < next_p && cost < threshold; next_b++)
            {
                cost += m_est.estimateFrameCost(frames, cur_p, next_p, next_b, 0);
            }
        }

        loc = next_p + 1;
        cur_p = next_p;
    }

    return cost;
}

void Lookahead::cuTree(Lowres **frames, int numframes, bool bIntra)
{
    int idx = !bIntra;
    int lastnonb, curnonb = 1;
    int bframes = 0;

    x265_emms();
    double totalDuration = 0.0;
    for (int j = 0; j <= numframes; j++)
        totalDuration += (double)m_param->fpsDenom / m_param->fpsNum;

    double averageDuration = totalDuration / (numframes + 1);

    int i = numframes;
    int cuCount = m_widthInCU * m_heightInCU;

    if (bIntra)
        m_est.estimateFrameCost(frames, 0, 0, 0, 0);

    while (i > 0 && frames[i]->sliceType == X265_TYPE_B)
        i--;

    lastnonb = i;

    /* Lookaheadless MB-tree is not a theoretically distinct case; the same extrapolation could
     * be applied to the end of a lookahead buffer of any size.  However, it's most needed when
     * lookahead=0, so that's what's currently implemented. */
    if (!m_param->lookaheadDepth)
    {
        if (bIntra)
        {
            memset(frames[0]->propagateCost, 0, cuCount * sizeof(uint16_t));
            memcpy(frames[0]->qpCuTreeOffset, frames[0]->qpAqOffset, cuCount * sizeof(double));
            return;
        }
        std::swap(frames[lastnonb]->propagateCost, frames[0]->propagateCost);
        memset(frames[0]->propagateCost, 0, cuCount * sizeof(uint16_t));
    }
    else
    {
        if (lastnonb < idx)
            return;
        memset(frames[lastnonb]->propagateCost, 0, cuCount * sizeof(uint16_t));
    }

    while (i-- > idx)
    {
        curnonb = i;
        while (frames[curnonb]->sliceType == X265_TYPE_B && curnonb > 0)
            curnonb--;

        if (curnonb < idx)
            break;

        m_est.estimateFrameCost(frames, curnonb, lastnonb, lastnonb, 0);
        memset(frames[curnonb]->propagateCost, 0, cuCount * sizeof(uint16_t));
        bframes = lastnonb - curnonb - 1;
        if (m_param->bBPyramid && bframes > 1)
        {
            int middle = (bframes + 1) / 2 + curnonb;
            m_est.estimateFrameCost(frames, curnonb, lastnonb, middle, 0);
            memset(frames[middle]->propagateCost, 0, cuCount * sizeof(uint16_t));
            while (i > curnonb)
            {
                int p0 = i > middle ? middle : curnonb;
                int p1 = i < middle ? middle : lastnonb;
                if (i != middle)
                {
                    m_est.estimateFrameCost(frames, p0, p1, i, 0);
                    estimateCUPropagate(frames, averageDuration, p0, p1, i, 0);
                }
                i--;
            }

            estimateCUPropagate(frames, averageDuration, curnonb, lastnonb, middle, 1);
        }
        else
        {
            while (i > curnonb)
            {
                m_est.estimateFrameCost(frames, curnonb, lastnonb, i, 0);
                estimateCUPropagate(frames, averageDuration, curnonb, lastnonb, i, 0);
                i--;
            }
        }
        estimateCUPropagate(frames, averageDuration, curnonb, lastnonb, lastnonb, 1);
        lastnonb = curnonb;
    }

    if (!m_param->lookaheadDepth)
    {
        m_est.estimateFrameCost(frames, 0, lastnonb, lastnonb, 0);
        estimateCUPropagate(frames, averageDuration, 0, lastnonb, lastnonb, 1);
        std::swap(frames[lastnonb]->propagateCost, frames[0]->propagateCost);
    }

    cuTreeFinish(frames[lastnonb], averageDuration, lastnonb);
    if (m_param->bBPyramid && bframes > 1 && !m_param->rc.vbvBufferSize)
        cuTreeFinish(frames[lastnonb + (bframes + 1) / 2], averageDuration, 0);
}

void Lookahead::estimateCUPropagate(Lowres **frames, double averageDuration, int p0, int p1, int b, int referenced)
{
    uint16_t *refCosts[2] = { frames[p0]->propagateCost, frames[p1]->propagateCost };
    int32_t distScaleFactor = (((b - p0) << 8) + ((p1 - p0) >> 1)) / (p1 - p0);
    int32_t bipredWeight = m_param->bEnableWeightedBiPred ? 64 - (distScaleFactor >> 2) : 32;
    MV *mvs[2] = { frames[b]->lowresMvs[0][b - p0 - 1], frames[b]->lowresMvs[1][p1 - b - 1] };
    int32_t bipredWeights[2] = { bipredWeight, 64 - bipredWeight };

    memset(m_scratch, 0, m_widthInCU * sizeof(int));

    uint16_t *propagateCost = frames[b]->propagateCost;

    x265_emms();
    double fpsFactor = CLIP_DURATION((double)m_param->fpsDenom / m_param->fpsNum) / CLIP_DURATION(averageDuration);

    /* For non-refferd frames the source costs are always zero, so just memset one row and re-use it. */
    if (!referenced)
        memset(frames[b]->propagateCost, 0, m_widthInCU * sizeof(uint16_t));

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
    fprintf(fp, "dist_scale_factor=%d bipred_weight=%d widthInCU=%d heightInCU=%d fps_factor="FLOAT_FORMAT"\n",
            distScaleFactor, bipredWeight, m_widthInCU, m_heightInCU, fpsFactor
           );
    fclose(fp);
  }
#endif
    int32_t StrideInCU = m_widthInCU;
    for (uint16_t blocky = 0; blocky < m_heightInCU; blocky++)
    {
        int cuIndex = blocky * StrideInCU;
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
      fprintf(fp, "cuIndex=%d\n", cuIndex);
      fclose(fp);
    }
#endif
        primitives.propagateCost(m_scratch, propagateCost,
                                 frames[b]->intraCost + cuIndex, frames[b]->lowresCosts[b - p0][p1 - b] + cuIndex,
                                 frames[b]->invQscaleFactor + cuIndex, &fpsFactor, m_widthInCU);

        if (referenced)
            propagateCost += m_widthInCU;
        for (uint16_t blockx = 0; blockx < m_widthInCU; blockx++, cuIndex++)
        {
            int32_t propagate_amount = m_scratch[blockx];
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
        fprintf(fp, "mb_x=%d cuIndex=%d propagate_amount=%d\n", blockx, cuIndex, propagate_amount);
        fclose(fp);
      }
#endif
            /* Don't propagate for an intra block. */
            if (propagate_amount > 0)
            {
                /* Access width-2 bitfield. */
                int32_t lists_used = frames[b]->lowresCosts[b - p0][p1 - b][cuIndex] >> LOWRES_COST_SHIFT;

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
          fprintf(fp, "lists_used=%d\n", lists_used);
          fclose(fp);
        }
#endif
                /* Follow the MVs to the previous frame(s). */
                for (uint16_t list = 0; list < 2; list++)
                {
                    if ((lists_used >> list) & 1)
                    {
#define CLIP_ADD(s, x) (s) = (uint16_t)X265_MIN((s) + (x), (1 << 16) - 1)
                        int32_t listamount = propagate_amount;

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
            {
              FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
              fprintf(fp, "list=%d lists_used=%d listamount=%d\n", list, lists_used, listamount);
              fclose(fp);
            }
#endif
                        /* Apply bipred weighting. */
                        if (lists_used == 3)
                            listamount = (listamount * bipredWeights[list] + 32) >> 6;

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
            {
              FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
              fprintf(fp, "listamount=%d mv=[%d,%d] ref_costs=%d\n", listamount, mvs[list][cuIndex].x, mvs[list][cuIndex].y, refCosts[list][cuIndex]);
              fclose(fp);
            }
#endif
                        /* Early termination for simple case of mv0. */
                        if (!mvs[list][cuIndex].word)
                        {
                            CLIP_ADD(refCosts[list][cuIndex], listamount);
                            continue;
                        }

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
            {
              FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
              fprintf(fp, "ref_costs=%d\n", refCosts[list][cuIndex]);
              fclose(fp);
            }
#endif
                        int32_t x = mvs[list][cuIndex].x;
                        int32_t y = mvs[list][cuIndex].y;
                        int32_t cux = (x >> 5) + blockx;
                        int32_t cuy = (y >> 5) + blocky;
                        int32_t idx0 = cux + cuy * StrideInCU;
                        int32_t idx1 = idx0 + 1;
                        int32_t idx2 = idx0 + StrideInCU;
                        int32_t idx3 = idx0 + StrideInCU + 1;

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
            {
              FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
              fprintf(fp, "x=%d y=%d mbx=%d mby=%d idx0=%d idx1=%d idx2=%d idx3=%d\n",
                      x, y, cux, cuy, idx0, idx1, idx2, idx3);
              fclose(fp);
            }
#endif
                        x &= 31;
                        y &= 31;
                        int32_t idx0weight = (32 - y) * (32 - x);
                        int32_t idx1weight = (32 - y) * x;
                        int32_t idx2weight = y * (32 - x);
                        int32_t idx3weight = y * x;

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
            {
              FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
              fprintf(fp, "x=%d y=%d idx0w=%d idx1w=%d idx2w=%d idx3w=%d\n",
                      x, y, idx0weight, idx1weight, idx2weight, idx3weight);
              fclose(fp);
            }
#endif
                        /* We could just clip the MVs, but pixels that lie outside the frame probably shouldn't
                         * be counted. */
                        if (cux < m_widthInCU - 1 && cuy < m_heightInCU - 1 && cux >= 0 && cuy >= 0)
                        {
                            CLIP_ADD(refCosts[list][idx0], (listamount * idx0weight + 512) >> 10);
                            CLIP_ADD(refCosts[list][idx1], (listamount * idx1weight + 512) >> 10);
                            CLIP_ADD(refCosts[list][idx2], (listamount * idx2weight + 512) >> 10);
                            CLIP_ADD(refCosts[list][idx3], (listamount * idx3weight + 512) >> 10);
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
              {
                FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
                fprintf(fp, "ref_costs0=%d ref_costs1=%d ref_costs2=%d ref_costs3=%d\n",
                        refCosts[list][idx0], refCosts[list][idx1], refCosts[list][idx2], refCosts[list][idx3]);
                fclose(fp);
              }
#endif
                        }
                        else /* Check offsets individually */
                        {
                            if (cux < m_widthInCU && cuy < m_heightInCU && cux >= 0 && cuy >= 0){
                                CLIP_ADD(refCosts[list][idx0], (listamount * idx0weight + 512) >> 10);
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
                {
                  FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
                  fprintf(fp, "ref_costs[%d][%d]=%d\n", list, idx0, refCosts[list][idx0]);
                  fclose(fp);
                }
#endif
              }
                            if (cux + 1 < m_widthInCU && cuy < m_heightInCU && cux + 1 >= 0 && cuy >= 0){
                                CLIP_ADD(refCosts[list][idx1], (listamount * idx1weight + 512) >> 10);
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
                {
                  FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
                  fprintf(fp, "ref_costs[%d][%d]=%d\n", list, idx1, refCosts[list][idx1]);
                  fclose(fp);
                }
#endif
              }
                            if (cux < m_widthInCU && cuy + 1 < m_heightInCU && cux >= 0 && cuy + 1 >= 0){
                                CLIP_ADD(refCosts[list][idx2], (listamount * idx2weight + 512) >> 10);
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
                {
                  FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
                  fprintf(fp, "ref_costs[%d][%d]=%d\n", list, idx2, refCosts[list][idx2]);
                  fclose(fp);
                }
#endif
              }
                            if (cux + 1 < m_widthInCU && cuy + 1 < m_heightInCU && cux + 1 >= 0 && cuy + 1 >= 0){
                                CLIP_ADD(refCosts[list][idx3], (listamount * idx3weight + 512) >> 10);
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
                {
                  FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
                  fprintf(fp, "ref_costs[%d][%d]=%d\n", list, idx3, refCosts[list][idx3]);
                  fclose(fp);
                }
#endif
              }
                        }
                    }
                }
            }
        }
    }

    if (m_param->rc.vbvBufferSize && m_param->lookaheadDepth && referenced)
        cuTreeFinish(frames[b], averageDuration, b == p1 ? b - p0 : 0);
}

void Lookahead::cuTreeFinish(Lowres *frame, double averageDuration, int ref0Distance)
{
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
    fprintf(fp, "average_durattion="FLOAT_FORMAT" m_frameDuration="FLOAT_FORMAT"\n", averageDuration, (double)m_param->fpsDenom / m_param->fpsNum);
    fclose(fp);
  }
#endif
    int fpsFactor = (int)(CLIP_DURATION(averageDuration) / CLIP_DURATION((double)m_param->fpsDenom / m_param->fpsNum) * 256);
    double weightdelta = 0.0;

    if (ref0Distance && frame->weightedCostDelta[ref0Distance - 1] > 0)
        weightdelta = (1.0 - frame->weightedCostDelta[ref0Distance - 1]);

    /* Allow the strength to be adjusted via qcompress, since the two
     * concepts are very similar. */

    int cuCount = m_widthInCU * m_heightInCU;
    double strength = 5.0 * (1.0 - m_param->rc.qCompress);

#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
    fprintf(fp, "fps_factor=%d weightdelta="FLOAT_FORMAT" strength="FLOAT_FORMAT" cuCount=%d\n",
            fpsFactor, weightdelta, strength, cuCount);
    fclose(fp);
  }
#endif
    for (int cuIndex = 0; cuIndex < cuCount; cuIndex++)
    {
        int intracost = (frame->intraCost[cuIndex] * frame->invQscaleFactor[cuIndex] + 128) >> 8;
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
      fprintf(fp, "index=%d intra_cost=%d inv_qscale_factor=%d intra_cost=%d\n",
              cuIndex, frame->intraCost[cuIndex], frame->invQscaleFactor[cuIndex], intracost);
      fclose(fp);
    }
#endif
        if (intracost)
        {
            int propagateCost = (frame->propagateCost[cuIndex] * fpsFactor + 128) >> 8;
            double log2_ratio = X265_LOG2(intracost + propagateCost) - X265_LOG2(intracost) + weightdelta;
            frame->qpCuTreeOffset[cuIndex] = frame->qpAqOffset[cuIndex] - strength * log2_ratio;
#if DEBUG_MBTREE_PROCESS&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_MBTREE_PROCESS), "a");
        fprintf(fp, "propagate_cost=%d log2_ratio="FLOAT_FORMAT" f_qp_offset="FLOAT_FORMAT" f_qp_offset_aq="FLOAT_FORMAT"\n",
                propagateCost, log2_ratio, frame->qpCuTreeOffset[cuIndex], frame->qpAqOffset[cuIndex]
               );
        fclose(fp);
      }
#endif
        }
    }
}

/* If MB-tree changes the quantizers, we need to recalculate the frame cost without
 * re-running lookahead. */
int64_t Lookahead::frameCostRecalculate(Lowres** frames, int p0, int p1, int b)
{
    int64_t score = 0;
    int *rowSatd = frames[b]->rowSatds[b - p0][p1 - b];
    double *qp_offset = (frames[b]->sliceType == X265_TYPE_B) ? frames[b]->qpAqOffset : frames[b]->qpCuTreeOffset;

    x265_emms();
    for (int cuy = m_heightInCU - 1; cuy >= 0; cuy--)
    {
        rowSatd[cuy] = 0;
        for (int cux = m_widthInCU - 1; cux >= 0; cux--)
        {
            int cuxy = cux + cuy * m_widthInCU;
            int cuCost = frames[b]->lowresCosts[b - p0][p1 - b][cuxy] & LOWRES_COST_MASK;
            double qp_adj = qp_offset[cuxy];
#if DEBUG_FRMCOST_RECALC&&DEBUG_FRMCOST_RECALC_REAL&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
        fprintf(fp, "poc=%d cuxy=%d cucost=%d qp_adj="FLOAT_FORMAT"\n",
                frames[b]->frameNum, cuxy, cuCost, qp_offset[cuxy]);
        fclose(fp);
      }
#endif
            cuCost = (cuCost * x265_exp2fix8(qp_adj) + 128) >> 8;
            rowSatd[cuy] += cuCost;
            if ((cuy > 0 && cuy < m_heightInCU - 1 &&
                 cux > 0 && cux < m_widthInCU - 1) ||
                m_widthInCU <= 2 || m_heightInCU <= 2)
            {
                score += cuCost;
            }
#if DEBUG_FRMCOST_RECALC&&DEBUG_FRMCOST_RECALC_REAL&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_FRMCOST_RECALC), "a");
        fprintf(fp, "cucost=%d score=%d\n", cuCost, score);
        fclose(fp);
      }
#endif
        }
    }

    return score;
}

CostEstimate::CostEstimate(ThreadPool *p)
    : WaveFront(p)
{
    m_param = NULL;
    m_curframes = NULL;
    m_wbuffer[0] = m_wbuffer[1] = m_wbuffer[2] = m_wbuffer[3] = 0;
    m_rows = NULL;
    m_paddedLines = m_widthInCU = m_heightInCU = 0;
    m_bDoSearch[0] = m_bDoSearch[1] = false;
    m_curb = m_curp0 = m_curp1 = 0;
    m_bFrameCompleted = false;
}

CostEstimate::~CostEstimate()
{
    for (int i = 0; i < 4; i++)
    {
        x265_free(m_wbuffer[i]);
    }

    delete[] m_rows;
}

void CostEstimate::init(x265_param *_param, Frame *curFrame)
{
    m_param = _param;
    m_widthInCU = ((m_param->sourceWidth / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;
    m_heightInCU = ((m_param->sourceHeight / 2) + X265_LOWRES_CU_SIZE - 1) >> X265_LOWRES_CU_BITS;

    m_rows = new EstimateRow[m_heightInCU];
    for (int i = 0; i < m_heightInCU; i++)
    {
        m_rows[i].m_widthInCU = m_widthInCU;
        m_rows[i].m_heightInCU = m_heightInCU;
        m_rows[i].m_param = m_param;
    }

    if (WaveFront::init(m_heightInCU))
        WaveFront::enableAllRows();
    else
        m_pool = NULL;

    if (m_param->bEnableWeightedPred)
    {
        PicYuv *orig = curFrame->m_fencPic;
        m_paddedLines = curFrame->m_lowres.lines + 2 * orig->m_lumaMarginY;
        intptr_t padoffset = curFrame->m_lowres.lumaStride * orig->m_lumaMarginY + orig->m_lumaMarginX;

        /* allocate weighted lowres buffers */
        for (int i = 0; i < 4; i++)
        {
            m_wbuffer[i] = (pixel*)x265_malloc(sizeof(pixel) * (curFrame->m_lowres.lumaStride * m_paddedLines));
            m_weightedRef.lowresPlane[i] = m_wbuffer[i] + padoffset;
        }

        m_weightedRef.fpelPlane[0] = m_weightedRef.lowresPlane[0];
        m_weightedRef.lumaStride = curFrame->m_lowres.lumaStride;
        m_weightedRef.isLowres = true;
        m_weightedRef.isWeighted = false;
    }
}

int64_t CostEstimate::estimateFrameCost(Lowres **frames, int p0, int p1, int b, bool bIntraPenalty)
{
    int64_t score = 0;
    Lowres *fenc = frames[b];

#if DEBUG_FRAME_COST_OUTPUT&&KEEP_AS265_SAME_WITH_X265
  {
    static int kk = 0;
    FILE* fp;
    if(kk == 0) {
      fp = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "w");
      kk = 1;
      fclose(fp);
    }
  }
#endif
    if (fenc->costEst[b - p0][p1 - b] >= 0 && fenc->rowSatds[b - p0][p1 - b][0] != -1)
        score = fenc->costEst[b - p0][p1 - b];
    else
    {
        m_weightedRef.isWeighted = false;
        if (m_param->bEnableWeightedPred && b == p1 && b != p0 && fenc->lowresMvs[0][b - p0 - 1][0].x == 0x7FFF)
        {
            if (!fenc->bIntraCalculated)
                estimateFrameCost(frames, b, b, b, 0);
            weightsAnalyse(frames, b, p0);
        }

        /* For each list, check to see whether we have lowres motion-searched this reference */
        m_bDoSearch[0] = b != p0 && fenc->lowresMvs[0][b - p0 - 1][0].x == 0x7FFF;
        m_bDoSearch[1] = b != p1 && fenc->lowresMvs[1][p1 - b - 1][0].x == 0x7FFF;

        if (m_bDoSearch[0]) fenc->lowresMvs[0][b - p0 - 1][0].x = 0;
        if (m_bDoSearch[1]) fenc->lowresMvs[1][p1 - b - 1][0].x = 0;

        m_curb = b;
        m_curp0 = p0;
        m_curp1 = p1;
        m_curframes = frames;
#if FIX_COSTEST_BUG&&KEEP_AS265_SAME_WITH_X265
    if(!fenc->bIntraCalculated) {
      fenc->costEst[0][0] = 0;
      fenc->costEstAq[0][0] = 0;
    }
    fenc->intraMbs[b - p0] = 0;
#endif
        fenc->costEst[b - p0][p1 - b] = 0;
        fenc->costEstAq[b - p0][p1 - b] = 0;

        for (int i = 0; i < m_heightInCU; i++)
        {
            m_rows[i].init();
            if (!fenc->bIntraCalculated)
                fenc->rowSatds[0][0][i] = 0;
            fenc->rowSatds[b - p0][p1 - b][i] = 0;
        }

        m_bFrameCompleted = false;

        if (m_pool)
        {
            WaveFront::enqueue();

            // enableAllRows must be already called
            enqueueRow(0);
            while (!m_bFrameCompleted)
                WaveFront::findJob(-1);

            WaveFront::dequeue();
        }
        else
        {
            for (int row = 0; row < m_heightInCU; row++)
                processRow(row, -1);

            x265_emms();
        }

        // Accumulate cost from each row
        for (int row = 0; row < m_heightInCU; row++)
        {
            fenc->costEst[b - p0][p1 - b] += m_rows[row].m_costEst;
            fenc->costEst[0][0] += m_rows[row].m_costIntra;
            if (m_param->rc.aqMode)
            {
                fenc->costEstAq[0][0] += m_rows[row].m_costIntraAq;
                fenc->costEstAq[b - p0][p1 - b] += m_rows[row].m_costEstAq;
            }
            fenc->intraMbs[b - p0] += m_rows[row].m_intraMbs;
        }

        fenc->bIntraCalculated = true;

        if (b != p1)
            score = (uint64_t)score * 100 / (130 + m_param->bFrameBias);
        if (b != p0 || b != p1) //Not Intra cost
            fenc->costEst[b - p0][p1 - b] = score;

#if DEBUG_FRAME_COST_OUTPUT&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
      int tmp_score = 0;
      int tmp_score2 = 0;
      fprintf(pf, "p0=%d b=%d p1=%d\n", p0, b, p1);
      for(int i = 0; i < m_heightInCU * m_widthInCU; i++) {
        int tmp_cost = fenc->lowresCosts[b - p0][p1 - b][i] & LOWRES_COST_MASK;
        int tmp_cost2 = fenc->lowresCosts[0][0][i] & LOWRES_COST_MASK;
        tmp_score += tmp_cost;
        tmp_score2 += tmp_cost2;
        fprintf(pf, "i=%d:tmp_cost=%d tmp_score=%d tmp_cost2=%d tmp_score2=%d\n",
                i, tmp_cost, tmp_score, tmp_cost2, tmp_score2);
      }
      fprintf(pf, "Frame cost(%d): %d %d\n", fenc->frameNum, fenc->costEst[b - p0][p1 - b], fenc->costEst[0][0]);
      fclose(pf);
    }

#endif
    }

    if (bIntraPenalty)
    {
        // arbitrary penalty for I-blocks after B-frames
        int ncu = NUM_CUS;
        score += (uint64_t)score * fenc->intraMbs[b - p0] / (ncu * 8);
    }
#if DEBUG_FRAME_COST_OUTPUT&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
    fprintf(pf, "Frame cost(%lld): %d intra_mbs=%d nmb=%d\n", fenc->frameNum, score, fenc->intraMbs[b - p0], NUM_CUS);
    fclose(pf);
  }

#endif
    return score;
}

uint32_t CostEstimate::weightCostLuma(Lowres **frames, int b, int p0, WeightParam *wp)
{
    Lowres *fenc = frames[b];
    Lowres *ref  = frames[p0];
    pixel *src = ref->fpelPlane[0];
    intptr_t stride = fenc->lumaStride;

    if (wp)
    {
        int offset = wp->inputOffset << (X265_DEPTH - 8);
        int scale = wp->inputWeight;
        int denom = wp->log2WeightDenom;
        int round = denom ? 1 << (denom - 1) : 0;
        int correction = IF_INTERNAL_PREC - X265_DEPTH; // intermediate interpolation depth
        int widthHeight = (int)stride;

        primitives.weight_pp(ref->buffer[0], m_wbuffer[0], stride, widthHeight, m_paddedLines,
                             scale, round << correction, denom + correction, offset);
        src = m_weightedRef.fpelPlane[0];
    }

    uint32_t cost = 0;
    intptr_t pixoff = 0;
    int mb = 0;

    for (int y = 0; y < fenc->lines; y += 8, pixoff = y * stride)
    {
        for (int x = 0; x < fenc->width; x += 8, mb++, pixoff += 8)
        {
            int satd = primitives.pu[LUMA_8x8].satd(src + pixoff, stride, fenc->fpelPlane[0] + pixoff, stride);
            cost += X265_MIN(satd, fenc->intraCost[mb]);
        }
    }

    return cost;
}

void CostEstimate::weightsAnalyse(Lowres **frames, int b, int p0)
{
    static const float epsilon = 1.f / 128.f;
    Lowres *fenc, *ref;

    fenc = frames[b];
    ref  = frames[p0];
    int deltaIndex = fenc->frameNum - ref->frameNum;

    /* epsilon is chosen to require at least a numerator of 127 (with denominator = 128) */
    float guessScale, fencMean, refMean;
    x265_emms();
    if (fenc->wp_ssd[0] && ref->wp_ssd[0])
        guessScale = sqrtf((float)fenc->wp_ssd[0] / ref->wp_ssd[0]);
    else
        guessScale = 1.0f;
    fencMean = (float)fenc->wp_sum[0] / (fenc->lines * fenc->width) / (1 << (X265_DEPTH - 8));
    refMean  = (float)ref->wp_sum[0] / (fenc->lines * fenc->width) / (1 << (X265_DEPTH - 8));

    /* Early termination */
    if (fabsf(refMean - fencMean) < 0.5f && fabsf(1.f - guessScale) < epsilon)
        return;

    int minoff = 0, minscale, mindenom;
    unsigned int minscore = 0, origscore = 1;
    int found = 0;

    m_w.setFromWeightAndOffset((int)(guessScale * 128 + 0.5f), 0, 7, true);
    mindenom = m_w.log2WeightDenom;
    minscale = m_w.inputWeight;

    origscore = minscore = weightCostLuma(frames, b, p0, NULL);

    if (!minscore)
        return;

    unsigned int s = 0;
    int curScale = minscale;
    int curOffset = (int)(fencMean - refMean * curScale / (1 << mindenom) + 0.5f);
    if (curOffset < -128 || curOffset > 127)
    {
        /* Rescale considering the constraints on curOffset. We do it in this order
         * because scale has a much wider range than offset (because of denom), so
         * it should almost never need to be clamped. */
        curOffset = x265_clip3(-128, 127, curOffset);
        curScale = (int)((1 << mindenom) * (fencMean - curOffset) / refMean + 0.5f);
        curScale = x265_clip3(0, 127, curScale);
    }
    SET_WEIGHT(m_w, 1, curScale, mindenom, curOffset);
    s = weightCostLuma(frames, b, p0, &m_w);
    COPY4_IF_LT(minscore, s, minscale, curScale, minoff, curOffset, found, 1);

    /* Use a smaller denominator if possible */
    while (mindenom > 0 && !(minscale & 1))
    {
        mindenom--;
        minscale >>= 1;
    }

    if (!found || (minscale == 1 << mindenom && minoff == 0) || (float)minscore / origscore > 0.998f)
        return;
    else
    {
        SET_WEIGHT(m_w, 1, minscale, mindenom, minoff);
        // set weighted delta cost
        fenc->weightedCostDelta[deltaIndex] = minscore / origscore;

        int offset = m_w.inputOffset << (X265_DEPTH - 8);
        int scale = m_w.inputWeight;
        int denom = m_w.log2WeightDenom;
        int round = denom ? 1 << (denom - 1) : 0;
        int correction = IF_INTERNAL_PREC - X265_DEPTH; // intermediate interpolation depth
        intptr_t stride = ref->lumaStride;
        int widthHeight = (int)stride;

        for (int i = 0; i < 4; i++)
            primitives.weight_pp(ref->buffer[i], m_wbuffer[i], stride, widthHeight, m_paddedLines,
                                 scale, round << correction, denom + correction, offset);

        m_weightedRef.isWeighted = true;
    }
}

void CostEstimate::processRow(int row, int /*threadId*/)
{
    ProfileScopeEvent(costEstimateRow);

    int realrow = m_heightInCU - 1 - row;
    Lowres **frames = m_curframes;
    ReferencePlanes *wfref0 = m_weightedRef.isWeighted ? &m_weightedRef : frames[m_curp0];

    /* Lowres lookahead goes backwards because the MVs are used as
     * predictors in the main encode.  This considerably improves MV
     * prediction overall. */
    for (int i = m_widthInCU - 1 - m_rows[row].m_completed; i >= 0; i--)
    {
        // TODO: use lowres MVs as motion candidates in full-res search
        m_rows[row].estimateCUCost(frames, wfref0, i, realrow, m_curp0, m_curp1, m_curb, m_bDoSearch);
        m_rows[row].m_completed++;

        if (m_rows[row].m_completed >= 2 && row < m_heightInCU - 1)
        {
            ScopedLock below(m_rows[row + 1].m_lock);
            if (m_rows[row + 1].m_active == false &&
                m_rows[row + 1].m_completed + 2 <= m_rows[row].m_completed)
            {
                m_rows[row + 1].m_active = true;
                enqueueRow(row + 1);
            }
        }

        ScopedLock self(m_rows[row].m_lock);
        if (row > 0 && (int32_t)m_rows[row].m_completed < m_widthInCU - 1 &&
            m_rows[row - 1].m_completed < m_rows[row].m_completed + 2)
        {
            m_rows[row].m_active = false;
            return;
        }
    }

    if (row == m_heightInCU - 1)
        m_bFrameCompleted = true;
}

void EstimateRow::init()
{
    m_costEst = 0;
    m_costEstAq = 0;
    m_costIntra = 0;
    m_costIntraAq = 0;
    m_intraMbs = 0;
    m_active = false;
    m_completed = 0;
}

#define OUTPUT_CU_PIXEL(frm,idx) do{\
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");\
    fprintf(pf, "%s:CU(%lld:%d,%d):\n",#frm,frm->frameNum, cux, cuy);\
    int i_stride=frm->lumaStride;\
    for(int y = -8; y < 8; y++) {\
      if(!y)\
        fprintf(pf, "-----------------------------------------------------------------\n");\
      for(int x = -8; x < 8; x++) {\
        if(!x)\
          fprintf(pf, "|%4d", (&frm->lowresPlane[idx][pelOffset])[x + y * i_stride]);\
        else\
          fprintf(pf, "%4d", (&frm->lowresPlane[idx][pelOffset])[x + y * i_stride]);\
      }\
      fprintf(pf, "\n");\
    }\
    fflush(pf);\
    fclose(pf);\
  }while(0)
void EstimateRow::estimateCUCost(Lowres **frames, ReferencePlanes *wfref0, int cux, int cuy, int p0, int p1, int b, bool bDoSearch[2])
{
    Lowres *fref1 = frames[p1];
    Lowres *fenc  = frames[b];

    const int bBidir = (b < p1);
    const int cuXY = cux + cuy * m_widthInCU;
    const int cuSize = X265_LOWRES_CU_SIZE;
    const intptr_t pelOffset = cuSize * cux + cuSize * cuy * fenc->lumaStride;

    // should this CU's cost contribute to the frame cost?
    const bool bFrameScoreCU = (cux > 0 && cux < m_widthInCU - 1 &&
                                cuy > 0 && cuy < m_heightInCU - 1) || m_widthInCU <= 2 || m_heightInCU <= 2;

    m_me.setSourcePU(fenc->lowresPlane[0], fenc->lumaStride, pelOffset, cuSize, cuSize);

    /* A small, arbitrary bias to avoid VBV problems caused by zero-residual lookahead blocks. */
    int lowresPenalty = 4;

    MV(*fenc_mvs[2]) = { &fenc->lowresMvs[0][b - p0 - 1][cuXY],
                         &fenc->lowresMvs[1][p1 - b - 1][cuXY] };
    int(*fenc_costs[2]) = { &fenc->lowresMvCosts[0][b - p0 - 1][cuXY],
                            &fenc->lowresMvCosts[1][p1 - b - 1][cuXY] };

    MV mvmin, mvmax;
    int bcost = m_me.COST_MAX;
#if KEEP_AS265_SAME_WITH_X265
  int i_inter_cost = 0, i_intra_cost = 0;
#endif
    int listused = 0;

    // establish search bounds that don't cross extended frame boundaries
    mvmin.x = (int16_t)(-cux * cuSize - 8);
    mvmin.y = (int16_t)(-cuy * cuSize - 8);
    mvmax.x = (int16_t)((m_widthInCU - cux - 1) * cuSize + 8);
    mvmax.y = (int16_t)((m_heightInCU - cuy - 1) * cuSize + 8);

#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_PIXEL_INFO&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
  OUTPUT_CU_PIXEL(fenc, 0);
#endif
    if (p0 != p1)
    {
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_PIXEL_INFO&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
    Lowres* fref0 = frames[p0];
    OUTPUT_CU_PIXEL(fref0, 0);
    OUTPUT_CU_PIXEL(fref0, 1);
    OUTPUT_CU_PIXEL(fref0, 2);
    OUTPUT_CU_PIXEL(fref0, 3);
    if(bBidir) {
      OUTPUT_CU_PIXEL(fref1, 0);
      OUTPUT_CU_PIXEL(fref1, 1);
      OUTPUT_CU_PIXEL(fref1, 2);
      OUTPUT_CU_PIXEL(fref1, 3);

    }
#endif
        for (int i = 0; i < 1 + bBidir; i++)
        {
            if (!bDoSearch[i])
            {
                /* Use previously calculated cost */
                COPY2_IF_LT(bcost, *fenc_costs[i], listused, i + 1);
                continue;
            }
            int numc = 0;
            MV mvc[4], mvp;
            MV *fenc_mv = fenc_mvs[i];

            /* Reverse-order MV prediction. */
            mvc[0] = 0;
            mvc[2] = 0;
#define MVC(mv) mvc[numc++] = mv;
            if (cux < m_widthInCU - 1)
                MVC(fenc_mv[1]);
            if (cuy < m_heightInCU - 1)
            {
                MVC(fenc_mv[m_widthInCU]);
                if (cux > 0)
                    MVC(fenc_mv[m_widthInCU - 1]);
                if (cux < m_widthInCU - 1)
                    MVC(fenc_mv[m_widthInCU + 1]);
            }
#undef MVC
            if (numc <= 1)
                mvp = mvc[0];
            else
            {
                median_mv(mvp, mvc[0], mvc[1], mvc[2]);
            }
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
        {
          FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
          for(int i = 0; i < numc; i++) {
            fprintf(pf, "\tmvc[%d]=(%d,%d) \n", i, mvc[i].x, mvc[i].y);
          }
          fclose(pf);
        }
#endif

            *fenc_costs[i] = m_me.motionEstimate(i ? fref1 : wfref0, mvmin, mvmax, mvp, numc, mvc, m_merange, *fenc_mvs[i]);
            COPY2_IF_LT(bcost, *fenc_costs[i], listused, i + 1);
        }

#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "\tfenc_mvs[%d]=(%d,%d) \n", cuXY, (*fenc_mvs[i]).x, (*fenc_mvs[i]).y);
        fclose(pf);
      }
#endif
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "Inter cost: %d\n", *fenc_costs[i]);
        fflush(pf);
        fclose(pf);
      }
#endif
        if (bBidir)
        {
            ALIGN_VAR_32(pixel, subpelbuf0[X265_LOWRES_CU_SIZE * X265_LOWRES_CU_SIZE]);
            ALIGN_VAR_32(pixel, subpelbuf1[X265_LOWRES_CU_SIZE * X265_LOWRES_CU_SIZE]);
            intptr_t stride0 = X265_LOWRES_CU_SIZE, stride1 = X265_LOWRES_CU_SIZE;
            pixel *src0 = wfref0->lowresMC(pelOffset, *fenc_mvs[0], subpelbuf0, stride0);
            pixel *src1 = fref1->lowresMC(pelOffset, *fenc_mvs[1], subpelbuf1, stride1);

            ALIGN_VAR_32(pixel, ref[X265_LOWRES_CU_SIZE * X265_LOWRES_CU_SIZE]);
            primitives.pu[LUMA_8x8].pixelavg_pp(ref, X265_LOWRES_CU_SIZE, src0, stride0, src1, stride1, 32);
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_ME_COST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        pixel* tmp_p1 = fenc->lowresPlane[0] + pelOffset;
        pixel* tmp_p2 = ref;
        fprintf(pf, "SRC:\t\t\t\tAVG:\n");
        for(int j = 0; j < 8; j++) {
          for(int i = 0; i < 8; i++) {
            fprintf(pf, "%4d", tmp_p1[i + j * fenc->lumaStride]);
          }
          fprintf(pf, "|");
          for(int i = 0; i < 8; i++) {
            fprintf(pf, "%4d", tmp_p2[i + j * X265_LOWRES_CU_SIZE]);
          }
          fprintf(pf, "\n");
        }
        fclose(pf);
      }
#endif
            int bicost = primitives.pu[LUMA_8x8].satd(fenc->lowresPlane[0] + pelOffset, fenc->lumaStride, ref, X265_LOWRES_CU_SIZE);
            COPY2_IF_LT(bcost, bicost, listused, 3);
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_ME_COST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "bidir: bcost=%d icost=%d\n", bcost, bicost);
        fclose(pf);
      }
#endif
            // Try 0,0 candidates
            src0 = wfref0->lowresPlane[0] + pelOffset;
            src1 = fref1->lowresPlane[0] + pelOffset;
            primitives.pu[LUMA_8x8].pixelavg_pp(ref, X265_LOWRES_CU_SIZE, src0, wfref0->lumaStride, src1, fref1->lumaStride, 32);
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_ME_COST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        pixel* tmp_p1 = fenc->lowresPlane[0] + pelOffset;
        pixel* tmp_p2 = ref;
        fprintf(pf, "SRC:\t\t\t\tAVG:\n");
        for(int j = 0; j < 8; j++) {
          for(int i = 0; i < 8; i++) {
            fprintf(pf, "%4d", tmp_p1[i + j * fenc->lumaStride]);
          }
          fprintf(pf, "|");
          for(int i = 0; i < 8; i++) {
            fprintf(pf, "%4d", tmp_p2[i + j * X265_LOWRES_CU_SIZE]);
          }
          fprintf(pf, "\n");
        }
        fclose(pf);
      }
#endif
            bicost = primitives.pu[LUMA_8x8].satd(fenc->lowresPlane[0] + pelOffset, fenc->lumaStride, ref, X265_LOWRES_CU_SIZE);
            COPY2_IF_LT(bcost, bicost, listused, 3);
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_ME_COST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "bidir: bcost=%d icost=%d\n", bcost, bicost);
        fclose(pf);
      }
#endif
        }
#if KEEP_AS265_SAME_WITH_X265
    i_inter_cost = bcost;
#endif
    }

    if (!fenc->bIntraCalculated)
    {
        ALIGN_VAR_32(pixel, prediction[X265_LOWRES_CU_SIZE * X265_LOWRES_CU_SIZE]);
        pixel neighbours[2][X265_LOWRES_CU_SIZE * 4 + 1];
        const int sizeIdx = X265_LOWRES_CU_BITS - 2; // partition size
        const int cuSize2 = cuSize << 1;

        pixel *pixCur = fenc->lowresPlane[0] + pelOffset;

        // Copy Above
        memcpy(neighbours[0], pixCur - 1 - fenc->lumaStride, (cuSize + 1) * sizeof(pixel));

        // Copy Left
        for (int i = 1; i < cuSize + 1; i++)
            neighbours[0][i + cuSize2] = pixCur[-1 - fenc->lumaStride + i * fenc->lumaStride];

        for (int i = 0; i < cuSize; i++)
        {
            // Copy above-last pixel
            neighbours[0][i + cuSize + 1] = neighbours[0][cuSize]; //neighbours[0][i + 9] = neighbours[0][8]
            // Copy left-last pixel
            neighbours[0][i + cuSize2 + cuSize + 1] = neighbours[0][cuSize2 + cuSize]; //neighbours[0][i + 25] = neighbours[0][24]
        }

        // Filter neighbour pixels with [1-2-1]
        neighbours[1][0]  = neighbours[0][0];  // Copy top-left pixel 
        neighbours[1][cuSize2] = neighbours[0][cuSize2]; //Copy top-right pixel
        neighbours[1][cuSize2 << 1] = neighbours[0][cuSize2 << 1]; // Bottom-left pixel

        neighbours[1][1]           = (neighbours[0][0] + (neighbours[0][1] << 1)           + neighbours[0][2] + 2)               >> 2;
        neighbours[1][cuSize2 + 1] = (neighbours[0][0] + (neighbours[0][cuSize2 + 1] << 1) + neighbours[0][cuSize2 + 1 + 1] + 2) >> 2;
        for (int i = 2; i < cuSize2; i++)
        {
            neighbours[1][i]           = (neighbours[0][i - 1]      + (neighbours[0][i] << 1)      + neighbours[0][i + 1]      + 2) >> 2;
            neighbours[1][cuSize2 + i] = (neighbours[0][cuSize2 + i - 1] + (neighbours[0][cuSize2 + i] << 1) + neighbours[0][cuSize2 + i + 1] + 2) >> 2;
        }

#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_PIXEL_INFO&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
      fprintf(pf, "TOP0 :");
      for(int i = 0; i < 2 * cuSize + 1; i++)
        fprintf(pf, "%4d", neighbours[0][i]);
      fflush(pf);
      fprintf(pf, "\nLeft0:");
      for(int i = 1; i < 2 * cuSize + 1; i++)
        fprintf(pf, "%4d", neighbours[0][i + cuSize2]);
      fflush(pf);

      fprintf(pf, "\nTOP1 :");
      for(int i = 0; i < 2 * cuSize + 1; i++)
        fprintf(pf, "%4d", neighbours[1][i]);
      fflush(pf);
      fprintf(pf, "\nLeft1:");
      for(int i = 1; i < 2 * cuSize + 1; i++)
        fprintf(pf, "%4d", neighbours[1][i + cuSize2]);
      fprintf(pf, "\n");
      fflush(pf);
      fclose(pf);
    }
#endif
#if DEBUG_INTRA_PREDICT&&KEEP_AS265_SAME_WITH_X265
#define PRINT_OUT_INTRA(INTRA_MODE,pix) do{\
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");\
    fprintf(pf, "Intra pred:%d \n", INTRA_MODE);\
    for(int j = 0; j < 8; j++) {\
      for(int i = 0; i < 8; i++) {\
        fprintf(pf, "%4d", pix[i + j * 8]);\
      }\
      fprintf(pf, "\n");\
    }\
    fclose(pf);\
  }while(0)
#define PRINT_OUT_INTRA_T(INTRA_MODE,pix) do{\
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");\
    fprintf(pf, "Intra pred:%d \n", INTRA_MODE);\
    for(int j = 0; j < 8; j++) {\
      for(int i = 0; i < 8; i++) {\
        fprintf(pf, "%4d", pix[j + i * 8]);\
      }\
      fprintf(pf, "\n");\
    }\
    fclose(pf);\
  }while(0)
#endif
        int icost = m_me.COST_MAX, ilowmode;
        primitives.cu[sizeIdx].intra_pred[DC_IDX](prediction, cuSize, neighbours[0], 0, (cuSize <= 16));
#if DEBUG_INTRA_PREDICT&&KEEP_AS265_SAME_WITH_X265
    PRINT_OUT_INTRA(DC_IDX, prediction);
#endif
        int cost = m_me.bufSATD(prediction, cuSize);
        COPY2_IF_LT(icost, cost, ilowmode, DC_IDX);

        pixel *planar = (cuSize >= 8) ? neighbours[1] : neighbours[0];
        primitives.cu[sizeIdx].intra_pred[PLANAR_IDX](prediction, cuSize, planar, 0, 0);
#if DEBUG_INTRA_PREDICT&&KEEP_AS265_SAME_WITH_X265
    PRINT_OUT_INTRA(PLANAR_IDX, prediction);
#endif
        cost = m_me.bufSATD(prediction, cuSize);
        COPY2_IF_LT(icost, cost, ilowmode, PLANAR_IDX);
#if 0//DEBUG_INTRA_PREDICT
    for(int i = 2; i < 35; i++) {
      if(i < 18)
        PRINT_OUT_INTRA_T(i, (m_predictions + i * predsize));
      else
        PRINT_OUT_INTRA(i, (m_predictions + i * predsize));
    }
#endif

        uint32_t mode, lowmode = 4;
        int acost = m_me.COST_MAX, filter;
#if USE_ALL_INTRA&&KEEP_AS265_SAME_WITH_X265
    for(mode = 2; mode < 35; mode += 1) {
      filter = !!(g_intraFilterFlags[mode] & cuSize);
      primitives.cu[sizeIdx].intra_pred[mode](prediction, cuSize, neighbours[filter], mode, cuSize <= 16);
      cost = m_me.bufSATD(prediction, cuSize);
      COPY2_IF_LT(acost, cost, lowmode, mode);
    }
#else
        for (mode = 5; mode < 35; mode += 5)
        {
            filter = !!(g_intraFilterFlags[mode] & cuSize);
            primitives.cu[sizeIdx].intra_pred[mode](prediction, cuSize, neighbours[filter], mode, cuSize <= 16);
            cost = m_me.bufSATD(prediction, cuSize);
            COPY2_IF_LT(acost, cost, lowmode, mode);
        }
        for (uint32_t dist = 2; dist >= 1; dist--)
        {
            int minusmode = lowmode - dist;
            int plusmode = lowmode + dist;

            mode = minusmode;
            filter = !!(g_intraFilterFlags[mode] & cuSize);
            primitives.cu[sizeIdx].intra_pred[mode](prediction, cuSize, neighbours[filter], mode, cuSize <= 16);
            cost = m_me.bufSATD(prediction, cuSize);
            COPY2_IF_LT(acost, cost, lowmode, mode);

            mode = plusmode;
            filter = !!(g_intraFilterFlags[mode] & cuSize);
            primitives.cu[sizeIdx].intra_pred[mode](prediction, cuSize, neighbours[filter], mode, cuSize <= 16);
            cost = m_me.bufSATD(prediction, cuSize);
            COPY2_IF_LT(acost, cost, lowmode, mode);
        }
#endif
        COPY2_IF_LT(icost, acost, ilowmode, lowmode);

        const int intraPenalty = 5 * m_lookAheadLambda;
        icost += intraPenalty + lowresPenalty; /* estimate intra signal cost */
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
    {
      FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
      fprintf(pf, "Intra cost:%d\n", icost);
      fclose(pf);
    }
#endif
#if FIX_INTRACOST_BUG&&KEEP_AS265_SAME_WITH_X265
    fenc->intra_cost[cuXY] = icost;
#endif
        fenc->intraCost[cuXY] = icost;
        fenc->intraMode[cuXY] = (uint8_t)ilowmode;

        int icostAq = icost;
        if (bFrameScoreCU)
        {
            m_costIntra += icost;
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_COST_EST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "output_intra[COST_EST]=%d\n", icost);
        fclose(pf);
      }
#endif
            if (fenc->invQscaleFactor)
            {
                icostAq = (icost * fenc->invQscaleFactor[cuXY] + 128) >> 8;
                m_costIntraAq += icostAq;
            }
        }
        fenc->rowSatds[0][0][cuy] += icostAq;
    }
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
    fprintf(pf, "bcost:%d lowres_penalty=%d\n", bcost, lowresPenalty);
    fclose(pf);
  }
#endif
#if KEEP_AS265_SAME_WITH_X265
  i_intra_cost = fenc->intraCost[cuXY];
#endif
    bcost += lowresPenalty;
    if (!bBidir)
    {  
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_INFO&&KEEP_AS265_SAME_WITH_X265
    {
      bool b_intra = fenc->intraCost[cuXY] < bcost;
      FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
      fprintf(pf, "b_intra=%d b_frame_score_cu=%d\n", b_intra, bFrameScoreCU);
      fclose(pf);
    }
#endif 
        if (fenc->intraCost[cuXY] < bcost)
        {
            if (bFrameScoreCU) m_intraMbs++;
            bcost = fenc->intraCost[cuXY];
            listused = 0;
        }

    }

    /* For I frames these costs were accumulated earlier */
    if (p0 != p1)
    {
        int bcostAq = bcost;
        if (bFrameScoreCU)
        {
            m_costEst += bcost;
#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_COST_EST&&KEEP_AS265_SAME_WITH_X265
      {
        FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
        fprintf(pf, "output_inter[COST_EST]=%d\n", bcost);
        fclose(pf);
      }
#endif
            if (fenc->invQscaleFactor)
            {
                bcostAq = (bcost * fenc->invQscaleFactor[cuXY] + 128) >> 8;
                m_costEstAq += bcostAq;
            }
        }
        fenc->rowSatds[b - p0][p1 - b][cuy] += bcostAq;
    }
    fenc->lowresCosts[b - p0][p1 - b][cuXY] = (uint16_t)(X265_MIN(bcost, LOWRES_COST_MASK) | (listused << LOWRES_COST_SHIFT));

#if DEBUG_FRAME_COST_OUTPUT&&DEBUG_CU_COST_OUTPUT&&KEEP_AS265_SAME_WITH_X265
  {
    FILE* pf = fopen(GET_FILENAME(DEBUG_FRAME_COST_OUTPUT), "a+");
    fprintf(pf, "Final MB cost(%lld:%d,%d):intra=%d inter=%d final=%d lowres_costs=%d\n",
            fenc->frameNum, cux, cuy, i_intra_cost, i_inter_cost, bcost,
            fenc->lowresCosts[b - p0][p1 - b][cuXY]);
    fclose(pf);
  }
#endif
}
