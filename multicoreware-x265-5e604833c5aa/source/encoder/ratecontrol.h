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

#ifndef X265_RATECONTROL_H
#define X265_RATECONTROL_H

#include "common.h"
#include "sei.h"

namespace x265 {
// encoder namespace

class Encoder;
class Frame;
struct SPS;
class SEIBufferingPeriod;
#define BASE_FRAME_DURATION 0.04

/* Arbitrary limitations as a sanity check. */
#define MAX_FRAME_DURATION 1.00
#define MIN_FRAME_DURATION 0.01

#define CLIP_DURATION(f) Clip3(MIN_FRAME_DURATION, MAX_FRAME_DURATION, f)

/* Current frame stats for 2 pass */
struct FrameStats
{
    int         mvBits;    /* MV bits (MV+Ref+Block Type) */
    int         coeffBits; /* Texture bits (DCT coefs) */
    int         miscBits;

    int         iCuCnt;
    int         pCuCnt;
    int         skipCuCnt;
    
    /* CU type counts stored as percentage */
    double      percentIntra;
    double      percentInter;
    double      percentSkip;
};

struct Predictor
{
    double coeff;
    double count;
    double decay;
    double offset;
};

struct HRDTiming
{
    double cpbInitialAT;
    double cpbFinalAT;
    double dpbOutputTime;
    double cpbRemovalTime;
};

struct RateControlEntry
{
    int64_t lastSatd; /* Contains the picture cost of the previous frame, required for resetAbr and VBV */
    int sliceType;
    int bframes;
    int poc;
    int encodeOrder;
    int64_t leadingNoBSatd;
    bool bLastMiniGopBFrame;
    double blurredComplexity;
    double qpaRc;
    double qpAq;
    double qRceq;
    double frameSizePlanned;  /* frame Size decided by RateCotrol before encoding the frame */
    double bufferRate;
    double movingAvgSum;
    double   rowCplxrSum;
    int64_t  rowTotalBits;  /* update cplxrsum and totalbits at the end of 2 rows */
    double qpNoVbv;
    double bufferFill;
    double frameDuration;
    double clippedDuration;
    Predictor rowPreds[3][2];
    Predictor* rowPred[2];
    double frameSizeEstimated;  /* hold frameSize, updated from cu level vbv rc */
    double frameSizeMaximum;  /* max frame Size according to minCR restrictions and level of the video */
    bool isActive;
    SEIPictureTiming *picTimingSEI;
    HRDTiming        *hrdTiming;
    /* Required in 2-pass rate control */
    double iCuCount;
    double pCuCount;
    double skipCuCount;
    bool keptAsRef;
    double expectedVbv;
    double qScale;
    double newQScale;
    double newQp;
    int mvBits;
    int miscBits;
    int coeffBits;
    uint64_t expectedBits; /* total expected bits up to the current frame (current one excluded) */
};

class RateControl
{
public:

    x265_param* m_param;
    Slice*      m_curSlice;      /* all info about the current frame */
    SliceType   m_sliceType;     /* Current frame type */
    int         m_ncu;           /* number of CUs in a frame */
    int         m_qp;            /* updated qp for current frame */

    bool   m_isAbr;
    bool   m_isVbv;
    bool   m_isCbr;
    bool   m_singleFrameVbv;

    bool   m_isAbrReset;
    int    m_lastAbrResetPoc;

    double m_frameDuration;     /* current frame duration in seconds */
    double m_bitrate;
    double m_rateFactorConstant;
    double m_bufferSize;
    double m_bufferFillFinal;  /* real buffer as of the last finished frame */
    double m_bufferFill;       /* planned buffer, if all in-progress frames hit their bit budget */
    double m_bufferRate;       /* # of bits added to buffer_fill after each frame */
    double m_vbvMaxRate;       /* in kbps */
    double m_rateFactorMaxIncrement; /* Don't allow RF above (CRF + this value). */
    double m_rateFactorMaxDecrement; /* don't allow RF below (this value). */

    Predictor m_pred[5];
    Predictor m_predBfromP;

    int       m_leadingBframes;
    int64_t   m_bframeBits;
    int64_t   m_currentSatd;
    int       m_qpConstant[3];
    double    m_ipOffset;
    double    m_pbOffset;

    int      m_lastNonBPictType;
    int64_t  m_leadingNoBSatd;

    double   m_cplxrSum;          /* sum of bits*qscale/rceq */
    double   m_wantedBitsWindow;  /* target bitrate * window */
    double   m_accumPQp;          /* for determining I-frame quant */
    double   m_accumPNorm;
    double   m_lastQScaleFor[3];  /* last qscale for a specific pict type, used for max_diff & ipb factor stuff */
    double   m_lstep;
    double   m_shortTermCplxSum;
    double   m_shortTermCplxCount;
    double   m_lastRceq;
    double   m_qCompress;
    int64_t  m_totalBits;        /* total bits used for already encoded frames */
    int      m_framesDone;       /* # of frames passed through RateCotrol already */
    double   m_fps;
    int64_t  m_satdCostWindow[50];
    int      m_sliderPos;
    int64_t  m_encodedBitsWindow[50];
    /* a common variable on which rateControlStart, rateControlEnd and rateControUpdateStats waits to
     * sync the calls to these functions. For example
     * -F2:
     * rceStart  10
     * rceUpdate 10
     * rceEnd    9
     * rceStart  11
     * rceUpdate 11
     * rceEnd    10
     * rceStart  12
     * rceUpdate 12
     * rceEnd    11 */
    ThreadSafeInteger m_startEndOrder;
    int      m_finalFrameCount;   /* set when encoder begins flushing */
    bool     m_bTerminated;       /* set true when encoder is closing */

    /* hrd stuff */
    SEIBufferingPeriod m_bufPeriodSEI;
    double   m_nominalRemovalTime;
    double   m_prevCpbFinalAT;

    /* 2 pass */
    bool     m_2pass;
    FILE*    m_statFileOut;
    FILE*    m_cutreeStatFileOut;
    FILE*    m_cutreeStatFileIn;
    int      m_numEntries;
    RateControlEntry *m_rce2Pass;
    double   m_lastAccumPNorm;
    int64_t  m_predictedBits;
    double   m_expectedBitsSum;   /* sum of qscale2bits after rceq, ratefactor, and overflow, only includes finished frames */
    struct
    {
        uint16_t *qpBuffer[2]; /* Global buffers for converting MB-tree quantizer data. */
        int qpBufPos;          /* In order to handle pyramid reordering, QP buffer acts as a stack.
                                * This value is the current position (0 or 1). */
    } m_cuTreeStats;

    RateControl(x265_param *p);
    void setFinalFrameCount(int count);
    void terminate();          /* un-block all waiting functions so encoder may close */
    void destroy();

    // to be called for each curFrame to process RateControl and set QP
    int rateControlStart(Frame* curFrame, RateControlEntry* rce, Encoder* enc);
    void calcAdaptiveQuantFrame(Frame *curFrame);
    void rateControlUpdateStats(RateControlEntry* rce);
    int rateControlEnd(Frame* curFrame, int64_t bits, RateControlEntry* rce, FrameStats* stats);
    int rowDiagonalVbvRateControl(Frame* curFrame, uint32_t row, RateControlEntry* rce, double& qpVbv);
    void hrdFullness(SEIBufferingPeriod* sei);
    bool init(const SPS* sps);
    void initHRD(SPS* sps);
    int rateControlSliceType(int frameNum);
    bool cuTreeReadFor2Pass(Frame* curFrame);

protected:

    static const double s_amortizeFraction;
    static const int    s_amortizeFrames;
    static const int    s_slidingWindowFrames;
    static const char  *s_defaultStatFileName;

    int m_residualFrames;
    int m_partialResidualFrames;
    int m_residualCost;
    int m_partialResidualCost;

    double getQScale(RateControlEntry *rce, double rateFactor);
    double rateEstimateQscale(Frame* pic, RateControlEntry *rce); // main logic for calculating QP based on ABR
    void accumPQpUpdate();
    uint32_t acEnergyCu(Frame* pic, uint32_t block_x, uint32_t block_y);

    void updateVbv(int64_t bits, RateControlEntry* rce);
    void updatePredictor(Predictor *p, double q, double var, double bits);
    double clipQscale(Frame* pic, RateControlEntry* rce, double q);
    void updateVbvPlan(Encoder* enc);
    double predictSize(Predictor *p, double q, double var);
    void checkAndResetABR(RateControlEntry* rce, bool isFrameDone);
    double predictRowsSizeSum(Frame* pic, RateControlEntry* rce, double qpm, int32_t& encodedBits);
    bool initPass2();
    double getDiffLimitedQScale(RateControlEntry *rce, double q);
    double countExpectedBits();
    bool vbv2Pass(uint64_t allAvailableBits);
    bool findUnderflow(double *fills, int *t0, int *t1, int over);
    bool fixUnderflow(int t0, int t1, double adjustment, double qscaleMin, double qscaleMax);
};
}
#endif // ifndef X265_RATECONTROL_H
