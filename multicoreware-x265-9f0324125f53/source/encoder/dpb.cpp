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

#include "common.h"
#include "frame.h"
#include "framedata.h"
#include "picyuv.h"
#include "slice.h"

#include "dpb.h"

using namespace x265;

DPB::~DPB()
{
    while (!m_freeList.empty())
    {
        Frame* curFrame = m_freeList.popFront();
        curFrame->destroy();
        delete curFrame;
    }

    while (!m_picList.empty())
    {
        Frame* curFrame = m_picList.popFront();
        curFrame->destroy();
        delete curFrame;
    }

    while (m_picSymFreeList)
    {
        FrameData* next = m_picSymFreeList->m_freeListNext;
        m_picSymFreeList->destroy();

        m_picSymFreeList->m_reconPic->destroy();
        delete m_picSymFreeList->m_reconPic;

        delete m_picSymFreeList;
        m_picSymFreeList = next;
    }
}

// move unreferenced pictures from picList to freeList for recycle
void DPB::recycleUnreferenced()
{
    Frame *iterFrame = m_picList.first();

    while (iterFrame)
    {
        Frame *curFrame = iterFrame;
        iterFrame = iterFrame->m_next;
        if (!curFrame->m_encData->m_bHasReferences && !curFrame->m_countRefEncoders)
        {
            curFrame->m_reconRowCount.set(0);
            curFrame->m_bChromaExtended = false;

            // iterator is invalidated by remove, restart scan
            m_picList.remove(*curFrame);
            iterFrame = m_picList.first();

            m_freeList.pushBack(*curFrame);
            curFrame->m_encData->m_freeListNext = m_picSymFreeList;
            m_picSymFreeList = curFrame->m_encData;
            curFrame->m_encData = NULL;
            curFrame->m_reconPic = NULL;
        }
    }
}

void DPB::prepareEncode(Frame *newFrame)
{
    Slice* slice = newFrame->m_encData->m_slice;
    slice->m_poc = newFrame->m_poc;

    int pocCurr = slice->m_poc;
    int type = newFrame->m_lowres.sliceType;
    bool bIsKeyFrame = newFrame->m_lowres.bKeyframe;

    slice->m_nalUnitType = getNalUnitType(pocCurr, bIsKeyFrame);

#if DEBUG_NALUTYPE
  char* x265_type2string[] = {
    "265_TYPE_AUTO",
    "265_TYPE_IDR",
    "265_TYPE_I",
    "265_TYPE_P",
    "265_TYPE_BREF",
    "265_TYPE_B",
    "265_TYPE_KEYFRAME",
  };
  fprintf(stdout, "NaluType:poc=%d(%d) type=%s lastIDR=%d\n", slice->m_poc, newFrame->m_lowres.frameNum, x265_type2string[type],
          m_lastIDR);
  fflush(stdout);
  //if(slice->m_nalUnitType==NAL_UNIT_CODED_SLICE_RADL_R){
  //  fprintf(stdout,"poc=%d(%d) type=%s\n",slice->m_poc,newFrame->m_lowres.frameNum,x265_type2string[type]);
  //  fflush(stdout);
  //}
#endif
    if (slice->m_nalUnitType == NAL_UNIT_CODED_SLICE_IDR_W_RADL)
        m_lastIDR = pocCurr;
    slice->m_lastIDR = m_lastIDR;
    slice->m_sliceType = IS_X265_TYPE_B(type) ? B_SLICE : (type == X265_TYPE_P) ? P_SLICE : I_SLICE;

    if (type == X265_TYPE_B)
    {
        // change from _R "referenced" to _N "non-referenced" NAL unit type
        switch (slice->m_nalUnitType)
        {
        case NAL_UNIT_CODED_SLICE_TRAIL_R:
            slice->m_nalUnitType = NAL_UNIT_CODED_SLICE_TRAIL_N;
            break;
        case NAL_UNIT_CODED_SLICE_RADL_R:
            slice->m_nalUnitType = NAL_UNIT_CODED_SLICE_RADL_N;
            break;
        case NAL_UNIT_CODED_SLICE_RASL_R:
            slice->m_nalUnitType = NAL_UNIT_CODED_SLICE_RASL_N;
            break;
        default:
            break;
        }
    }

#if DEBUG_NALUTYPE
  char* nalutype2string[] = {
    "NAL_TRAIL_N",
    "NAL_TRAIL_R",
    "NAL_TSA_N",
    "NAL_TSA_R",
    "NAL_STSA_N",
    "NAL_STSA_R",
    "NAL_RADL_N",
    "NAL_RADL_R",
    "NAL_RASL_N",
    "NAL_RASL_R",
    "NAL_RSV_VCL_N10",
    "NAL_RSV_VCL_N12",
    "NAL_RSV_VCL_N14",
    "NAL_RSV_VCL_R11",
    "NAL_RSV_VCL_R13",
    "NAL_RSV_VCL_R15",
    "NAL_BLA_W_LP",
    "NAL_BLA_W_RADL",
    "NAL_BLA_N_LP",
    "NAL_IDR_W_RADL",
    "NAL_IDR_N_LP",
    "NAL_CRA_NUT",
  };
  fprintf(stdout, "NaluType:\tnalu=%s\n", nalutype2string[slice->m_nalUnitType]);
  fflush(stdout);
#endif

    /* m_bHasReferences starts out as true for non-B pictures, and is set to false
     * once no more pictures reference it */
    newFrame->m_encData->m_bHasReferences = IS_REFERENCED(newFrame);

    m_picList.pushFront(*newFrame);

    // Do decoding refresh marking if any
    decodingRefreshMarking(pocCurr, slice->m_nalUnitType);

    computeRPS(pocCurr, slice->isIRAP(), &slice->m_rps, slice->m_sps->maxDecPicBuffering);

    // Mark pictures in m_piclist as unreferenced if they are not included in RPS
    applyReferencePictureSet(&slice->m_rps, pocCurr);

    slice->m_numRefIdx[0] = X265_MIN(m_maxRefL0, slice->m_rps.numberOfNegativePictures); // Ensuring L0 contains just the -ve POC
    slice->m_numRefIdx[1] = X265_MIN(m_maxRefL1, slice->m_rps.numberOfPositivePictures);
    slice->setRefPicList(m_picList);

    X265_CHECK(slice->m_sliceType != B_SLICE || slice->m_numRefIdx[1], "B slice without L1 references (non-fatal)\n");

    if (slice->m_sliceType == B_SLICE)
    {
        /* TODO: the lookahead should be able to tell which reference picture
         * had the least motion residual.  We should be able to use that here to
         * select a colocation reference list and index */
        slice->m_colFromL0Flag = false;
        slice->m_colRefIdx = 0;
        slice->m_bCheckLDC = false;
    }
    else
    {
        slice->m_bCheckLDC = true;
        slice->m_colFromL0Flag = true;
        slice->m_colRefIdx = 0;
    }
    slice->m_sLFaseFlag = (SLFASE_CONSTANT & (1 << (pocCurr % 31))) > 0;

    /* Increment reference count of all motion-referenced frames to prevent them
     * from being recycled. These counts are decremented at the end of
     * compressFrame() */
    int numPredDir = slice->isInterP() ? 1 : slice->isInterB() ? 2 : 0;
    for (int l = 0; l < numPredDir; l++)
    {
        for (int ref = 0; ref < slice->m_numRefIdx[l]; ref++)
        {
            Frame *refpic = slice->m_refPicList[l][ref];
            ATOMIC_INC(&refpic->m_countRefEncoders);
        }
    }
}

void DPB::computeRPS(int curPoc, bool isRAP, RPS * rps, unsigned int maxDecPicBuffer)
{
    unsigned int poci = 0, numNeg = 0, numPos = 0;

    Frame* iterPic = m_picList.first();

#if DEBUG_CRA
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_CRA), "a");
    fprintf(fp, "curPoc=%lld isRAP=%d maxDecPicBuffer=%d\n",
            curPoc, isRAP, maxDecPicBuffer);
    fclose(fp);
  }
#endif
    while (iterPic && (poci < maxDecPicBuffer - 1))
    {
        if ((iterPic->m_poc != curPoc) && iterPic->m_encData->m_bHasReferences)
        {
            rps->poc[poci] = iterPic->m_poc;
            rps->deltaPOC[poci] = rps->poc[poci] - curPoc;
            (rps->deltaPOC[poci] < 0) ? numNeg++ : numPos++;
            rps->bUsed[poci] = !isRAP;
#if DEBUG_CRA
      {
        FILE* fp = fopen(GET_FILENAME(DEBUG_CRA), "a");
        fprintf(fp, "poc=%lld refidc=%d pocd_refpic=%d b_used_by_currpic=%d poci=%d numNeg=%d numPos=%d\n",
                iterPic->m_poc, rps->poc[poci], rps->deltaPOC[poci], rps->bUsed[poci], poci, numNeg, numPos);
        fclose(fp);
      }
#endif
            poci++;
        }
        iterPic = iterPic->m_next;
    }

    rps->numberOfPictures = poci;
    rps->numberOfPositivePictures = numPos;
    rps->numberOfNegativePictures = numNeg;

#if DEBUG_CRA
  {
    FILE* fp = fopen(GET_FILENAME(DEBUG_CRA), "a");
    fprintf(fp, "curPoc=%lld poci=%d numPos=%d numNeg=%d\n",
            curPoc, poci, numPos, numNeg);
    fclose(fp);
  }
#endif
    rps->sortDeltaPOC();
}

/* Marking reference pictures when an IDR/CRA is encountered. */
void DPB::decodingRefreshMarking(int pocCurr, NalUnitType nalUnitType)
{
    if (nalUnitType == NAL_UNIT_CODED_SLICE_IDR_W_RADL)
    {
        /* If the nal_unit_type is IDR, all pictures in the reference picture
         * list are marked as "unused for reference" */
        Frame* iterFrame = m_picList.first();
        while (iterFrame)
        {
            if (iterFrame->m_poc != pocCurr)
                iterFrame->m_encData->m_bHasReferences = false;
            iterFrame = iterFrame->m_next;
        }
    }
    else // CRA or No DR
    {
        if (m_bRefreshPending && pocCurr > m_pocCRA)
        {
            /* If the bRefreshPending flag is true (a deferred decoding refresh
             * is pending) and the current temporal reference is greater than
             * the temporal reference of the latest CRA picture (pocCRA), mark
             * all reference pictures except the latest CRA picture as "unused
             * for reference" and set the bRefreshPending flag to false */
            Frame* iterFrame = m_picList.first();
            while (iterFrame)
            {
                if (iterFrame->m_poc != pocCurr && iterFrame->m_poc != m_pocCRA)
                    iterFrame->m_encData->m_bHasReferences = false;
                iterFrame = iterFrame->m_next;
            }

            m_bRefreshPending = false;
        }
        if (nalUnitType == NAL_UNIT_CODED_SLICE_CRA)
        {
            /* If the nal_unit_type is CRA, set the bRefreshPending flag to true
             * and pocCRA to the temporal reference of the current picture */
            m_bRefreshPending = true;
            m_pocCRA = pocCurr;
        }
    }

    /* Note that the current picture is already placed in the reference list and
     * its marking is not changed.  If the current picture has a nal_ref_idc
     * that is not 0, it will remain marked as "used for reference" */
}

/** Function for applying picture marking based on the Reference Picture Set */
void DPB::applyReferencePictureSet(RPS *rps, int curPoc)
{
    // loop through all pictures in the reference picture buffer
    Frame* iterFrame = m_picList.first();
    while (iterFrame)
    {
        if (iterFrame->m_poc != curPoc && iterFrame->m_encData->m_bHasReferences)
        {
            // loop through all pictures in the Reference Picture Set
            // to see if the picture should be kept as reference picture
            bool referenced = false;
            for (int i = 0; i < rps->numberOfPositivePictures + rps->numberOfNegativePictures; i++)
            {
                if (iterFrame->m_poc == curPoc + rps->deltaPOC[i])
                {
                    referenced = true;
                    break;
                }
            }
            if (!referenced)
                iterFrame->m_encData->m_bHasReferences = false;
        }
        iterFrame = iterFrame->m_next;
    }
}

/* deciding the nal_unit_type */
NalUnitType DPB::getNalUnitType(int curPOC, bool bIsKeyFrame)
{
    if (!curPOC)
        return NAL_UNIT_CODED_SLICE_IDR_W_RADL;

    if (bIsKeyFrame)
        return m_bOpenGOP ? NAL_UNIT_CODED_SLICE_CRA : NAL_UNIT_CODED_SLICE_IDR_W_RADL;

    if (m_pocCRA && curPOC < m_pocCRA)
        // All leading pictures are being marked as TFD pictures here since
        // current encoder uses all reference pictures while encoding leading
        // pictures. An encoder can ensure that a leading picture can be still
        // decodable when random accessing to a CRA/CRANT/BLA/BLANT picture by
        // controlling the reference pictures used for encoding that leading
        // picture. Such a leading picture need not be marked as a TFD picture.
        return NAL_UNIT_CODED_SLICE_RASL_R;

    if (m_lastIDR && curPOC < m_lastIDR)
        return NAL_UNIT_CODED_SLICE_RADL_R;

    return NAL_UNIT_CODED_SLICE_TRAIL_R;
}
