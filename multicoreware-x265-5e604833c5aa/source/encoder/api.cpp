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
#include "bitstream.h"
#include "param.h"

#include "encoder.h"
#include "entropy.h"
#include "level.h"
#include "nal.h"
#include "bitcost.h"

using namespace x265;

extern "C"
x265_encoder *x265_encoder_open(x265_param *p)
{
    if (!p)
        return NULL;

    x265_param *param = X265_MALLOC(x265_param, 1);
    if (!param)
        return NULL;

    memcpy(param, p, sizeof(x265_param));
    x265_log(param, X265_LOG_INFO, "HEVC encoder version %s\n", x265_version_str);
    x265_log(param, X265_LOG_INFO, "build info %s\n", x265_build_info_str);

    x265_setup_primitives(param, param->cpuid);

    if (x265_check_params(param))
        return NULL;

    if (x265_set_globals(param))
        return NULL;

    Encoder *encoder = new Encoder;
    if (!param->rc.bEnableSlowFirstPass)
        x265_param_apply_fastfirstpass(param);

    // may change params for auto-detect, etc
    encoder->configure(param);
    
    // may change rate control and CPB params
    if (!enforceLevel(*param, encoder->m_vps))
    {
        delete encoder;
        return NULL;
    }

    // will detect and set profile/tier/level in VPS
    determineLevel(*param, encoder->m_vps);

    encoder->create();
    encoder->init();

    x265_print_params(param);

    return encoder;
}

extern "C"
int x265_encoder_headers(x265_encoder *enc, x265_nal **pp_nal, uint32_t *pi_nal)
{
    if (pp_nal && enc)
    {
        Encoder *encoder = static_cast<Encoder*>(enc);
        Entropy sbacCoder;
        Bitstream bs;
        encoder->getStreamHeaders(encoder->m_nalList, sbacCoder, bs);
        *pp_nal = &encoder->m_nalList.m_nal[0];
        if (pi_nal) *pi_nal = encoder->m_nalList.m_numNal;
        return encoder->m_nalList.m_occupancy;
    }

    return -1;
}

extern "C"
void x265_encoder_parameters(x265_encoder *enc, x265_param *out)
{
    if (enc && out)
    {
        Encoder *encoder = static_cast<Encoder*>(enc);
        memcpy(out, encoder->m_param, sizeof(x265_param));
    }
}

extern "C"
int x265_encoder_encode(x265_encoder *enc, x265_nal **pp_nal, uint32_t *pi_nal, x265_picture *pic_in, x265_picture *pic_out)
{
    if (!enc)
        return -1;

    Encoder *encoder = static_cast<Encoder*>(enc);
    int numEncoded;

    // While flushing, we cannot return 0 until the entire stream is flushed
    do
    {
        numEncoded = encoder->encode(pic_in, pic_out);
    }
    while (numEncoded == 0 && !pic_in && encoder->m_numDelayedPic);

    // do not allow reuse of these buffers for more than one picture. The
    // encoder now owns these analysisData buffers.
    if (pic_in)
    {
        pic_in->analysisData.intraData = NULL;
        pic_in->analysisData.interData = NULL;
    }

    if (pp_nal && numEncoded > 0)
    {
        *pp_nal = &encoder->m_nalList.m_nal[0];
        if (pi_nal) *pi_nal = encoder->m_nalList.m_numNal;
    }
    else if (pi_nal)
        *pi_nal = 0;

    return numEncoded;
}

extern "C"
void x265_encoder_get_stats(x265_encoder *enc, x265_stats *outputStats, uint32_t statsSizeBytes)
{
    if (enc && outputStats)
    {
        Encoder *encoder = static_cast<Encoder*>(enc);
        encoder->fetchStats(outputStats, statsSizeBytes);
    }
}

extern "C"
void x265_encoder_log(x265_encoder* enc, int argc, char **argv)
{
    if (enc)
    {
        Encoder *encoder = static_cast<Encoder*>(enc);
        encoder->writeLog(argc, argv);
    }
}

extern "C"
void x265_encoder_close(x265_encoder *enc)
{
    if (enc)
    {
        Encoder *encoder = static_cast<Encoder*>(enc);

        encoder->printSummary();
        encoder->destroy();
        delete encoder;
    }
}

extern "C"
void x265_cleanup(void)
{
    destroyROM();
    BitCost::destroy();
}

extern "C"
x265_picture *x265_picture_alloc()
{
    return (x265_picture*)x265_malloc(sizeof(x265_picture));
}

extern "C"
void x265_picture_init(x265_param *param, x265_picture *pic)
{
    memset(pic, 0, sizeof(x265_picture));

    pic->bitDepth = param->internalBitDepth;
    pic->colorSpace = param->internalCsp;
    pic->forceqp = X265_QP_AUTO;
    if (param->analysisMode)
    {
        uint32_t numPartitions   = 1 << (g_maxFullDepth * 2);
        uint32_t widthInCU       = (param->sourceWidth  + g_maxCUSize - 1) >> g_maxLog2CUSize;
        uint32_t heightInCU      = (param->sourceHeight + g_maxCUSize - 1) >> g_maxLog2CUSize;

        uint32_t numCUsInFrame   = widthInCU * heightInCU;
        pic->analysisData.numCUsInFrame = numCUsInFrame;
        pic->analysisData.numPartitions = numPartitions;
    }
}

extern "C"
void x265_picture_free(x265_picture *p)
{
    return x265_free(p);
}

int x265_alloc_analysis_data(x265_picture* pic)
{
    CHECKED_MALLOC(pic->analysisData.interData, x265_inter_data, pic->analysisData.numCUsInFrame * 85);
    CHECKED_MALLOC(pic->analysisData.intraData, x265_intra_data, 1);
    pic->analysisData.intraData->cuAddr     = NULL;
    pic->analysisData.intraData->depth      = NULL;
    pic->analysisData.intraData->modes      = NULL;
    pic->analysisData.intraData->partSizes  = NULL;
    pic->analysisData.intraData->poc        = NULL;
    CHECKED_MALLOC(pic->analysisData.intraData->depth, uint8_t, pic->analysisData.numPartitions * pic->analysisData.numCUsInFrame);
    CHECKED_MALLOC(pic->analysisData.intraData->modes, uint8_t, pic->analysisData.numPartitions * pic->analysisData.numCUsInFrame);
    CHECKED_MALLOC(pic->analysisData.intraData->partSizes, char, pic->analysisData.numPartitions * pic->analysisData.numCUsInFrame);
    CHECKED_MALLOC(pic->analysisData.intraData->cuAddr, uint32_t, pic->analysisData.numCUsInFrame);
    CHECKED_MALLOC(pic->analysisData.intraData->poc, int, pic->analysisData.numCUsInFrame);
    return 0;

fail:
    x265_free_analysis_data(pic);
    return -1;
}

void x265_free_analysis_data(x265_picture* pic)
{
    X265_FREE(pic->analysisData.interData);
    pic->analysisData.interData = NULL;
    X265_FREE(pic->analysisData.intraData->depth);
    X265_FREE(pic->analysisData.intraData->modes);
    X265_FREE(pic->analysisData.intraData->partSizes);
    X265_FREE(pic->analysisData.intraData->cuAddr);
    X265_FREE(pic->analysisData.intraData->poc);
    X265_FREE(pic->analysisData.intraData);
    pic->analysisData.intraData = NULL;
}
