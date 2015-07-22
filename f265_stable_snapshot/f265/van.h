// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.


// Used for HM tests, included in all files of the encoder.

#ifndef VAN_H
#define VAN_H

#include <stdio.h>
#include "van_cfg.h"

// Permanently enable the use of HM ME.
#define VAN_USE_HM_ME

#ifdef VAN_LOAD_CTB_ANALYSIS
#define VAN_DUMP_CTB_ANALYSIS_PATH "/tmp/hevc_ctb_analysis.dump"
#endif

// Print the block pixels specified.
static inline void print_block_pix(uint8_t *p, int stride, int width, int height)
{
    for (uint32_t i = 0; i < height; i++)
    {
        printf("\t");
        for (uint32_t j = 0; j < width; j++)
        {
            uint32_t pix = p[i*stride + j];
            printf("%02x ", pix);
        }
        printf("\n");
    }
}

static inline void print_block_pix_s16(int16_t *p, int stride, int width, int height)
{
    for (uint32_t i = 0; i < height; i++)
    {
        printf("\t");
        for (uint32_t j = 0; j < width; j++)
        {
            int32_t pix = p[i*stride + j];
            char c = ' ';
            if (pix < 0)
            {
                pix = -pix;
                c = '-';
            }
            printf("%c%04x ", c, pix);
        }
        printf("\n");
    }
}

// Print the macroblock pixels specified. FIXME: remove this.
static inline void print_mb_pix(uint8_t *p, uint32_t stride, uint32_t mb_size)
{
    for (uint32_t i = 0; i < mb_size; i++)
    {
        printf("\t");
        for (uint32_t j = 0; j < mb_size; j++)
        {
            uint32_t pix = p[i*stride + j];
            printf("%02x ", pix);
        }
        printf("\n");
    }
}

// Print the macroblock pixels specified with two bytes per pixel, signed.
// FIXME: remove this.
static inline void print_mb_pix_s16(int16_t *p, uint32_t stride, uint32_t mb_size)
{
    for (uint32_t i = 0; i < mb_size; i++)
    {
        printf("\t");
        for (uint32_t j = 0; j < mb_size; j++)
        {
            int32_t pix = p[i*stride + j];
            char c = ' ';
            if (pix < 0)
            {
                pix = -pix;
                c = '-';
            }
            printf("%c%04x ", c, pix);
        }
        printf("\n");
    }
}

#ifdef VAN_USE_HM_ME
// Custom assert (and break point) for debug.
#if 1
#define MC_ASSERT(cond) {if (!(cond)) {printf("Assert failed in: %s: %s: %d\n", __FILE__, __FUNCTION__, __LINE__);\
                         fflush(0); asm("int $3"); }}
#else
#define MC_ASSERT(cond) { }
#endif
#endif

// Custom log function.
#if 1
#define MC_PRINTF(...) { printf("[F265] "); printf(__VA_ARGS__); printf("\n"); fflush(0);}
#else
#define MC_PRINTF(...) do {} while (0)
#endif

#endif
