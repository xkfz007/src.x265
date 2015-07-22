// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.


///////////////////////////////////////////////////////////////////////////////
// Documentation.
//
// This file contains bit depth independent code.


///////////////////////////////////////////////////////////////////////////////
// Includes.

#ifndef F265_BDI_H
#define F265_BDI_H

// Public interface.
#include "f265/f265.h"

// Benchmarking support.
#ifdef F265_PERF_BENCHMARK
#include "f265/perf.h"
#endif

// Standard headers.
#ifdef F265_HAVE_STDIO
#ifndef _LARGEFILE64_SOURCE     // Needed for YUV dump features.
#define _LARGEFILE64_SOURCE
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <stddef.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <stdlib.h>             // Needed for exit.
#include <assert.h>
#include <string.h>
#include <math.h>
#include <limits.h>             // Needed for HM compatibility.
#include <float.h>              // Needed for HM compatibility.
#endif

// Threading.
#ifdef F265_HAVE_PTHREAD
#include <pthread.h>
#endif


///////////////////////////////////////////////////////////////////////////////
// Standard macros.

// Branch prediction hints.
#ifndef likely
#define likely(x)   (__builtin_expect((x),1))
#endif
#ifndef unlikely
#define unlikely(x) (__builtin_expect((x),0))
#endif

// Prevent and force inlining.
#undef noinline
#undef finline

#ifdef F265_NO_INLINING
#define noinline
#define finline
#else
#define noinline __attribute__ ((noinline))
// GCC randomly fails to inline functions after trivial modifications, with a
// message like "sorry, unimplemented, function not inlinable". Uncomment when
// gcc gets repaired or we switch to a better compiler.
// #define finline __attribute__ ((always_inline))
#define finline inline
#endif

// Data alignment specifications.
// NOTE This prefix screws up Doxygen's parser, disable it when producing the
// documentation to get all symbols.
#ifndef DOXY
#define F265_ALIGN16 __attribute__((aligned(16)))
#define F265_ALIGN64 __attribute__((aligned(64)))
#else
#define F265_ALIGN16
#define F265_ALIGN64
#endif

// Prevent unused warning on static function.
#define F265_UNUSED __attribute__ ((unused))

// Shared constant data. FIXME, find out a definition that works everywhere.
#ifdef __cplusplus
#define F265_SHARED_CONST extern __attribute__((weak)) const
#else
#define F265_SHARED_CONST __attribute__((weak)) const
#endif

// Breakpoint macro.
#define F265_BREAK __asm__ __volatile__ ("int $03")

// Align X to boundary A. The ~(A-1) part is the alignment mask, the X+A-1 part
// allows the mask to be applied without the aligned value becoming smaller than
// X.
#define F265_ALIGN_VAL(X, A) (((X)+((A)-1))&~((A)-1))

// Math utilities.
#define F265_ABS(a) (((a) < 0) ? (-(a)) : (a))
#define F265_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define F265_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define F265_CLAMP(v, min, max) ((v) < (min) ? (min) : (v) > (max) ? (max) : (v))

// Swap two values of type X.
#define F265_SWAP(X, A, B) { X f265_swap_tmp = (A); (A) = (B); (B) = f265_swap_tmp; }

// Swap bytes in a 32 bits variable (equivalent to bswap asm instruction).
#define F265_BYTESWAP32(a) ((((a) >> 24)) +  (((a) >> 8)  & 0x0000FF00) + (((a) << 8)  & 0x00FF0000) + (((a) << 24)))

// Set or clear the flag specified in the bitfield specified according to the
// value specified. Value must be 0 or 1.
#define F265_SET_FLAG(bitfield, flag, value) (bitfield) = (bitfield)&(~(flag))|((~((value)-1))&(flag))

// Return a boolean value reflecting whether the flag specified is set in the
// bitfield specified.
#define F265_GET_FLAG(bitfield, flag) (!!((bitfield)&(flag)))

// Threading primitives.
#if defined F265_HAVE_PTHREAD
    #define F265_HAVE_MT
    #define f265_mutex pthread_mutex_t
    #define f265_cond_var pthread_cond_t
    #define f265_thread_handle pthread_t
    #define f265_mutex_init(M) pthread_mutex_init(M, NULL)
    #define f265_mutex_deinit(M) pthread_mutex_destroy(M)
    #define f265_cond_init(C) pthread_cond_init(C, NULL)
    #define f265_cond_deinit(C) pthread_cond_destroy(C)
    #define f265_thread_create(T, F, C) pthread_create(T, NULL, F, C)
    #define f265_thread_join(T) pthread_join(T, NULL)
    #define f265_mutex_lock(M) pthread_mutex_lock(M)
    #define f265_mutex_unlock(M) pthread_mutex_unlock(M)
    #define f265_cond_wait(C, M) pthread_cond_wait(C, M)
    #define f265_cond_signal(C) pthread_cond_signal(C)
#else
    #undef F265_HAVE_MT
    #define f265_mutex_lock(M) do { } while (0)
    #define f265_mutex_unlock(M) do { } while (0)
    #define f265_cond_wait(C, M) do { } while (0)
    #define f265_cond_signal(C) do { } while (0)
#endif


///////////////////////////////////////////////////////////////////////////////
// Constants.

// Internal operation return codes.
#define F265_RET_RETRY                  5       // The operation must be retried.

// Padding outside the luma frame.
#define F265_LUMA_PLANE_PADDING         72

// Maximum fullpel motion vector range.
#define F265_MAX_MV_FPEL                127

// PMV out-of-bounds tolerance for cost computation.
#define F265_PMV_OOB                    32

// Search out-of-bounds tolerance.
#define F265_SEARCH_OOB                 6

// Maximum SAD value (64*64*3*2^14*2). This leaves enough room for the
// distortion plus the syntax bits. Fits in 32-bits.
#define F265_MAX_SAD                    0x18000000

// Maximum SSD value (64*64*3*(2^14)^2*2*32). Same as above for the SSD, plus
// some scaling for precision. Fits in 64-bits.
#define F265_MAX_SSD                   0xc00000000000ll

// Number of costs in the quartelpel motion vector cost table. The actual and
// the predicted motion vectors can be opposite. The table is indexed by the
// predicted motion vector, which can be negative. There is an extra cost for
// the zero motion vector.
//
//     |MAX|MAX+OOB|MAX+OOB|MAX|
//     <---+------>C<------+--->
//         ^       ^       ^
//       -PMV   Center   +PMV
#define F265_NB_MV_COSTS                ((F265_MAX_MV_FPEL*4+F265_SEARCH_OOB*2)*4+1)

// The default weighted prediction scaling denominator is 64. The scaling
// numerator can range between 0 and 127, allowing pixel scaling from 0% to
// 198%. We do not use negative scaling numerators.
#define F265_WEIGHTP_LOG2_DENOM         6


// Thread status.
#define F265_THREAD_BUSY                1       // Processing a job.
#define F265_THREAD_IDLE                2       // Waiting for a job.
#define F265_THREAD_EXIT                3       // Exiting.


// Frame flags visible from the main thread only.
#define F265_FF_ENC                     (1<<0)  // Encoded frame.
#define F265_FF_LAST_ENC                (1<<1)  // Last encoded frame.
#define F265_FF_DPB_REF                 (1<<2)  // The frame is used for reference in the current DPB.
#define F265_FF_LA_REF                  (1<<3)  // The frame is used for reference by the lookahead.

// Frame flags visible from all threads.
#define F265_FF_KEY                     (1<<4)  // Key frame (IDR/CRA/periodic intra refresh/etc).
#define F265_FF_CRA                     (1<<5)  // Clean random access frame.
#define F265_FF_IDR                     (1<<6)  // Instantaneous decoding refresh frame.
#define F265_FF_REORDER                 (1<<7)  // The reference lists were reordered.
#define F265_FF_REF                     (1<<8)  // The frame is intended to be used by reference
                                                // (e.g. not a disposable B frame).
#define F265_FF_DEBLOCK                 (1<<9)  // The frame must be deblocked.


// Parameter flags. At encoder level, a flag enables a feature globally. At
// slice or CTB level, a flag enables or disables a feature for the duration of
// the frame or the CTB.

// Use the deblocking filter.
#define F265_PF_DEBLOCK                 (1<<1)

// Deblock frames not used for reference or SAO.
#define F265_PF_DEBLOCK_ALL             (1<<2)

// Use the SAO filter.
#define F265_PF_SAO                     (1<<3)

// Use scaling lists.
#define F265_PF_SCALING                 (1<<4)

// Use pulse code modulation.
#define F265_PF_PCM                     (1<<5)

// Use transform & quantization bypass.
#define F265_PF_TRANSQUANT_BYPASS       (1<<6)

// Use sign data hiding.
#define F265_PF_SIGN_HIDING             (1<<7)

// Use transform skip.
#define F265_PF_TRANSFORM_SKIP          (1<<8)

// Use strong intra smoothing.
#define F265_PF_SMOOTH_INTRA            (1<<9)

// Use assymetric motion partitions.
#define F265_PF_AMP                     (1<<10)

// Use temporal motion vectors.
#define F265_PF_TMV                     (1<<11)

// Use chroma motion estimation.
#define F265_PF_CHROMA_ME               (1<<12)

// True if the subpel interpolation is used.
#define F265_PF_SUBPEL_ME               (1<<13)

// True if the lookahead is enabled.
#define F265_PF_LA                      (1<<14)

// Allow the frame to be reencoded.
#define F265_PF_ALLOW_FRAME_RECODE      (1<<15)

// Allow the CTB to be reencoded.
#define F265_PF_ALLOW_CTB_RECODE        (1<<16)

// Use RDO for quantization.
#define F265_PF_RDOQ                    (1<<17)

// Use wavefront parallel processing.
#define F265_PF_WPP                     (1<<18)

// Number of CABAC contexts.
#define F265_NB_CABAC_CTX               154

// CABAC context offset by syntax element.
#define F265_CO_SAO_MERGE               0x0
#define F265_CO_SAO_TYPE                0x1
#define F265_CO_SPLIT_CU                0x2
#define F265_CO_TRANSQUANT_BYPASS       0x5
#define F265_CO_CU_SKIP                 0x6
#define F265_CO_PRED_MODE               0x9
#define F265_CO_PART_MODE               0xa
#define F265_CO_INTRA_LUMA_PRED         0xe
#define F265_CO_INTRA_CHROMA_PRED       0xf
#define F265_CO_MERGE_FLAG              0x10
#define F265_CO_MERGE_IDX               0x11
#define F265_CO_INTER_PRED              0x12
#define F265_CO_REF_IDX                 0x17
#define F265_CO_MVP                     0x19
#define F265_CO_SPLIT_TRANSFORM         0x1a
#define F265_CO_ROOT_CBF                0x1d
#define F265_CO_CBF_LUMA                0x1e
#define F265_CO_CBF_CHROMA              0x20
#define F265_CO_MVD_GREATER0            0x24
#define F265_CO_MVD_GREATER1            0x25
#define F265_CO_CU_QP_DELTA             0x26
#define F265_CO_TRANSFORM_SKIP          0x28
#define F265_CO_LAST_SIG_COEFF          0x2a
#define F265_CO_CODED_SUB_BLOCK         0x4e
#define F265_CO_SIG_COEFF               0x52
#define F265_CO_COEFF_GREATER1          0x7c
#define F265_CO_COEFF_GREATER2          0x94

// CABAC encoding modes:
// - raw: regular CABAC encoding.
// - sao: store CABAC bins for delayed encoding.
// - rdo: measure the entropy loss.
#define F265_CABAC_RAW                  1
#define F265_CABAC_SAO                  2
#define F265_CABAC_RDO                  3

// Syntax element cost approximation categories. There is an entry for each
// possible binarization of a syntax element. Notes:
// - "Split CU", "skip CU" and "split transform" have an entry for each CABAC
//   context (3).
// - "Intra luma mode" uses the MPM idx, or 3 for non-MPM.
// - "Intra chroma mode" is 0 for the luma mode, 1 for the other modes.
// - "Inter part" uses 3 entries for min-BS, 8 entries for non-min-BS.
// - "Merge idx" uses 5 entries for the merge indices, 1 for non-merge.
// - "Ref idx" has 16 entries per list.
// - Bypass is the cost of a bypass bin.
// - SE_SIZE is the table size.
#define F265_SE_SPLIT_CU                0
#define F265_SE_CU_SKIP                 6
#define F265_SE_SPLIT_TRANSFORM         12
#define F265_SE_PRED_MODE               18
#define F265_SE_INTRA_PART              20
#define F265_SE_INTER_PART              22
#define F265_SE_INTRA_LUMA_MODE         33
#define F265_SE_INTRA_CHROMA_MODE       37
#define F265_SE_INTER_PRED              39
#define F265_SE_MERGE_IDX               42
#define F265_SE_MVP                     48
#define F265_SE_REF_IDX                 50
#define F265_SE_BYPASS                  82
#define F265_SE_SIZE                    83

// Segmentation control flags.

// End the current segment after encoding the CTB.
#define F265_SC_SEGMENT_END             (1<<0)

// End the current chunk after encoding the CTB.
#define F265_SC_CHUNK_END               (1<<1)

// Initialize the CABAC contexts before encoding this CTB.
#define F265_SC_CABAC_INIT              (1<<2)

// Load the WPP CABAC contexts before encoding this CTB.
#define F265_SC_CABAC_LOAD              (1<<3)

// Save the WPP CABAC contexts after encoding this CTB.
#define F265_SC_CABAC_SAVE              (1<<4)

// Create a dependent segment (the segmentation layer knows when to create a
// segment, this merely indicates the segment type).
#define F265_SC_DEPENDENT               (1<<5)

// Unused object categories in main data thread.
#define F265_UNUSED_FRM_OBJ 0
#define F265_UNUSED_SRC_OBJ 1
#define F265_UNUSED_LA_OBJ  2
#define F265_UNUSED_ENC_OBJ 3

// Coding block flags.

// The CB is intra predicted (this includes PCM).
#define F265_CB_INTRA                   (1<<0)

// The CB is not completely outside the frame.
#define F265_CB_PRESENT                 (1<<1)

// The CB split is mandatory. This implies F265_CB_PRESENT and F265_CB_SPLIT.
#define F265_CB_FORBIDDEN               (1<<2)

// The CB is split in four child CBs.
#define F265_CB_SPLIT                   (1<<3)

// The CB uses transquant bypass.
#define F265_CB_TRANSQUANT_BYPASS       (1<<4)

// The CB uses skip (merge mode with no residual).
#define F265_CB_SKIP                    (1<<5)

// Coding block partitions.
#define F265_PART_UN                    0
#define F265_PART_H2                    1
#define F265_PART_V2                    2
#define F265_PART_HV                    3
#define F265_PART_H1                    4
#define F265_PART_H3                    5
#define F265_PART_V1                    6
#define F265_PART_V3                    7

// Index of the coding block that represents an unavailable out-of-CTB
// neighbour.
#define F265_UNAVAIL_CB_IDX             103

// Hadamard matrix dimension.
#define F265_H_DIM                      8

// Rate control modes.

// Constant QP.
#define F265_RCM_CQP                    0

// Average bitrate.
#define F265_RCM_ABR                    1

// Complexity estimation methods.

// No estimation is performed.
#define F265_EST_METHOD_BLIND           0

// The frames are subsampled. The halfpel planes are computed by filtering the
// fullpel planes. The motion estimation proceeds as usual. Weighted prediction
// can be performed inexpensively in this mode.
#define F265_EST_METHOD_SUB             1


///////////////////////////////////////////////////////////////////////////////
// Globals.

extern const int8_t f265_h8x8[F265_H_DIM][F265_H_DIM];
extern const int8_t f265_if_luma[4][8];
extern const int8_t f265_if_chroma[8][4];
extern const uint8_t f265_cabac_ctx_init_table[3][F265_NB_CABAC_CTX];
extern const uint8_t f265_cabac_range_table[64][4];
extern const uint8_t f265_cabac_transit_table[128][2];
extern const uint32_t f265_cabac_entropy_table[128];
extern const uint32_t f265_cabac_unskewed_entropy_table[128];
extern const uint8_t f265_vlc_table[64][2];
extern const int16_t f265_quant_mult[6];
extern const int16_t f265_dequant_mult[6];
extern const int8_t f265_dst_mat[4][4];
extern const int8_t f265_dct_mat_4[4][4];
extern const int8_t f265_dct_mat_8[8][8];
extern const int8_t f265_dct_mat_16[16][16];
extern const int8_t f265_dct_mat_32[32][32];
extern const uint8_t f265_scan_map_idx[4][3];
extern const uint8_t f265_scan_map_data[2*256];
extern const uint8_t f265_last_coeff_table[32];
extern const uint8_t f265_coeff_nz_ctx_table[5*16];
extern const uint8_t f265_gt1_ctx_counter_table[8];
extern const uint8_t f265_mpm_bin_table[3];
extern const uint8_t f265_chroma_mode_table[4];
extern const int8_t f265_intra_angle_table[17];
extern const int16_t f265_intra_inv_angle_table[8];
extern const uint8_t f265_intra_mode_dist[35];
extern const uint8_t f265_intra_dist_thresholds[4];
extern const int8_t f265_ref_bits[3][16];
extern const uint8_t f265_part_table[8][4][2];
extern const uint8_t f265_nb_parts[8];
extern const uint8_t f265_tc_table[54];
extern const uint8_t f265_width_to_awi[65];
extern const int16_t f265_lambdas[52];
extern uint16_t f265_mv_costs[52][F265_NB_MV_COSTS];
extern const int8_t f265_hpel_src0[16];
extern const int8_t f265_hpel_src1[16];


///////////////////////////////////////////////////////////////////////////////
// Include the auto-generated dispatch code.

#include "f265/asm.h"

#endif

