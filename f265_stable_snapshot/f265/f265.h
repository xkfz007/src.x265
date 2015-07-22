// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.


///////////////////////////////////////////////////////////////////////////////
// Documentation.
//
// This file defines the public interface of the f265 codec.
//
// There are three contexts implemented in the codec.
//
// The parameter parsing context converts a string to a set of encoder
// parameters. It is independent of the other contexts. It is initialized when
// the library is loaded. The context uses malloc() internally. The context
// aborts the process if memory is exhausted.
//
// The global context contains the global tables used by the encoder instances.
// It is initialized when the library is loaded. The context does not use
// malloc().
//
// The encoder context contains the data used by an encoder instance. The
// context is initialized when an instance is created and deinitialized when the
// instance is destroyed. The context uses malloc() and file I/O if they are
// available. The context returns an error if a resource is exhausted or an
// internal problem occurs.


///////////////////////////////////////////////////////////////////////////////
// Headers.

#ifndef F265_F265_H
#define F265_F265_H

// Configuration.
#include "f265/f265_config.h"

// Standard headers.
#ifdef F265_HAVE_STDINT
#include <stdint.h>
#endif

// Needed for HM tests. To be removed eventually.
#include "van.h"

#ifdef __cplusplus
extern "C" {
#endif


///////////////////////////////////////////////////////////////////////////////
// Constants.

// Operation return codes.
#define F265_RET_VALID                  0       // The operation completed successfully.
#define F265_RET_EMPTY                  1       // The operation completed successfully without producing a frame.
#define F265_RET_NORES                  2       // The operation aborted because a resource is not available.
#define F265_RET_ABORT                  3       // The operation aborted due to an external exit request.
#define F265_RET_ERROR                  4       // The operation aborted because an internal error occurred.

// Bitstream format.
#define F265_FORMAT_BYTE                0
#define F265_FORMAT_AVC                 1

// Multithreading mode for encoding.
#define F265_MT_ENC_NONE                0
#define F265_MT_ENC_TILE                1
#define F265_MT_ENC_WPP                 2
#define F265_MT_ENC_FRAME               3

// Frame types.
#define F265_FRAME_I                    0
#define F265_FRAME_P                    1
#define F265_FRAME_B                    2

// Search algorithm identifiers.
#define F265_ME_ALGO_DIA                0
#define F265_ME_ALGO_XDIA               1
#define F265_ME_ALGO_HEX                2


///////////////////////////////////////////////////////////////////////////////
// Data structures.

struct f265_enc; typedef struct f265_enc f265_enc;

// Encoder parameters.
typedef struct f265_enc_params
{
    // Memory requirements.

    // Memory required for the encoder instance.
    int64_t enc_mem;

    // Memory required for the bitstream of a frame.
    int32_t bs_mem;

    // True if the encoding uses high bit depths. Computed during analysis.
    int8_t hbd_flag;


    // General.

    // Clipped frame width and height (not a multiple of CB size).
    int32_t clip_dim[2];

    // Number of encoding and lookahead worker threads.
    int8_t nb_workers[2];

    // Multithreading mode for encoding.
    int8_t mt_mode;

    // Chroma format (0 => 4:0:0, 1 => 4:2:0, 2 => 4:2:2, 3 => 4:4:4).
    int8_t chroma_format;

    // Luma, chroma, PCM luma, PCM chroma bit depths.
    int8_t bit_depth[4];

    // CB minimum/maximum log size (3..6).
    int8_t cb_range[2];

    // Quantization group log size (3..6, -1: disabled).
    int8_t pcm_range[2];

    // TB minimum/maximum log size (2..5).
    int8_t tb_range[2];

    // Maximum transform depths for intra/inter.
    int8_t tb_depth[2];

    // Quantization group log size (3..6, -1: disabled).
    int8_t qg_log;

    // Number of reference frames.
    int8_t nb_refs;

    // Number of B frames.
    int8_t nb_b_frames;

    // Profile indicator. Set to 0 for automatic management.
    int8_t profile_idc;

    // Level indicator.
    int8_t level_idc;

    // Chroma QP index offset.
    int8_t chroma_qp_idx_off;

    // Deblocking filter alpha/beta offsets.
    int8_t deblock_off[2];

    // Number of merge candidates.
    uint8_t merge_cand;

    // Log 2 of the parallel merge level.
    uint8_t parallel_merge_level;

    // Specify the default bitstream format. The encoder does not use this
    // field. It is set for the convenience of the user.
    int8_t format;

    // True if wavefront parallel processing is enabled.
    int8_t wpp_flag;

    // True if the deblocking filter is enabled.
    int8_t deblock_flag;

    // True if the SAO filter is enabled.
    int8_t sao_flag;

    // Use scaling lists.
    int scaling_flag;

    // Use pulse code modulation.
    int8_t pcm_flag;

    // Use transform & quantization bypass.
    int8_t transquant_bypass_flag;

    // Use sign data hiding.
    int8_t sign_hiding_flag;

    // Use RDO for quantization.
    int8_t rdoq_flag;

    // Use transform skip.
    int8_t transform_skip_flag;

    // Use strong intra smoothing.
    int smooth_intra_flag;

    // Use assymetric motion partitions.
    int8_t amp_flag;

    // Use temporal motion vectors.
    int8_t tmv_flag;

    // True if weighted prediction is used.
    int8_t weight_flag;

    // True if chroma motion estimation is enabled.
    int8_t chroma_me_flag;

    // Temporary quality settings. Rate distortion optimization level, use HM's
    // search algorithm, test all intra modes, nullify inter TBs.
    int8_t rdo_level;
    int8_t hm_me_flag;
    int8_t all_intra_flag;
    int8_t nullify_inter_tb_flag;

    // True if the lookahead is enabled.
    int8_t la_flag;

    // True if the frame bitstream can be reformatted multiple times.
    int8_t reformat_flag;


    // Lookahead.

    // Decision and process delays.
    int16_t la_decision_delay;
    int16_t la_process_delay;


    // Rate control.

    // Bit rate in kbps. Target, min and max.
    int64_t bitrate[3];

    // Total number of frames in the stream. 0 if unknown.
    int64_t stream_nb_frames;

    // Interval between two key frames.
    int32_t key_interval;

    // Key frame type (0 => CRA, 1 => IDR).
    int32_t key_type;

    // Frame rate numerator and denominator.
    int32_t frame_rate_num;
    int32_t frame_rate_den;

    // Rate control method.
    int8_t rc_method;

    // Constant QP.
    int8_t qp;

    // Minimum frames and exponent used to compute convergence window
    // for ABR RC method.
    int32_t lt_conv_min;
    float lt_conv_exp;

    // Minimum and maximum QP.
    int8_t qp_bounds[2];

    // Allow frame/macroblock recoding.
    int8_t allow_recode_flag[2];


    // Inter.

    // Algorithm, iteration, distortion metric ID (fpel, hpel, qpel, refine).
    int8_t me_algo[4];
    int16_t me_iter[4];
    int8_t me_dist[4];

    // Candidate evaluation levels (non-weighted, weighted).
    int8_t cand_level[2];

    // Path to the reconstructed YUV frame file. NULL if not used.
    char *yuv_dump_path;

    // Path to the HM GOP file. NULL if not used.
    char *hm_gop_path;

    // Bitfield used to selectively enable experimental algorithms. The bitfield
    // flags are hardcoded and not documented on purpose.
    uint64_t algo;

} f265_enc_params;

// Input frame data.
typedef struct f265_input_frame
{
    // Frame timestamp and duration, in nano seconds.
    int64_t timestamp;
    int64_t duration;

    // Pointers to the beginning of each plane in a raw YUV frame buffer.
    void *planes[3];

    // Width in bytes of each plane.
    int32_t stride[3];

    // Bit depth of the YUV pixels.
    int8_t bit_depth;

} f265_input_frame;

// Output frame data.
typedef struct f265_output_frame
{
    // Bytestream buffer (output). The buffer must be large enough for the worst
    // case.
    uint8_t *bs;

    // Bytestream size (output).
    int32_t bs_size;

    // Frame timestamp and duration, in nano seconds (output).
    int64_t timestamp;
    int64_t duration;

    // Frame type (output).
    int8_t frame_type;

    // Bytestream format (input).
    int8_t format;

    // Number of bytes used for the NAL size (in AVC format) (input).
    int8_t nal_size_len;

    // Specify if the encoder must write an access unit delimiter (input).
    int8_t aud_flag;

    // Specify if the encoder must write a VPS NAL (input).
    int8_t vps_flag;

    // Specify if the encoder must write a PPS NAL (input).
    int8_t pps_flag;

    // Specify if the encoder must write a SPS NAL (input).
    int8_t sps_flag;

} f265_output_frame;

// An encoder request instructs the encoder what to do and how to present it.
typedef struct f265_enc_req
{
    // Input frame, owned by caller, if NULL triggers a flush of the encoder.
    f265_input_frame *input;

    // Output frame array, owned by caller. Each entry corresponds to a
    // formatting request of the current frame being outputted, if any.
    f265_output_frame *outputs;

    // Output frame array size.
    int32_t nb_outputs;

    // Force the frame to be a (IDR) key frame.
    int8_t force_key_flag;
    int8_t force_idr_flag;

    // Error string. Non-null if the return code is NORES or ERROR, null
    // otherwise.
    char *error_string;

} f265_enc_req;


///////////////////////////////////////////////////////////////////////////////
// Functions.

// Initialize the parameter parsing context.
void f265_init_parsing();

// Initialize the global context. If no_asm_flag is true, disable the assembly
// support even if it is available.
void f265_init_global(int no_asm_flag);

// Parse the encoder parameter string specified. A malloc'ed error string is set
// if an error occurs during parsing. The error string must be freed by the
// caller. This function calls f265_set_default_params() on entry. This function
// calls f265_free_params() on failure.
void f265_parse_params(f265_enc_params *params, const char *param_string, char **error_handle);

// Free the parameter strings allocated during parsing. This function can be
// called multiple times as long as the parameters have been zeroed initially,
// which f265_parse_params() does on entry.
void f265_free_params(f265_enc_params *params);

// Set the default parameter values.
void f265_set_default_params(f265_enc_params *params);

// Normalize the parameters of the encoder. This must be called unless it is
// certain that the codec parameters are valid.
void f265_normalize_params(f265_enc_params *params);

// Compute the memory requirements of the encoder instance.
void f265_analyze_params(f265_enc_params *params);

// Initialize an encoder instance with the parameters specified. The buffer
// provided must be large enough to store the instance data. If a resource is
// exhausted, set the error string handle to a static string, otherwise set the
// encoder handle.
void f265_init_enc(f265_enc_params *params, uint8_t *buf, f265_enc **enc_handle, char **error_handle);

// Deinitialize an encoder instance if it is non-null.
void f265_deinit_enc(f265_enc *enc);

// Process an encoder request. Return one of VALID, EMPTY, NORES, ABORT, ERROR.
int f265_process_enc_req(f265_enc *enc, f265_enc_req *req);

// Request an encoder instance to abort its processing as soon as possible. This
// can be called from any thread as long as the encoder instance is guaranteed
// to be initialized (even if another thread is currently processing a request).
void f265_abort_enc(f265_enc *enc);

#ifdef __cplusplus
}
#endif

#endif

