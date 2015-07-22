// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include "f265/enc.h"

// Check the encoding result and adjust the rate control accordingly. Return
// VALID, RETRY or ERROR.
int venc_rc_check_frame(f265_enc_thread *t)
{
    f265_enc *enc = t->enc;
    f265_rate_control *rc = &t->enc->md.rc;
    f265_frame *f = t->src_frame;
    int res = F265_RET_VALID;

    // Get output frame.
    int of_idx = 0;
    f265_output_frame *of = enc->md.req->outputs + of_idx;

    // Update consumed bits.
    f->actual_bits = of->bs_size << 3;

    // Update ABR rate control model.
    if (rc->method == F265_RCM_ABR) venc_rc_frame_end(t, f->actual_bits, f->qp);

    #ifdef VAN_DEBUG_RC
    // FIXME. Debug code
    double avrg = (double)rc->past_enc_bits_sum / enc->md.nb_past_frames;
    avrg *= rc->frame_rate;

    MC_PRINTF(" Encoded frames %ld", enc->md.nb_past_frames);
    MC_PRINTF(" Frame rate %f", rc->frame_rate);
    MC_PRINTF(" Frame qp %d", f->qp);
    MC_PRINTF(" Frame size in bits %d  Total %ld (%ld Bytes)", f->actual_bits, rc->past_enc_bits_sum, rc->past_enc_bits_sum >> 3);
    MC_PRINTF(" Delta bits %ld", rc->past_enc_bits_sum - rc->past_exp_bits_sum);
    MC_PRINTF(" Average bitrate %.2f Kb/s - Target %.2f Kb/s", avrg/1000.0, rc->bitrate/1000.0);
    // End of debug.
    #endif

    return res;
}

// The complexity doubles every 6 QP. The QS (quantizer scale factor) represents
// the complexity factor associated to a QP. The complexity of a macroblock can
// be expressed in qbits (quantizer-scaled bits) given the number of bits and
// the QP used to encode it: qbits = QS * bits. This is the basis used to
// compare the complexity of macroblocks encoded with different QPs.
static inline float venc_rc_qp_to_qs(float qp)
{
    return 0.85 * powf(2.0, (qp-12.0)/6.0);
}

static inline float venc_rc_qs_to_qp(float qs)
{
    return 12.0 + 6.0*log2f(qs/0.85);
}

// FIXME. Straight import from v264. Check if it needs adaptation for HEVC context.
static inline int8_t venc_rc_init_qp(f265_enc *enc, f265_enc_params *params)
{
    f265_rate_control *rc = &enc->md.rc;

    double x0 = 0, y0 = 1.19, x1 = 1.75, y1 = 1.75;
    double fs;

    // Init with luma size (width * height).
    fs = enc->gd.clip_dim[0] * enc->gd.clip_dim[1];

    // Add chroma.
    fs += fs / 2;

    int8_t qp = 1. / 1.2 * pow(10., (log10(fs * 2. / 3. * rc->frame_rate / (double)rc->bitrate) - x0)
              * (y1 - y0) / (x1 - x0) + y0) + 0.5;
    qp = F265_CLAMP(qp, rc->qp_bounds[0], rc->qp_bounds[1]);

    return qp;
}

int64_t venc_rc_get_bitrate(f265_enc *enc)
{
    return enc->md.rc.bitrate;
}

// Set bit rate for dynamic bit rate adaptation.
int64_t venc_rc_set_bitrate(f265_enc *enc, int64_t bitrate)
{
    f265_rate_control *rc = &enc->md.rc;
    bitrate = (bitrate > rc->dbra_range[0] ? bitrate : rc->dbra_range[0]);
    bitrate = (bitrate < rc->dbra_range[1] ? bitrate : rc->dbra_range[1]);
    rc->bitrate = bitrate;
    return bitrate;
}

// Stream level initialization.
void venc_rc_init_stream(f265_enc *enc, f265_enc_params *params)
{
    #ifdef VAN_DEBUG_RC
    char *rcm []= {"cqp", "abr"};
    MC_PRINTF("RC method: %s", rcm[params->rc_method]);
    #endif

    f265_rate_control *rc = &enc->md.rc;

    // Set the RC method.
    rc->method = params->rc_method;

    // Set the frame rate.
    rc->frame_rate = (double)params->frame_rate_num / params->frame_rate_den;

    // init_qp will be used as the constant QP.
    if (rc->method == F265_RCM_CQP)
    {
        enc->gd.init_qp = params->qp;
        return;
    }

    // FIXME. Set the estimation method to blind.
    rc->est_method = F265_EST_METHOD_BLIND;

    rc->past_enc_bits_sum = 0;
    rc->past_exp_bits_sum = 0;
    for (int i = 0; i < 3; i++)
    {
        rc->est_cmpl_sum[i] = 0;
        rc->est_cmpl_count[i] = 0;
    }
    rc->est_cmpl_decay = 0;
    rc->dec_qs_sum = 0;
    rc->dec_qs_count = 0;
    rc->dec_qs_decay = 0;
    rc->qs_sum = 0;
    rc->qs_count = 0;
    rc->next_frames = 1;

    rc->stream_nb_frames = params->stream_nb_frames;
    rc->key_interval = params->key_interval;
    rc->key_countdown = params->key_interval;
    rc->key_idr_flag = 1;
    rc->bitrate = params->bitrate[0]*1000;

    // Adjust dbra parameters.
    {
        int64_t min = params->bitrate[1]*1000;
        int64_t max = params->bitrate[2]*1000;

        if (min < 0)
            min = (rc->bitrate >> 1);
        else if ((min == 0) | (min > rc->bitrate))
            min = rc->bitrate;

        if (max < 0)
           max = (rc->bitrate << 1);
        else if ((max == 0) | (max < rc->bitrate))
            max = rc->bitrate;

        rc->dbra_range[0] = min;
        rc->dbra_range[1] = max;
    }

    rc->lt_conv_exp = params->lt_conv_exp;
    rc->lt_conv_min = params->lt_conv_min;
    rc->qp_bounds[0] = params->qp_bounds[0];
    rc->qp_bounds[1] = params->qp_bounds[1];

    uint32_t frame_ctb = enc->gd.nb_ctb;

    // FIXME. Training data, untweaked. Need ranges for different bits/MB
    // ratios.
    uint32_t test_bits = 14388200;
    uint32_t test_frames = 501;
    uint32_t test_ctb_per_frame = 40 * 30;
    uint32_t test_qp = 30;
    uint32_t test_ctb = test_frames * test_ctb_per_frame;
    float test_qs = venc_rc_qp_to_qs(test_qp);
    uint32_t test_qbits = test_qs * test_bits;
    double test_qbits_per_ctb = (double)test_qbits / (double)test_ctb;

    // Map the training data to the actual video data.
    double frame_bits = params->bitrate[0]*1000 / rc->frame_rate;
    double frame_qbits = test_qbits_per_ctb * (double)frame_ctb;
    double frame_qs = frame_qbits / frame_bits;

    // Initialize the rate-control variables.
    // FIXME: untweaked.
    for (int i = 0; i < 3; i++)
    {
        // Set the expected QS.
        uint64_t init_sample_count = 2;
        rc->qs_count = init_sample_count;
        rc->qs_sum = frame_qs * rc->qs_count;

        // Blind mode: enc_qbits.
        // Pre-estimation mode: enc_qbits/pre_est_cost.
        rc->est_cmpl_count[i] = init_sample_count;
        rc->est_cmpl_sum[i] = frame_qbits * rc->est_cmpl_count[i];
    }

    // FIXME : for reduction mode ?
    rc->est_cmpl_decay = 0.5;
    rc->qs_factor_sum = 1.0 / rc->reduction_mult;
    rc->qs_factor_count = 1;
    rc->qs_factor_decay = 0.96;

    // Set the initial QP.
    enc->gd.init_qp = venc_rc_init_qp(enc, params);

    return;
}

static inline void venc_rc_set_frame_qs(f265_frame *frame, float qs)
{
    frame->qs = F265_CLAMP(qs, 0.2, 76.9);
}

// Get the long-term convergence window size in bits and frames, and the number
// of future frames in the long-term window.
static void venc_rc_get_long_term_window(f265_enc_thread *t, int64_t *win_bits, int32_t *win_frames, int32_t *future_frames)
{
    f265_enc *enc = t->enc;
    f265_main_data *md = &enc->md;
    f265_rate_control *rc = &md->rc;

    int64_t min_frames = rc->lt_conv_min;
    double expo = rc->lt_conv_exp;
    int64_t stream_nb_frames = rc->stream_nb_frames;

    // FIXME. Check these.
    int64_t nb_past_frames = md->nb_past_frames;
    int32_t thread_frames = 0;

    // Compute the number of frames in the long-term convergence window. The
    // length is smaller on each end of the video sequence and higher in the
    // middle.
    int64_t base_frames = nb_past_frames + thread_frames;
    if (stream_nb_frames && nb_past_frames > stream_nb_frames/2 && nb_past_frames <= stream_nb_frames)
        base_frames = rc->stream_nb_frames - nb_past_frames;
    *win_frames = min_frames + pow(base_frames, expo);
    *win_bits = *win_frames * rc->bitrate / rc->frame_rate;

    // Cap the convergence window length to the number of future frames we have
    // to respect the convergence constraints.
    *future_frames = F265_MIN(rc->next_frames, *win_frames);
}

// Get the actual and expected committed bits in the long term window.
static void venc_rc_get_long_term_committed_bits(f265_enc_thread *t, int64_t *actual, int64_t *expected)
{
    f265_enc *enc = t->enc;
    f265_main_data *md = &enc->md;
    f265_rate_control *rc = &md->rc;

    // Add the bits from the past frames.
    *actual = rc->past_enc_bits_sum;
    *expected = rc->past_exp_bits_sum;
}

// Allocate bits to the future frames using the long-term window.
static void venc_rc_handle_long_term(f265_enc_thread *t)
{
    f265_enc *enc = t->enc;
    f265_main_data *md = &enc->md;
    f265_rate_control *rc = &md->rc;

    // Get the long-term window information.
    int64_t lt_win_bits;
    int32_t lt_win_frames;
    int32_t lt_future_frames;
    int64_t lt_committed_actual;
    int64_t lt_committed_expected;
    venc_rc_get_long_term_window(t,  &lt_win_bits, &lt_win_frames, &lt_future_frames);
    venc_rc_get_long_term_committed_bits(t,  &lt_committed_actual, &lt_committed_expected);

    // Update the future frame statistics.
    int64_t lt_future_actual = 0;
    int64_t lt_future_expected = 0;

    // Iterate on future frames that was already analyzed (look ahead).
    lt_future_frames = 1;
    for (int32_t i = 0; i < lt_future_frames; i++)
    {
        // FIXME when lookahead is functional, frame should point to every analyzed frame,
        // for now, it points on the current frame and lt_future_frames is forced to 1.

        f265_frame *frame = t->src_frame;

        // On average bitrate, use the average QS for that frame type.
        double base_qs = rc->qs_sum / (double)rc->qs_count;
        double est_cmpl_sum = rc->est_cmpl_sum[frame->frame_type];
        double est_cmpl_count = rc->est_cmpl_count[frame->frame_type];
        venc_rc_set_frame_qs(frame, base_qs);
        frame->qbits = est_cmpl_sum/est_cmpl_count * frame->est_cost;
        frame->actual_bits = frame->qbits / frame->qs;

        lt_future_actual += frame->actual_bits;
        lt_future_expected += frame->expected_bits;
    }

    // Compute how far off we are from the target bitrate according to the past
    // and the estimates for the in-flight and future frames. We increase the
    // weight of the future frames so that it has more effect on the QS issued.
    // This affects the trade-off between the average quality and the variance
    // of the quality in the case of average bitrate.
    int64_t lt_committed_delta = lt_committed_actual - lt_committed_expected;
    int64_t lt_future_delta = lt_future_actual - lt_future_expected;
    int64_t lt_bits_delta = lt_committed_delta + lt_future_delta*((double)lt_win_frames/4.0);

    // Compute the long-term correction factor.
    double lt_correction_factor = F265_CLAMP(1.0 + (double)lt_bits_delta / (double)lt_win_bits, 0.01, 100.0);

    // Apply the long-term correction.
    for (int32_t i = 0; i < lt_future_frames; i++)
    {
        f265_frame *frame = t->src_frame;
        venc_rc_set_frame_qs(frame, frame->qs*lt_correction_factor);
        frame->actual_bits = frame->qbits / frame->qs;
    }
}

// Compute the frame QP.
int8_t venc_rc_frame_start(f265_enc_thread *t, f265_frame *prev)
{
    f265_rate_control *rc = &t->enc->md.rc;
    f265_frame *frame = t->src_frame;

    // Set the expected number of bits for the current frame.
    frame->expected_bits = (double)t->src_frame->duration / 1000000000 * rc->bitrate;

    if (rc->method == F265_RCM_CQP)
    {
        return t->enc->gd.init_qp;
    }

    if (rc->method == F265_RCM_ABR)
    {
        // If it's the first encoded frame, use the initial QP.
        if (!t->enc->md.nb_past_frames) return t->enc->gd.init_qp;

        // FIXME. Look-ahead based complexity estimation not implemented yet. Blind is assumed.
        frame->est_cost = 1;

        // Allocate bits to the future frames.
        venc_rc_handle_long_term(t);

        // Get the QP corresponding to the QS.
        int8_t qp = venc_rc_qs_to_qp(frame->qs) + 0.5;

        // Clamp the QP.
        qp = F265_CLAMP(qp, rc->qp_bounds[0], rc->qp_bounds[1]);

        // Temporary fix until intra frame support improves: use the QP of the last
        // frame for key frames.
        if (t->enc->md.nb_past_frames && frame->frame_type == F265_FRAME_I) qp = prev->qp;

        return qp;
    }

    // Invalid rate control method. Should not get here.
    return -1;
}

void venc_rc_frame_end(f265_enc_thread *t, int32_t actual_bits, float avg_qp)
{
    // Get the current frame object.
    f265_enc *enc = t->enc;
    f265_rate_control *rc = &enc->md.rc;
    f265_frame *frame = t->src_frame;

    // Update the frame statistics.
    venc_rc_set_frame_qs(frame, venc_rc_qp_to_qs(avg_qp));
    frame->actual_bits = actual_bits;

    frame->qbits = frame->qs * frame->actual_bits;

    // Update the past.
    rc->past_enc_bits_sum += frame->actual_bits;
    rc->past_exp_bits_sum += frame->expected_bits;
    rc->qs_sum += frame->qs;
    rc->qs_count++;

    // On average bitrate, update the estimation cost model.
    rc->est_cmpl_sum[frame->frame_type] *= rc->est_cmpl_decay;
    rc->est_cmpl_count[frame->frame_type] *= rc->est_cmpl_decay;
    rc->est_cmpl_sum[frame->frame_type] += (double)frame->qbits / frame->est_cost;
    rc->est_cmpl_count[frame->frame_type]++;

    // FIXME is it the right place to increment this ?
    enc->md.nb_past_frames++;
}
