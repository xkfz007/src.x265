// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// Lookahead frame analysis.

#include "enc.h"

// Pop a frame from a frame queue.
static inline f265_frame* venc_frame_queue_pop(f265_frame **queue, int16_t *queue_size)
{
    assert(*queue_size);
    f265_frame *f = queue[0];
    memmove(queue, queue + 1, --(*queue_size) * sizeof(void*));
    return f;
}

// Prepare the frame for insertion in the lookahead.
f265_frame* venc_la_make_frame(f265_enc *e)
{
    f265_gen_data *gd = &e->gd;
    f265_main_data *md = &e->md;
    f265_enc_req *req = md->req;

    // Link the frame, source and lookahead objects.
    f265_frame *f = (f265_frame*)venc_get_unused_obj(e, F265_UNUSED_FRM_OBJ);
    f->mem_buf[0] = (uint8_t*)f;
    venc_link_src_obj(e, f);
    f->mem_buf[3] = NULL;

    // Reset the frame flags.
    f->main_flags = f->gen_flags = 0;

    // The frame is used for reference in the lookahead.
    f->main_flags |= F265_FF_LA_REF;

    // Set the absolute picture order count.
    f->abs_poc = md->abs_poc_counter++;

    // Set the timestamp and the duration.
    f->timestamp = req->input->timestamp;
    f->duration = req->input->duration;

    // Force key frame if requested.
    F265_SET_FLAG(f->gen_flags, F265_FF_KEY, req->force_key_flag|req->force_idr_flag);
    F265_SET_FLAG(f->gen_flags, F265_FF_IDR, req->force_idr_flag);

    // Force IDR if this is the first frame or we're overflowing the picture
    // order count.
    if (!f->abs_poc || md->enc_frames[0] && md->enc_frames[0]->h265_poc >= 2000000000)
    {
        F265_SET_FLAG(f->gen_flags, F265_FF_KEY|F265_FF_IDR, 1);
    }

    // FIXME: assuming 4:2:0 data for now. We need to figure out how we want to
    // handle the format conversions. In particular, we may need to increase the
    // bit depth of the input.

    // Get the plane information. The luma plane needs to be padded with an
    // extra pixel if the lookahead is used (subsampled planes + AQ intra
    // prediction).
    int32_t stride = gd->stride;
    f265_pix **planes = (f265_pix**)f->src_planes;
    int32_t *pix_dim = gd->pix_dim;
    int32_t *clip_dim = gd->clip_dim;
    int32_t pad_luma[4] = { 0, pix_dim[0] - clip_dim[0], 0, pix_dim[1] - clip_dim[1] };
    int32_t pad_chroma[2] = { pad_luma[1]>>1, pad_luma[3]>>1 };
    for (int i = 0; i < 4; i++) pad_luma[i] += F265_GET_FLAG(gd->eflags, F265_PF_LA);
    for (int i = 0; i < 2; i++) pad_luma[i] = F265_ALIGN_VAL(pad_luma[i], 8);
    pad_chroma[0] = F265_ALIGN_VAL(pad_chroma[0], 8);

    // Copy the plane data.
    for (int i = 0; i < 3; i++)
    {
        int32_t in_stride = req->input->stride[i];
        uint8_t *in_plane = req->input->planes[i];
        f265_pix *p = planes[i];
        for (int32_t j = 0; j < clip_dim[1]>>!!i; j++, in_plane += in_stride, p += stride)
            memcpy(p, in_plane, clip_dim[0]>>!!i);
    }

    // Pad the planes.
    venc_pad_plane(planes[0], stride, 0, clip_dim[0] - 1, 0, clip_dim[1] - 1,
                   pad_luma[0], pad_luma[1], pad_luma[2], pad_luma[3]);
    for (int i = 0; i < 2; i++)
        venc_pad_plane(planes[1+i], stride, 0, (clip_dim[0] - 1)>>1, 0, (clip_dim[1] - 1)>>1,
                       0, pad_chroma[0], 0, pad_chroma[1]);

    return f;
}

// Process the received frame in regular lookahead. Return the output frame, if
// any. We process the buffered frames as soon as possible. For real-time
// processing this can help to distribute the CPU load evenly.
f265_frame* venc_la_process_frame_regular(f265_enc *enc, f265_frame *in)
{
    f265_lookahead *la = &enc->la;

    // Buffer the input frame.
    if (in) la->display[1 + la->nb_undecided++] = in;

    // Process the buffered frames.
    if (la->nb_undecided)
    {
        // Current frame to decide.
        f265_frame *fc = la->display[1];

        // Minimum/maximum key frame interval. Arbitrary.
        int keyint_max = la->key_interval;
        int keyint_min = F265_MAX(F265_MIN(keyint_max>>3, 30), 1);

        // Number of frames between fc and the next key frame.
        int next_key_dist = la->key_countdown - 1;

        // Number of frames between fc and the previous key frame.
        int prev_key_dist = la->key_interval - la->key_countdown + 1;

        // Size of the current group of pictures (P or I frame followed by 0 or
        // more B frames).
        int nb_gop = 0;

        // True if the leading frame in the current GOP is an intra/key frame.
        int i_flag = 0, key_flag = 0;

        // Key frame type, if the leading frame is used as a key frame.
        int key_type = la->key_type;

        // Import the key frame status and type (if it is already set) from the
        // external interface.
        if (F265_GET_FLAG(fc->gen_flags, F265_FF_KEY))
        {
            i_flag = key_flag = 1;
            if (F265_GET_FLAG(fc->gen_flags, F265_FF_CRA)) key_type = 0;
            else if (F265_GET_FLAG(fc->gen_flags, F265_FF_IDR)) key_type = 1;
        }

        // Enforce the key frame interval for the current frame.
        if (!next_key_dist) i_flag = key_flag = 1;

        // Enforce use of intra frames.
        i_flag |= !enc->gd.nb_refs;

        // Use the lookahead (future support).

        // Convert I to key frame if the GOP size is higher than the minimum.
        // FIXME: enable when we don't need HM compatibility.
        #if 0
        key_flag |= (i_flag && prev_key_dist >= keyint_min);
        #else
        if (prev_key_dist + keyint_min) {}
        #endif

        // Compute the number of B frames naively. Enforce the B frame limit,
        // the key frame interval limit, the frame type and the key frame flags.
        // We include the key frame if the GOP is open and there are no
        // interfering key frames. Should probably be done earlier eventually.
        int open_gop_flag = !key_flag && key_type == 0;
        int nb_seq_b = (!i_flag)*F265_MIN(F265_MIN(la->nb_undecided - 1, la->nb_b_frames),
                                          next_key_dist + open_gop_flag - 1);
        for (int i = 0; i < nb_seq_b; i++)
        {
            if (la->display[1+i]->gen_flags&F265_FF_KEY)
            {
                nb_seq_b = i;
                break;
            }
        }
        nb_gop = nb_seq_b + 1;
        if (nb_seq_b == next_key_dist) key_flag = 1;

        // Write the GOP.
        f265_frame **p = la->coded + la->nb_committed;

        // Leading GOP frame.
        p[0] = la->display[1+nb_seq_b];
        p[0]->la_frame_type = i_flag ? F265_FRAME_I : F265_FRAME_P;
        F265_SET_FLAG(p[0]->gen_flags, F265_FF_REF, !!enc->gd.nb_refs);

        // B frames.
        for (int i = 0; i < nb_seq_b; i++)
        {
            p[1+i] = la->display[1+i];
            p[1+i]->la_frame_type = F265_FRAME_B;
            F265_SET_FLAG(p[1+i]->gen_flags, F265_FF_REF, 0);
        }

        // Commit the current GOP if there are no committed frames and the
        // lookahead is full or flushing.
        if (!la->nb_committed && (la->nb_undecided == la->decision_delay + 1 || !in))
        {
            // Reorder the HM GOP. This is a kludge that can force the last
            // frame of the GOP to become a key frame.
            if (enc->gd.hm_gop_size)
            {
                venc_hm_reorder_gop(enc, &nb_gop, &key_flag);
                i_flag |= key_flag;
            }

            // Commit the key frame status of the leading GOP frame.
            if (key_flag)
            {
                f265_frame *f = p[0];
                f->la_frame_type = F265_FRAME_I;
                F265_SET_FLAG(f->gen_flags, F265_FF_KEY, 1);
                if (key_type == 0) F265_SET_FLAG(f->gen_flags, F265_FF_CRA, 1);
                else if (key_type == 1) F265_SET_FLAG(f->gen_flags, F265_FF_IDR, 1);
            }

            // Update the key frame counter.
            if (key_flag) la->key_countdown = la->key_interval;
            else la->key_countdown -= nb_gop;

            // Remove the committed frames from the display array.
            la->nb_committed = nb_gop;
            la->nb_undecided -= nb_gop;
            memmove(la->display + 1, la->display + 1 + la->nb_committed, la->nb_undecided*sizeof(void*));

            // Update the unification frame.
            la->display[0] = p[0];
        }
    }

    // Output frame.
    f265_frame *out = NULL;

    // Retrieve the output frame from the coded array if there are committed
    // frames.
    if (la->nb_committed)
    {
        out = venc_frame_queue_pop(la->coded, &la->nb_committed);
    }

    return out;
}

// This function is executed by the primary lookahead thread to analyze frames.
// The mutex is locked on entry and on exit.
void venc_la_process_frames(f265_enc *enc)
{
    f265_lookahead *la = &enc->la;

    // Loop until all the frames were processed.
    while (la->in_queue_len)
    {
        // Pop the frame from the input queue.
        f265_frame *in = venc_frame_queue_pop(la->in_queue, &la->in_queue_len);

        // Compute the subsampled planes. FIXME.
        #if 0
        if (in && (gd->eflags&F265_PF_LA))
        {
            // Interpolate the planes.
            venc_compute_subsampled(in->sub_planes, in->src_planes[0], gd->stride, gd->mb_size[0], gd->mb_size[1]);

            // Pad the subsampled planes.
            int32_t br = (gd->pix_size[0]>>1) - 1, bb = (gd->pix_size[1]>>1) - 1;
            int32_t pad = F265_LUMA_PLANE_PADDING>>1;
            for (int i = 0; i < 4; i++) venc_pad_plane(in->sub_planes[i], gd->stride, 0, br, 0, bb, pad, pad, pad, pad);
        }
        #endif

        // Process the frame and get the output frame, which can be null.
        f265_frame *out = venc_la_process_frame_regular(enc, in);

        // Add the frame in the output queue and notify the main thread.
        la->out_queue[la->out_queue_len++] = out;
    }
}

// Add a frame to the lookahead and return the next frame to encode. The input
// frame is null when the lookahead is being flushed. The output frame is null
// when the frames are being buffered.
f265_frame* venc_la_add_frame(f265_enc *e, f265_frame *in)
{
    f265_gen_data *gd = &e->gd;
    f265_main_data *md = &e->md;
    f265_lookahead *la = &e->la;
    int mt_flag = e->gd.nb_workers[1];

    // Add the frame to the input queue.
    la->in_queue[la->in_queue_len++] = in;
    if (in) md->nb_la_frames++;

    // Run the lookahead.
    if (!mt_flag) venc_la_process_frames(e);

    // Read a frame from the output queue. The frame can be null.
    f265_frame *out = venc_frame_queue_pop(la->out_queue, &la->out_queue_len);
    if (out) md->nb_la_frames--;

    // No output frame.
    if (!out) return NULL;

    // If the frame is not a B frame, then it might still be used for reference
    // in the lookahead as the last frame of the last group of pictures.
    // FIXME: revisit that assumption, set flag status in lookahead.
    if ((gd->eflags&F265_PF_LA) && out->la_frame_type != F265_FRAME_B)
    {
        // Clear the previous frame used for lookahead reference. The memory is
        // only released if the frame was encoded.
        if (md->la_ref_frame)
        {
            F265_SET_FLAG(md->la_ref_frame->main_flags, F265_FF_LA_REF, 0);
            if (md->la_ref_frame->main_flags&F265_FF_ENC) venc_release_frame_mem(e, md->la_ref_frame);
        }
        md->la_ref_frame = out;
    }

    // Mark the frame as no longer used for reference by the lookahead.
    else
    {
        F265_SET_FLAG(out->main_flags, F265_FF_LA_REF, 0);
    }

    return out;
}

