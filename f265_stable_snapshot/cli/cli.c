// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <getopt.h>
#include <ctype.h>
#include "f265/f265.h"

#ifdef F265_PERF_BENCHMARK
#include "f265/perf.h"
#endif

#ifndef F265_NO_LIBAV
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#endif

/* Represent the data of an input file. */
typedef struct input_file_data
{
    /* Path to the file. */
    char *path;

    /* True if the file is a YUV file. */
    int yuv_flag;

    /* For a raw YUV file, the file handle. */
    FILE *file_handle;

    /* Size in bytes of the YUV input */
    int64_t file_size;

    /* Number of frames in the YUV input based on file size and frame dimensions */
    int64_t frame_count;

    #ifndef F265_NO_LIBAV
    /* For an non-YUV file, the AV library frame data. */
    AVFormatContext *format_ctx;
    AVCodecContext *codec_ctx;
    int stream_id;
    AVStream *stream;
    AVCodec *codec;
    AVFrame *frame;
    #endif

} input_file_data;

/* Represent the data of an output file. */
typedef struct output_file_data
{
    /* Path to the file. */
    char *path;

    /* File handle. */
    FILE *file_handle;

} output_file_data;

/* Represent the data of the entire video. */
typedef struct video_data
{
    /* Disable printing various stats */
    int no_stats_flag;

    /* True if we're verbose. */
    int verbose_flag;

    /* True to disable assembly support even if it is available. */
    int no_asm_flag;

    /* Frame dimensions per component. */
    uint32_t width[3];
    uint32_t height[3];

    /* Frame rate. */
    float frame_rate;

    /* Input file data. */
    input_file_data ifd;

    /* Pointer to the planes of the current frame of the input file. */
    uint8_t *planes[3];

    /* Encoder object. */
    f265_enc* enc;
    uint8_t* enc_buf;

    /* Encoder parameters. */
    f265_enc_params enc_params;

    /* Encoder input and output frames. */
    f265_input_frame enc_input_frame;
    f265_output_frame enc_output_frame;

    /* Encoder command context */
    f265_enc_req enc_req;

    /* Output file data. */
    output_file_data ofd;

    /* Number of frames to encode. 0 for the end of file. */
    uint32_t nb_frame_to_encode;

    /* Codec parameters. */
    char *codec_params;

    /* Quality and Performance Statistics */

    /* Accumulator for the time passed in the encoder */
    int64_t total_encoding_time;

    /* Longest encoding time for a single frame */
    int64_t max_encoding_time;

    /* Accumulator for the number of frames encoded */
    int64_t total_frame_count;

    /* Accumulator for all output buffer sizes */
    uint64_t total_data_size;

    /* Accumulator for all filler used in OTFMP4 */
    uint64_t total_filler_size;

    /* Timestamp of last update */
    uint64_t last_progress_update;

} video_data;

/* Malloc a block of zeroed bytes. */
void* zalloc(size_t size)
{
    void *data = calloc(1, size);
    if (!data)
    {
        printf("OOM\n");
        exit(1);
    }

    return data;
}

/* Malloc a block of non-zeroes bytes. */
void* cli_malloc(size_t size)
{
    void *data = malloc(size);
    if (!data)
    {
        printf("OOM\n");
        exit(1);
    }

    return data;
}

/* Close an input video file. */
void close_input_file(input_file_data *ifd)
{
   if (ifd->file_handle) { fclose(ifd->file_handle); ifd->file_handle = NULL; }
   #ifndef F265_NO_LIBAV
   if (ifd->frame) { av_free(ifd->frame); ifd->frame = NULL; }
   if (ifd->codec_ctx) { avcodec_close(ifd->codec_ctx); ifd->codec_ctx = NULL; }
   #ifdef F265_USE_NEW_LIBAV
   if (ifd->format_ctx) { avformat_close_input(&ifd->format_ctx); ifd->format_ctx = NULL; }
   #else
   if (ifd->format_ctx) { av_close_input_file(ifd->format_ctx); ifd->format_ctx = NULL; }
   #endif
   #endif
}

/* Close an output video file. */
void close_output_file(output_file_data *ofd)
{
   if (ofd->file_handle) { fclose(ofd->file_handle); ofd->file_handle = NULL; }
}

/* Initialize the video data. */
void video_data_init(video_data *v)
{
    memset(v, 0, sizeof (video_data));
    v->codec_params = "";

    f265_output_frame *out = v->enc_req.outputs = &v->enc_output_frame;
    v->enc_req.nb_outputs = 1;
    out->aud_flag = 1;
    out->nal_size_len = 4;
    out->vps_flag = 1;
    out->sps_flag = 1;
    out->pps_flag = 1;
}

/* Clear the video data. */
void video_data_clear(video_data *v)
{
    f265_free_params(&v->enc_params);
    if (v->enc) { f265_deinit_enc(v->enc); v->enc = NULL; }
    free(v->enc_buf); v->enc_buf = NULL;
    free(v->enc_output_frame.bs); v->enc_output_frame.bs = NULL;
    close_input_file(&v->ifd);
    for (int i = 0; i < 3; i++) { free(v->planes[i]); v->planes[i] = NULL; }
    close_output_file(&v->ofd);
}

/* Compute the execution time and add it to collector */
static void enc_timer (video_data *v, struct timeval *t1, struct timeval *t2)
{
    int64_t encoding_time;
    encoding_time  = (int64_t) ((t2->tv_sec  - t1->tv_sec) * 1000000LL);
    encoding_time += (int64_t) ((t2->tv_usec - t1->tv_usec));

    v->total_encoding_time += encoding_time;

    if (v->max_encoding_time < encoding_time)
        v->max_encoding_time = encoding_time;
}

/* Process an input frame, if any, and return true if a frame was outputted. */
int process_frame(video_data *v, uint32_t nb_frame, int input_flag)
{
    f265_input_frame *in = NULL;
    f265_output_frame *out = &v->enc_output_frame;
    out->bs_size = 0;

    struct timeval t1;
    struct timeval t2;

    if (input_flag)
    {
        double dur = (1.0 / (double)v->frame_rate)*1000000000.0;
        in = &v->enc_input_frame;
        in->timestamp = (double)nb_frame * dur;
        in->duration = dur;
    }

    gettimeofday (&t1, 0);
    v->enc_req.input = in;
    int err = f265_process_enc_req(v->enc, &v->enc_req);
    int out_flag = 0;
    gettimeofday (&t2, 0);

    enc_timer(v, &t1, &t2);

    switch (err)
    {
        case F265_RET_VALID :

            v->total_data_size += out->bs_size;
            v->total_frame_count ++;

            fwrite(out->bs, out->bs_size, 1, v->ofd.file_handle);
            out_flag = 1;
            break;
        case F265_RET_EMPTY :
            break;
        case F265_RET_ABORT :
            printf("The codec returned ABORT.\n");
            video_data_clear(v);
            exit(1);
            break;
        case F265_RET_NORES :
        case F265_RET_ERROR :
            printf("Codec error: %s.\n", v->enc_req.error_string);
            video_data_clear(v);
            exit(1);
            break;
    };

    return out_flag;
}

/* Initialize the encoder. */
void init_encoder(video_data *v)
{
    uint32_t i;
    char* err;

    f265_enc_params* params = &v->enc_params;

    f265_init_global(v->no_asm_flag);
    f265_init_parsing();

    f265_parse_params(params, v->codec_params, &err);

    if (err)
    {
        printf("Parameter error: %s.\n", err);
        video_data_clear(v);
        exit(1);
    }

    params->clip_dim[0] = v->width[0];
    params->clip_dim[1] = v->height[0];

    f265_normalize_params(params);
    f265_analyze_params(params);

    v->enc_output_frame.bs = (uint8_t*)cli_malloc(params->bs_mem);
    v->enc_buf = (uint8_t*)zalloc(params->enc_mem);
    v->enc_output_frame.format = params->format;

    f265_init_enc(params, v->enc_buf, &v->enc, &err);

    if (err)
    {
        printf("Initialization error: %s.\n", err);
        video_data_clear(v);
        exit(1);
    }

    v->frame_rate = (float)params->frame_rate_num/(float)params->frame_rate_den;

    for (i = 0; i < 3; i++)
    {
        v->enc_input_frame.stride[i] = v->width[i];
        v->enc_input_frame.planes[i] = v->planes[i];
    }
}

#ifndef F265_NO_LIBAV
/* Read the next non-YUV frame data. Return true on success. */
int read_nonyuv_frame(video_data *v)
{
    input_file_data *ifd = &v->ifd;
    int frame_flag = 0;

    /* Read the packets until we get a frame, if any. */
    while (1)
    {
        AVPacket packet;
        int done_flag = 0;

        /* Read the next frame, if any. */
        int read_success = av_read_frame(ifd->format_ctx, &packet);

        /* Video frame. */
        if (read_success < 0 || packet.stream_index == ifd->stream_id)
        {
            if (read_success < 0)
            {
                packet.size = 0;
                packet.data = 0;
            }

            /* We decoded a frame. */
            #ifdef F265_USE_OLD_LIBAV
            int byte_used = avcodec_decode_video(ifd->codec_ctx, ifd->frame, &frame_flag, packet.data, packet.size);
            #elif defined F265_USE_NEW_LIBAV
            int byte_used = avcodec_decode_video2(ifd->codec_ctx, ifd->frame, &frame_flag, &packet);
            #endif
            /* We're done if the codec outputted a frame. */
            if ((byte_used > 0) | frame_flag) done_flag = !!frame_flag;

            /* Error decoding the frame. */
            else
            {
                frame_flag = 0;
                done_flag = 1;
            }
        }

        /* Free the packet. */
        av_free_packet(&packet);

        if (done_flag) break;
    }

    /* Read the frame data. */
    if (frame_flag)
    {
        /* Read each component. The YUV samples are ordered per line. */
        for (uint32_t i = 0; i < 3; i++)
        {
            uint8_t *read_plane = ifd->frame->data[i];
            uint8_t *write_plane = v->planes[i];
            uint32_t read_pitch = ifd->frame->linesize[i];
            uint32_t nb = v->width[i];

            for (uint32_t line = 0; line < v->height[i]; line++)
            {
                memcpy(write_plane, read_plane, nb);
                read_plane += read_pitch;
                write_plane += nb;
            }
        }
    }

    return frame_flag;
}
#endif

/* Read the next YUV frame data. Return true on success. */
int read_yuv_frame(video_data *v)
{
    input_file_data *ifd = &v->ifd;

    for (uint32_t i = 0; i < 3; i++)
    {
        int nb = fread(v->planes[i], v->width[i]*v->height[i], 1, ifd->file_handle);
        if (nb != 1) return 0;
    }

    return 1;
}

/* Open a YUV video file. */
void open_yuv_file(video_data *v, input_file_data *ifd)
{
    FILE *f = fopen(ifd->path, "rb");
    if (!f)
    {
        printf("Cannot open %s: %s.\n", ifd->path, strerror(errno));
        video_data_clear(v);
        exit(1);
    }
    ifd->file_handle = f;

    fseek(f, 0, SEEK_END);
    ifd->file_size = ftell(f);
    fseek(f, 0, SEEK_SET);

    ifd->frame_count = ifd->file_size / (v->width[0] * v->height[0] * 3 / 2);

    if (!v->nb_frame_to_encode)
        v->nb_frame_to_encode = ifd->frame_count;
}

#ifndef F265_NO_LIBAV
/* Open a non-YUV video file. */
void open_nonyuv_file(video_data *v, input_file_data *ifd)
{
    const char *fail = NULL;

    /* Try. */
    do
    {
        /* Open the video file. */
        #ifdef F265_USE_OLD_LIBAV
        if (av_open_input_file(&ifd->format_ctx, ifd->path, NULL, 0, NULL) != 0)
        #elif defined F265_USE_NEW_LIBAV
        if (avformat_open_input(&ifd->format_ctx, ifd->path, NULL, NULL) != 0)
        #endif
        {
            fail = "av_open_input_file() failed"; break;
        }

        /* Find the stream information. */
        #ifdef F265_USE_OLD_LIBAV
        if (av_find_stream_info(ifd->format_ctx) < 0)
        #else
        if (avformat_find_stream_info(ifd->format_ctx, NULL) < 0)
        #endif
        {
            fail = "av_find_stream_info() failed"; break;
        }

        /* Find the first video stream. */
        ifd->stream_id = -1;
        for (uint32_t i = 0; i < ifd->format_ctx->nb_streams; i++)
        {
            ifd->stream = ifd->format_ctx->streams[i];
            #ifdef VAN_USE_OLD_LIBAV
            if (ifd->stream->codec->codec_type == CODEC_TYPE_VIDEO)
            #elif defined F265_USE_NEW_LIBAV
            if (ifd->stream->codec->codec_type == AVMEDIA_TYPE_VIDEO)
            #endif
            {
                ifd->stream_id = i;
                break;
            }
        }

        if (ifd->stream_id == -1)
        {
            fail = "cannot find video stream"; break;
        }

        /* Get a pointer to the codec context and set the pixel format. */
        ifd->codec_ctx = ifd->stream->codec;
        ifd->codec_ctx->pix_fmt = PIX_FMT_YUV420P;

        /* Find the decoder for the video stream. */
        ifd->codec = avcodec_find_decoder(ifd->codec_ctx->codec_id);
        if (ifd->codec == NULL)
        {
            fail = "avcodec_find_decoder() failed"; break;
        }

        /* Open the codec. */
        #ifdef F265_USE_OLD_LIBAV
        if (avcodec_open(ifd->codec_ctx, ifd->codec) < 0)
        #elif defined F265_USE_NEW_LIBAV
        if (avcodec_open2(ifd->codec_ctx, ifd->codec, NULL) < 0)
        #endif
        {
            fail = "avcodec_open() failed"; break;
        }

        /* Allocate a video frame. */
        ifd->frame = avcodec_alloc_frame();

        /* Set the luma width and height. */
        v->width[0] = ifd->codec_ctx->width;
        v->height[0] = ifd->codec_ctx->height;
        if (!v->width[0] || !v->height[0])
        {
            fail = "bad frame dimensions"; break;
        }

    } while (0);

    if (fail)
    {
        printf("Cannot open %s: %s.\n", ifd->path, fail);
        video_data_clear(v);
        exit(1);
    }
}
#endif

/* Open an input video file. */
void open_input_file(video_data *v, input_file_data *ifd)
{
    /* Determine whether the file is an YUV file based on the extension, if -w
     * has not been specified.
     */
    if (!ifd->yuv_flag)
    {
        char ext[4] = { 0, 0, 0, 0 };
        int path_len = strlen(ifd->path);
        if (path_len >= 4) memcpy(ext, ifd->path + path_len - 4, 4);
        for (int i = 0; i < 4; i++) ext[i] = tolower(ext[i]);
        ifd->yuv_flag = !memcmp(ext, ".yuv", 4);
    }

    /* Open the file. */
    #ifndef F265_NO_LIBAV
    ifd->yuv_flag ? open_yuv_file(v, ifd) : open_nonyuv_file(v, ifd);
    #else
    if (!ifd->yuv_flag)
    {
        fprintf(stderr, "Sorry, non-YUV support not compiled in.\n");
        exit(1);
    }
    open_yuv_file(v, ifd);
    #endif
}

/* Open an output video file. */
void open_output_file(output_file_data *ofd)
{
    FILE *f = fopen(ofd->path, "wb");
    if (!f)
    {
        printf("Cannot open %s: %s.\n", ofd->path, strerror(errno));
        exit(1);
    }
    ofd->file_handle = f;
}

/* Split two strings separated by a 'x' or a colon. Return true on success. */
int split_str(char *src, char *first, char *second)
{
    char *colon = strstr(src, ":");
    if (colon == NULL) colon = strstr(src, "x");
    if (colon == NULL) return 0;
    int first_len = colon - src;
    int second_len = strlen(src) - first_len - 1;
    if (first_len == 0 || second_len == 0) return 0;
    memcpy(first, src, first_len);
    memcpy(second, colon + 1, second_len);
    first[first_len] = 0;
    second[second_len] = 0;
    return 1;
}

/* Print the usage and exit. */
void print_usage(uint32_t exit_code)
{
    printf("Usage: f265cli <input> <output>\n");
    printf("Options:\n");
    printf("  -h: print help and exit.\n");
    printf("  -V: print version and exit.\n");
    printf("  -v: show progress.\n");
    printf("  -n: disable stats output (for profiling and benchmarking).\n");
    printf("  -a: disable assembly support even if it is available.\n");
    printf("  -w <width:height>: set the YUV frame size. Force the input file to be\n");
    printf("                     processed as YUV.\n");
    printf("  -c <count>: set the number of frames to encode. Loopover on underflow.\n");
    printf("              Default to 0 for the end of file.\n");
    printf("  -p <params>: set the codec parameters.\n");
    printf("\n");
    exit(exit_code);
}

/* Set the frame dimensions. */
void set_frame_dimensions(video_data *v)
{
    /* The frame dimensions of a YUV file must be set explicitly. */
    if (!v->width[0] || !v->height[0])
    {
        printf("Please specify the YUV frame dimensions (-w).\n");
        print_usage(1);
    }

    for (uint32_t i = 0; i < 2; i++)
    {
        v->width[1+i] = v->width[0] / 2;
        v->height[1+i] = v->height[0] / 2;
    }

    /* Allocate the planes. */
    for (uint32_t i = 0; i < 3; i++)
    {
        v->planes[i] = (uint8_t*)zalloc(v->height[i]*v->width[i]);
    }
}

/* Parse the command line. */
void parse_cmd_line(int argc, char **argv, video_data *v)
{
    while (1)
    {
        int ret = getopt(argc, argv, "hVvnaw:c:p:");
        if (ret == -1) break;
        else if (ret == '?' || ret == ':') print_usage(1);
        else if (ret == 'h') print_usage(0);

        else if (ret == 'V')
        {
            printf("Version %s\n", F265_VERSION);
            exit(0);
        }

        else if (ret == 'v') v->verbose_flag = 1;
        else if (ret == 'n') v->no_stats_flag = 1;
        else if (ret == 'a') v->no_asm_flag = 1;
        else if (ret == 'w')
        {
            char width_str[100], height_str[100];
            if (!split_str(optarg, width_str, height_str))
            {
                printf("Invalid YUV frame size format.\n");
                print_usage(1);
            }
            v->ifd.yuv_flag = 1;
            v->width[0] = atoi(width_str);
            v->height[0] = atoi(height_str);
        }
        else if (ret == 'c')
        {
            v->nb_frame_to_encode = atoi(optarg);
            if (v->nb_frame_to_encode < 0)
            {
                printf("Invalid number of frames to encode.\n");
                print_usage(1);
            }
        }
        else if (ret == 'p') v->codec_params = optarg;
    }

    int nb_arg_left = argc - optind;
    if (nb_arg_left != 2) print_usage(1);

    v->ifd.path = argv[optind + 0];
    v->ofd.path = argv[optind + 1];
}

void dump_perf(video_data *v)
{
    if (v->total_frame_count && v->total_data_size && v->enc_params.frame_rate_num)
    {
        double fps;               /* encoding framerate */
        double min;               /* lowest encoding rate */
        double dur;               /* clock time elapsed while encoding */
        double br;                /* effective bitrate */
        double eff;               /* efficiency ratio */
        double len;               /* length of media in seconds */

        dur  = (double) v->total_encoding_time / 1000000.0;
        fps  = (double) v->total_frame_count / dur;
        min  = (double) 1000000.0 / v->max_encoding_time;
        br   = (double) v->total_data_size - v->total_filler_size;
        eff  = br * 100.0 / (double) v->total_data_size;
        len  = (double) v->total_frame_count * (double) v->enc_params.frame_rate_den;
        len /= (double) v->enc_params.frame_rate_num;
        br  /= len;
        br  /= 125.0;

        printf("%-15s : %ld frames, %.2f seconds\n", "Performance",
                v->total_frame_count, dur);

        printf("%-15s : avg %.1f fps, min %.1f fps\n", "Throughput",
                fps, min);

        printf("%-15s : %.2f kbps, %.0f%% useful\n", "Efficiency",
                br, eff);
    }
}

int main(int argc, char **argv)
{
    /* Initialize. */
    video_data v;
    video_data_init(&v);
    #ifndef F265_NO_LIBAV
    av_register_all();
    #endif

    /* Parse the command line. */
    parse_cmd_line(argc, argv, &v);

    /* Open the files. */
    open_input_file(&v, &v.ifd);
    open_output_file(&v.ofd);

    /* Set the frame dimensions. */
    set_frame_dimensions(&v);

    /* Start benchmarking. */
    #ifdef F265_PERF_BENCHMARK
    PERF_PROG_ENTER();
    #endif

    /* Initialize the encoder. */
    init_encoder(&v);

    if (v.verbose_flag)
        printf("\n\n");

    /* Loop until all the frames have been processed or an error occurs. */
    int read_one_flag = 0;
    for (uint32_t nb_frame = 0; nb_frame < v.nb_frame_to_encode || !v.nb_frame_to_encode;)
    {
        /* Read the current frame. */
        int done_flag = 0;
        while (1)
        {
            #ifndef F265_NO_LIBAV
            int frame_flag = v.ifd.yuv_flag ? read_yuv_frame(&v) : read_nonyuv_frame(&v);
            #else
            int frame_flag = read_yuv_frame(&v);
            #endif

            /* We got a frame. */
            if (frame_flag)
            {
                read_one_flag = 1;
                break;
            }

            /* We processed up to the end of the file. We're done. */
            else if (!v.nb_frame_to_encode)
            {
                done_flag = 1;
                break;
            }

            /* We did not read a frame yet. We're failing badly. */
            else if (!read_one_flag)
            {
                printf("Failed to read a frame from the input file.\n");
                done_flag = 1;
                break;
            }

            /* We read at least one frame, so we can loopover safely. */
            else
            {
                close_input_file(&v.ifd);
                open_input_file(&v, &v.ifd);
                read_one_flag = 0;
            }
        }
        if (done_flag) break;

        /* Process the frame. */
        process_frame(&v, nb_frame++, 1);

        if (v.verbose_flag && (v.total_encoding_time - v.last_progress_update > 250000))
        {
            double spf = ((double) v.total_encoding_time / 1000000.0) / (double) nb_frame;
            v.last_progress_update = v.total_encoding_time;

            if (v.nb_frame_to_encode)
            {
                double progress = (double) nb_frame / (double) v.nb_frame_to_encode;
                double eta = spf * (double) (v.nb_frame_to_encode - nb_frame);

                printf("\rEncoded : %d/%d (%.1f%%) @ %.ffps ETA %.2fsec  ",
                    nb_frame, v.nb_frame_to_encode, progress * 100.0, 1.0 / spf, eta);
                fflush(stdout);
            }
            else
            {
                printf("\rEncoded : %d @ %.ffps ", nb_frame, 1.0 / spf);
                fflush(stdout);
            }
        }
    }

    /* Flush the encoder output. */
    while (process_frame(&v, 0, 0)) {}

    if (v.verbose_flag)
        printf("\n\n");

    /* Finish benchmarking. */
    #ifdef F265_PERF_BENCHMARK
    PERF_PROG_LEAVE();
    #endif

    if (!v.no_stats_flag)
    {
        // FIXME: dump stats.
    }

    /* Clear the video data. */
    video_data_clear(&v);
    return 0;
}

