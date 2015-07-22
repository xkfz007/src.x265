// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

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
#include <sys/stat.h>
#include <f265/f265_config.h>

#ifndef F265_NO_LIBAV
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#endif

/* Math utilities. */
#define ABS(a) (((a) < 0) ? (-(a)) : (a))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define CLAMP(v, min, max) ((v) < (min) ? (min) : (v) > (max) ? (max) : (v))

/* Metrics. */
#define METRIC_MATCH        0
#define METRIC_SIZE         1
#define METRIC_DURATION     2
#define METRIC_PSNR         3
#define METRIC_SSIM         4
#define NB_METRICS          5

/* Represent a macroblock mismatch in an image component. */
typedef struct mismatch_data
{
    /* Macroblock number and location (in terms of macroblock). */
    uint32_t uMb;
    uint32_t xMb;
    uint32_t yMb;

    /* Component where the mismatch was detected. */
    uint32_t comp;

    /* Location of the first mismatching pixel in the macroblock (X, Y). */
    uint32_t px;
    uint32_t py;
} mismatch_data;

/* Represent the data of a frame pair. */
typedef struct frame_data
{
    /* Pointer to the components of the frame in file1 and file2. */
    uint8_t *p[3][2];

    /* Pointer to the first mismatch in the pair. */
    mismatch_data *first_mismatch;

    /* True if there is a mismatch in the frame component. */
    int comp_mismatch[3];

    /* Size of the frame. */
    uint32_t size[2];

    /* Duration of the frame. */
    double duration[2];

} frame_data;

/* Metric statistics. */
typedef struct metric_data
{
    /* True if this metric is used. */
    int active_flag;

    /* True if this metric is computed once per component. */
    int comp_flag;

    /* Name of the metric. */
    char *name;

    /* Padding before the metric name and digits. */
    int pad_name;
    int pad_digits;

    /* Format to print/dump the metric. */
    char *print_fmt;
    char *dump_fmt;

    /* File in which the metrics are dumped, if any. */
    FILE *file_handle;

    /* Values for the current frame, per component. */
    double cur[3];

    /* Running sums of values, per component. */
    double sums[3];

    /* Running sums of squared values, per component. */
    double ssums[3];

    /* Handler function for computing the metric, per frame or frame component. */
    union
    {
        double (*compute_comp)(frame_data *f, uint32_t comp, uint8_t *p[2], uint32_t width, uint32_t height,
                               uint32_t pitch);
        double (*compute_frame)(frame_data *f);
        void *compute_func;
    };

} metric_data;

/* Represent the data of an input file. */
typedef struct input_file_data
{
    /* Path to the file. */
    char *path;

    /* True if the file is a YUV file, false if it is an non-yuv file. */
    int yuv_flag;

    /* True if the file is a raw HEVC file (h265) */
    int hevc_flag;

    /* Frame dimensions, in pixels, of the luma plane. */
    uint32_t luma_dim[2];

    /* For a raw YUV file, the file handle. */
    FILE *file_handle;

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

/* Represent the data of the entire video. */
typedef struct video_data
{
    /* True if the chroma components are used. */
    int chroma_flag;

    /* True if the frame metrics are printed. */
    int print_frame_metric_flag;

    /* True if the metric average is printed. */
    int print_average_flag;

    /* True if the metric standard deviation is printed. */
    int print_std_deviation_flag;

    /* True if chroma is considered for macroblock mismatches. */
    int chroma_mismatch_flag;

    /* True if the text mode is active. */
    int text_mode_flag;

    /* True if the macroblock mismatches are printed. */
    int print_mismatch_flag;

    /* True if only the specified macroblock in the specified frame is assumed
     * to be mismatching.
     */
    int fake_mismatch_flag;

    /* Input file data. */
    input_file_data file_data[2];

    /* Directory in which the metrics are dumped. */
    char *dump_dir;

    /* Frame display dimensions (not necessarily a multiple of 16), in pixels, per component. */
    uint32_t clip_height[3];
    uint32_t clip_width[3];

    /* Frame dimensions, in pixels, per component. */
    uint32_t width[3];
    uint32_t height[3];

    /* Macroblock size, per component. */
    uint32_t mb_size[3];

    /* Faked mismatch macroblock frame and number. */
    uint32_t fake_frame;
    uint32_t fake_mb;

    /* Number of components in the frame. */
    uint32_t nb_comp;

    /* Number of frames processed. */
    uint32_t nb_frame;

    /* Field size information. */
    int total_line_size;
    int field_min_size;
    int field_sep;
    int metric_digits;

    /* True if there is a mismatch in the video component. */
    int comp_mismatch[3];

    /* Metric data array. */
    metric_data mda[NB_METRICS];

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

/* Allocate a macroblock mismatch. */
mismatch_data* mismatch_data_new()
{
    mismatch_data *m = zalloc(sizeof(mismatch_data));
    return m;
}

/* Allocate a frame. */
frame_data* frame_data_new(video_data *v)
{
    frame_data *f = zalloc(sizeof(frame_data));
    for (uint32_t i = 0; i < v->nb_comp; i++)
        for (uint32_t j = 0; j < 2; j++)
            f->p[i][j] = zalloc(v->width[i]*v->height[i]);
    return f;
}

/* Free a frame. */
void frame_data_free(frame_data *f, video_data *v)
{
    if (!f) return;
    for (uint32_t i = 0; i < v->nb_comp; i++)
        for (uint32_t j = 0; j < 2; j++)
            free(f->p[i][j]);
    free(f->first_mismatch);
    free(f);
}

/* Print the macroblock pixels specified. */
void print_mb_pix(uint8_t *p, uint32_t stride, uint32_t mb_size)
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

/* Print the macroblock pixels specified assuming terminal output. */
void print_mb_pix_nice(uint8_t *r, uint8_t *c, uint32_t stride, uint32_t mb_size)
{
    for (uint32_t i = 0; i < mb_size; i++)
    {
        uint32_t j;
        if ((i & 3) == 0)
        {
            printf ("  \e[1;30m+-------------+-------------+-------------+-------------+\e[0m");
            printf ("  \e[1;30m+-------------+-------------+-------------+-------------+\e[0m");
            printf ("\n");
        }

        printf ("  \e[1;30m|\e[0m ");
        for (j = 0; j < mb_size; j++)
        {
            uint32_t ref = r[i*stride + j];
            uint32_t pix = c[i*stride + j];

            if (ref != pix)
                printf("\e[1;37;44m%02x\e[0m ", ref);
            else
                printf("%02x ", ref);

            if ((j & 3) == 3)
            {
                printf ("\e[1;30m|\e[0m ");
            }
        }

        printf (" \e[1;30m|\e[0m ");
        for (j = 0; j < mb_size; j++)
        {
            uint32_t ref = r[i*stride + j];
            uint32_t pix = c[i*stride + j];

            if (ref != pix)
                printf("\e[1;37;41m%02x\e[0m ", pix);
            else
                printf("%02x ", pix);

            if ((j & 3) == 3)
            {
                printf ("\e[1;30m|\e[0m ");
            }
        }

        printf("\n");
    }

    printf ("  \e[1;30m+-------------+-------------+-------------+-------------+\e[0m");
    printf ("  \e[1;30m+-------------+-------------+-------------+-------------+\e[0m");
    printf ("\n");
}

/* Print the number of spaces specified. */
void print_spaces(int size)
{
    for (int i = 0; i < size; i++) printf(" ");
}

/* Print the string specified in the field specified. */
void print_string_pad(char *s, int pad_before, int pad_after)
{
    print_spaces(pad_before);
    printf("%s", s);
    print_spaces(pad_after);
}

/* Print the program header. */
void print_header(video_data *v)
{
    if (!v->print_frame_metric_flag && !v->print_average_flag && !v->print_std_deviation_flag) return;
    for (int i = 0; i < v->total_line_size; i++) printf("=");
    printf("\n");
    print_string_pad("frame", v->field_min_size - strlen("frame"), v->field_sep);
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->active_flag) continue;
        print_string_pad(md->name, md->pad_name, v->field_sep);
    }
    printf("\n");
}

/* Print the frame metrics. */
void print_frame_metrics(frame_data *f, video_data *v)
{
    if (!v->print_frame_metric_flag) return;
    print_spaces(v->field_min_size - v->metric_digits);
    printf("% 6d", v->nb_frame);
    print_spaces(v->field_sep);
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->active_flag) continue;
        int nb_comp = md->comp_flag ? v->nb_comp : 1;
        for (int comp = 0; comp < nb_comp; comp++)
        {
            print_spaces(md->pad_digits);
            if (!md->comp_flag || f->comp_mismatch[comp]) printf(md->print_fmt, md->cur[comp]);
            else printf("******");
            if (comp < nb_comp - 1) printf(" ");
        }
        print_spaces(v->field_sep);
    }
    printf("\n");
}

/* Print the program footer. */
void print_footer(video_data *v)
{
    if (!v->print_average_flag && !v->print_std_deviation_flag) return;
    for (int i = 0; i < v->total_line_size; i++) printf("=");
    printf("\n");

    if (v->print_average_flag)
    {
        printf("Average");
        print_spaces(v->field_min_size + v->field_sep - strlen("Average"));
        for (int i = 0; i < NB_METRICS; i++)
        {
            metric_data *md = v->mda + i;
            if (!md->active_flag) continue;
            int nb_comp = md->comp_flag ? v->nb_comp : 1;
            for (int comp = 0; comp < nb_comp; comp++)
            {
                print_spaces(md->pad_digits);
                double avg = 0;
                if (v->nb_frame) avg = md->sums[comp] / (double)v->nb_frame;
                if (!md->comp_flag || v->comp_mismatch[comp]) printf(md->print_fmt, avg);
                else printf("******");
                if (comp < nb_comp - 1) printf(" ");
            }
            print_spaces(v->field_sep);
        }
        printf("\n");
    }

    if (v->print_std_deviation_flag)
    {
        printf("Std dev");
        print_spaces(v->field_min_size + v->field_sep - strlen("Std dev"));
        for (int i = 0; i < NB_METRICS; i++)
        {
            metric_data *md = v->mda + i;
            if (!md->active_flag) continue;
            int nb_comp = md->comp_flag ? v->nb_comp : 1;
            for (int comp = 0; comp < nb_comp; comp++)
            {
                print_spaces(md->pad_digits);
                double dev = 0;
                if (v->nb_frame > 1)
                {
                    double avg = md->sums[comp] / (double)v->nb_frame;
                    double var = (md->ssums[comp] - md->sums[comp]*avg) / (double)(v->nb_frame - 1);
                    dev = sqrt(var);
                }
                if (!md->comp_flag || v->comp_mismatch[comp]) printf(md->print_fmt, dev);
                else printf("******");
                if (comp < nb_comp - 1) printf(" ");
            }
            print_spaces(v->field_sep);
        }
        printf("\n");
    }
}

/* Dump the frame metrics. */
void dump_frame_metrics(frame_data *f, video_data *v)
{
    if (!v->dump_dir) return;
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->active_flag) continue;
        int nb_comp = md->comp_flag ? v->nb_comp : 1;
        for (int comp = 0; comp < nb_comp; comp++)
        {
            fprintf(md->file_handle, md->dump_fmt, md->cur[comp]);
            if (comp < nb_comp - 1) fprintf(md->file_handle, " ");
        }
        fprintf(md->file_handle, "\n");
    }
}

/* Print the mismatch specified. */
void print_mb_mismatch(frame_data *f, video_data *v, mismatch_data *m)
{
    uint32_t uMb = m->uMb;
    uint32_t xMb = m->xMb;
    uint32_t yMb = m->yMb;
    uint32_t comp = m->comp;
    uint32_t px = m->px;
    uint32_t py = m->py;
    uint32_t width = v->width[comp];
    uint32_t mb_size = v->mb_size[comp];
    uint8_t *file_pix[2];
    char *comp_names[] = { "luma", "blue (U)", "red (V)" };
    char *comp_name = comp_names[comp];
    for (uint32_t i = 0; i < 2; i++) file_pix[i] = f->p[comp][i] + yMb*mb_size*width + xMb*mb_size;

    if (v->print_frame_metric_flag) printf("\n");
    printf("Macroblock mismatch in frame %d, %s component.\n", v->nb_frame, comp_name);
    printf("uMb %u, xMb %u, yMb %u, mb pixel (%d,%d):\n", uMb, xMb, yMb, px, py);

    if (v->text_mode_flag)
    {
        for (uint32_t i = 0; i < 2; i++)
        {
            printf("\n");
            printf("\tMacroblock in %s:\n", v->file_data[i].path);
            print_mb_pix(file_pix[i], width, mb_size);
        }
    }

    else
    {
        printf ("\n  %-57s  %-57s\n", v->file_data[0].path, v->file_data[1].path);
        print_mb_pix_nice(file_pix[0], file_pix[1], width, mb_size);
    }

    printf("\n");
}

/* Handle a macroblock mismatch that has just been detected. */
void handle_mb_mismatch(frame_data *f, video_data *v, uint32_t uMb, uint32_t xMb, uint32_t yMb, uint32_t comp,
                        uint32_t px, uint32_t py)
{
    /* Register the mismatch in the component. */
    f->comp_mismatch[comp] = 1;
    v->comp_mismatch[comp] = 1;

    /* Return if we're not further interested in the mismatch. */
    if (!v->print_mismatch_flag ||
        (!v->chroma_mismatch_flag && comp > 0) ||
        f->first_mismatch ||
        (v->fake_mismatch_flag && (v->fake_frame != v->nb_frame || v->fake_mb != uMb)))
        return;

    /* Remember the mismatch data. */
    mismatch_data *m = f->first_mismatch = mismatch_data_new();
    m->uMb = uMb;
    m->xMb = xMb;
    m->yMb = yMb;
    m->comp = comp;
    m->px = px;
    m->py = py;
}

/* Analyze a macroblock pair component for mismatches. */
void process_mb_pair_comp(frame_data *f, video_data *v, uint32_t uMb, uint32_t xMb, uint32_t yMb, uint32_t comp)
{
    uint32_t width = v->width[comp];
    uint32_t mb_size = v->mb_size[comp];
    uint32_t reported_flag = 0;
    uint8_t *p[2];
    for (uint32_t i = 0; i < 2; i++) p[i] = f->p[comp][i] + yMb*mb_size*width + xMb*mb_size;

    for (uint32_t py = 0; py < mb_size; py++)
    {
        for (uint32_t px = 0; px < mb_size; px++)
        {
            int p0 = p[0][py*width + px];
            int p1 = p[1][py*width + px];
            if (p0 != p1)
            {
                handle_mb_mismatch(f, v, uMb, xMb, yMb, comp, px, py);
                reported_flag = 1;
                break;
            }
        }
        if (reported_flag) break;
    }

    /* Fake a mismatch if required. */
    if (!reported_flag && v->fake_mismatch_flag && v->fake_frame == v->nb_frame && v->fake_mb == uMb)
    {
        handle_mb_mismatch(f, v, uMb, xMb, yMb, comp, 16, 16);
        reported_flag = 1;
    }
}

/* Match metric. */
double compute_comp_match(frame_data *f, uint32_t comp, uint8_t *p[2], uint32_t width, uint32_t height, uint32_t pitch)
{
    return !f->comp_mismatch[comp];
}

/* Compute the SSIM of the 8x8 block specified. */
float block_ssim(uint8_t *a, uint8_t *b, uint32_t stride)
{
    /* Constant values. */
    float c1 = (0.01*255.0)*(0.01*255.0)*64.0;
    float c2 = (0.03*255.0)*(0.03*255.0)*64.0*63.0;

    /* Sum of each source. */
    uint32_t sa = 0, sb = 0;

    /* Sum of a*b. */
    uint32_t sab = 0;

    /* Sum of a^2 + b^2. */
    uint32_t sa2b2 = 0;

    /* Add up the pixels. */
    for (uint32_t y = 0; y < 8; y++, a += stride, b += stride)
    {
        for (uint32_t x = 0; x < 8; x++)
        {
            uint32_t pa = a[x];
            uint32_t pb = b[x];
            sa += pa;
            sb += pb;
            sab += pa*pb;
            sa2b2 += pa*pa + pb*pb;
        }
    }

    /* Convert the variables to float. */
    float fsa = sa;
    float fsb = sb;
    float fsab = sab;
    float fsa2b2 = sa2b2;

    /* ssim = (2*E(a)*E(b) + c1)*(2*(E(a*b) - E(a)*E(b)) + c2) /
     *        ((E(a)^2 + E(b)^2 + c1)*(E(a*a) - E(a)^2 + E(b*b) - E(b)^2 + c2)).
     */
    return (2.0*fsa*fsb + c1)*(2*(64.0*fsab - fsa*fsb) + c2) /
        ((fsa*fsa + fsb*fsb + c1)*(64.0*fsa2b2 - fsa*fsa - fsb*fsb + c2));
}

/* SSIM metric. */
double compute_comp_ssim(frame_data *f, uint32_t comp, uint8_t *p[2], uint32_t width, uint32_t height, uint32_t pitch)
{
    /* We process 8x8 blocks, offset by 2 to avoid alignment on macroblock
     * boundaries, and with 4 pixels overlap between the blocks.
     */
    uint32_t nbX = (width-6)/4;
    uint32_t nbY = (height-6)/4;
    uint8_t *fa = p[0] + 2*pitch + 2;
    uint8_t *fb = p[1] + 2*pitch + 2;
    double ssim = 0;
    for (uint32_t y = 0; y < nbY; y++, fa += 4*pitch, fb += 4*pitch)
        for (uint32_t x = 0; x < nbX; x++)
            ssim += block_ssim(fa + 4*x, fb + 4*x, pitch);
    ssim /= (float)(nbX*nbY);

    double max_ssim_db = 50;
    double ssim_db = max_ssim_db;
    if (ssim < 1) ssim_db = -10.0 * log10(1.0 - ssim);
    ssim_db = MIN(ssim_db, max_ssim_db);
    return ssim_db;
}

/* PSNR metric. */
double compute_comp_psnr(frame_data *f, uint32_t comp, uint8_t *p[2], uint32_t width, uint32_t height, uint32_t pitch)
{
    /* Compute the squared pixel sum. */
    uint64_t ssum = 0;
    for (uint32_t y = 0; y < height; y++)
    {
        for (uint32_t x = 0; x < width; x++)
        {
            int p0 = p[0][y*pitch + x];
            int p1 = p[1][y*pitch + x];
            int64_t diff = (p0-p1)*(p0-p1);
            ssum += diff;
        }
    }

    /* Compute and cap the PSNR. */
    double max_psnr = 99;
    double psnr = max_psnr;
    double mse = (double)ssum / (double)(width*height);
    if (mse) psnr = 10.0*log(255.0*255.0/mse)/log(10);
    psnr = MIN(psnr, max_psnr);
    return psnr;
}

/* Frame duration metric. */
double compute_frame_duration(frame_data *f)
{
    return f->duration[1];
}

/* Frame size metric. */
double compute_frame_size(frame_data *f)
{
    return f->size[1];
}

/* Analyze a frame pair. */
void process_frame_pair(frame_data *f, video_data *v)
{
    /* Process the frame metrics. */
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->active_flag || md->comp_flag) continue;
        double val = md->compute_frame(f);
        md->cur[0] = val;
        md->sums[0] += val;
        md->ssums[0] += val*val;
    }

    /* Process the components. */
    for (uint32_t comp = 0; comp < v->nb_comp; comp++)
    {
        uint32_t uMb = 0;
        uint32_t width = v->width[comp];
        uint32_t height = v->height[comp];
        uint32_t clip_width = v->clip_width[comp];
        uint32_t clip_height = v->clip_height[comp];
        uint32_t mb_size = v->mb_size[comp];

        /* Process the macroblock pairs. */
        for (uint32_t yMb = 0; yMb < height/mb_size; yMb++)
            for (uint32_t xMb = 0; xMb < width/mb_size; xMb++, uMb++)
                process_mb_pair_comp(f, v, uMb, xMb, yMb, comp);

        /* Process the frame component metrics. */
        for (int i = 0; i < NB_METRICS; i++)
        {
            metric_data *md = v->mda + i;
            if (!md->active_flag || !md->comp_flag) continue;
            uint8_t *p[2];
            for (uint32_t i = 0; i < 2; i++) p[i] = f->p[comp][i];
            double val = md->compute_comp(f, comp, p, clip_width, clip_height, width);
            md->cur[comp] = val;
            md->sums[comp] += val;
            md->ssums[comp] += val*val;
        }
    }

    /* Dump the frame metrics. */
    dump_frame_metrics(f, v);

    /* Print the frame metrics. */
    print_frame_metrics(f, v);

    /* Print the first mismatch. */
    if (f->first_mismatch)
        print_mb_mismatch(f, v, f->first_mismatch);
}

#ifndef F265_NO_LIBAV
/* Read the next non-YUV frame data. Return true on success. */
int read_nonyuv_frame(video_data *v, frame_data *f, int frame_idx)
{
    input_file_data *ifd = v->file_data + frame_idx;
    int frame_flag = 0;

    /* Read the packets until we get a frame, if any. */
    while (1)
    {
        AVPacket packet;
        int done_flag = 0;

        /* Read the next frame, if any. */
        if (av_read_frame(ifd->format_ctx, &packet) < 0)
            break;

        /* Video frame. */
        if (packet.stream_index == ifd->stream_id)
        {
            /* We decoded a frame. */
            #ifdef F265_USE_NEW_LIBAV
            if (avcodec_decode_video2(ifd->codec_ctx, ifd->frame, &frame_flag, &packet /*.data, packet.size*/) >= 0)
            #else
            if (avcodec_decode_video(ifd->codec_ctx, ifd->frame, &frame_flag, packet.data, packet.size) >= 0)
            #endif
            {
                /* We're done if the codec outputted a frame. */
                done_flag = frame_flag;

                /* Get the size and the duration. */
                f->size[frame_idx] = packet.size;
                f->duration[frame_idx] = 0;
                if (packet.duration)
                    f->duration[frame_idx] = (double)packet.duration *
                                             ((double)ifd->stream->time_base.num / (double)ifd->stream->time_base.den);
            }

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
        for (uint32_t i = 0; i < v->nb_comp; i++)
        {
            uint8_t *read_plane = ifd->frame->data[i];
            uint8_t *write_plane = f->p[i][frame_idx];
            uint32_t read_pitch = ifd->frame->linesize[i];
            uint32_t write_pitch = v->width[i];
            uint32_t width = v->clip_width[i];
            uint32_t height = v->clip_height[i];

            for (uint32_t line = 0; line < height; line++)
            {
                memcpy(write_plane, read_plane, width);
                read_plane += read_pitch;
                write_plane += write_pitch;
            }
        }
    }

    return frame_flag;
}
#endif

/* Read the next YUV frame data. Return true on success. */
int read_yuv_frame(video_data *v, frame_data *f, int frame_idx)
{
    input_file_data *ifd = v->file_data + frame_idx;

    for (uint32_t i = 0; i < v->nb_comp; i++)
    {
        uint8_t *write_plane = f->p[i][frame_idx];
        uint32_t write_pitch = v->width[i];
        uint32_t width = v->clip_width[i];
        uint32_t height = v->clip_height[i];

        for (uint32_t line = 0; line < height; line++)
        {
            int nb = fread(write_plane, width, 1, ifd->file_handle);
            if (nb != 1) return 0;
            write_plane += write_pitch;
        }
    }

    return 1;
}

/* Pad the frame if the dimensions are not a multiple of the macroblock size. */
void pad_frame(video_data *v, frame_data *f, int frame_idx)
{
    /* Pad the plane width. */
    if (v->clip_width[0] != v->width[0])
    {
        for (int i = 0; i < 3; i++)
        {
            uint8_t *plane = f->p[i][frame_idx];
            int32_t width = v->width[i];
            int32_t clip_width = v->clip_width[i];
            int32_t clip_height = v->clip_height[i];
            for (int32_t j = 0; j < clip_height; j++)
            {
                uint8_t *src = plane + width*j;
                memset(src + clip_width, src[clip_width-1], width - clip_width);
            }
        }
    }

    /* Pad the plane height. */
    if (v->clip_height[0] != v->height[0])
    {
        for (int i = 0; i < 3; i++)
        {
            uint8_t *plane = f->p[i][frame_idx];
            int32_t width = v->width[i];
            int32_t height = v->height[i];
            int32_t clip_height = v->clip_height[i];
            uint8_t *src = plane + width*(clip_height-1);
            for (int32_t j = clip_height; j < height; j++)
                memcpy(plane + width*j, src, width);
        }
    }
}

/* Read the next frame pair. Return NULL on failure. */
frame_data* read_next_frame_pair(video_data *v)
{
    int frame_flag = 0;
    frame_data *f = frame_data_new(v);

    for (int i = 0; i < 2; i++)
    {
        /* Read the frame. */
        #ifndef F265_NO_LIBAV
        frame_flag = v->file_data[i].yuv_flag ? read_yuv_frame(v, f, i) : read_nonyuv_frame(v, f, i);
        #else
        frame_flag = read_yuv_frame(v, f, i);
        #endif

        if (!frame_flag) break;

        /* Pad the frame. */
        pad_frame(v, f, i);
    }

    if (!frame_flag)
    {
        frame_data_free(f, v);
        f = NULL;
    }

    return f;
}

/* Close an input video file. */
void close_video_file(input_file_data *ifd)
{
   if (ifd->file_handle) fclose(ifd->file_handle);
   #ifndef F265_NO_LIBAV
   if (ifd->frame) av_free(ifd->frame);
   if (ifd->codec_ctx) avcodec_close(ifd->codec_ctx);
   if (ifd->format_ctx)
   #ifdef F265_USE_NEW_LIBAV
   avformat_close_input(&ifd->format_ctx);
   #else
   av_close_input_file(ifd->format_ctx);
   #endif
   #endif
   memset(ifd, 0, sizeof(*ifd));
}

/* Open a YUV video file. */
void open_yuv_file(input_file_data *ifd)
{
    FILE *f = fopen(ifd->path, "rb");
    if (!f)
    {
        printf("Cannot open %s: %s.\n", ifd->path, strerror(errno));
        exit(1);
    }
    ifd->file_handle = f;
}

/* Open a non-yuv video file. */
#ifndef F265_NO_LIBAV
void open_nonyuv_file(input_file_data *ifd)
{
    char *fail = NULL;

    /* Try. */
    do
    {
        /* Open the video file. */
        #ifdef F265_USE_NEW_LIBAV
        AVInputFormat* avinputFormat = NULL;
        if (ifd->hevc_flag)
        {
           avinputFormat = av_find_input_format("hevc");
        }

        if (avformat_open_input(&ifd->format_ctx, ifd->path, avinputFormat, NULL) != 0)
        #else
        if (av_open_input_file(&ifd->format_ctx, ifd->path, NULL, 0, NULL) != 0)
        #endif
        {
            fail = "av_open_input_file() failed"; break;
        }

        /* Find the stream information. */
        #ifdef F265_USE_NEW_LIBAV
        if (avformat_find_stream_info(ifd->format_ctx, NULL) < 0)
        #else
        if (av_find_stream_info(ifd->format_ctx) < 0)
        #endif
        {
            fail = "av_find_stream_info() failed"; break;
        }

        /* Find the first video stream. */
        ifd->stream_id = -1;
        for (int i = 0; i < ifd->format_ctx->nb_streams; i++)
        {
            ifd->stream = ifd->format_ctx->streams[i];
            #ifdef F265_USE_NEW_LIBAV
            if (ifd->stream->codec->codec_type == AVMEDIA_TYPE_VIDEO)
            #else
            if (ifd->stream->codec->codec_type == CODEC_TYPE_VIDEO)
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
        #ifdef F265_USE_NEW_LIBAV
        if (avcodec_open2(ifd->codec_ctx, ifd->codec, NULL) < 0)
        #else
        if (avcodec_open(ifd->codec_ctx, ifd->codec) < 0)
        #endif
        {
            fail = "avcodec_open() failed"; break;
        }

        /* Allocate a video frame. */
        ifd->frame = avcodec_alloc_frame();

        /* Set the luma width and height. */
        ifd->luma_dim[0] = ifd->codec_ctx->width;
        ifd->luma_dim[1] = ifd->codec_ctx->height;

    } while (0);

    if (fail)
    {
        printf("Cannot open %s: %s.\n", ifd->path, fail);
        exit(1);
    }
}
#endif

/* Open an input video file. */
void open_video_file(input_file_data *ifd)
{
    /* Determine whether the file is an YUV file based on the extension. */
    char ext[4] = { 0, 0, 0, 0 };
    int path_len = strlen(ifd->path);
    if (path_len >= 4) memcpy(ext, ifd->path + path_len - 4, 4);
    for (int i = 0; i < 4; i++) ext[i] = tolower(ext[i]);
    ifd->yuv_flag = !memcmp(ext, ".yuv", 4);

    if ((memcmp(ext, ".265", 4) == 0) || (memcmp(ext, ".hevc", 5) == 0) || (memcmp(ext, ".h265", 5) == 0))
    {
        #ifdef F265_USE_NEW_LIBAV
        ifd->hevc_flag = 1;
        #else
        printf("You need to use libav-hm10.0 to have hevc (h265) decoding support\n");
        exit(-1);
        #endif
    }
    else
    {
       ifd->hevc_flag = 0;
    }
    /* Open the file. */
    #ifndef F265_NO_LIBAV
    ifd->yuv_flag ? open_yuv_file(ifd) : open_nonyuv_file(ifd);
    #else
    open_yuv_file(ifd);
    #endif
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
    printf("Usage: yuvdiff <file1> <file2>\n");
    printf("Options:\n");
    printf("  -h: print help and exit.\n");
    printf("  -c: assume files contain no chroma data (luma only).\n");
    printf("  -C: do not consider chroma for macroblock mismatch.\n");
    printf("  -w <width:height>: set the YUV frame size.\n");
    printf("  -r: report the first macroblock mismatch.\n");
    printf("  -t: assume text output (not terminal).\n");
    printf("  -f <frame:mb>: pretend that the only mismatch is macroblock 'mb' in frame\n");
    printf("                 'frame'.\n");
    printf("  -m <metric,...>: set the list of quality metrics to report.\n");
    printf("  -p: print the metric values for each frame.\n");
    printf("  -a: print the average of the metric values.\n");
    printf("  -d: print the standard deviation of the metric values.\n");
    printf("  -D <dir>: dump the metric values in a directory.\n");
    printf("Metric list:\n");
    printf("  size, duration, match, psnr, ssim.\n");
    printf("  The size and duration metrics are unavailable if file2 is a YUV file.\n");
    #ifdef F265_USE_NEW_LIBAV
    printf("This version support HEVC raw bytestream file (extension .hevc, .h265 or .265)\n");
    #endif
    printf("\n");
    exit(exit_code);
}

/* Set the frame dimensions. */
void set_frame_dimensions(video_data *v)
{
    input_file_data *ifd = v->file_data;

    /* If both files are YUV files, then the frame dimensions must be set
     * explicitly.
     */
    if (ifd[0].yuv_flag && ifd[1].yuv_flag)
    {
        if (!ifd[0].luma_dim[0] || !ifd[0].luma_dim[1])
        {
            printf("Please specify the YUV frame dimensions (-w).\n");
            print_usage(1);
        }
    }

    /* Import the YUV frame dimensions from the non-yuv file, if required. */
    else
    {
        uint32_t *dim = NULL;
        for (int i = 0; i < 2; i++)
        {
            if (!ifd[i].yuv_flag)
            {
                dim = ifd[i].luma_dim;
                break;
            }
        }
        for (int i = 0; i < 2; i++)
        {
            if (ifd[i].yuv_flag && !ifd[i].luma_dim[0])
            {
                ifd[i].luma_dim[0] = dim[0];
                ifd[i].luma_dim[1] = dim[1];
            }
        }
    }

    /* Check for luma dimension mismatch. */
    if (ifd[0].luma_dim[0] != ifd[1].luma_dim[0] || ifd[0].luma_dim[1] != ifd[1].luma_dim[1])
    {
        printf("Frame dimension mismatch (%dx%d vs %dx%d).\n",
               ifd[0].luma_dim[0], ifd[0].luma_dim[1], ifd[1].luma_dim[0], ifd[1].luma_dim[1]);
        print_usage(1);
    }

    /* Set the plane dimensions and validate. */
    v->clip_width[0] = ifd[0].luma_dim[0];
    v->clip_height[0] = ifd[0].luma_dim[1];

    if (!v->clip_width[0] || !v->clip_height[0])
    {
        printf("Bad frame dimensions.\n");
        print_usage(1);
    }

    v->width[0] = (v->clip_width[0] + 15) & ~15;
    v->height[0] = (v->clip_height[0] + 15) & ~15;

    for (uint32_t i = 0; i < 2; i++)
    {
        v->clip_width[1+i] = v->clip_width[0] / 2;
        v->clip_height[1+i] = v->clip_height[0] / 2;
        v->width[1+i] = v->width[0] / 2;
        v->height[1+i] = v->height[0] / 2;
    }

    /* Set the field sizes. */
    v->field_min_size = 8;
    v->field_sep = 3;
    v->metric_digits = 6;
    v->total_line_size = v->field_min_size;
    int sep_needed = 0;
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->active_flag) continue;
        int nb_comp = md->comp_flag ? v->nb_comp : 1;
        int digits_size = v->metric_digits*nb_comp + 1*(nb_comp-1);
        int field_length = MAX(digits_size, v->field_min_size);
        md->pad_name = field_length - strlen(md->name);
        md->pad_digits = field_length - digits_size;
        v->total_line_size += field_length;
        sep_needed++;
    }
    v->total_line_size += v->field_sep*sep_needed;
}

/* Parse the command line. */
void parse_cmd_line(int argc, char **argv, video_data *v)
{
    while (1)
    {
        uint32_t ret = getopt(argc, argv, "hcCw:rtf:m:padD:");
        if (ret == -1) break;
        else if (ret == '?' || ret == ':') print_usage(1);
        else if (ret == 'h') print_usage(0);
        else if (ret == 'c') { v->chroma_flag = 0; v->nb_comp = 1; }
        else if (ret == 'C') { v->chroma_mismatch_flag = 0; }
        else if (ret == 'w')
        {
            char width_str[100], height_str[100];
            if (!split_str(optarg, width_str, height_str))
            {
                printf("Invalid YUV frame size format.\n");
                print_usage(1);
            }
            for (int i = 0; i < 2; i++)
            {
                v->file_data[i].luma_dim[0] = atoi(width_str);
                v->file_data[i].luma_dim[1] = atoi(height_str);
            }
        }
        else if (ret == 't') v->text_mode_flag = 1;
        else if (ret == 'r') v->print_mismatch_flag = 1;
        else if (ret == 'f')
        {
            char frame_str[100], mb_str[100];
            if (!split_str(optarg, frame_str, mb_str))
            {
                printf("Invalid fake macroblock mismatch string format.\n");
                print_usage(1);
            }
            v->print_mismatch_flag = 1;
            v->fake_mismatch_flag = 1;
            v->fake_frame = atoi(frame_str);
            v->fake_mb = atoi(mb_str);
        }
        else if (ret == 'm')
        {
            char *start = optarg;
            char *end = start;
            while (1)
            {
                if (*end == 0 || *end == ',')
                {
                    int len = end-start;
                    if (len)
                    {
                        char str[100];
                        memcpy(str, start, len);
                        str[len] = 0;
                        int found_flag = 0;
                        for (int i = 0; i < NB_METRICS; i++)
                        {
                            if (!strcmp(v->mda[i].name, str))
                            {
                                v->mda[i].active_flag = 1;
                                found_flag = 1;
                                break;
                            }
                        }
                        if (!found_flag)
                        {
                            printf("Invalid metric name %s.\n", str);
                            print_usage(1);
                        }
                    }
                    if (!*end) break;
                    end++;
                    start = end;
                }
                else end++;
            }
        }
        else if (ret == 'p') v->print_frame_metric_flag = 1;
        else if (ret == 'a') v->print_average_flag = 1;
        else if (ret == 'd') v->print_std_deviation_flag = 1;
        else if (ret == 'D') v->dump_dir = optarg;
    }

    int nb_arg_left = argc - optind;
    if (nb_arg_left != 2) print_usage(1);

    v->file_data[0].path = argv[optind + 0];
    v->file_data[1].path = argv[optind + 1];
}

/* Initialize the video data. */
void init_video(video_data *v)
{
    memset(v, 0, sizeof(video_data));

    v->chroma_flag = 1;
    v->chroma_mismatch_flag = 1;
    v->mb_size[0] = 16;
    v->mb_size[1] = 8;
    v->mb_size[2] = 8;
    v->nb_comp = 3;

    static char* metric_names[NB_METRICS] = { "match", "size", "duration", "psnr", "ssim" };
    static char* metric_print_fmts[NB_METRICS] = { "%6.4f", "% 6.0f", "% 6.3f", "% 6.2f", "% 6.2f" };
    static char* metric_dump_fmts[NB_METRICS] = { "%1.0f", "%06.0f", "%06.3f", "%05.2f", "%05.2f" };
    static int metric_comp_flags[NB_METRICS] = { 1, 0, 0, 1, 1 };
    static void* metric_funcs[NB_METRICS] = { compute_comp_match, compute_frame_size, compute_frame_duration,
                                              compute_comp_psnr, compute_comp_ssim };
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        md->name = metric_names[i];
        md->print_fmt = metric_print_fmts[i];
        md->dump_fmt = metric_dump_fmts[i];
        md->comp_flag = metric_comp_flags[i];
        md->compute_func = metric_funcs[i];
    }

    #ifndef F265_NO_LIBAV
    av_register_all();
    #endif
}

/* Open the video files. */
void open_files(video_data *v)
{
    /* Open the input files. */
    for (int i = 0; i < 2; i++)
        open_video_file(v->file_data + i);

    /* Disable the size and the duration metrics if the second file is a YUV
     * file, since the data is not available in that case.
     */
    if (v->file_data[1].yuv_flag)
    {
        v->mda[METRIC_SIZE].active_flag = 0;
        v->mda[METRIC_DURATION].active_flag = 0;
    }

    /* Open the metric files. */
    if (v->dump_dir)
    {
        if (access(v->dump_dir, F_OK))
        {
            #ifdef _WIN32
            if (mkdir(v->dump_dir))
            #else
            if (mkdir(v->dump_dir, 0777))
            #endif
            {
                printf("Cannot create directory %s: %s.\n", v->dump_dir, strerror(errno));
                exit(1);
            }
        }

        for (int i = 0; i < NB_METRICS; i++)
        {
            metric_data *md = v->mda + i;
            if (!md->active_flag) continue;
            char path[1000] = { 0 };
            strcat(path, v->dump_dir);
            strcat(path, "/");
            strcat(path, md->name);
            md->file_handle = fopen(path, "wb");
            if (!md->file_handle)
            {
                printf("Cannot open %s: %s.\n", path, strerror(errno));
                exit(1);
            }
        }
    }
}

/* Close the video files. */
void close_files(video_data *v)
{
    /* Close the input files. */
    for (int i = 0; i < 2; i++)
        close_video_file(v->file_data + i);

    /* Close the metric files. */
    for (int i = 0; i < NB_METRICS; i++)
    {
        metric_data *md = v->mda + i;
        if (!md->file_handle) continue;
        if (fclose(md->file_handle))
        {
            printf("Cannot close metric file: %s.\n", strerror(errno));
            exit(1);
        }
        md->file_handle = NULL;
    }
}

int main(int argc, char **argv)
{
    video_data v;

    /* Initialize. */
    init_video(&v);

    /* Parse the command line. */
    parse_cmd_line(argc, argv, &v);

    /* Open the files. */
    open_files(&v);

    /* Set the frame dimensions. */
    set_frame_dimensions(&v);

    /* Print the header. */
    print_header(&v);

    /* Read and process the frame pairs until requested to stop. */
    while (1)
    {
        frame_data *f = read_next_frame_pair(&v);
        if (f == NULL) break;
        process_frame_pair(f, &v);

        if (v.print_mismatch_flag)
            if (f->comp_mismatch[0] + f->comp_mismatch[1] + f->comp_mismatch[2])
                break;

        frame_data_free(f, &v);
        v.nb_frame++;
    }

    /* Print the footer. */
    print_footer(&v);

    /* Close the files. */
    close_files(&v);

    return 0;
}
