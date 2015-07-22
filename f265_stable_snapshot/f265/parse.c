// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include "f265/bdi.h"
#include "ktools/kstr.h"
#include "ktools/karray.h"
#include "ktools/khash.h"
#include "ktools/kutil.h"

typedef struct f265_parse_ctx f265_parse_ctx;
typedef struct f265_parse_entry f265_parse_entry;
typedef union f265_parse_arg f265_parse_arg;

// Context used by the encoder's parameter-parsing code.
struct f265_parse_ctx
{
    // Entry being processed, if any.
    f265_parse_entry *entry;

    // Error string.
    kstr err;
};

// Encoder parameter table entry.
struct f265_parse_entry
{
    // Name of the key.
    const char *key;

    // Parameter handler function. Must call f265_handle_parse_error() on error.
    void (*handler)(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *args, int32_t nb_args);

    // Number of arguments expected. Zero means a variable number of arguments
    // is allowed.
    int32_t nb_args;

    // Argument type: 0=int, 1=float, 2=string. All arguments have the same
    // type. The string type can be converted to other types as needed.
    int32_t arg_type;
};

// Encoder parameter argument (int, float or string).
union f265_parse_arg
{
    int64_t i;
    float f;
    char *s;
};

// Abort the program when a memory request cannot be fulfilled.
static void f265_out_of_memory()
{
    printf("f265 encoder out of memory.\n");
    abort();
}

// Malloc wrappers that abort the program if the memory request cannot be
// fulfilled.
static void* f265_malloc(size_t count)
{
    void *ptr = malloc(count);
    if (ptr == NULL) f265_out_of_memory();
    return ptr;
}
static void* f265_calloc(size_t count)
{
    void *ptr = calloc(1, count);
    if (ptr == NULL) f265_out_of_memory();
    return ptr;
}
static void* f265_realloc(void *ptr, size_t count)
{
    ptr = realloc(ptr, count);
    if (ptr == NULL) f265_out_of_memory();
    return ptr;
}

// Handle an error encountered during parameter parsing. Only keep the first
// error.
static void f265_handle_parse_error(f265_parse_ctx *ctx, const char *format, ...)
{
    kstr *err = &ctx->err;
    if (err->slen) return;
    if (ctx->entry) kstr_sf(err, "%s: ", ctx->entry->key);
    va_list arg;
    va_start(arg, format);
    kstr_append_sfv(err, format, arg);
    va_end(arg);
}

// Convert the argument specified to integer. Return 0 on success. Set the
// context error on error.
static int f265_get_int_param(f265_parse_ctx *ctx, char *s, int64_t *i)
{
    char *end;
    *i = strtol(s, &end, 10);
    if (strlen(s) && !*end) return 0;
    f265_handle_parse_error(ctx, "invalid integer value '%s'", s);
    return 1;
}

// Same as above for float.
static int f265_get_float_param(f265_parse_ctx *ctx, char *s, float *f)
{
    char *end;
    *f = strtof(s, &end);
    if (strlen(s) && !*end) return 0;
    f265_handle_parse_error(ctx, "invalid floating point value '%s'", s);
    return 1;
}

// Parameter handler functions.
static void handle_param_format(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    if (!strcmp(a->s, "byte")) p->format = F265_FORMAT_BYTE;
    else if (!strcmp(a->s, "avc")) p->format = F265_FORMAT_AVC;
    else f265_handle_parse_error(ctx, "invalid format '%s'", a->s);
}

static void handle_param_chroma_format(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->chroma_format = a->i;
}

static void handle_param_bit_depth(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    for (int i = 0; i < 4; i++) p->bit_depth[i] = a[i].i;
}

static void handle_param_cb_range(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->cb_range[0] = a[0].i;
    p->cb_range[1] = a[1].i;
}

static void handle_param_pcm_range(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->pcm_range[0] = a[0].i;
    p->pcm_range[1] = a[1].i;
}

static void handle_param_tb_range(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->tb_range[0] = a[0].i;
    p->tb_range[1] = a[1].i;
}

static void handle_param_tb_depth(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->tb_depth[0] = a[0].i;
    p->tb_depth[1] = a[1].i;
}

static void handle_param_qg(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->qg_log = a->i;
}

static void handle_param_ref(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->nb_refs = a->i;
}

static void handle_param_bframes(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->nb_b_frames = a->i;
}

static void handle_param_wpp(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->wpp_flag = a->i;
}

static void handle_param_deblock(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->deblock_flag = a->i;
}

static void handle_param_sao(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->sao_flag = a->i;
}

static void handle_param_scaling(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->scaling_flag = a->i;
}

static void handle_param_transquant_bypass(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->transquant_bypass_flag = a->i;
}

static void handle_param_sign_hiding(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->sign_hiding_flag = a->i;
}

static void handle_param_rdoq(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->rdoq_flag = a->i;
}

static void handle_param_transform_skip(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->transform_skip_flag = a->i;
}

static void handle_param_smooth_intra(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->smooth_intra_flag = a->i;
}

static void handle_param_amp(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->amp_flag = a->i;
}

static void handle_param_tmv(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->tmv_flag = a->i;
}

static void handle_param_merge_candidate(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->merge_cand = a->i;
}

static void handle_param_parallel_merge_level(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->parallel_merge_level = a->i;
}

static void handle_param_weight(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->weight_flag = a->i;
}

static void handle_param_chroma_me(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->chroma_me_flag = a->i;
}

static void handle_param_quality(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    int quality = a->i;

    // Temporary.
    if (quality <= 0)
    {
        p->cb_range[0] = 4;
        p->cb_range[1] = 4;
        p->tb_range[0] = 3;
        p->tb_range[1] = 4;
        p->tb_depth[0] = 0;
        p->tb_depth[1] = 0;
        p->merge_cand = 1;
        for (int i = 0; i < 3; i++) p->me_algo[i] = 0;
        p->algo = (1<<10);
    }

    else if (quality <= 25)
    {
        p->cb_range[0] = 3;
        p->cb_range[1] = 4;
        p->tb_range[0] = 2;
        p->tb_range[1] = 4;
        p->tb_depth[0] = 1;
        p->tb_depth[1] = 1;
        p->rdo_level = 1;
        p->amp_flag = 0;
        p->tmv_flag = 1;
        p->me_dist[2] = 1;
        p->algo = (1<<10) + 32 + 7;
    }

    else if (quality <= 50)
    {
        p->cb_range[0] = 3;
        p->cb_range[1] = 6;
        p->tb_range[0] = 2;
        p->tb_range[1] = 5;
        p->tb_depth[0] = 4;
        p->tb_depth[1] = 4;
        p->nb_refs = 3;
        p->rdo_level = 1;
        p->amp_flag = 1;
        p->tmv_flag = 1;
        p->me_dist[2] = 1;
        p->hm_me_flag = 1;
        p->all_intra_flag = 1;
        p->nullify_inter_tb_flag = 1;
    }
}

static void handle_param_rdo(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->rdo_level = a->i;
}

static void handle_param_hm_me(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->hm_me_flag = a->i;
}

static void handle_param_all_intra(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->all_intra_flag = a->i;
}

static void handle_param_nullify_inter_tb(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->nullify_inter_tb_flag = a->i;
}

static void handle_param_qp(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->qp = a->i;
}

static void handle_param_me_algo(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int index)
{
    int64_t algo = 0, iter = 0, dist = 0;

    if (!strcmp(a[0].s, "dia")) algo = 0;
    else if (!strcmp(a[0].s, "xdia")) algo = 1;
    else if (!strcmp(a[0].s, "hex")) algo = 2;
    else if (!strcmp(a[0].s, "square")) algo = 3;
    else f265_handle_parse_error(ctx, "invalid algorithm '%s'", a[0].s);

    f265_get_int_param(ctx, a[1].s, &iter);

    if (!strcmp(a[2].s, "sad")) dist = 0;
    else if (!strcmp(a[2].s, "satd")) dist = 1;
    else f265_handle_parse_error(ctx, "invalid distortion metric '%s'", a[2].s);

    p->me_algo[index] = algo;
    p->me_iter[index] = iter;
    p->me_dist[index] = dist;
}

static void handle_param_fpel(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
    { handle_param_me_algo(ctx, p, a, 0); }
static void handle_param_hpel(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
    { handle_param_me_algo(ctx, p, a, 1); }
static void handle_param_qpel(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
    { handle_param_me_algo(ctx, p, a, 2); }

static void handle_param_key_frame_spacing(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->key_interval = a->i;
}

static void handle_param_key_frame_type(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->key_type = a->i;
}

static void handle_param_yuv_dump(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    free(p->yuv_dump_path);
    p->yuv_dump_path = strdup(a->s);
    if (!p->yuv_dump_path) f265_out_of_memory();
}

static void handle_param_hm_gop(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    free(p->hm_gop_path);
    p->hm_gop_path = strdup(a->s);
    if (!p->hm_gop_path) f265_out_of_memory();
}

static void handle_param_rc(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    if (!strcmp(a->s, "cqp")) p->rc_method = F265_RCM_CQP;
    else if (!strcmp(a->s, "abr")) p->rc_method = F265_RCM_ABR;
    else f265_handle_parse_error(ctx, "invalid method '%s'", a->s);
}

static void handle_param_bitrate(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->bitrate[0] = a->i;
}

static void handle_param_bitrate_range(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->bitrate[1] = a[0].i;
    p->bitrate[2] = a[1].i;
}

static void handle_param_qp_bounds(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->qp_bounds[0] = a[0].i;
    p->qp_bounds[1] = a[1].i;
}

static void handle_param_fps(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    if (nb_args < 1 || nb_args > 2)
        f265_handle_parse_error(ctx, "at least the numerator has to be provided for the frame rate e.g. 30");

    p->frame_rate_num = a[0].i;

    p->frame_rate_den = 1;
    if (nb_args == 2)
        p->frame_rate_den = a[1].i;
}

static void handle_param_lt_conv_min(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->lt_conv_min = a[0].i;
}

static void handle_param_lt_conv_exp(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    p->lt_conv_exp = a[1].f;
}

static void handle_param_algo(f265_parse_ctx *ctx, f265_enc_params *p, f265_parse_arg *a, int32_t nb_args)
{
    int nb_flags = F265_MIN(strlen(a->s), 64);
    for (int i = 0; i < nb_flags; i++)
        if (a->s[i] == '1')
            p->algo |= (1<<(nb_flags-i-1));
}

// Parameter dispatch table.
static const f265_parse_entry f265_enc_params_table[] =
{
    // Name, handler, argument count, argument type (0=int, 1=float, 2=string).
    { "format", handle_param_format, 1, 2 },
    { "chroma-format", handle_param_chroma_format, 1, 0 },
    { "bit-depth", handle_param_bit_depth, 4, 0 },
    { "cb-range", handle_param_cb_range, 2, 0 },
    { "pcm-range", handle_param_pcm_range, 2, 0 },
    { "tb-range", handle_param_tb_range , 2, 0 },
    { "tb-depth", handle_param_tb_depth, 2, 0 },
    { "qg", handle_param_qg, 1, 0 },
    { "ref", handle_param_ref, 1, 0 },
    { "bframes", handle_param_bframes, 1, 0 },
    { "wpp", handle_param_wpp, 1, 0 },
    { "deblock", handle_param_deblock, 1, 0 },
    { "sao", handle_param_sao, 1, 0 },
    { "scaling", handle_param_scaling, 1, 0 },
    { "transquant-bypass", handle_param_transquant_bypass, 1, 0 },
    { "rdoq", handle_param_rdoq, 1, 0 },
    { "sign-hiding", handle_param_sign_hiding, 1, 0 },
    { "transform-skip", handle_param_transform_skip, 1, 0 },
    { "smooth-intra", handle_param_smooth_intra, 1, 0 },
    { "amp", handle_param_amp, 1, 0 },
    { "tmv", handle_param_tmv, 1, 0 },
    { "nb-merge", handle_param_merge_candidate, 1, 0 },
    { "pml", handle_param_parallel_merge_level, 1, 0 },
    { "weight", handle_param_weight, 1, 0 },
    { "fpel", handle_param_fpel, 3, 2 },
    { "hpel", handle_param_hpel, 3, 2 },
    { "qpel", handle_param_qpel, 3, 2 },
    { "chroma-me", handle_param_chroma_me, 1, 0 },
    { "quality", handle_param_quality, 1, 0 },
    { "rdo", handle_param_rdo, 1, 0 },
    { "hm-me", handle_param_hm_me, 1, 0 },
    { "all-intra", handle_param_all_intra, 1, 0 },
    { "nullify-inter-tb", handle_param_nullify_inter_tb, 1, 0, },
    { "qp", handle_param_qp, 1, 0 },
    { "key-frame-spacing", handle_param_key_frame_spacing, 1, 0 },
    { "key-frame-type", handle_param_key_frame_type, 1, 0 },
    { "yuv-dump", handle_param_yuv_dump, 1, 2 },
    { "hm-gop", handle_param_hm_gop, 1, 2 },
    { "rc", handle_param_rc, 1, 2 },
    { "bitrate", handle_param_bitrate, 1, 0 },
    { "bitrate-range", handle_param_bitrate_range, 2, 0 },
    { "qp-bounds", handle_param_qp_bounds, 2, 0 },
    { "fps", handle_param_fps, 0, 0 },
    { "lt-conv-min", handle_param_lt_conv_min, 1, 0 },
    { "lt-conv-exp", handle_param_lt_conv_exp, 1, 1 },
    { "algo", handle_param_algo, 1, 2 },
};

void f265_parse_params(f265_enc_params *params, const char *param_string, char **error_handle)
{
    // Parameter context.
    f265_parse_ctx ctx;

    // Work string.
    kstr str;

    // Array of temporary string values.
    karray str_array;

    // Argument array.
    int nb_args = 0;
    int arg_array_size = 0;
    f265_parse_arg *arg_array = NULL;

    // key->entry hash.
    khash key_hash;

    // Current mode: 0 (skip mode), 1 (key mode), 2 (value mode).
    int mode = 0;

    // Start-of-string pointer.
    const char *start = NULL;

    ctx.entry = NULL;
    kstr_init(&ctx.err);
    kstr_init(&str);
    karray_init(&str_array);
    khash_init_func(&key_hash, kutil_cstr_hc, kutil_cstr_eq);

    // Fill the key hash.
    for (int i = 0; i < sizeof(f265_enc_params_table)/sizeof(f265_parse_entry); i++)
        khash_add(&key_hash, (char*)(f265_enc_params_table[i].key), (void*)&f265_enc_params_table[i]);

    // Set the default parameter values.
    f265_set_default_params(params);

    // Kludge. Apply quality=25 if no parameter is specified.
    if (!strlen(param_string))
    {
        f265_parse_arg quality;
        quality.i = 25;
        handle_param_quality(&ctx, params, &quality, 1);
    }

    // Loop until an exit condition is reached.
    for (const char *cur = param_string; ; cur++)
    {
        char c = *cur;

        // Skip mode.
        if (mode == 0)
        {
            if (!c) break;
            else if (c == ';' || c == ' ') {}
            else if (c >= 'a' && c <= 'z' || c == '-')
            {
                mode = 1;
                start = cur;
            }
            else
            {
                f265_handle_parse_error(&ctx, "invalid character '%c' in parameter name", c);
                break;
            }
        }

        // key mode.
        else if (mode == 1)
        {
            if (c >= 'a' && c <= 'z' || c == '-') {}
            else if (c == '=' || c == ':')
            {
                kstr_assign_buf(&str, start, cur-start);
                ctx.entry = (f265_parse_entry*)khash_get(&key_hash, str.data);
                if (!ctx.entry)
                {
                    f265_handle_parse_error(&ctx, "invalid parameter '%s'", str.data);
                    break;
                }
                mode = 2;
                start = cur+1;
            }
            else if (!c)
            {
                f265_handle_parse_error(&ctx, "unterminated parameter name");
                break;
            }
            else
            {
                f265_handle_parse_error(&ctx, "invalid character '%c' in parameter name", c);
                break;
            }
        }

        // Value mode.
        else
        {
            if (c == '=' || c == ':')
            {
                f265_handle_parse_error(&ctx, "invalid character '%c' in value", c);
                break;
            }

            else if (!c || c == ',' || c == ';' || c == ' ')
            {
                if (cur == start)
                {
                    f265_handle_parse_error(&ctx, "no value specified");
                    break;
                }
                kstr_assign_buf(&str, start, cur-start);

                // Ignore 'none' arguments.
                if (strcmp(str.data, "none"))
                {
                    // Map true and false to numeric arguments.
                    if (!strcmp(str.data, "true")) kstr_assign_cstr(&str, "1");
                    if (!strcmp(str.data, "false")) kstr_assign_cstr(&str, "0");

                    if (nb_args == arg_array_size)
                    {
                        arg_array_size += 10;
                        arg_array = (f265_parse_arg*)f265_realloc(arg_array, arg_array_size*sizeof(f265_parse_arg));
                    }
                    f265_parse_arg *arg = arg_array + nb_args++;

                    if (ctx.entry->arg_type == 0)
                    {
                        if (f265_get_int_param(&ctx, str.data, &arg->i)) break;
                    }

                    else if (ctx.entry->arg_type == 1)
                    {
                        if (f265_get_float_param(&ctx, str.data, &arg->f)) break;
                    }

                    else
                    {
                        kstr *s = kstr_new();
                        kstr_assign_kstr(s, &str);
                        karray_push(&str_array, s);
                        arg->s = s->data;
                    }
                }

                if (c == ',')
                {
                    start = cur + 1;
                }

                else
                {
                    if (ctx.entry->nb_args && nb_args != ctx.entry->nb_args)
                    {
                        f265_handle_parse_error(&ctx, "parameter requires %d arguments", ctx.entry->nb_args);
                        break;
                    }

                    ctx.entry->handler(&ctx, params, arg_array, nb_args);
                    karray_cleanse(&str_array, (kcleaner)kstr_destroy);
                    nb_args = 0;
                    ctx.entry = NULL;
                    mode = 0;
                    if (!c || ctx.err.slen) break;
                }
            }
        }
    }

    // Error.
    if (ctx.err.slen)
    {
        f265_free_params(params);
        *error_handle = (char*)f265_malloc(ctx.err.slen+1);
        strcpy(*error_handle, ctx.err.data);
    }

    // Success.
    else *error_handle = NULL;

    kstr_clean(&ctx.err);
    kstr_clean(&str);
    karray_cleanse(&str_array, (kcleaner)kstr_destroy);
    karray_clean(&str_array);
    free(arg_array);
    khash_clean(&key_hash);
}

void f265_free_params(f265_enc_params *params)
{
    if (!params) return;
    free(params->yuv_dump_path); params->yuv_dump_path = NULL;
    free(params->hm_gop_path); params->hm_gop_path = NULL;
}

void f265_init_parsing()
{
    ktools_init(f265_malloc, f265_calloc, f265_realloc, free);
}

