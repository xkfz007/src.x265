// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

// HM-related tests. Put ugly code here.

#include "f265/enc.h"

// Open an HM dump file for reading or writing.
void venc_open_dump_file(FILE **f, const char *path, const char *mode)
{
    if (*f) return;
    *f = fopen(path, mode);
    if (!*f) { printf("Cannot open %s.\n", path); exit(1); }
}

// Show the current block location.
void venc_show_loc(f265_enc_thread *t, int ct_ox, int ct_oy, int size)
{
    printf("Frame %d CTB %d (%d,%d) offsets (%d,%d) size %d pix (%d,%d).\n",
           (int)t->src_frame->abs_poc, t->ctb_xy, t->ctb_x, t->ctb_y, ct_ox, ct_oy, size,
           t->ctb_off[0]+ct_ox, t->ctb_off[1]+ct_oy);
}

// Print the (filtered) neighbours of a block.
void venc_print_intra_neighbours(f265_pix nbuf[129], int bs, const char *prefix)
{
    printf("%sH: ", prefix);
    for (int i = 0; i < 2*bs+1; i++) printf("%02x ", nbuf[i]);
    printf("\n");
    printf("%sV:    ", prefix);
    for (int i = 0; i < 2*bs; i++) printf("%02x ", nbuf[65+i]);
    printf("\n");
}

// Convert the pixel offsets of a block from the CTB origin to the uiAbsPartIdx
// variable used by HM.
int venc_hm_part_idx(f265_gen_data *gd, int b_ox, int b_oy)
{
    int shift = gd->tb_range[0];
    return venc_get_depth_scan_idx(b_ox>>shift, b_oy>>shift, 1<<(gd->cb_range[1]-gd->tb_range[0]));
}

#ifdef VAN_LOAD_CTB_ANALYSIS
FILE *van_an_dump_file;

// Read from the analysis file.
static void van_an_read(void *dst, int size, FILE *file)
{
    if (fread(dst, size, 1, file) != 1)
    {
        printf("Analysis file end-of-file.\n");
        exit(1);
    }
}

// Load the transform tree recursively.
static void van_load_an_tt(f265_enc_thread *t)
{
    // Load the node.
    uint8_t tn;
    van_an_read(&tn, 1, van_an_dump_file);
    *t->tt.tn++ = tn;

    // Handle split.
    if (tn&8)
        for (int i = 0; i < 4; i++)
            van_load_an_tt(t);
}

// Load the CB quadtree recursively. We assume the CB flags are clean.
static void van_load_an_cb_qt(f265_enc_thread *t, f265_cb *cb)
{
    if (!(cb->flags&F265_CB_PRESENT)) return;

    uint8_t split_flag;
    van_an_read(&split_flag, 1, van_an_dump_file);

    // Handle split.
    F265_SET_FLAG(cb->flags, F265_CB_SPLIT, split_flag);
    if (split_flag)
    {
        for (int i = 0; i < 4; i++) van_load_an_cb_qt(t, t->cb + cb->child_idx + i);
        return;
    }

    // CB mode (intra 0, PCM 1, inter 2, skip 3).
    uint8_t cb_mode;
    van_an_read(&cb_mode, 1, van_an_dump_file);
    if (cb_mode <= 1) F265_SET_FLAG(cb->flags, F265_CB_INTRA, 1);

    // Set PCM, update the prediction map and return.
    if (cb_mode == 1)
    {
        cb->intra_luma_mode[0] = 35;
        venc_update_pmap_unsplit_cb(t, cb);
        return;
    }

    // Set skip.
    if (cb_mode == 3) F265_SET_FLAG(cb->flags, F265_CB_SKIP, 1);

    // Set the partitioning mode.
    uint8_t part_mode;
    van_an_read(&part_mode, 1, van_an_dump_file);
    cb->intra_luma_mode[1] = -1;
    cb->inter_part = part_mode;

    // Validate that we're not importing inter HV.
    if (cb_mode == 2 && part_mode == F265_PART_HV)
    {
        printf("Cannot import inter HV mode.\n");
        exit(1);
    }

    // Pass each partition.
    int nb_parts = f265_nb_parts[part_mode];
    for (int part_idx = 0; part_idx < nb_parts; part_idx++)
    {
        // Intra luma modes.
        if (cb_mode == 0)
        {
            uint8_t intra_mode;
            van_an_read(&intra_mode, 1, van_an_dump_file);
            cb->intra_luma_mode[part_idx] = intra_mode;
        }

        // Inter or skip.
        else
        {
            // Merge index and list flags (1: L0, 2: L1, 3: Bi).
            uint8_t merge_idx, list_flags;
            van_an_read(&merge_idx, 1, van_an_dump_file);
            van_an_read(&list_flags, 1, van_an_dump_file);
            uint8_t inter_bf = merge_idx;

            // Pass each list.
            for (int list = 0; list < 2; list++)
            {
                cb->ref_idx[part_idx][list] = -1;
                if (!(list_flags&(1+list))) continue;

                // Reference index, MV and PMV index.
                uint8_t ref_idx, pmv_idx;
                f265_mv mv;
                van_an_read(&ref_idx, 1, van_an_dump_file);
                van_an_read(&mv.x, 2, van_an_dump_file);
                van_an_read(&mv.y, 2, van_an_dump_file);
                van_an_read(&pmv_idx, 1, van_an_dump_file);
                cb->mv[part_idx][list] = mv;
                cb->ref_idx[part_idx][list] = ref_idx;
                inter_bf |= pmv_idx<<(3+list);
            }

            cb->inter_bf[part_idx] = inter_bf;

            #ifdef VAN_LOAD_HM_ME
            // Load inter search data from HM dump file.
            cb->hm_ref_list[part_idx] = (list_flags-1)&3;
            cb->hm_cb_mode = cb_mode;

            for (int list = 0; list < 2; list++)
            {
                for (int ref=0; ref < 16; ref++)
                {
                    f265_mv mvp_list[2];

                    van_an_read(&mvp_list[0].x, 2, van_an_dump_file);
                    van_an_read(&mvp_list[0].y, 2, van_an_dump_file);
                    van_an_read(&mvp_list[1].x, 2, van_an_dump_file);
                    van_an_read(&mvp_list[1].y, 2, van_an_dump_file);

                    cb->hm_mvp_list[part_idx][list][ref][0].p = mvp_list[0].p;
                    cb->hm_mvp_list[part_idx][list][ref][1].p = mvp_list[1].p;
                }

                int lambda_me;
                int pmode_bits, num_of_ref;
                int8_t mvd_l1_zero_flag;

                van_an_read(&lambda_me, sizeof(int), van_an_dump_file);
                van_an_read(&mvd_l1_zero_flag, sizeof(int8_t), van_an_dump_file);
                van_an_read(&pmode_bits, sizeof(int), van_an_dump_file);
                van_an_read(&num_of_ref, sizeof(int), van_an_dump_file);

                cb->hm_lambda_me = lambda_me;
                cb->hm_mvd_l1_zero_flag = mvd_l1_zero_flag;
                cb->hm_pmode_bits[part_idx][list] = pmode_bits;
                cb->hm_num_of_ref[part_idx][list] = num_of_ref;

                // Read only from list 0 (same data as list 1).
                if (!list)
                    for (int l=0; l<3; l++)
                    {
                        int val;
                        van_an_read(&val, sizeof(int), van_an_dump_file);
                        cb->hm_cost[part_idx][l] = val;
                        van_an_read(&val, sizeof(int), van_an_dump_file);
                        cb->hm_bits[part_idx][l] = val;
                        van_an_read(&val, sizeof(int), van_an_dump_file);
                        cb->hm_ref_idx[part_idx][l] = val;
                    }
            }
            #endif
        }
    }

    // Intra chroma mode.
    if (cb_mode == 0)
    {
        uint8_t intra_mode;
        van_an_read(&intra_mode, 1, van_an_dump_file);
        cb->intra_chroma_mode = intra_mode;
    }

    // Load the transform tree.
    uint8_t *cb_tn = t->tt.tn;
    van_load_an_tt(t);

    // Correct the HM TT node if the 64x64 TB is unsplit so that we don't end up
    // trying to reconstruct 64x64 TBs. This happens when the YUV residual is 0.
    if (cb->lg_bs == 6 && !*cb_tn)
    {
        t->tt.tn = cb_tn;
        *t->tt.tn++ = 8;
        for (int i = 0; i < 4; i++) *t->tt.tn++ = 0;
    }

    // Update the prediction map.
    venc_update_pmap_unsplit_cb(t, cb);
}

// Load the CTB analysis.
void venc_load_hm_ctb_analysis(f265_enc_thread *t)
{
    venc_open_dump_file(&van_an_dump_file, VAN_DUMP_CTB_ANALYSIS_PATH, "rb");
    venc_init_transform_tree(t);
    van_load_an_cb_qt(t, t->cb);
    uint32_t sanity, expected = 0x42ffff42;
    van_an_read(&sanity, 4, van_an_dump_file);
    if (sanity != expected)
    {
        printf("Sanity mismatch, expected %x got %x.\n", expected, sanity);
        exit(1);
    }
}
#endif

// Parse the HM gop-like file.
void venc_parse_hm_gop_file(f265_enc *enc, char *path)
{
    enc->gd.hm_gop_compat_flag = 1;

    // Track the most frequent nb_active value below.
    uint8_t nb_active_freqs[16];
    memset(nb_active_freqs, 0, 16);

    FILE *file = fopen(path, "rb");
    if (!file)
    {
        printf("Failed to open %s: %s.\n", path, strerror(errno));
        exit(1);
    }

    int *gop_size = &enc->gd.hm_gop_size;
    while (*gop_size < 16 && !feof(file))
    {
        // Parse the current line. Skip it silently if the pattern doesn't
        // match.
        f265_hm_gop_entry *entry = enc->gd.hm_gop + *gop_size;
        int num, poc, qp_off, tc, beta, id, nb_active, nb_refs;
        char type;
        double qp_fac;
        if (fscanf(file, "Frame%d: %c %d %d %lf %d %d %d %d %d",
                   &num, &type, &poc, &qp_off, &qp_fac, &tc, &beta, &id, &nb_active, &nb_refs) == 10)
        {
            if (num != *gop_size + 1) { printf("Invalid GOP frame number %d.\n", num); exit(1); }
            if (type == 'I') entry->type = F265_FRAME_I;
            else if (type == 'P') entry->type = F265_FRAME_P;
            else if (type == 'B') entry->type = F265_FRAME_B;
            else { printf("Invalid frame type %c.\n", type); exit(1); }
            entry->poc_off = poc;
            entry->qp_offset = qp_off;
            entry->qp_factor = qp_fac;
            entry->nb_active = nb_active;
            nb_active_freqs[nb_active]++;
            entry->nb_refs = nb_refs;
            // Kludge: we assume P frames are used for reference. Other cases
            // are handled below.
            entry->ref_flag = type == 'P';
            for (int i = 0; i < nb_refs; i++)
                if (fscanf(file, "%d", entry->refs + i) != 1) { printf("Invalid reference POC.\n"); exit(1); }
            (*gop_size)++;
        }

        // Skip the rest of the line.
        while (1)
        {
            int c = fgetc(file);
            if (c == '\n' || c == EOF) break;
        }
    }

    fclose(file);

    // Update default_nb_ref_idx.
    int max_freq = -1;
    for (int i = 0; i < 16; i++)
    {
        if (nb_active_freqs[i] > max_freq)
        {
            max_freq = nb_active_freqs[i];
            enc->gd.default_nb_ref_idx[0] = enc->gd.default_nb_ref_idx[1] = i;
        }
    }

    // Update nb_reordered_frames. We compute the delay before outputting each
    // frame. We also update ref_flag.
    int reorder = 0;
    for (int i = 0; i < *gop_size; i++)
    {
        f265_hm_gop_entry *entry = enc->gd.hm_gop + i;
        reorder = F265_MAX(i+1 - entry->poc_off, reorder);
        for (int j = 0; j < entry->nb_refs; j++)
        {
            int ref_poc = entry->poc_off + entry->refs[j];
            while (ref_poc < 0) ref_poc += *gop_size;
            for (int k = 0; k < *gop_size; k++)
            {
                f265_hm_gop_entry *e2 = enc->gd.hm_gop + k;
                if (e2->poc_off == ref_poc) e2->ref_flag = 1;
            }
        }
    }
    enc->gd.nb_reordered_frames = reorder;
}

#ifdef VAN_VALIDATE_MODE_PRED
// Import and validate inter merge candidates extracted with HM decoder.
static void venc_hm_validate_merge(f265_enc_thread *t, f265_cb *cb, f265_inter_neighbour_mv *merge_cand, int part_idx)
{
    static FILE *merge_file = 0;
    venc_open_dump_file(&merge_file, "/tmp/merge.bin", "rb");

    uint16_t val16;

    van_an_read(&val16,2, merge_file);
    if (val16 != t->src_frame->abs_poc)
    {
        printf("[ERR] MCD - Invalid POC (f265:%li hm:%i)\n", t->src_frame->abs_poc, val16);
        exit(1);
    }

    van_an_read(&val16,2, merge_file);
    if (val16 != (t->ctb_off[0] + cb->cb_off[0]))
    {
        printf("[ERR] MCD - Invalid position X (f265:%i hm:%i)\n", t->ctb_off[0] + cb->cb_off[0], val16);
        exit(1);
    }

    van_an_read(&val16,2, merge_file);
    if (val16 != (t->ctb_off[1] + cb->cb_off[1]))
    {
        printf("[ERR] MCD - Invalid position Y (f265:%i hm:%i)\n", t->ctb_off[1] + cb->cb_off[1], val16);
        exit(1);
    }

    uint8_t val;
    van_an_read(&val,1, merge_file);
    if (val != part_idx)
    {
        printf("[ERR] MCD - Invalid partition index (f265:%i hm:%i)\n", part_idx, val);
        exit(1);
    }

    uint8_t nb_cand;
    van_an_read(&nb_cand, 1, merge_file);
    if (nb_cand == 0)
        return;

    if (nb_cand != t->enc->gd.merge_cand)
    {
        printf("[ERR] Invalid number of merge candidate. HM: %i. F265: %i\n", nb_cand, t->enc->gd.merge_cand);
        exit(1);
    }

    for (int i = 0; i < nb_cand; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            int8_t ref_idx;
            van_an_read(&ref_idx, 1, merge_file);

            if (ref_idx >= t->src_frame->nb_ref_idx[j])
            {
                printf("[ERR] Invalid merge index (%i, max ref=%i)\n", ref_idx, t->src_frame->nb_ref_idx[j]);
                exit(1);
            }

            f265_mv mv;
            van_an_read(&mv.x, 2, merge_file);
            van_an_read(&mv.y, 2, merge_file);

            if (ref_idx != merge_cand[i].ref_idx[j] || (ref_idx > -1  && (mv.p != merge_cand[i].mv[j].p)))
            {
                printf("[ERR] Invalid merge candidate %i L%i. f265 %i %i:%i. hm %i %i:%i\n",
                       i, j,
                       merge_cand[i].ref_idx[j], merge_cand[i].mv[j].x, merge_cand[i].mv[j].y,
                       ref_idx, mv.x, mv.y);
                exit(1);
            }
        }
    }

    int expected_sanity = 0x42ffff42;
    int sanity;
    van_an_read(&sanity, 4, merge_file);
    if (expected_sanity != sanity)
    {
        printf("[ERR] MCD - Sanity check failed: %i\n", sanity);
        exit(1);
    }
}

// Import and validate inter PMV extracted with HM Encoder.
static void venc_hm_validate_pmv(f265_enc_thread *t, f265_cb *cb, int part_idx, int list_idx, int ref_idx, f265_mv pmv[2])
{
    static FILE *pmv_file = 0;
    venc_open_dump_file(&pmv_file, "/tmp/pmv.bin", "rb");

    uint16_t val16;

    van_an_read(&val16, 2, pmv_file);
    if (val16 != t->src_frame->abs_poc)
    {
        printf("[ERR] PMV - Invalid POC (f265:%li hm:%i)\n", t->src_frame->abs_poc, val16);
        exit(1);
    }

    van_an_read(&val16, 2, pmv_file);
    if (val16 != (t->ctb_off[0] + cb->cb_off[0]))
    {
        printf("[ERR] PMV - Invalid position X (f265:%i hm:%i)\n", t->ctb_off[0] + cb->cb_off[0], val16);
        exit(1);
    }

    van_an_read(&val16,2, pmv_file);
    if (val16 != (t->ctb_off[1] + cb->cb_off[1]))
    {
        printf("[ERR] PMV - Invalid position Y (f265:%i hm:%i)\n", t->ctb_off[1] + cb->cb_off[1], val16);
        exit(1);
    }

    uint8_t val;
    van_an_read(&val,1, pmv_file);
    if (val != part_idx)
    {
        printf("[ERR] PMV - Invalid partition index (f265:%i hm:%i)\n", part_idx, val);
        exit(1);
    }

    uint8_t nb_cand;
    van_an_read(&nb_cand, 1, pmv_file);
    if (nb_cand > 0)
    {
        if (nb_cand != 2)
        {
            printf("[ERR] Invalid number of PMV. hm %i. f265 2\n", nb_cand);
            exit(1);
        }

        for (int j = 0; j < nb_cand; j++)
        {
            f265_mv mv;
            van_an_read(&mv.x, 2, pmv_file);
            van_an_read(&mv.y, 2, pmv_file);

            if (mv.p != pmv[j].p)
            {
                printf("[ERR] Invalid PMV (%i): f265 %i:%i. HM %i:%i\n", j, pmv[j].x, pmv[j].y, mv.x, mv.y);
                exit(1);
            }
        }
    }

    int sanity;
    int expected_sanity = 0x42ffff42;

    van_an_read(&sanity, 4, pmv_file);

    if (sanity != expected_sanity)
    {
        printf("[ERR] PMV - Failed PMV sanity check %x\n", sanity);
        exit(1);
    }
}

// Import and validate intra MPM extracted with HM Decoder.
static void venc_hm_validate_intra(f265_enc_thread *t, f265_cb *cb, int part_idx, uint8_t *pred_mode)
{
    static FILE *pmv_file = 0;
    venc_open_dump_file(&pmv_file, "/tmp/intra.bin", "rb");

    uint8_t intra_mode[3];
    for (int i = 0; i < 3; i++)
        van_an_read(intra_mode + i, 1, pmv_file);

    for (int i = 0; i < 3; i++)
    {
        if (pred_mode[i] != intra_mode[i])
        {
            printf("[Err] Invalid prediction %i: f265 %u, %u, %u, hm %u, %u, %u\n",
                    i,
                    pred_mode[0],   pred_mode[1],   pred_mode[2],
                    intra_mode[0], intra_mode[1], intra_mode[2]);
            exit(1);
        }
    }

    int sanity;
    int expected_sanity = 0x42ffff42;
    van_an_read(&sanity, 4, pmv_file);
    if (sanity != expected_sanity)
    {
        printf(" [ERR] Failed PMV sanity check %x\n", sanity);
        exit(1);
    }
}

// Generate and test intra MPM.
static void venc_hm_test_cb_intra_pred(f265_enc_thread *t, f265_cb *cb)
{
    int partition_count = f265_nb_parts[cb->inter_part];

    venc_update_pmap_unsplit_cb(t, cb);
    for (int part_idx = 0; part_idx < partition_count; part_idx++)
    {
        uint8_t pred_mode[3];
        venc_get_intra_pred_mode(t, cb, part_idx, pred_mode);
        venc_hm_validate_intra(t, cb, part_idx, pred_mode);
    }
}

// Generate and test inter PMV.
static void venc_hm_test_pmv(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    f265_inter_neighbour_mv neighbour[6];
    venc_get_neighbour_mvs(t, cb, part_idx, neighbour);

    for (int list_idx = 0; list_idx < t->src_frame->nb_lists; list_idx++)
    {
        int ref_idx = cb->ref_idx[part_idx][list_idx];

        if (ref_idx == -1) continue;

        f265_mv pmv[2];
        venc_get_pmv(t, part_idx, neighbour, ref_idx, list_idx, pmv);
        venc_hm_validate_pmv(t, cb, part_idx, list_idx, ref_idx, pmv);
    }
}

// Generate and test inter merge candidates.
static void venc_hm_test_merge(f265_enc_thread *t, f265_cb *cb, int part_idx)
{
    f265_inter_neighbour_mv neighbour[6];
    venc_get_neighbour_mvs(t, cb, part_idx, neighbour);

    f265_inter_neighbour_mv merge_cand[5];
    venc_get_merge_candidate(t, cb, part_idx, neighbour, merge_cand);
    venc_hm_validate_merge(t, cb, merge_cand, part_idx);
}

// Validate the prediction on the CB. Split the CB if required.
static void venc_hm_validate_internal(f265_enc_thread *t, f265_cb *cb)
{
    if (cb->flags & F265_CB_SPLIT)
    {
        // Recursively call ourselves if the block is splitted.
        for (int i = 0; i < 4; i++) venc_hm_validate_internal(t, t->cb + cb->child_idx + i);
        return;
    }

    if (cb->flags & F265_CB_FORBIDDEN) return;

    // Test intra CB.
    if (cb->flags & F265_CB_INTRA)
    {
        // Skip PCM.
        if (cb->intra_luma_mode[0] != 35)
            venc_hm_test_cb_intra_pred(t, cb);

        venc_update_pmap_unsplit_cb(t, cb);
    }

    // Test inter CB.
    else
    {
        int partition_count = f265_nb_parts[cb->inter_part];
        for (int part_idx = 0; part_idx < partition_count; ++part_idx)
        {
            // Test the PMV.
            if ((cb->inter_bf[part_idx] & 7) == 5)
                venc_hm_test_pmv(t, cb, part_idx);

            // Test the merge candidate.
            // They are always generated, even if the merge mode is not used.
            venc_hm_test_merge(t, cb, part_idx);

            // Make sure the current partion is available to the next partition.
            venc_update_pmap_unsplit_cb(t, cb);
        }
    }

    venc_save_cb(t, cb);
}

// Validate inter predicted mode, inter merge candidate and intre MPM on the CB.
void venc_hm_validate(f265_enc_thread *t)
{
    if (t->src_frame->frame_type != F265_FRAME_I)
    {
        venc_frame_set_up_inter_pred(t);
        venc_ctb_set_up_inter_pred(t);
    }

    // Reset the transform node pointer to save the deblocking information.
    t->tt.tn = t->tmap;

    venc_hm_validate_internal(t, t->cb);
}
#endif

// Reorder the HM GOP.
void venc_hm_reorder_gop(f265_enc *enc, int *nb_gop, int *key_flag)
{
    f265_lookahead *la = &enc->la;
    f265_frame **display = la->display;
    f265_frame **p = la->coded + la->nb_committed;
    int nb_frames = la->nb_undecided;
    int write_pos = 0;

    // Special case for the first frame.
    if (!display[1]->abs_poc) return;

    // Pass every GOP entry in order.
    for (int i = 0; i < enc->gd.hm_gop_size; i++)
    {
        f265_hm_gop_entry *entry = enc->gd.hm_gop + i;

        // Skip the entry if it's beyond the end of the stream.
        if (entry->poc_off > nb_frames) continue;

        // Select the frame.
        f265_frame *f = display[entry->poc_off];
        f->gop_entry = entry;
        f->la_frame_type = entry->type;
        F265_SET_FLAG(f->gen_flags, F265_FF_REF, 1);
        p[write_pos++] = f;
    }

    // All frames become part of the current GOP.
    *nb_gop = nb_frames;

    // Force the first GOP frame to become a key frame as necessary.
    *key_flag = la->key_countdown == nb_frames;
}

// Test intra bit-exactitude with HM.
#ifdef VAN_DUMP_INTRA
// Set the intra information for the CTB.
void venc_hm_set_intra(f265_enc_thread *t)
{
    f265_cb *cb = t->cb;
    F265_SET_FLAG(cb->flags, F265_CB_INTRA, 1);
    cb->intra_luma_mode[0] = cb->intra_chroma_mode = 26;
    cb->intra_luma_mode[1] = -1;
}
#endif

// Future implementation notes.

// Assembly for intra angular:
// - Consider sharing functions between block sizes.
// - Consider doing a function that will do -135, -90, -45, 0, 45 degrees in one
//   shot (fast estimation). May also do planar + DC. Handle or ignore filtering
//   issues.
// - Consider computing the SAD instead of storing rows (fast estimation).
// - Consider using pshufb based on the inverse angle to project neighbours.
// - If horizontal, translate the horizontal case to the vertical case (flip
//                  the neighbours).
// - Call the appropriate dispatch function for the projection angle (3 or 4
//   cases).
//   - If the projection is vertical:
//     - Broadcast the pixels.
//     - When filtering:
//       - Compute column filtering in a register.
//       - Replace the first byte of every row using palignr and vpblendvb.
//   - Else if the projection is -45 or 45 degrees:
//     - Find a way to use palignr to avoid loads.
//   - Else:
//     - Go load heavy. Load neigbhours from cache for every row at the computed
//       offset to avoid branches.
//     - Pack neighbours as ABBCCDDE. Use pmaddubsw to multiply by the fractions
//       (32-iFact, iFact). Shift and pack two rows at a time.
//     - Assume the fractional case even when iFact==0 for a row to avoid
//       branches.
// - If horizontal, flip using punpck.

// Reconst:8 -> unfiltered/filtered:tmp@16,out@8 -> prediction:tmp@16,out@8 ->
// (src-pred):tmp@16,out@16 -> DCT1D&clip:tmp@32,out@16 -> DCT1D&clip:tmp@32,out@16 ->
// quantization:tmp@32,out@16 -> dequantization:same ->
// DCT1D&clip:tmp@32,out@16 -> DCT1D&add&clip:tmp@32,out@8
//
// Assuming high bit depth:
// Same, but:
// - unfiltered/filtered/prediction: double "tmp" sizes.
// - input/output: double size.

// Coefficient encoding:
//
// One quantization function per transform block size.
// Special case:
// - May use a function to process 4 4x4 blocks together.
// - May split between a function that does pure quant and a function that does
//   quant + nz_flags stuff.
// Quantization function assembly:
// - Quantize every row in raster scan.
// - PACKSS columns together to get 1 byte per 4-coeff group.
// - OR rows to get 1 byte per subblock.
// - PCMPEQ to get 1 flag per subblock.
// - PSHUFB/reorder to get the raster order and the encode order (nz_flags[0..1]).
// - PMOVMSKB to get the flags in general purpose.
// - Compute tb->nz_flags[2..3]:
//   - Use nz_flags[0] as base.
//   - Below: PDEP all but the bottom row.
//   - Right: PEXT to remove first column, PDEP to add zeroes to right column.
//
// One preprocessing function for all block sizes.
// Preprocessing function assembly:
// - Get the coefficient stride from lg_bs.
// - Set up vpgatherqq index register.
// - Set up function pointer to reorder the coefficients.
// - BSF over the non-zero subblocks:
//   - first_coeff_pos = coeff_off[pos].
//   - Load coefficients using vpgatherqq (clone the index register first).
//   - Call func to reorder the 16 coefficients in YMM in encode order:
//     - PSHUFB to reorder each lane.
//     - Swap lanes in new register.
//     - PALIGNR to the merge the high-low lanes together.
//     - PSHUFB to reorder each lane.
//   - (Copy coeffs so we don't lose their content as needed here).
//   - PACKSSWD to have one byte per coeff.
//   - PCMPGTB -1 to get mask of signs.
//   - PABSB to get absolute coeff values (bytes).
//   - PCMPGTB 0 to get non-zero coeff flags.
//   - PCMPGTB 1 to get greater-than-1 flags.
//   - PCMPGTB 2 to get greater-than-2 flags.
//   - PMOVMSKB to nz_flags.
//   - PMOVMSKB to signs.
//   - PMOVMSKB to gt1.
//   - PMOVMSKB to gt2.
//   - PEXT signs using nz_flags.
//   - PEXT gt1 using nz_flags.
//   - PEXT gt2 using nz_flags.
//   - AND gt1, 0xff (keep only first 8 bits).
//   - BLSI TMP, gt1 (extract first gt1 bit set, if any).
//   - AND gt2, TMP (extract first gt2 bit).
//   - SETNE gt2 (set gt2 bit value, 0 if it doesn't exist).
//   - POPCNT nb_nz, nz_flags.
//     (There are remaining levels if nb_gt1 + gt2 > 1 || nb_nz > 8).
//     (Alternatively, check if remain_flags is non-zero).
//   - POPCNT nb_gt1, gt1
//   - ADD nb_gt1, gt2
//   - CMP nb_gt1, 1
//   - SETGT gt1_overflow
//   - CMP nb_nz, 8
//   - SETGT nb_nz_overflow
//   - OR remain_flag, gt1_overflow, nb_nz_overflow
//   - Possibly branch to compute coefficient base levels below.
//   - Else, just store the absolute coeff in order:
//     - PABSW 16-bit coefficients.
//     - MOVDQA 16-bit coefficients in store location.
//     - remain_flag << 5 (update store location if the store was required).
//     - ADD store pointer, remain_flag.
//   - SB pointer += sizeof(SB).
//   - How to reverse the 16 signs:
//     - BSWAP.
//     - PDEP to insert 4 zeros per group of 4 bits.
//     - MOV to YMM.
//     - PHUFB.
//     - MOV to GPR.
//     - PEXT to extract the reordered bits.
//     - Shift right to align with bit 0.
//     - Other possibility:
//       - Broadcast in YMM.
//       - AND with bit identification mask.
//       - Mov high lane to low.
//       - PACKSSWB.
//       - PSHUFB.
//       - PCMPEQ
//       - PMOVMSKB.
//     - Other possibility:
//       - Process in group of four in C loop (map).
//    - How to put back the gt1 flags in place:
//       - PDEP using the coefficient non-zero flags.
//       - Can be done for gt2 directly.
//   // Branch to compute exactly the base levels, slow so don't do unless
//   // required. NOT VERY USEFUL, the base levels must be computed anyway for the
//   // rice param update.
//   - PEXT A, 0xff, nz_flags (remove 1 for the 8 first non-zero coefficients).
//   - VBROADCASTW A       (move coefficient bits back to ymm).
//   - PAND A, [word_mask] (keep only the bit corresponding to the coeff).
//   - PCMPEQW A, 0        (-1 if the bit is not set).
//   - PCMPEQW A, 0        (-1 if the bit is set).
//   - PADDW coeffs, A     (remove 1 for the 8 first non-zero coefficients).
//   - PEXT  A, gt2, nz_flags (put back gt2 bit at its coeff location).
//   - (same instructions as above to remove 1 for the gt2 coefficient).
//   - PACKSSWD B, A to have one byte per coeff.
//   - PCMPEQB B, 0 to get mask of remaining coeffs (2x).
//   - PMOVMSKB B to remain_flags.
//   - PADDW A, -1                (remove 1 for each non-zero coefficient).
//   - MOVDQA A, [remain]         (store remaining coefficient level).
//   - Update store pointer, etc.

