// Generate the CABAC context initialization tables and context offsets.

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>


///////////////////////////////////////////////////////////////////////////////
// Code imported from HM. It was modified to remove their idiosyncrasies.

#define MAX_NUM_CTX_MOD             512       ///< maximum number of supported contexts
#define NUM_SPLIT_FLAG_CTX            3       ///< number of context models for split flag
#define NUM_SKIP_FLAG_CTX             3       ///< number of context models for skip flag
#define NUM_MERGE_FLAG_EXT_CTX        1       ///< number of context models for merge flag of merge extended
#define NUM_MERGE_IDX_EXT_CTX         1       ///< number of context models for merge index of merge extended
#define NUM_PART_MODE_CTX             4       // ==> ADDED THIS.
#define NUM_PART_SIZE_CTX             4       ///< number of context models for partition size ==> REMOVED THIS.
#define NUM_CU_AMP_CTX                1       ///< number of context models for partition size (AMP) ==> REMOVED THIS.
#define NUM_PRED_MODE_CTX             1       ///< number of context models for prediction mode
#define NUM_ADI_CTX                   1       ///< number of context models for intra prediction
#define NUM_CHROMA_PRED_CTX_FIXED     1       // ==> ADDED THIS.
#define NUM_CHROMA_PRED_CTX           2       ///< number of context models for intra prediction (chroma) ==> REMOVED THIS.
#define NUM_INTER_DIR_CTX             5       ///< number of context models for inter prediction direction
#define NUM_MV_RES_CTX                2       ///< number of context models for motion vector difference ==> REMOVED THIS.
#define NUM_REF_NO_CTX                2       ///< number of context models for reference index
#define NUM_TRANS_SUBDIV_FLAG_CTX     3       ///< number of context models for transform subdivision flags
#define NUM_QT_CBF_LUMA_CTX           2       // ==> ADDED THIS.
#define NUM_QT_CBF_CHROMA_CTX         4       // ==> ADDED THIS.
#define NUM_QT_CBF_CTX                5       ///< number of context models for QT CBF ==> REMOVED THIS.
#define NUM_QT_ROOT_CBF_CTX           1       ///< number of context models for QT ROOT CBF
#define NUM_DELTA_QP_CTX_FIXED        2       // ==> ADDED THIS.
#define NUM_DELTA_QP_CTX              3       ///< number of context models for dQP ==> REMOVED THIS.
#define NUM_SIG_CG_FLAG_CTX           4       ///< this is coded_sub_block ==> FIXED NUMBER.
#define NUM_SIG_FLAG_CTX              42      ///< number of context models for sig flag
#define NUM_SIG_FLAG_CTX_LUMA         27      ///< number of context models for luma sig flag ==> UNUSED.
#define NUM_SIG_FLAG_CTX_CHROMA       15      ///< number of context models for chroma sig flag ==> UNUSED.
#define NUM_CTX_LAST_FLAG_XY_FIXED    18      // ==> ADDED THIS.
#define NUM_CTX_LAST_FLAG_XY          15      ///< number of context models for last coefficient position ==> REMOVED THIS.
#define NUM_ONE_FLAG_CTX              24      ///< number of context models for greater than 1 flag
#define NUM_ONE_FLAG_CTX_LUMA         16      ///< number of context models for greater than 1 flag of luma ==> UNUSED.
#define NUM_ONE_FLAG_CTX_CHROMA        8      ///< number of context models for greater than 1 flag of chroma ==> UNUSED.
#define NUM_ABS_FLAG_CTX               6      ///< number of context models for greater than 2 flag
#define NUM_ABS_FLAG_CTX_LUMA          4      ///< number of context models for greater than 2 flag of luma ==> UNUSED.
#define NUM_ABS_FLAG_CTX_CHROMA        2      ///< number of context models for greater than 2 flag of chroma ==> UNUSED.
#define NUM_MVP_IDX_CTX_FIXED         1       // ==> ADDED THIS.
#define NUM_MVP_IDX_CTX               2       ///< number of context models for MVP index ==> REMOVED THIS.
#define NUM_SAO_MERGE_FLAG_CTX        1       ///< number of context models for SAO merge flags
#define NUM_SAO_TYPE_IDX_CTX          1       ///< number of context models for SAO type index
#define NUM_TRANSFORMSKIP_FLAG_CTX    2       ///< number of context models for transform skipping  ==> FIXED NUMBER.
#define NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX  1
#define CNU                          154      ///< dummy initialization value for unused context models 'Context model Not Used'

#define UChar uint8_t

// HM uses the initialization order B, P, I (the reverse order of the spec).

static const UChar
INIT_CU_TRANSQUANT_BYPASS_FLAG[3][NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX] =
{
  { 154 },
  { 154 },
  { 154 },
};

static const UChar
INIT_SPLIT_FLAG[3][NUM_SPLIT_FLAG_CTX] =
{
  { 107,  139,  126, },
  { 107,  139,  126, },
  { 139,  141,  157, },
};

static const UChar
INIT_SKIP_FLAG[3][NUM_SKIP_FLAG_CTX] =
{
  { 197,  185,  201, },
  { 197,  185,  201, },
  { CNU,  CNU,  CNU, },
};

static const UChar
INIT_MERGE_FLAG_EXT[3][NUM_MERGE_FLAG_EXT_CTX] =
{
  { 154, },
  { 110, },
  { CNU, },
};

static const UChar
INIT_MERGE_IDX_EXT[3][NUM_MERGE_IDX_EXT_CTX] =
{
  { 137, },
  { 122, },
  { CNU, },
};

// ==> ADDED. NOTE: swapping contexts 2 and 3 here (the initialization values
// are the same so the initialization is HM-compatible.
static const UChar
INIT_PART_MODE[3][NUM_PART_MODE_CTX] =
{
  { 154,  139,  154,  154, },
  { 154,  139,  154,  154, },
  { 184,  CNU,  CNU,  CNU, },
};

// ==> REMOVED.
static const UChar
INIT_PART_SIZE[3][NUM_PART_SIZE_CTX] =
{
  { 154,  139,  CNU,  CNU, },
  { 154,  139,  CNU,  CNU, },
  { 184,  CNU,  CNU,  CNU, },
};

// ==> REMOVED.
static const UChar
INIT_CU_AMP_POS[3][NUM_CU_AMP_CTX] =
{
  { 154, },
  { 154, },
  { CNU, },
};

static const UChar
INIT_PRED_MODE[3][NUM_PRED_MODE_CTX] =
{
  { 134, },
  { 149, },
  { CNU, },
};

static const UChar
INIT_INTRA_PRED_MODE[3][NUM_ADI_CTX] =
{
  { 183, },
  { 154, },
  { 184, },
};

// ==> ADDED.
static const UChar
INIT_CHROMA_PRED_MODE_FIXED[3][NUM_CHROMA_PRED_CTX_FIXED] =
{
  { 152 },
  { 152 },
  {  63 },
};

// ==> REMOVED.
static const UChar
INIT_CHROMA_PRED_MODE[3][NUM_CHROMA_PRED_CTX] =
{
  { 152,  139, },
  { 152,  139, },
  {  63,  139, },
};

static const UChar
INIT_INTER_DIR[3][NUM_INTER_DIR_CTX] =
{
  {  95,   79,   63,   31,  31, },
  {  95,   79,   63,   31,  31, },
  { CNU,  CNU,  CNU,  CNU, CNU, },
};

// ==> REMOVED.
static const UChar
INIT_MVD[3][NUM_MV_RES_CTX] =
{
  { 169,  198, },
  { 140,  198, },
  { CNU,  CNU, },
};

// ==> ADDED.
static const UChar
INIT_MVD_GT0[3][1] =
{
  { 169 },
  { 140 },
  { CNU },
};

// ==> ADDED.
static const UChar
INIT_MVD_GT1[3][1] =
{
  { 198 },
  { 198 },
  { CNU },
};

static const UChar
INIT_REF_PIC[3][NUM_REF_NO_CTX] =
{
  { 153,  153 },
  { 153,  153 },
  { CNU,  CNU },
};

// ==> ADDED.
static const UChar
INIT_DQP_FIXED[3][NUM_DELTA_QP_CTX_FIXED] =
{
  { 154,  154 },
  { 154,  154 },
  { 154,  154 },
};

// ==> REMOVED.
static const UChar
INIT_DQP[3][NUM_DELTA_QP_CTX] =
{
  { 154,  154,  154, },
  { 154,  154,  154, },
  { 154,  154,  154, },
};

// ==> ADDED.
static const UChar
INIT_QT_LUMA_CBF[3][NUM_QT_CBF_LUMA_CTX] =
{
  { 153, 111 },
  { 153, 111 },
  { 111, 141 },
};

// ==> ADDED.
static const UChar
INIT_QT_CHROMA_CBF[3][NUM_QT_CBF_CHROMA_CTX] =
{
  { 149, 92, 167, 154 },
  { 149, 107, 167, 154 },
  { 94, 138, 182, 154 },
};

// ==> REMOVED.
static const UChar
INIT_QT_CBF[3][2*NUM_QT_CBF_CTX] =
{
  { 153,  111,  CNU,  CNU,  CNU,  149,   92,  167,  CNU,  CNU, },
  { 153,  111,  CNU,  CNU,  CNU,  149,  107,  167,  CNU,  CNU, },
  { 111,  141,  CNU,  CNU,  CNU,   94,  138,  182,  CNU,  CNU, },
};

static const UChar
INIT_QT_ROOT_CBF[3][NUM_QT_ROOT_CBF_CTX] =
{
  {  79, },
  {  79, },
  { CNU, },
};

// ==> ADDED.
static const UChar
INIT_LAST_FIXED[3][NUM_CTX_LAST_FLAG_XY_FIXED] =
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79, 108,  123,   93 },
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94, 108,  123,  108 },
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79, 108,  123,   63 },
};

// ==> REMOVED.
static const UChar
INIT_LAST[3][2*NUM_CTX_LAST_FLAG_XY] =
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79,
    108,  123,   93,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,
  },
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94,
    108,  123,  108,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,
  },
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79,
    108,  123,   63,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,
  },
};

// ==> FIXED NUMBER.
static const UChar
INIT_SIG_CG_FLAG[3][NUM_SIG_CG_FLAG_CTX] =
{
  { 121,  140,
    61,  154,
  },
  { 121,  140,
    61,  154,
  },
  {  91,  171,
    134,  141,
  },
};

static const UChar
INIT_SIG_FLAG[3][NUM_SIG_FLAG_CTX] =
{
  { 170,  154,  139,  153,  139,  123,  123,   63,  124,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  170,  153,  138,  138,  122,  121,  122,  121,  167,  151,  183,  140,  151,  183,  140,  },
  { 155,  154,  139,  153,  139,  123,  123,   63,  153,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  166,  183,  140,  136,  153,  154,  170,  153,  123,  123,  107,  121,  107,  121,  167,  151,  183,  140,  151,  183,  140,  },
  { 111,  111,  125,  110,  110,   94,  124,  108,  124,  107,  125,  141,  179,  153,  125,  107,  125,  141,  179,  153,  125,  107,  125,  141,  179,  153,  125,  140,  139,  182,  182,  152,  136,  152,  136,  153,  136,  139,  111,  136,  139,  111,  },
};

static const UChar
INIT_ONE_FLAG[3][NUM_ONE_FLAG_CTX] =
{
  { 154,  196,  167,  167,  154,  152,  167,  182,  182,  134,  149,  136,  153,  121,  136,  122,  169,  208,  166,  167,  154,  152,  167,  182, },
  { 154,  196,  196,  167,  154,  152,  167,  182,  182,  134,  149,  136,  153,  121,  136,  137,  169,  194,  166,  167,  154,  167,  137,  182, },
  { 140,   92,  137,  138,  140,  152,  138,  139,  153,   74,  149,   92,  139,  107,  122,  152,  140,  179,  166,  182,  140,  227,  122,  197, },
};

static const UChar
INIT_ABS_FLAG[3][NUM_ABS_FLAG_CTX] =
{
  { 107,  167,   91,  107,  107,  167, },
  { 107,  167,   91,  122,  107,  167, },
  { 138,  153,  136,  167,  152,  152, },
};

// ADDED.
static const UChar
INIT_MVP_IDX_FIXED[3][NUM_MVP_IDX_CTX_FIXED] =
{
  { 168 },
  { 168 },
  { CNU },
};

// ==> REMOVED.
static const UChar
INIT_MVP_IDX[3][NUM_MVP_IDX_CTX] =
{
  { 168,  CNU, },
  { 168,  CNU, },
  { CNU,  CNU, },
};

static const UChar
INIT_SAO_MERGE_FLAG[3][NUM_SAO_MERGE_FLAG_CTX] =
{
  { 153,  },
  { 153,  },
  { 153,  },
};

// ==> REMOVED DEFINE.
static const UChar
INIT_SAO_TYPE_IDX[3][NUM_SAO_TYPE_IDX_CTX] =
{
  { 160, },
  { 185, },
  { 200, },
};

// ==> REMOVED DEFINE.
static const UChar
INIT_TRANS_SUBDIV_FLAG[3][NUM_TRANS_SUBDIV_FLAG_CTX] =
{
  { 224,  167,  122, },
  { 124,  138,   94, },
  { 153,  138,  138, },
};

// ==> FIXED NUMBER.
static const UChar
INIT_TRANSFORMSKIP_FLAG[3][NUM_TRANSFORMSKIP_FLAG_CTX] =
{
  { 139,  139},
  { 139,  139},
  { 139,  139},
};


///////////////////////////////////////////////////////////////////////////////
// Dump code.

// Represent a table in the spec.
struct table_info
{
    // Name of the syntax element. Null if the offset need not be printed.
    const char *name;

    // Pointer to the initialization table.
    uint8_t *init;

    // Size of the initialization table per frame type.
    int size;
};

// Represent a syntax category.
struct cat_info
{
    // Name of the syntax element.
    const char *name;

    // Number of entries.
    int size;
};

int main()
{
    // Reorder the tables in the spec order. Unnamed offsets indicate chroma /
    // cr / L1 / Y, i.e. the logical extension.
    #define NB_TABLES 28
    struct table_info tables[NB_TABLES] =
    {
        { "sao_merge", (uint8_t*)INIT_SAO_MERGE_FLAG, NUM_SAO_MERGE_FLAG_CTX },
        { "sao_type", (uint8_t*)INIT_SAO_TYPE_IDX, NUM_SAO_TYPE_IDX_CTX },

        { "split_cu", (uint8_t*)INIT_SPLIT_FLAG, NUM_SPLIT_FLAG_CTX },

        { "transquant_bypass", (uint8_t*)INIT_CU_TRANSQUANT_BYPASS_FLAG, NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX },
        { "cu_skip", (uint8_t*)INIT_SKIP_FLAG, NUM_SKIP_FLAG_CTX },
        { "pred_mode", (uint8_t*)INIT_PRED_MODE, NUM_PRED_MODE_CTX },
        { "part_mode", (uint8_t*)INIT_PART_MODE, NUM_PART_MODE_CTX },
        { "intra_luma_pred", (uint8_t*)INIT_INTRA_PRED_MODE, NUM_ADI_CTX },
        { "intra_chroma_pred", (uint8_t*)INIT_CHROMA_PRED_MODE_FIXED, NUM_CHROMA_PRED_CTX_FIXED },

        { "merge_flag", (uint8_t*)INIT_MERGE_FLAG_EXT, NUM_MERGE_FLAG_EXT_CTX },
        { "merge_idx", (uint8_t*)INIT_MERGE_IDX_EXT, NUM_MERGE_IDX_EXT_CTX },
        { "inter_pred", (uint8_t*)INIT_INTER_DIR, NUM_INTER_DIR_CTX },
        { "ref_idx", (uint8_t*)INIT_REF_PIC, NUM_REF_NO_CTX },
        { "mvp", (uint8_t*)INIT_MVP_IDX_FIXED, NUM_MVP_IDX_CTX_FIXED },

        { "split_transform", (uint8_t*)INIT_TRANS_SUBDIV_FLAG, NUM_TRANS_SUBDIV_FLAG_CTX },
        { "root_cbf", (uint8_t*)INIT_QT_ROOT_CBF, NUM_QT_ROOT_CBF_CTX },
        { "cbf_luma", (uint8_t*)INIT_QT_LUMA_CBF, NUM_QT_CBF_LUMA_CTX },
        { "cbf_chroma", (uint8_t*)INIT_QT_CHROMA_CBF, NUM_QT_CBF_CHROMA_CTX },

        { "mvd_greater0", (uint8_t*)INIT_MVD_GT0, 1 },
        { "mvd_greater1", (uint8_t*)INIT_MVD_GT1, 1 },

        { "cu_qp_delta", (uint8_t*)INIT_DQP_FIXED, NUM_DELTA_QP_CTX_FIXED },

        { "transform_skip", (uint8_t*)INIT_TRANSFORMSKIP_FLAG, NUM_TRANSFORMSKIP_FLAG_CTX },
        { "last_sig_coeff", (uint8_t*)INIT_LAST_FIXED, NUM_CTX_LAST_FLAG_XY_FIXED },
        { "", (uint8_t*)INIT_LAST_FIXED, NUM_CTX_LAST_FLAG_XY_FIXED },
        { "coded_sub_block", (uint8_t*)INIT_SIG_CG_FLAG, NUM_SIG_CG_FLAG_CTX },
        { "sig_coeff", (uint8_t*)INIT_SIG_FLAG, NUM_SIG_FLAG_CTX },
        { "coeff_greater1", (uint8_t*)INIT_ONE_FLAG, NUM_ONE_FLAG_CTX },
        { "coeff_greater2", (uint8_t*)INIT_ABS_FLAG, NUM_ABS_FLAG_CTX },
    };

    // Output the context offsets.
    int nb_offs = 0;
    for (int i = 0; i < NB_TABLES; i++)
    {
        struct table_info *t = tables + i;
        if (*t->name)
        {
            printf("#define F265_CO_");
            for (int j = 0; j < strlen(t->name); j++) printf("%c", toupper(t->name[j]));
            for (int j = 0; j < 24 - strlen(t->name); j++) printf(" ");
            printf("0x%x\n", nb_offs);
        }
        nb_offs += t->size;
    }

    // Output the number of offsets.
    printf("\n");
    printf("#define F265_NB_CABAC_CTX %d\n", nb_offs);
    printf("\n");

    // Output the context table for each frame type.
    printf("Table values follow:\n");
    for (int ft = 0; ft < 3; ft++)
    {
        // Get HM order.
        int hm_ft = 2 - ft;

        // Populate the table entries.
        int entries[1024];
        for (int i = 0, off = 0; i < NB_TABLES; i++)
        {
            struct table_info *t = tables + i;
            for (int j = 0; j < t->size; j++) entries[off + j] = t->init[t->size*hm_ft + j];
            off += t->size;
        }

        // Dump the entries, 16 per line.
        printf("    {\n");
        for (int i = 0; i < nb_offs; i++)
        {
            if (i % 16 == 0)
            {
                if (i) printf("\n");
                printf("        ");
            }
            if (i % 16 != 0) printf(" ");
            printf("%3d,", entries[i]);
        }
        printf("\n");
        printf("    },\n");
    }

    // Output the syntax element approximation categories.
    #define NB_SYN_CAT 13
    struct cat_info cat_tables[NB_SYN_CAT] =
    {
        { "SPLIT_CU", 3*2 },
        { "CU_SKIP", 3*2 },
        { "SPLIT_TRANSFORM", 3*2 },
        { "PRED_MODE", 2 },
        { "INTRA_PART", 2 },
        { "INTER_PART", 3+8 },
        { "INTRA_LUMA_MODE", 4 },
        { "INTRA_CHROMA_MODE", 2 },
        { "INTER_PRED", 3 },
        { "MERGE_IDX", 6 },
        { "MVP", 2 },
        { "REF_IDX", 2*16 },
        { "BYPASS", 1 },
    };

    printf("\n");
    nb_offs = 0;
    for (int i = 0; i < NB_SYN_CAT; i++)
    {
        struct cat_info *t = cat_tables + i;
        printf("#define F265_SE_");
        for (int j = 0; j < strlen(t->name); j++) printf("%c", toupper(t->name[j]));
        for (int j = 0; j < 24 - strlen(t->name); j++) printf(" ");
        printf("%d\n", nb_offs);
        nb_offs += t->size;
    }

    // Output the number of offsets.
    printf("#define F265_SE_SIZE                    %d\n", nb_offs);
    printf("\n");

    return 0;
}

