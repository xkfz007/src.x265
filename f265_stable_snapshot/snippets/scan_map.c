// Generate the scanning maps for each transform size and encoding order.
// Generate the non-zero coefficient flag context index derivation table.

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

// Block count: 1 + 4 + 16 + 64 = 85.
// Size required: 85 * 3 = 255. Round to 256 for the coefficient map.
uint8_t venc_scan_map_idx[4][3];
uint8_t venc_scan_map_data[2*256];
int write_pos = 0;

// Write the current map value.
void update_map(int lbs, int pos)
{
    int tb_bs = 4*(1<<lbs);
    int y = pos>>lbs;
    int x = pos&((1<<lbs)-1);
    int off = y*tb_bs + x;
    venc_scan_map_data[write_pos] = pos;
    venc_scan_map_data[256 + write_pos] = off;
    write_pos++;
}

void gen_scan_map()
{
    // Zero the unused entries.
    venc_scan_map_data[255] = venc_scan_map_data[256+255] = 0;

    // Iterate on the block sizes in reverse order for alignment purposes.
    for (int lbs = 3; lbs >= 0; lbs--)
    {
        int bs = 1<<lbs;
        int bs2 = bs*bs;

        // Set the positions.
        for (int i = 0; i < 3; i++) venc_scan_map_idx[lbs][i] = write_pos + i*bs2;

        // Up-right. We do diagonals starting from bottom-right and discard the
        // values outside the block.
        //          / <= dia=3
        //         // <= dia=2
        //        /// <= dia=1
        for (int dia = 1; dia < 2*bs; dia++)
        {
            // Pass the values on the diagonal.
            for (int i = 0; i < dia; i++)
            {
                int y = bs - dia + i;
                int x = bs - 1 - i;
                if (x < 0 || y < 0) continue;
                int pos = y*bs + x;
                update_map(lbs, pos);
            }
        }

        // Horizontal.
        for (int i = bs2 - 1; i >= 0; i--) update_map(lbs, i);

        // Vertical.
        for (int x = bs - 1; x >= 0; x--)
        {
            for (int y = bs - 1; y >= 0; y--)
            {
                int pos = y*bs + x;
                update_map(lbs, pos);
            }
        }
    }

    printf("const uint8_t venc_scan_map_idx[4][3] =\n{\n");
    for (int i = 0; i < 4; i++)
    {
        printf("    { 0x%02x, 0x%02x, 0x%02x },", venc_scan_map_idx[i][0], venc_scan_map_idx[i][1], venc_scan_map_idx[i][2]);
        printf("\n");
    }
    printf("};\n");

    printf("const uint8_t venc_scan_map_data[2*256] =\n{\n");
    for (int i = 0; i < 2*256; i++)
    {
        if (i % 16 == 0)
        {
            if (i) printf("\n");
            printf("    ");
        }
        if (i % 16 != 0) printf(" ");
        printf("0x%02x,", venc_scan_map_data[i]);
    }
    printf("\n};\n");
}

void gen_last_coeff_table()
{
    int t[5*16];
    for (int sb_ctx = 0; sb_ctx < 4; sb_ctx++)
        for (int y = 0; y < 4; y++)
            for (int x = 0; x < 4; x++)
            {
                int pos = y*4 + x;
                int val = 0;
                if (sb_ctx == 0) val = (x+y == 0) ? 2 : (x+y < 3) ? 1 : 0;
                else if (sb_ctx == 1) val = (y == 0) ? 2 : (y == 1) ? 1 : 0;
                else if (sb_ctx == 2) val = (x == 0) ? 2 : (x == 1) ? 1 : 0;
                else if (sb_ctx == 3) val = 2;
                t[sb_ctx*16 + pos] = val;
            }
    int tb_4x4[16] = { 0, 1, 4, 5, 2, 3, 4, 5, 6, 6, 8, 8, 7, 7, 8, 0 };
    memcpy(t + 4*16, tb_4x4, 4*16);

    printf("const uint8_t coeff_nz_ctx_table[5*16] =\n{\n");
    for (int i = 0; i < 5*16; i++)
    {
        if (i % 16 == 0)
        {
            if (i) printf("\n");
            printf("    ");
        }
        if (i % 16 != 0) printf(" ");
        printf("%d,", t[i]);
    }
    printf("\n};\n");
}

int main()
{
    gen_scan_map();
    gen_last_coeff_table();
    return 0;
}

