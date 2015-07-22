; Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
; for the full license text.

%include "x86inc.asm"

section .data

align 4
pat_w_32:           dw 32,32                ; We duplicate the factors to broadcast them faster.
pat_dw_2048:        dd 2048

align 16
pat_if_luma_8bit:   db -1,4,-1,4, -10,58,-10,58, 17,-5,17,-5, 1,0,1,0
                    db -1,4,-1,4, -11,40,-11,40, 40,-11,40,-11, 4,-1,4,-1
                    db 0,1,0,1, -5,17,-5,17, 58,-10,58,-10, 4,-1,4,-1
pat_if_luma_16bit:  dw -1,4, -10,58, 17,-5, 1,0
                    dw -1,4, -11,40, 40,-11, 4,-1
                    dw 0,1, -5,17, 58,-10, 4,-1
pat_avg_pix_12:     dd 5, 6, 0, 0, 0, 1, 2, 4

section .text

; int sad(f265_pix *src, int src_stride, f265_pix *ref, int ref_stride, int packed_dims)
; Input parameters:
; - g0:     source.
; - g1:     source stride.
; - g2:     reference.
; - g3:     reference stride.
; - g4:     packed_dims.
%macro DEFINE_SAD 1                         ; %1: block width.

; Declare the function. We do just enough work per loop iteration to amortize the
; loop overhead (either two or four rows per iteration).
%if %1 == 32
DEFFUN f265_lbd_fsad_%1_avx2, ia=5, at=84844, ti=2, tv=3, ym=1
%else
DEFFUN f265_lbd_fsad_%1_avx2, ia=5, at=84844, ti=0, tv=3, ym=1
%endif
    ; Initialization.
    vpxor           y0, y0, y0              ; Initialize the SAD accumulator.
    and             g4, 0x7f                ; Get the height.
    %if %1 == 32
    lea             g5, [3*g1]              ; Load 3*stride.
    lea             g6, [3*g3]
    %endif

    ; Loop body before pointer update.
    .loop:
    %if %1 == 64
    vmovdqu         y1, [g0]                ; Load source.
    vpsadbw         y1, y1, [g2]            ; SAD with reference.
    vpaddd          y0, y0, y1              ; Accumulate.

    vmovdqu         y1, [g0 + 32]
    vpsadbw         y1, y1, [g2 + 32]
    vpaddd          y0, y0, y1

    vmovdqu         y1, [g0 + g1]
    vpsadbw         y1, y1, [g2 + g3]
    vpaddd          y0, y0, y1

    vmovdqu         y1, [g0 + g1 + 32]
    vpsadbw         y1, y1, [g2 + g3 + 32]
    vpaddd          y0, y0, y1

    %elif %1 == 32
    vmovdqu         y1, [g0]
    vpsadbw         y1, y1, [g2]
    vpaddd          y0, y0, y1

    vmovdqu         y1, [g0 + g1]
    vpsadbw         y1, y1, [g2 + g3]
    vpaddd          y0, y0, y1

    vmovdqu         y1, [g0 + 2*g1]
    vpsadbw         y1, y1, [g2 + 2*g3]
    vpaddd          y0, y0, y1

    vmovdqu         y1, [g0 + g5]
    vpsadbw         y1, y1, [g2 + g6]
    vpaddd          y0, y0, y1

    %elif %1 == 16
    vmovdqu         y1, [g0]
    vinserti128     y1, y1, [g0 + g1], 1
    vmovdqu         y2, [g2]
    vinserti128     y2, y2, [g2 + g3], 1
    vpsadbw         y1, y1, y2
    vpaddd          y0, y0, y1

    %elif %1 == 8
    vmovq           x1, [g0]
    vmovhps         x1, [g0 + g1]
    vmovq           x2, [g2]
    vmovhps         x2, [g2 + g3]
    vpsadbw         y1, y1, y2
    vpaddd          y0, y0, y1

    %elif %1 == 4
    vmovd           x1, [g0]
    vpunpckldq      x1, x1, [g0 + g1]
    vmovd           x2, [g2]
    vpunpckldq      x2, x2, [g2 + g3]
    vpsadbw         y1, y1, y2
    vpaddd          y0, y0, y1
    %endif

    ; Pointer update.
    %if %1 == 32
    lea             g0, [g0 + 4*g1]
    lea             g2, [g2 + 4*g3]
    sub             g4, 4
    %else
    lea             g0, [g0 + 2*g1]
    lea             g2, [g2 + 2*g3]
    sub             g4, 2
    %endif
    jnz             .loop

    ; Final combination.
    %if %1 >= 16
    vextracti128    x1, y0, 1
    vpaddd          y0, y0, y1
    %endif
    %if %1 >= 8
    vpshufd         y1, y0, 2
    vpaddd          y0, y0, y1
    %endif
    vmovd           gad, x0
    RET
%endmacro
DEFINE_SAD 4
DEFINE_SAD 8
DEFINE_SAD 16
DEFINE_SAD 32
DEFINE_SAD 64
%unmacro DEFINE_SAD 1


; void sadx(int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims).
;
; Input parameters:
; - g0:     costs.
; - g1:     source.
; - g2:     source stride.
; - g3:     references.
; - g4:     reference stride.
; - g5:     packed_dims.
;
; Register usage:
; - g0:     costs.
; - g1:     source.
; - g2:     source stride.
; - g3:     even stride accumulator (0, 2*stride, 4*stride, ...)
;           OR single stride accumulator (block sizes 32 and 64).
; - g4:     odd stride accumulator (stride, 3*stride, 5*stride, ...)
;           OR first reference pointer (block sizes 32 and 64).
; - ga:     2*stride OR 1*stride.
; - g5:     height.
; - g6-9:   reference pointers.
; - y0:     scratch register.
; - y1:     source pixels.
; - y2-5:   SAD accumulators.
%macro DEFINE_SAD 2                         ; %1: block width, %2: number of references.

; We process single rows for widths 32 and 64.
; We need one extra vector register for SAD4.
%assign nbti 2
%assign nbtv 5
%if %1 <= 16
%assign nbti nbti+1
%endif
%if %2 == 4
%assign nbti nbti+1
%assign nbtv nbtv+1
%endif
DEFFUN f265_lbd_sad%2_%1_avx2, ia=6, at=884844, ti=nbti, tv=nbtv, ym=1
    ; 16 and below.
    %if %1 <= 16
    mov             g6, [g3]                ; Load the reference pointers.
    mov             g7, [g3+8]
    mov             g8, [g3+16]
    %if %2 == 4
    mov             g9, [g3+24]
    %endif
    lea             ga, [2*g4]              ; Load the stride increment.

    ; 32 and above.
    %else
    mov             ga, g4                  ; Load the stride increment.
    mov             g4, [g3]
    mov             g6, [g3+8]
    mov             g7, [g3+16]
    %if %2 == 4
    mov             g8, [g3+24]
    %endif
    %endif

    xor             g3, g3                  ; Initialize the stride counter.
    and             g5, 0xff                ; Extract the height.
    vpxor           y2, y2, y2              ; Initialize the SAD accumulators.
    vpxor           y3, y3, y3
    vpxor           y4, y4, y4
    %if %2 == 4
    vpxor           y5, y5, y5
    %endif

    .loop:                                  ; Process 1 or 2 rows per loop iteration.

    ; Width 64.
    %if %1 == 64
    vmovdqu         y1, [g1]
    vpsadbw         y0, y1, [g4 + g3]
    vpaddd          y2, y2, y0
    vpsadbw         y0, y1, [g6 + g3]
    vpaddd          y3, y3, y0
    vpsadbw         y0, y1, [g7 + g3]
    vpaddd          y4, y4, y0
    %if %2 == 4
    vpsadbw         y0, y1, [g8 + g3]
    vpaddd          y5, y5, y0
    %endif

    vmovdqu         y1, [g1 + 32]
    vpsadbw         y0, y1, [g4 + g3 + 32]
    vpaddd          y2, y2, y0
    vpsadbw         y0, y1, [g6 + g3 + 32]
    vpaddd          y3, y3, y0
    vpsadbw         y0, y1, [g7 + g3 + 32]
    vpaddd          y4, y4, y0
    %if %2 == 4
    vpsadbw         y0, y1, [g8 + g3 + 32]
    vpaddd          y5, y5, y0
    %endif

    ; Width 32.
    %elif %1 == 32
    vmovdqu         y1, [g1]

    vpsadbw         y0, y1, [g4 + g3]
    vpaddd          y2, y2, y0

    vpsadbw         y0, y1, [g6 + g3]
    vpaddd          y3, y3, y0

    vpsadbw         y0, y1, [g7 + g3]
    vpaddd          y4, y4, y0

    %if %2 == 4
    vpsadbw         y0, y1, [g8 + g3]
    vpaddd          y5, y5, y0
    %endif

    ; Width 16.
    %elif %1 == 16
    vmovdqu         y1, [g1]
    vinserti128     y1, y1, [g1 + g2], 1

    vmovdqu         y0, [g6 + g3]
    vinserti128     y0, y0, [g6 + g4], 1
    vpsadbw         y0, y0, y1
    vpaddd          y2, y2, y0

    vmovdqu         y0, [g7 + g3]
    vinserti128     y0, y0, [g7 + g4], 1
    vpsadbw         y0, y0, y1
    vpaddd          y3, y3, y0

    vmovdqu         y0, [g8 + g3]
    vinserti128     y0, y0, [g8 + g4], 1
    vpsadbw         y0, y0, y1
    vpaddd          y4, y4, y0

    %if %2 == 4
    vmovdqu         y0, [g9 + g3]
    vinserti128     y0, y0, [g9 + g4], 1
    vpsadbw         y0, y0, y1
    vpaddd          y5, y5, y0
    %endif

    ; Width 8.
    %elif %1 == 8
    vmovq           x1, [g1]
    vmovhps         x1, [g1 + g2]

    vmovq           x0, [g6 + g3]
    vmovhps         x0, [g6 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y2, y2, y0

    vmovq           x0, [g7 + g3]
    vmovhps         x0, [g7 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y3, y3, y0

    vmovq           x0, [g8 + g3]
    vmovhps         x0, [g8 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y4, y4, y0

    %if %2 == 4
    vmovq           x0, [g9 + g3]
    vmovhps         x0, [g9 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y5, y5, y0
    %endif

    ; Width 4.
    %elif %1 == 4

    vmovd           x1, [g1]
    vpunpckldq      x1, x1, [g1 + g2]

    vmovd           x0, [g6 + g3]
    vpunpckldq      x0, x0, [g6 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y2, y2, y0

    vmovd           x0, [g7 + g3]
    vpunpckldq      x0, x0, [g7 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y3, y3, y0

    vmovd           x0, [g8 + g3]
    vpunpckldq      x0, x0, [g8 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y4, y4, y0

    %if %2 == 4
    vmovd           x0, [g9 + g3]
    vpunpckldq      x0, x0, [g9 + g4]
    vpsadbw         y0, y0, y1
    vpaddd          y5, y5, y0
    %endif

    %endif

    %if %1 <= 16
    lea             g1, [g1 + 2*g2]         ; Pointer update.
    add             g3, ga
    add             g4, ga
    sub             g5, 2
    %else
    add             g1, g2
    add             g3, ga
    sub             g5, 1
    %endif
    jnz             .loop

    ; Width 8 and above.
    %if %1 >= 8
    vpsllq          y3, y3, 32              ; Interleave costs 0 and 1.
    vpor            y2, y2, y3
    %if %2 == 4
    vpsllq          y5, y5, 32              ; Interleave costs 2 and 3.
    vpor            y4, y4, y5
    %endif

    %if %1 >= 16
    vperm2i128      y1, y2, y4, 0x20        ; Combine lanes.
    vperm2i128      y2, y2, y4, 0x31
    vpaddd          y2, y1, y2
    %else
    vinserti128     y2, y2, x4, 1
    %endif
    vpshufd         y1, y2, 14              ; Combine low and high (xxxx 1110 == 14).
    vpaddd          y2, y1, y2
    vpermq          y2, y2, 8               ; Put quadwords 0 and 2 together (xxxx 1000 == 8).

    ; Width 4.
    %else
    vpunpckldq      y2, y2, y3
    %if %2 == 4
    vpunpckldq      y4, y4, y5
    %endif
    vpunpcklqdq     y2, y2, y4
    %endif

    vmovdqu         [g0], x2                ; Store.
    RET
%endmacro
DEFINE_SAD 4, 3
DEFINE_SAD 4, 4
DEFINE_SAD 8, 3
DEFINE_SAD 8, 4
DEFINE_SAD 16, 3
DEFINE_SAD 16, 4
DEFINE_SAD 32, 3
DEFINE_SAD 32, 4
DEFINE_SAD 64, 3
DEFINE_SAD 64, 4
%unmacro DEFINE_SAD 2


; void venc_avg_pix(f265_pix *dst, f265_pix *src0, int src0_stride, f265_pix *src1, int src1_stride, int packed_dims)
; Input parameters:
; - g0:     destination.
; - g1:     source 0.
; - g2:     source 0 stride.
; - g3:     source 1.
; - g4:     source 1 stride.
; - g5:     packed_dims.
%macro DEFINE_AVG_PIX 1                     ; %1: width.

%assign nbti 0
%if %1 <= 16 || %1 == 32
%assign nbti 2
%endif

DEFFUN f265_lbd_avg_pix_%1_avx2, ia=6, at=884844, ti=nbti, tv=3, ym=1
    ; Initialize.
    and             g5, 0x7f                ; Height.
    %if %1 <= 16 || %1 == 32
    lea             g6, [3*g2]              ; 3*stride.
    lea             g7, [3*g4]
    %endif
    %if %1 == 12
    vmovdqu          y2, [pat_avg_pix_12]
    %endif

    ; Compute the average.
    .loop:

    ; Width 4.
    %if %1 == 4
    vmovdqu         y0, [g1]
    vpunpckldq      y0, y0, [g1+g2]
    vmovdqu         y1, [g3]
    vpunpckldq      y1, y1, [g3+g4]
    vpavgb          y0, y0, y1
    vmovq           [g0], x0

    vmovdqu         y0, [g1+2*g2]
    vpunpckldq      y0, y0, [g1+g6]
    vmovdqu         y1, [g3+2*g4]
    vpunpckldq      y1, y1, [g3+g7]
    vpavgb          y0, y0, y1
    vmovq           [g0+8], x0

    ; Width 8.
    %elif %1 == 8
    vmovdqu         x0, [g1]
    vmovhps         x0, [g1+g2]
    vmovdqu         x1, [g3]
    vmovhps         x1, [g3+g4]
    vpavgb          y0, y0, y1
    vmovdqu         [g0], x0

    vmovdqu         x0, [g1+2*g2]
    vmovhps         x0, [g1+g6]
    vmovdqu         x1, [g3+2*g4]
    vmovhps         x1, [g3+g7]
    vpavgb          y0, y0, y1
    vmovdqu         [g0+16], x0

    ; Width 16.
    %elif %1 == 16
    vmovdqu         y0, [g1]
    vinserti128     y0, y0, [g1+g2], 1
    vmovdqu         y1, [g3]
    vinserti128     y1, y1, [g3+g4], 1
    vpavgb          y0, y0, y1
    vmovdqu         [g0], y0

    vmovdqu         y0, [g1+2*g2]
    vinserti128     y0, y0, [g1+g6], 1
    vmovdqu         y1, [g3+2*g4]
    vinserti128     y1, y1, [g3+g7], 1
    vpavgb          y0, y0, y1
    vmovdqu         [g0+32], y0

    ; Width 32.
    %elif %1 == 32
    vmovdqu         y0, [g1]
    vpavgb          y0, y0, [g3]
    vmovdqu         [g0], y0

    vmovdqu         y0, [g1+g2]
    vpavgb          y0, y0, [g3+g4]
    vmovdqu         [g0+32], y0

    vmovdqu         y0, [g1+2*g2]
    vpavgb          y0, y0, [g3+2*g4]
    vmovdqu         [g0+64], y0

    vmovdqu         y0, [g1+g6]
    vpavgb          y0, y0, [g3+g7]
    vmovdqu         [g0+96], y0

    ; Width 64.
    %elif %1 == 64
    vmovdqu         y0, [g1]
    vpavgb          y0, y0, [g3]
    vmovdqu         [g0], y0

    vmovdqu         y0, [g1+32]
    vpavgb          y0, y0, [g3+32]
    vmovdqu         [g0+32], y0

    vmovdqu         y0, [g1+g2]
    vpavgb          y0, y0, [g3+g4]
    vmovdqu         [g0+64], y0

    vmovdqu         y0, [g1+g2+32]
    vpavgb          y0, y0, [g3+g4+32]
    vmovdqu         [g0+96], y0

    ; Width 12.
    %elif %1 == 12
    vmovdqu         y0, [g1]
    vinserti128     y0, y0, [g1+g2], 1
    vmovdqu         y1, [g3]
    vinserti128     y1, y1, [g3+g4], 1
    vpavgb          y0, y0, y1
    vpermd          y0, y2, y0              ; Pack 16 bytes high, 8 bytes low.
    vextracti128    [g0], y0, 1
    vmovq           [g0+16], x0

    vmovdqu         y0, [g1+2*g2]
    vinserti128     y0, y0, [g1+g6], 1
    vmovdqu         y1, [g3+2*g4]
    vinserti128     y1, y1, [g3+g7], 1
    vpavgb          y0, y0, y1
    vpermd          y0, y2, y0
    vextracti128    [g0+24], y0, 1
    vmovq           [g0+40], x0

    ; Width 24.
    %elif %1 == 24
    vmovdqu         y0, [g1]
    vpavgb          y0, y0, [g3]
    vmovdqu         [g0], x0
    vextracti128    x0, y0, 1
    vmovq           [g0+16], x0

    vmovdqu         y0, [g1+g2]
    vpavgb          y0, y0, [g3+g4]
    vmovdqu         [g0+24], x0
    vextracti128    x0, y0, 1
    vmovq           [g0+40], x0

    ; Width 48.
    %elif %1 == 48
    vmovdqu         y0, [g1]
    vpavgb          y0, y0, [g3]
    vmovdqu         [g0+0*16], y0

    vmovdqu         y0, [g1+32]
    vpavgb          y0, y0, [g3+32]
    vmovdqu         [g0+2*16], x0

    vmovdqu         y0, [g1+g2]
    vpavgb          y0, y0, [g3+g4]
    vmovdqu         [g0+3*16], y0

    vmovdqu         y0, [g1+g2+32]
    vpavgb          y0, y0, [g3+g4+32]
    vmovdqu         [g0+5*16], x0

    %endif

    ; Update the pointer and loop.
    %if %1 <= 16 || %1 == 32
    add             g0, %1*4
    lea             g1, [g1+4*g2]
    lea             g3, [g3+4*g4]
    sub             g5, 4
    %else
    add             g0, %1*2
    lea             g1, [g1+2*g2]
    lea             g3, [g3+2*g4]
    sub             g5, 2
    %endif
    jnz             .loop
    RET
%endmacro
DEFINE_AVG_PIX 4
DEFINE_AVG_PIX 8
DEFINE_AVG_PIX 16
DEFINE_AVG_PIX 32
DEFINE_AVG_PIX 64
DEFINE_AVG_PIX 12
DEFINE_AVG_PIX 24
DEFINE_AVG_PIX 48
%unmacro DEFINE_AVG_PIX 1


; void venc_interpol_luma_qpel_pix(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,
;                                  int packed_dims, uint8_t *spill)
; Input parameters:
; - g0:     destination.
; - g1:     destination stride.
; - g2:     source.
; - g3:     source stride.
; - g4:     interpolation fraction.
; - g5:     packed dimensions.
; - g6:     spill.
;
; We want an interpolation function for 3 block sizes (luma 16,8,4, chroma 8,4,2), fractions h,v,d,
; in clipping and unclipping versions. We probably want to use macros to generate those 36 functions.
DEFFUN f265_lbd_interpol_luma_qpel_pix_8_h_avx2, ia=7, at=8484448, ti=0, tv=9, ym=1
    ; Load the factors.
    and             g4, 3                   ; Get the factor index.
    shl             g4, 4
    lea             ga, [pat_if_luma_8bit - 16] ; Load the factor base.
    lea             ga, [ga + g4]           ; Add the factor index.
    vpbroadcastd    y5, [ga + 0]            ; Broadcast the factors.
    vpbroadcastd    y6, [ga + 4]
    vpbroadcastd    y7, [ga + 8]
    vpbroadcastd    y8, [ga + 12]
    vpbroadcastd    y9, [pat_w_32]          ; Load the bias.

    ; Initialize.
    and             g5, 0x7f                ; Get the height.
    sub             g2, 3                   ; Offset the source at -3.
    lea             ga, [3*g1]              ; 3 * destination stride.
    lea             g4, [3*g3]              ; 3 * source stride.

    ; Align the factors, multiply, and add.
    %macro ADD_FACTORS 3                    ; %1: accumulator, %2: alignment, %3: factor.
    vpalignr        y0, y4, y3, %2
    vpmaddubsw      y0, y0, %3
    vpaddw          %1, %1, y0
    %endmacro

    ; Load a row in each lane and process them.
    %macro PROCESS_PAIR 3                   ; %1: accumulator, %2: first source, %3: second source.
    vmovdqu         y0, [%2]                ; Load the rows at offset -3.
    vinserti128     y0, y0, [%3], 1
    vpalignr        %1, y0, y0, 1           ; Align each pixel with the next.
    vpunpcklbw      y3, y0, %1              ; Interleave 0,1, 1,2, ... low and high.
    vpunpckhbw      y4, y0, %1
    vpmaddubsw      %1, y3, y5              ; Factors -3,-2.
    ADD_FACTORS     %1, 4, y6               ; Factors -1,0.
    ADD_FACTORS     %1, 8, y7               ; Factors 1,2.
    ADD_FACTORS     %1, 12, y8              ; Factors 3,4.
    vpaddw          %1, %1, y9              ; Add bias.
    vpsraw          %1, %1, 6               ; Shift.
    %endmacro

    ; Process four rows at a time.
    .loop:
    PROCESS_PAIR    y1, g2, g2+2*g3         ; Rows 0 and 2.
    PROCESS_PAIR    y2, g2+g3, g2+g4        ; Rows 1 and 3.
    vpackuswb       y0, y1, y2              ; Convert to 8-bit and pack the rows together.
    vextracti128    x1, y0, 1               ; Store the four rows.
    vmovq           [g0], x0
    vmovhps         [g0 + g1], x0
    vmovq           [g0 + 2*g1], x1
    vmovhps         [g0 + ga], x1
    lea             g0, [g0 + 4*g1]         ; Loop.
    lea             g2, [g2 + 4*g3]
    sub             g5, 4
    jnz             .loop
    RET
    %unmacro ADD_FACTORS 3
    %unmacro PROCESS_PAIR 3


DEFFUN f265_lbd_interpol_luma_qpel_pix_8_v_avx2, ia=7, at=8484448, ti=0, tv=15, ym=1
    ; We process 4 rows per iteration. Two of these rows are from the bottom half
    ; of the block. This reduces cross-lane moves and the pressure on port 5.

    ; Load one row and optionally interleave with the previous row.
    %macro LOAD_ROW 5                       ; %1: row dst, %2-3: src, %4: interleave dst (0 if none), %5: prev row.
    vmovdqu         %1, [%2]                ; Load the two sources for the destination row.
    vinserti128     %1, %1, [%3], 1
    %if %4 != 0                             ; Interleave with the previous row.
    vpunpcklbw      %4, %5, %1
    %endif
    %endmacro

    ; Process an interleaved pair of rows.
    %macro PROCESS_PAIR 4                   ; %1: work, %2: accumulator, %3: source, %4: factor.
    vpmaddubsw      %1, %3, %4              ; Multiply.
    %if %1 != %2
    vpaddw          %2, %2, %1              ; Add.
    %endif
    %endmacro

    ; Load the factors.
    shr             g4, 2                   ; Get the factor index.
    shl             g4, 4
    lea             ga, [pat_if_luma_8bit - 16] ; Load the factor base.
    lea             ga, [ga + g4]           ; Add the factor index.
    vpbroadcastd    y3, [ga + 0]            ; Broadcast the factors.
    vpbroadcastd    y4, [ga + 4]
    vpbroadcastd    y5, [ga + 8]
    vpbroadcastd    y6, [ga + 12]
    vpbroadcastd    y13, [pat_w_32]         ; Load the bias.

    ; Initialize.
    and             g5, 0x7f                ; Get the height.
    shr             g5, 1                   ; Height/2.
    lea             g6, [3*g3]              ; 3 * source stride.
    sub             g2, g6                  ; Offset the source at -3*stride.
    mov             ga, g5                  ; Source -3*stride + (height/2)*stride.
    imul            ga, g3
    add             ga, g2
    mov             g4, g5                  ; Destination + (height/2)*stride.
    imul            g4, g1
    add             g4, g0

    ; Interleave rows (-3,-2), (-2,-1), (-1,0), (0,1), (1,2), (2,3).
    LOAD_ROW        y0, g2, ga, 0, 0        ; Load the first three pairs.
    LOAD_ROW        y1, g2+g3, ga+g3, y7, y0
    LOAD_ROW        y0, g2+2*g3, ga+2*g3, y8, y1
    LOAD_ROW        y1, g2+g6, ga+g6, y9, y0
    lea             g2, [g2 + 4*g3]         ; Offset the sources.
    lea             ga, [ga + 4*g3]
    LOAD_ROW        y0, g2, ga, y10, y1     ; Load the last three pairs.
    LOAD_ROW        y1, g2+g3, ga+g3, y11, y0
    LOAD_ROW        y14, g2+2*g3, ga+2*g3, y12, y1 ; Keep the data of the last row.
    add             g2, g6                  ; Offset the sources.
    add             ga, g6

    ; Process four rows at a time.
    .loop:

    ; Process the 6 pairs currently available.
    PROCESS_PAIR    y1, y1, y7, y3
    PROCESS_PAIR    y2, y2, y8, y3
    PROCESS_PAIR    y0, y1, y9, y4
    PROCESS_PAIR    y0, y2, y10, y4
    PROCESS_PAIR    y0, y1, y11, y5
    PROCESS_PAIR    y0, y2, y12, y5

    ; Shift the rows down.
    vmovdqu         y7, y9
    vmovdqu         y8, y10
    vmovdqu         y9, y11
    vmovdqu         y10, y12

    ; Load rows (3,4) and (4,5) and offset the sources.
    LOAD_ROW        y0, g2, ga, y11, y14
    LOAD_ROW        y14, g2+g3, ga+g3, y12, y0
    lea             g2, [g2 + 2*g3]
    lea             ga, [ga + 2*g3]

    ; Process the last 2 pairs.
    PROCESS_PAIR    y0, y1, y11, y6
    PROCESS_PAIR    y0, y2, y12, y6

    ; Add the bias, shift and pack.
    vpaddw          y1, y13, y1
    vpaddw          y2, y13, y2
    vpsraw          y1, y1, 6
    vpsraw          y2, y2, 6
    vpackuswb       y0, y1, y2

    ; Store the four rows and offset the destination.
    vextracti128    x1, y0, 1
    vmovq           [g0], x0
    vmovhps         [g0 + g1], x0
    vmovq           [g4], x1
    vmovhps         [g4 + g1], x1
    lea             g0, [g0 + 2*g1]
    lea             g4, [g4 + 2*g1]

    ; Loop.
    sub             g5, 2
    jnz             .loop
    RET
    %unmacro LOAD_ROW 5
    %unmacro PROCESS_PAIR 4


DEFFUN f265_lbd_interpol_luma_qpel_pix_8_d_avx2, ia=7, at=8484448, ti=0, tv=16, ym=1

    ; The first pass is similar to the horizontal case above, except that we
    ; process two rows at a time.

    ; Load the factors.
    vmovq           x15, g4                 ; Save the fraction (ymm15).
    and             g4, 3                   ; Get the factor index.
    shl             g4, 4
    lea             ga, [pat_if_luma_8bit - 16] ; Load the factor base.
    lea             ga, [ga + g4]           ; Add the factor index.
    vpbroadcastd    y4, [ga + 0]            ; Broadcast the factors.
    vpbroadcastd    y5, [ga + 4]
    vpbroadcastd    y6, [ga + 8]
    vpbroadcastd    y7, [ga + 12]

    ; Initialize.
    and             g5, 0x7f                ; Get the height.
    mov             g4, g5                  ; Save the height.
    add             g5, 8                   ; Process height+8 rows.
    lea             ga, [3*g3 + 3]          ; Offset the source at -3 - 3*stride.
    sub             g2, ga
    mov             ga, g6                  ; Spill destination.

    ; Align the factors, multiply, and add.
    %macro ADD_FACTORS 3                    ; %1: accumulator, %2: alignment, %3: factor.
    vpalignr        y0, y3, y2, %2
    vpmaddubsw      y0, y0, %3
    vpaddw          %1, %1, y0
    %endmacro

    ; Load a row in each lane and process them.
    %macro PROCESS_PAIR 3                   ; %1: accumulator, %2: first source, %3: second source.
    vmovdqu         y0, [%2]                ; Load the rows at offset -3.
    vinserti128     y0, y0, [%3], 1
    vpalignr        %1, y0, y0, 1           ; Align each pixel with the next.
    vpunpcklbw      y2, y0, %1              ; Interleave 0,1, 1,2, ... low and high.
    vpunpckhbw      y3, y0, %1
    vpmaddubsw      %1, y2, y4              ; Factors -3,-2.
    ADD_FACTORS     %1, 4, y5               ; Factors -1,0.
    ADD_FACTORS     %1, 8, y6               ; Factors 1,2.
    ADD_FACTORS     %1, 12, y7              ; Factors 3,4.
    %endmacro

    ; Process two rows at a time.
    .loop:
    PROCESS_PAIR    y1, g2, g2+g3           ; Rows 0 and 1.
    vmovdqu         [ga], y1                ; Store.
    add             ga, 32                  ; Loop.
    lea             g2, [g2 + 2*g3]
    sub             g5, 2
    jnz             .loop
    %unmacro ADD_FACTORS 3
    %unmacro PROCESS_PAIR 3

    ; The second pass is different from the pure vertical case. We process two
    ; rows per loop interation, in order. The even/odd rows at processed in the
    ; low/high lane. We pack as shown below, then we interleave rows for pmadd.
    ;   1 | 0
    ;   2 | 1
    ;   3 | 2
    ;   4 | 3

    ; Load one row and optionally interleave with the previous row.
    %macro LOAD_ROW 3                       ; %1: row dst, %2: src, %3: prev row (0 if none).
    vbroadcasti128  %1, [%2]                ; Broadcast the row.
    %if %3 != 0                             ; Interleave with the previous row.
    vpblendd        %3, %3, %1, 0xf0        ; Put the current row in the high lane of the previous row.
    %endif
    %endmacro

    ; Process an interleaved pair of rows.
    %macro PROCESS_PAIR 4                   ; %1-2: sources, %3: factor, %4: true to store directly in acc.
    %if %4
    vpunpcklwd      y0, %1, %2              ; Interleave.
    vpmaddwd        y1, y0, %3              ; Multiply.
    vpunpckhwd      y0, %1, %2
    vpmaddwd        y2, y0, %3
    %else
    vpunpcklwd      y0, %1, %2
    vpmaddwd        y0, y0, %3
    vpaddd          y1, y1, y0
    vpunpckhwd      y0, %1, %2
    vpmaddwd        y0, y0, %3
    vpaddd          y2, y2, y0
    %endif
    %endmacro

    ; Load the factors.
    movq            g5, x15                 ; Get the factor index.
    shr             g5, 2
    shl             g5, 4
    lea             ga, [pat_if_luma_16bit - 16] ; Load the factor base.
    lea             ga, [ga + g5]           ; Add the factor index.
    vpbroadcastd    y3, [ga + 0]            ; Broadcast the factors.
    vpbroadcastd    y4, [ga + 4]
    vpbroadcastd    y5, [ga + 8]
    vpbroadcastd    y6, [ga + 12]
    vpbroadcastd    y13, [pat_dw_2048]      ; Load the bias.

    ; Initialize.
    mov             ga, g6                  ; Set the source.

    ; Interleave rows (-3,-2), (-2,-1), (-1,0), (0,1), (1,2), (2,3).
    LOAD_ROW        y7,  ga+0*16, 0
    LOAD_ROW        y8,  ga+1*16, y7
    LOAD_ROW        y9,  ga+2*16, y8
    LOAD_ROW        y10, ga+3*16, y9
    LOAD_ROW        y11, ga+4*16, y10
    LOAD_ROW        y12, ga+5*16, y11
    LOAD_ROW        y14, ga+6*16, y12
    add             ga, 7*16

    ; Process four rows at a time.
    .loop2:

    ; Process the 6 pairs currently available.
    PROCESS_PAIR    y7,  y8,  y3, 1
    PROCESS_PAIR    y9,  y10, y4, 0
    PROCESS_PAIR    y11, y12, y5, 0

    ; Shift the rows down.
    vmovdqu         y7, y9
    vmovdqu         y8, y10
    vmovdqu         y9, y11
    vmovdqu         y10, y12
    vmovdqu         y11, y14

    ; Load rows (3,4) and (4,5) and offset the source.
    LOAD_ROW        y12, ga, y11
    LOAD_ROW        y14, ga+16, y12
    add             ga, 32

    ; Process the last 2 pairs.
    PROCESS_PAIR    y11, y12, y6, 0

    ; Add the bias, shift and pack twice.
    vpaddd          y1, y13, y1
    vpaddd          y2, y13, y2
    vpsrad          y1, y1, 12
    vpsrad          y2, y2, 12
    vpackssdw       y0, y1, y2
    vpackuswb       y0, y0, y0

    ; Store the two rows and offset the destination.
    vextracti128    x1, y0, 1
    vmovq           [g0], x0
    vmovq           [g0 + g1], x1
    lea             g0, [g0 + 2*g1]

    ; Loop.
    sub             g4, 2
    jnz             .loop2
    RET
    %unmacro LOAD_ROW 3
    %unmacro PROCESS_PAIR 4

