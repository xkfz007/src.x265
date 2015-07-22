; Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
; for the full license text.

%include "x86inc.asm"

section .data
align 4
pat_quant_dw_1:     dd 1
pat_dequant_dw_1:   dw 1,1

section .text
; ---------------------- QUANT/DEQUANT macros ---------------------
%macro MULTIPLY 3                           ; %1: tmp register, %2: src 0, %3: src 1.
    ; For QUANT: Interleaving gives - input, Add >> (Shift - 12) | input <n>, Add >> (Shift - 12).
    ; For DEQUANT: Interleaving gives - 1, input | 1, input.
    vpunpckhwd      y5, %2, %3
    vpunpcklwd      %1, %2, %3
    ; For QUANT: Multiplication gives- (input*mult) + ((Add >> (Shift - 12))*(1 << (Shift - 12))).
    ; For DEQUANT: Multiplication gives- (add + input*mult).
    vpmaddwd        y5, y5, y1              ; Multiply.
    vpmaddwd        %1, %1, y1
    vpsrad          y5, y5, x3              ; Shift.
    vpsrad          %1, %1, x3
    vpackssdw       %1, %1, y5              ; Pack to 16-bits.
%endmacro

; QUANT 4x4.
;
; Input parameters:
; - g0:     destination.
; - g1:     source.
; - g2:     bs.
; - g3:     multiplication value.
; - g4:     addition value.
; - g5:     shift value.
DEFFUN f265_lbd_quant_4_avx2, ia=6, at=884444, ti=0, tv=6, ym=1

    vmovdqu         y0, [g1]                ; Input.
    vmovd           x1, g3d                 ; Multiplication Factor.
    vmovd           x2, g4d                 ; Addition Factor.
    vmovd           x3, g5d                 ; Shift Factor.

    vpabsw          y4, y0                  ; Absolute value of input.
    vpbroadcastd    y1, x1
    vpbroadcastd    y2, x2
    ; We need to do the following multiplication: mult*input. However, since we use "vpmaddwd", the instruction does
    ; mult*input + x*y. Either x or y should be 0, and the other (y or x respectively) can be any value.
    ; In the instruction "vpbroadcastd  y1, x1", alternate 16-bit words of y1 are 0. Hence, we require that alternate
    ; 16-bit words of the input are redundant (to do mult*input + 0*y), and hence, we interleave with any register.
    vpunpckhwd      y5, y4, y0              ; Interleave to have the required input values on every alternate word.
    vpunpcklwd      y4, y4, y0
    vpmaddwd        y4, y4, y1              ; Multiply. The multiplication factor fits in 16 bits.
    vpmaddwd        y5, y5, y1
    vpaddd          y4, y4, y2              ; Add bias.
    vpaddd          y5, y5, y2
    vpsrad          y4, y4, x3              ; Shift.
    vpsrad          y5, y5, x3
    vpackssdw       y4, y4, y5              ; Pack to 16-bits.
    vpsignw         y5, y4, y0              ; Restore the sign.
    vmovdqu         [g0], y5                ; Store.
    ; Now, compute the non-zero flag: return 1 if there are coefficients, and 0 if no coefficients.
    mov             ga, 0                   ; Initialize return value to 0.
    vpcmpeqd        y1, y1                  ; Load '1' in all bits.
    vptest          y5, y1
    setnz           gab                     ; Coefficients present.
    RET

; DEQUANT 4x4.
;
; Input parameters:
; - g0:     destination.
; - g1:     source.
; - g2:     bs.
; - g3:     multiplication value.
; - g4:     addition value.
; - g5:     shift value.
DEFFUN f265_lbd_dequant_4_avx2, ia=6, at=884444, ti=0, tv=6, ym=1

    vmovd           x1, g3d                 ; Multiplication Factor.
    vmovd           x2, g4d                 ; Addition Factor.
    vmovdqu         y0, [g1]                ; Input.
    vpbroadcastd    y4, [pat_dequant_dw_1]

    vpbroadcastw    y1, x1                  ; Multiplier and bias both fit in 16 bits.
    vpbroadcastw    y2, x2
    vmovd           x3, g5d                 ; Shift Factor.
    vpblendw        y1, y1, y2, 0xAA        ; Interleave: add, mult | add, mult.
    MULTIPLY        y4, y0, y4              ; Interleave:  1, input | 1, input.  Multiply: (add + input*mult).
    vmovdqu         [g0], y4                ; Store.
    RET

; QUANT 32x32, 16x16 and 8x8.
;
; Input parameters:
; - g0:     destination.
; - g1:     source.
; - g2:     bs.
; - g3:     multiplication value.
; - g4:     addition value.
; - g5:     shift value.
DECLARE_FUNCTION f265_lbd_quant_32_avx2
DECLARE_FUNCTION f265_lbd_quant_16_avx2
DEFFUN f265_lbd_quant_8_avx2, ia=6, at=884444, ti=0, tv=7, ym=1

    vmovd           x3, g5d                 ; Shift Factor.
    vmovd           x1, g3d                 ; Multiplication Factor.
    sub             g5, 12                  ; Shift - 12.
    vmovd           x2, g4d                 ; Addition Factor.
    imul            g2, g2                  ; bs*bs.
    vpbroadcastd    y4, [pat_quant_dw_1]    ; Constant 1.
    vpbroadcastw    y1, x1
    shl             g2, 1
    vmovd           x5, g5d
    vpxor           y6, y6
    vpslld          y4, y4, x5              ; 1 << (Shift - 12).
    vpsrad          y2, y2, x5              ; Add >> (Shift - 12).
    xor             ga, ga                  ; Initialize return value to 0.
    xor             g3, g3
    vpblendw        y1, y1, y4, 0x55        ; Mult, 1 << (Shift - 12) | Mult, 1 << (Shift - 12).
    vpbroadcastw    y2, x2

    .loop_quant:
    vmovdqu         y0, [g1 + g3]           ; Input:  'n+1' | 'n'.
    vpabsw          y4, y0                  ; Absolute value of input.
    ; Interleaving:
    ; Input, Add >> (Shift - 12) | Input, Add >> (Shift - 12).
    ; Multiply: (input*mult) + (Add >> (Shift - 12)*1 << (Shift - 12)).
    MULTIPLY        y4, y2, y4
    vpsignw         y5, y4, y0              ; Restore the sign.
    vmovdqu         [g0 + g3], y5           ; Store.
    add             g3, 32
    vpor            y6, y6, y5
    cmp             g3, g2
    jnz             .loop_quant
    ; Now, compute the non-zero flag: return 1 if there are coefficients, and 0 if no coefficients.
    vpcmpeqd        y4, y4                  ; Load '1' in all bits.
    vptest          y6, y4
    setnz           gab                     ; Coefficients present.
    RET

; DEQUANT 16x16 and 8x8.
;
; Input parameters:
; - g0:     destination.
; - g1:     source.
; - g2:     bs.
; - g3:     multiplication value.
; - g4:     addition value.
; - g5:     shift value.
DECLARE_FUNCTION f265_lbd_dequant_32_avx2
DECLARE_FUNCTION f265_lbd_dequant_16_avx2
DEFFUN f265_lbd_dequant_8_avx2, ia=6, at=884444, ti=0, tv=6, ym=1

    vmovd           x1, g3d                 ; Multiplication Factor.
    vpbroadcastd    y4, [pat_dequant_dw_1]
    vmovd           x2, g4d                 ; Addition Factor.
    vpbroadcastw    y1, x1                  ; Multiplier fits in 16 bits.
    imul            g2, g2                  ; bs*bs.
    vpbroadcastw    y2, x2                  ; Bias fits in 16 bits.
    vpblendw        y1, y1, y2, 0xAA        ; Interleave: add, mult | add, mult.
    vmovd           x3, g5d                 ; Shift Factor.
    shl             g2, 1
    xor             g3, g3

    .loop_dequant:
    vmovdqu         y0, [g1 + g3]           ; Input.
    MULTIPLY        y2, y0, y4              ; Interleave:  1, input | 1, input. Multiply: (add + input*mult).
    vmovdqu         [g0 + g3], y2           ; Store.
    add             g3, 32
    cmp             g3, g2
    jnz             .loop_dequant
    RET

%unmacro CALC_LOOP_ITER 1
%unmacro MULTIPLY 3

