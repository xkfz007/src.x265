;*****************************************************************************
;* Copyright (C) 2013 x265 project
;*
;* Authors: Min Chen <chenm003@163.com>
;*          Nabajit Deka <nabajit@multicorewareinc.com>
;*          Praveen Kumar Tiwari <praveen@multicorewareinc.com>
;*
;* This program is free software; you can redistribute it and/or modify
;* it under the terms of the GNU General Public License as published by
;* the Free Software Foundation; either version 2 of the License, or
;* (at your option) any later version.
;*
;* This program is distributed in the hope that it will be useful,
;* but WITHOUT ANY WARRANTY; without even the implied warranty of
;* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;* GNU General Public License for more details.
;*
;* You should have received a copy of the GNU General Public License
;* along with this program; if not, write to the Free Software
;* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
;*
;* This program is also available under a commercial proprietary license.
;* For more information, contact us at license @ x265.com.
;*****************************************************************************/

%include "x86inc.asm"
%include "x86util.asm"

SECTION_RODATA 32
tab_Tm:    db 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6
           db 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10
           db 8, 9,10,11, 9,10,11,12,10,11,12,13,11,12,13, 14

ALIGN 32
const interp4_vpp_shuf, times 2 db 0, 4, 1, 5, 2, 6, 3, 7, 8, 12, 9, 13, 10, 14, 11, 15

ALIGN 32
const interp_vert_shuf, times 2 db 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9
                        times 2 db 4, 6, 5, 7, 6, 8, 7, 9, 8, 10, 9, 11, 10, 12, 11, 13

ALIGN 32
const interp4_vpp_shuf1, dd 0, 1, 1, 2, 2, 3, 3, 4
                         dd 2, 3, 3, 4, 4, 5, 5, 6

ALIGN 32
const pb_8tap_hps_0, times 2 db 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8
                     times 2 db 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10
                     times 2 db 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12
                     times 2 db 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14

ALIGN 32
tab_Lm:    db 0, 1, 2, 3, 4,  5,  6,  7,  1, 2, 3, 4,  5,  6,  7,  8
           db 2, 3, 4, 5, 6,  7,  8,  9,  3, 4, 5, 6,  7,  8,  9,  10
           db 4, 5, 6, 7, 8,  9,  10, 11, 5, 6, 7, 8,  9,  10, 11, 12
           db 6, 7, 8, 9, 10, 11, 12, 13, 7, 8, 9, 10, 11, 12, 13, 14

tab_Vm:    db 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1
           db 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3

tab_Cm:    db 0, 2, 1, 3, 0, 2, 1, 3, 0, 2, 1, 3, 0, 2, 1, 3

tab_c_526336:   times 4 dd 8192*64+2048

pd_526336:      times 8 dd 8192*64+2048

tab_ChromaCoeff: db  0, 64,  0,  0
                 db -2, 58, 10, -2
                 db -4, 54, 16, -2
                 db -6, 46, 28, -4
                 db -4, 36, 36, -4
                 db -4, 28, 46, -6
                 db -2, 16, 54, -4
                 db -2, 10, 58, -2
ALIGN 32
tab_ChromaCoeff_V: times 8 db 0, 64
                   times 8 db 0,  0

                   times 8 db -2, 58
                   times 8 db 10, -2

                   times 8 db -4, 54
                   times 8 db 16, -2

                   times 8 db -6, 46
                   times 8 db 28, -4

                   times 8 db -4, 36
                   times 8 db 36, -4

                   times 8 db -4, 28
                   times 8 db 46, -6

                   times 8 db -2, 16
                   times 8 db 54, -4

                   times 8 db -2, 10
                   times 8 db 58, -2

tab_ChromaCoeffV: times 4 dw 0, 64
                  times 4 dw 0, 0

                  times 4 dw -2, 58
                  times 4 dw 10, -2

                  times 4 dw -4, 54
                  times 4 dw 16, -2

                  times 4 dw -6, 46 
                  times 4 dw 28, -4

                  times 4 dw -4, 36
                  times 4 dw 36, -4

                  times 4 dw -4, 28
                  times 4 dw 46, -6

                  times 4 dw -2, 16
                  times 4 dw 54, -4

                  times 4 dw -2, 10
                  times 4 dw 58, -2

ALIGN 32
pw_ChromaCoeffV:  times 8 dw 0, 64
                  times 8 dw 0, 0

                  times 8 dw -2, 58
                  times 8 dw 10, -2

                  times 8 dw -4, 54
                  times 8 dw 16, -2

                  times 8 dw -6, 46 
                  times 8 dw 28, -4

                  times 8 dw -4, 36
                  times 8 dw 36, -4

                  times 8 dw -4, 28
                  times 8 dw 46, -6

                  times 8 dw -2, 16
                  times 8 dw 54, -4

                  times 8 dw -2, 10
                  times 8 dw 58, -2

tab_LumaCoeff:   db   0, 0,  0,  64,  0,   0,  0,  0
                 db  -1, 4, -10, 58,  17, -5,  1,  0
                 db  -1, 4, -11, 40,  40, -11, 4, -1
                 db   0, 1, -5,  17,  58, -10, 4, -1

tab_LumaCoeffV: times 4 dw 0, 0
                times 4 dw 0, 64
                times 4 dw 0, 0
                times 4 dw 0, 0

                times 4 dw -1, 4
                times 4 dw -10, 58
                times 4 dw 17, -5
                times 4 dw 1, 0

                times 4 dw -1, 4
                times 4 dw -11, 40
                times 4 dw 40, -11
                times 4 dw 4, -1

                times 4 dw 0, 1
                times 4 dw -5, 17
                times 4 dw 58, -10
                times 4 dw 4, -1

ALIGN 32
pw_LumaCoeffVer: times 8 dw 0, 0
                 times 8 dw 0, 64
                 times 8 dw 0, 0
                 times 8 dw 0, 0

                 times 8 dw -1, 4
                 times 8 dw -10, 58
                 times 8 dw 17, -5
                 times 8 dw 1, 0

                 times 8 dw -1, 4
                 times 8 dw -11, 40
                 times 8 dw 40, -11
                 times 8 dw 4, -1

                 times 8 dw 0, 1
                 times 8 dw -5, 17
                 times 8 dw 58, -10
                 times 8 dw 4, -1

pb_LumaCoeffVer: times 16 db 0, 0
                 times 16 db 0, 64
                 times 16 db 0, 0
                 times 16 db 0, 0

                 times 16 db -1, 4
                 times 16 db -10, 58
                 times 16 db 17, -5
                 times 16 db 1, 0

                 times 16 db -1, 4
                 times 16 db -11, 40
                 times 16 db 40, -11
                 times 16 db 4, -1

                 times 16 db 0, 1
                 times 16 db -5, 17
                 times 16 db 58, -10
                 times 16 db 4, -1

tab_LumaCoeffVer: times 8 db 0, 0
                  times 8 db 0, 64
                  times 8 db 0, 0
                  times 8 db 0, 0

                  times 8 db -1, 4
                  times 8 db -10, 58
                  times 8 db 17, -5
                  times 8 db 1, 0

                  times 8 db -1, 4
                  times 8 db -11, 40
                  times 8 db 40, -11
                  times 8 db 4, -1

                  times 8 db 0, 1
                  times 8 db -5, 17
                  times 8 db 58, -10
                  times 8 db 4, -1

ALIGN 32
tab_LumaCoeffVer_32: times 16 db 0, 0
                     times 16 db 0, 64
                     times 16 db 0, 0
                     times 16 db 0, 0

                     times 16 db -1, 4
                     times 16 db -10, 58
                     times 16 db 17, -5
                     times 16 db 1, 0

                     times 16 db -1, 4
                     times 16 db -11, 40
                     times 16 db 40, -11
                     times 16 db 4, -1

                     times 16 db 0, 1
                     times 16 db -5, 17
                     times 16 db 58, -10
                     times 16 db 4, -1

ALIGN 32
tab_ChromaCoeffVer_32: times 16 db 0, 64
                       times 16 db 0, 0

                       times 16 db -2, 58
                       times 16 db 10, -2

                       times 16 db -4, 54
                       times 16 db 16, -2

                       times 16 db -6, 46
                       times 16 db 28, -4

                       times 16 db -4, 36
                       times 16 db 36, -4

                       times 16 db -4, 28
                       times 16 db 46, -6

                       times 16 db -2, 16
                       times 16 db 54, -4

                       times 16 db -2, 10
                       times 16 db 58, -2

tab_c_64_n64:   times 8 db 64, -64

const interp4_shuf, times 2 db 0, 1, 8, 9, 4, 5, 12, 13, 2, 3, 10, 11, 6, 7, 14, 15

ALIGN 32
interp4_horiz_shuf1:    db 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6
                        db 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14

ALIGN 32
interp4_hpp_shuf: times 2 db 0, 1, 2, 3, 1, 2, 3, 4, 8, 9, 10, 11, 9, 10, 11, 12

ALIGN 32
interp8_hps_shuf: dd 0, 4, 1, 5, 2, 6, 3, 7

ALIGN 32
interp4_hps_shuf: times 2 db 0, 1, 2, 3, 1, 2, 3, 4, 8, 9, 10, 11, 9, 10, 11, 12

SECTION .text

cextern pb_128
cextern pw_1
cextern pw_512
cextern pw_2000

%macro FILTER_H4_w2_2 3
    movh        %2, [srcq - 1]
    pshufb      %2, %2, Tm0
    movh        %1, [srcq + srcstrideq - 1]
    pshufb      %1, %1, Tm0
    punpcklqdq  %2, %1
    pmaddubsw   %2, coef2
    phaddw      %2, %2
    pmulhrsw    %2, %3
    packuswb    %2, %2
    movd        r4, %2
    mov         [dstq], r4w
    shr         r4, 16
    mov         [dstq + dststrideq], r4w
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_2x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_2x4, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

%rep 2
FILTER_H4_w2_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
%endrep

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_2x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_2x8, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

%rep 4
FILTER_H4_w2_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
%endrep

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_2x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_2x16, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

mov         r5d,        16/2

.loop:
FILTER_H4_w2_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
dec         r5d
jnz         .loop

RET

%macro FILTER_H4_w4_2 3
    movh        %2, [srcq - 1]
    pshufb      %2, %2, Tm0
    pmaddubsw   %2, coef2
    movh        %1, [srcq + srcstrideq - 1]
    pshufb      %1, %1, Tm0
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    pmulhrsw    %2, %3
    packuswb    %2, %2
    movd        [dstq], %2
    palignr     %2, %2, 4
    movd        [dstq + dststrideq], %2
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4x2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_4x2, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

FILTER_H4_w4_2   t0, t1, t2

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_4x4, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

%rep 2
FILTER_H4_w4_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
%endrep

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_4x8, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

%rep 4
FILTER_H4_w4_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
%endrep

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_4x16, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

%rep 8
FILTER_H4_w4_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
%endrep

RET

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4x32(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_4x32, 4, 6, 5, src, srcstride, dst, dststride
%define coef2       m4
%define Tm0         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]

mov         r5d,        32/2

.loop:
FILTER_H4_w4_2   t0, t1, t2
lea         srcq,       [srcq + srcstrideq * 2]
lea         dstq,       [dstq + dststrideq * 2]
dec         r5d
jnz         .loop

RET

ALIGN 32
const interp_4tap_8x8_horiz_shuf,   dd 0, 4, 1, 5, 2, 6, 3, 7


%macro FILTER_H4_w6 3
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    pmulhrsw    %2, %3
    packuswb    %2, %2
    movd        [dstq],      %2
    pextrw      [dstq + 4], %2, 2
%endmacro

%macro FILTER_H4_w8 3
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    pmulhrsw    %2, %3
    packuswb    %2, %2
    movh        [dstq],      %2
%endmacro

%macro FILTER_H4_w12 3
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    pmulhrsw    %2, %3
    movu        %1, [srcq - 1 + 8]
    pshufb      %1, %1, Tm0
    pmaddubsw   %1, coef2
    phaddw      %1, %1
    pmulhrsw    %1, %3
    packuswb    %2, %1
    movh        [dstq],      %2
    pextrd      [dstq + 8], %2, 2
%endmacro

%macro FILTER_H4_w16 4
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq - 1 + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    pmulhrsw    %2, %3
    pmulhrsw    %4, %3
    packuswb    %2, %4
    movu        [dstq],      %2
%endmacro

%macro FILTER_H4_w24 4
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq - 1 + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    pmulhrsw    %2, %3
    pmulhrsw    %4, %3
    packuswb    %2, %4
    movu        [dstq],          %2
    movu        %1, [srcq - 1 + 16]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    pmulhrsw    %2, %3
    packuswb    %2, %2
    movh        [dstq + 16],     %2
%endmacro

%macro FILTER_H4_w32 4
    movu        %1, [srcq - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq - 1 + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    pmulhrsw    %2, %3
    pmulhrsw    %4, %3
    packuswb    %2, %4
    movu        [dstq],      %2
    movu        %1, [srcq - 1 + 16]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq - 1 + 24]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    pmulhrsw    %2, %3
    pmulhrsw    %4, %3
    packuswb    %2, %4
    movu        [dstq + 16],      %2
%endmacro

%macro FILTER_H4_w16o 5
    movu        %1, [srcq + %5 - 1]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + %5 - 1 + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    pmulhrsw    %2, %3
    pmulhrsw    %4, %3
    packuswb    %2, %4
    movu        [dstq + %5],      %2
%endmacro

%macro FILTER_H4_w48 4
    FILTER_H4_w16o %1, %2, %3, %4, 0
    FILTER_H4_w16o %1, %2, %3, %4, 16
    FILTER_H4_w16o %1, %2, %3, %4, 32
%endmacro

%macro FILTER_H4_w64 4
    FILTER_H4_w16o %1, %2, %3, %4, 0
    FILTER_H4_w16o %1, %2, %3, %4, 16
    FILTER_H4_w16o %1, %2, %3, %4, 32
    FILTER_H4_w16o %1, %2, %3, %4, 48
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro IPFILTER_CHROMA 2
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_%1x%2, 4, 6, 6, src, srcstride, dst, dststride
%define coef2       m5
%define Tm0         m4
%define Tm1         m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,        r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

mov           r5d,       %2

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]
mova        Tm1,         [tab_Tm + 16]

.loop:
FILTER_H4_w%1   t0, t1, t2
add         srcq,        srcstrideq
add         dstq,        dststrideq

dec         r5d
jnz        .loop

RET
%endmacro


IPFILTER_CHROMA 6,   8
IPFILTER_CHROMA 8,   2
IPFILTER_CHROMA 8,   4
IPFILTER_CHROMA 8,   6
IPFILTER_CHROMA 8,   8
IPFILTER_CHROMA 8,  16
IPFILTER_CHROMA 8,  32
IPFILTER_CHROMA 12, 16

IPFILTER_CHROMA 6,  16
IPFILTER_CHROMA 8,  12
IPFILTER_CHROMA 8,  64
IPFILTER_CHROMA 12, 32

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro IPFILTER_CHROMA_W 2
INIT_XMM sse4
cglobal interp_4tap_horiz_pp_%1x%2, 4, 6, 7, src, srcstride, dst, dststride
%define coef2       m6
%define Tm0         m5
%define Tm1         m4
%define t3          m3
%define t2          m2
%define t1          m1
%define t0          m0

mov         r4d,         r4m

%ifdef PIC
lea         r5,          [tab_ChromaCoeff]
movd        coef2,       [r5 + r4 * 4]
%else
movd        coef2,       [tab_ChromaCoeff + r4 * 4]
%endif

mov         r5d,          %2

pshufd      coef2,       coef2,      0
mova        t2,          [pw_512]
mova        Tm0,         [tab_Tm]
mova        Tm1,         [tab_Tm + 16]

.loop:
FILTER_H4_w%1   t0, t1, t2, t3
add         srcq,        srcstrideq
add         dstq,        dststrideq

dec         r5d
jnz        .loop

RET
%endmacro

IPFILTER_CHROMA_W 16,  4
IPFILTER_CHROMA_W 16,  8
IPFILTER_CHROMA_W 16, 12
IPFILTER_CHROMA_W 16, 16
IPFILTER_CHROMA_W 16, 32
IPFILTER_CHROMA_W 32,  8
IPFILTER_CHROMA_W 32, 16
IPFILTER_CHROMA_W 32, 24
IPFILTER_CHROMA_W 24, 32
IPFILTER_CHROMA_W 32, 32

IPFILTER_CHROMA_W 16, 24
IPFILTER_CHROMA_W 16, 64
IPFILTER_CHROMA_W 32, 48
IPFILTER_CHROMA_W 24, 64
IPFILTER_CHROMA_W 32, 64

IPFILTER_CHROMA_W 64, 64
IPFILTER_CHROMA_W 64, 32
IPFILTER_CHROMA_W 64, 48
IPFILTER_CHROMA_W 48, 64
IPFILTER_CHROMA_W 64, 16


%macro FILTER_H8_W8 7-8   ; t0, t1, t2, t3, coef, c512, src, dst
    movu        %1, %7
    pshufb      %2, %1, [tab_Lm +  0]
    pmaddubsw   %2, %5
    pshufb      %3, %1, [tab_Lm + 16]
    pmaddubsw   %3, %5
    phaddw      %2, %3
    pshufb      %4, %1, [tab_Lm + 32]
    pmaddubsw   %4, %5
    pshufb      %1, %1, [tab_Lm + 48]
    pmaddubsw   %1, %5
    phaddw      %4, %1
    phaddw      %2, %4
  %if %0 == 8
    pmulhrsw    %2, %6
    packuswb    %2, %2
    movh        %8, %2
  %endif
%endmacro

%macro FILTER_H8_W4 2
    movu        %1, [r0 - 3 + r5]
    pshufb      %2, %1, [tab_Lm]
    pmaddubsw   %2, m3
    pshufb      m7, %1, [tab_Lm + 16]
    pmaddubsw   m7, m3
    phaddw      %2, m7
    phaddw      %2, %2
%endmacro

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
%macro IPFILTER_LUMA 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4,7,8

    mov       r4d, r4m

%ifdef PIC
    lea       r6, [tab_LumaCoeff]
    movh      m3, [r6 + r4 * 8]
%else
    movh      m3, [tab_LumaCoeff + r4 * 8]
%endif
    punpcklqdq  m3, m3

%ifidn %3, pp 
    mova      m2, [pw_512]
%else
    mova      m2, [pw_2000]
%endif

    mov       r4d, %2
%ifidn %3, ps
    add       r3, r3
    cmp       r5m, byte 0
    je        .loopH
    lea       r6, [r1 + 2 * r1]
    sub       r0, r6
    add       r4d, 7
%endif

.loopH:
    xor       r5, r5
%rep %1 / 8
  %ifidn %3, pp 
    FILTER_H8_W8  m0, m1, m4, m5, m3, m2, [r0 - 3 + r5], [r2 + r5]
  %else
    FILTER_H8_W8  m0, m1, m4, m5, m3, UNUSED, [r0 - 3 + r5]
    psubw     m1, m2
    movu      [r2 + 2 * r5], m1
  %endif
    add       r5, 8
%endrep

%rep (%1 % 8) / 4
    FILTER_H8_W4  m0, m1
  %ifidn %3, pp 
    pmulhrsw  m1, m2
    packuswb  m1, m1
    movd      [r2 + r5], m1
  %else
    psubw     m1, m2
    movh      [r2 + 2 * r5], m1
  %endif
%endrep

    add       r0, r1
    add       r2, r3

    dec       r4d
    jnz       .loopH
    RET
%endmacro


INIT_YMM avx2
cglobal interp_8tap_horiz_pp_4x4, 4,6,6
    mov             r4d, r4m

%ifdef PIC
    lea             r5, [tab_LumaCoeff]
    vpbroadcastq    m0, [r5 + r4 * 8]
%else
    vpbroadcastq    m0, [tab_LumaCoeff + r4 * 8]
%endif

    mova            m1, [tab_Lm]
    vpbroadcastd    m2, [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    sub             r0, 3
    ; Row 0-1
    vbroadcasti128  m3, [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m3, m1
    pmaddubsw       m3, m0
    pmaddwd         m3, m2
    vbroadcasti128  m4, [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    phaddd          m3, m4                          ; DWORD [R1D R1C R0D R0C R1B R1A R0B R0A]

    ; Row 2-3
    lea             r0, [r0 + r1 * 2]
    vbroadcasti128  m4, [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    vbroadcasti128  m5, [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m1
    pmaddubsw       m5, m0
    pmaddwd         m5, m2
    phaddd          m4, m5                          ; DWORD [R3D R3C R2D R2C R3B R3A R2B R2A]

    packssdw        m3, m4                          ; WORD [R3D R3C R2D R2C R1D R1C R0D R0C R3B R3A R2B R2A R1B R1A R0B R0A]
    pmulhrsw        m3, [pw_512]
    vextracti128    xm4, m3, 1
    packuswb        xm3, xm4                        ; BYTE [R3D R3C R2D R2C R1D R1C R0D R0C R3B R3A R2B R2A R1B R1A R0B R0A]
    pshufb          xm3, [interp4_shuf]             ; [row3 row1 row2 row0]

    lea             r0, [r3 * 3]
    movd            [r2], xm3
    pextrd          [r2+r3], xm3, 2
    pextrd          [r2+r3*2], xm3, 1
    pextrd          [r2+r0], xm3, 3
    RET

%macro FILTER_HORIZ_LUMA_AVX2_4xN 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_horiz_pp_4x%1, 4, 6, 9
    mov             r4d, r4m

%ifdef PIC
    lea             r5, [tab_LumaCoeff]
    vpbroadcastq    m0, [r5 + r4 * 8]
%else
    vpbroadcastq    m0, [tab_LumaCoeff + r4 * 8]
%endif

    mova            m1, [tab_Lm]
    mova            m2, [pw_1]
    mova            m7, [interp8_hps_shuf]
    mova            m8, [pw_512]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1
    lea             r4, [r1 * 3]
    lea             r5, [r3 * 3]
    sub             r0, 3
%rep %1 / 8
    ; Row 0-1
    vbroadcasti128  m3, [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m3, m1
    pmaddubsw       m3, m0
    pmaddwd         m3, m2
    vbroadcasti128  m4, [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    phaddd          m3, m4                          ; DWORD [R1D R1C R0D R0C R1B R1A R0B R0A]

    ; Row 2-3
    vbroadcasti128  m4, [r0 + r1 * 2]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    vbroadcasti128  m5, [r0 + r4]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m1
    pmaddubsw       m5, m0
    pmaddwd         m5, m2
    phaddd          m4, m5                          ; DWORD [R3D R3C R2D R2C R3B R3A R2B R2A]

    packssdw        m3, m4                          ; WORD [R3D R3C R2D R2C R1D R1C R0D R0C R3B R3A R2B R2A R1B R1A R0B R0A]
    lea             r0, [r0 + r1 * 4]
    ; Row 4-5
    vbroadcasti128  m5, [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m1
    pmaddubsw       m5, m0
    pmaddwd         m5, m2
    vbroadcasti128  m4, [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    phaddd          m5, m4                          ; DWORD [R5D R5C R4D R4C R5B R5A R4B R4A]

    ; Row 6-7
    vbroadcasti128  m4, [r0 + r1 * 2]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddwd         m4, m2
    vbroadcasti128  m6, [r0 + r4]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m6, m1
    pmaddubsw       m6, m0
    pmaddwd         m6, m2
    phaddd          m4, m6                          ; DWORD [R7D R7C R6D R6C R7B R7A R6B R6A]

    packssdw        m5, m4                          ; WORD [R7D R7C R6D R6C R5D R5C R4D R4C R7B R7A R6B R6A R5B R5A R4B R4A]
    vpermd          m3, m7, m3
    vpermd          m5, m7, m5
    pmulhrsw        m3, m8
    pmulhrsw        m5, m8
    packuswb        m3, m5
    vextracti128    xm5, m3, 1

    movd            [r2], xm3
    pextrd          [r2 + r3], xm3, 1
    movd            [r2 + r3 * 2], xm5
    pextrd          [r2 + r5], xm5, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm3, 2
    pextrd          [r2 + r3], xm3, 3
    pextrd          [r2 + r3 * 2], xm5, 2
    pextrd          [r2 + r5], xm5, 3
    lea             r0, [r0 + r1 * 4]
    lea             r2, [r2 + r3 * 4]
%endrep
    RET
%endif
%endmacro

FILTER_HORIZ_LUMA_AVX2_4xN 8
FILTER_HORIZ_LUMA_AVX2_4xN 16

INIT_YMM avx2
cglobal interp_8tap_horiz_pp_8x4, 4, 6, 7
    mov             r4d, r4m

%ifdef PIC
    lea             r5, [tab_LumaCoeff]
    vpbroadcastq    m0, [r5 + r4 * 8]
%else
    vpbroadcastq    m0, [tab_LumaCoeff + r4 * 8]
%endif

    mova            m1, [tab_Lm]
    mova            m2, [tab_Lm + 32]

    ; register map
    ; m0     - interpolate coeff
    ; m1, m2 - shuffle order table

    sub             r0, 3
    lea             r5, [r1 * 3]
    lea             r4, [r3 * 3]

    ; Row 0
    vbroadcasti128  m3, [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m3, m2
    pshufb          m3, m1
    pmaddubsw       m3, m0
    pmaddubsw       m4, m0
    phaddw          m3, m4
    ; Row 1
    vbroadcasti128  m4, [r0 + r1]                   ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m4, m2
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddubsw       m5, m0
    phaddw          m4, m5

    phaddw          m3, m4                          ; WORD [R1H R1G R1D R1C R0H R0G R0D R0C R1F R1E R1B R1A R0F R0E R0B R0A]
    pmulhrsw        m3, [pw_512]

    ; Row 2
    vbroadcasti128  m4, [r0 + r1 * 2]               ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m4, m2
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddubsw       m5, m0
    phaddw          m4, m5
    ; Row 3
    vbroadcasti128  m5, [r0 + r5]                   ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m6, m5, m2
    pshufb          m5, m1
    pmaddubsw       m5, m0
    pmaddubsw       m6, m0
    phaddw          m5, m6

    phaddw          m4, m5                          ; WORD [R3H R3G R3D R3C R2H R2G R2D R2C R3F R3E R3B R3A R2F R2E R2B R2A]
    pmulhrsw        m4, [pw_512]

    packuswb        m3, m4
    vextracti128    xm4, m3, 1
    punpcklwd       xm5, xm3, xm4

    movq            [r2], xm5
    movhps          [r2 + r3], xm5

    punpckhwd       xm5, xm3, xm4
    movq            [r2 + r3 * 2], xm5
    movhps          [r2 + r4], xm5
    RET

%macro IPFILTER_LUMA_AVX2_8xN 2
INIT_YMM avx2
cglobal interp_8tap_horiz_pp_%1x%2, 4, 7, 7
    mov             r4d, r4m

%ifdef PIC
    lea             r5, [tab_LumaCoeff]
    vpbroadcastq    m0, [r5 + r4 * 8]
%else
    vpbroadcastq    m0, [tab_LumaCoeff + r4 * 8]
%endif

    mova            m1, [tab_Lm]
    mova            m2, [tab_Lm + 32]

    ; register map
    ; m0     - interpolate coeff
    ; m1, m2 - shuffle order table

    sub             r0, 3
    lea             r5, [r1 * 3]
    lea             r6, [r3 * 3]
    mov             r4d, %2 / 4
.loop:
    ; Row 0
    vbroadcasti128  m3, [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m4, m3, m2
    pshufb          m3, m1
    pmaddubsw       m3, m0
    pmaddubsw       m4, m0
    phaddw          m3, m4
    ; Row 1
    vbroadcasti128  m4, [r0 + r1]                   ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m4, m2
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddubsw       m5, m0
    phaddw          m4, m5

    phaddw          m3, m4                          ; WORD [R1H R1G R1D R1C R0H R0G R0D R0C R1F R1E R1B R1A R0F R0E R0B R0A]
    pmulhrsw        m3, [pw_512]

    ; Row 2
    vbroadcasti128  m4, [r0 + r1 * 2]               ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m5, m4, m2
    pshufb          m4, m1
    pmaddubsw       m4, m0
    pmaddubsw       m5, m0
    phaddw          m4, m5
    ; Row 3
    vbroadcasti128  m5, [r0 + r5]                   ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb          m6, m5, m2
    pshufb          m5, m1
    pmaddubsw       m5, m0
    pmaddubsw       m6, m0
    phaddw          m5, m6

    phaddw          m4, m5                          ; WORD [R3H R3G R3D R3C R2H R2G R2D R2C R3F R3E R3B R3A R2F R2E R2B R2A]
    pmulhrsw        m4, [pw_512]

    packuswb        m3, m4
    vextracti128    xm4, m3, 1
    punpcklwd       xm5, xm3, xm4

    movq            [r2], xm5
    movhps          [r2 + r3], xm5

    punpckhwd       xm5, xm3, xm4
    movq            [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm5

    lea             r0, [r0 + r1 * 4]
    lea             r2, [r2 + r3 * 4]
    dec             r4d
    jnz             .loop
    RET
%endmacro

IPFILTER_LUMA_AVX2_8xN 8, 8
IPFILTER_LUMA_AVX2_8xN 8, 16
IPFILTER_LUMA_AVX2_8xN 8, 32

%macro IPFILTER_LUMA_AVX2 2
INIT_YMM avx2
cglobal interp_8tap_horiz_pp_%1x%2, 4,6,8
    sub               r0,        3
    mov               r4d,       r4m
%ifdef PIC
    lea               r5,        [tab_LumaCoeff]
    vpbroadcastd      m0,        [r5 + r4 * 8]
    vpbroadcastd      m1,        [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,         [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,         [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,         [tab_Tm + 16]
    vpbroadcastd      m7,         [pw_1]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,        %2/2
.loop:
    ; Row 0
    vbroadcasti128    m4,         [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,         m4,     m3
    pshufb            m4,         [tab_Tm]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7
    vbroadcasti128    m5,         [r0 + 8]                    ; second 8 elements in Row0 
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m4,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,         [pw_512]
    vbroadcasti128    m2,         [r0 + r1]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,         m2,     m3
    pshufb            m2,         [tab_Tm]
    pmaddubsw         m2,         m0
    pmaddubsw         m5,         m1
    paddw             m2,         m5
    pmaddwd           m2,         m7
    vbroadcasti128    m5,         [r0 + r1 + 8]                    ; second 8 elements in Row0 
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m2,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m2,         [pw_512]
    packuswb          m4,         m2
    vpermq            m4,         m4,     11011000b
    vextracti128      xm5,        m4,     1
    pshufd            xm4,        xm4,    11011000b
    pshufd            xm5,        xm5,    11011000b
    movu              [r2],       xm4
    movu              [r2+r3],    xm5
    lea               r0,         [r0 + r1 * 2]
    lea               r2,         [r2 + r3 * 2]
    dec               r4d
    jnz              .loop
    RET
%endmacro

%macro IPFILTER_LUMA_32x_avx2 2
INIT_YMM avx2
cglobal interp_8tap_horiz_pp_%1x%2, 4,6,8
    sub               r0,         3
    mov               r4d,        r4m
%ifdef PIC
    lea               r5,         [tab_LumaCoeff]
    vpbroadcastd      m0,         [r5 + r4 * 8]
    vpbroadcastd      m1,         [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,         [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,         [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,         [tab_Tm + 16]
    vpbroadcastd      m7,         [pw_1]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,        %2
.loop:
    ; Row 0
    vbroadcasti128    m4,         [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,         m4,     m3
    pshufb            m4,         [tab_Tm]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7
    vbroadcasti128    m5,         [r0 + 8]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m4,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,         [pw_512]
    vbroadcasti128    m2,         [r0 + 16]
    pshufb            m5,         m2,     m3
    pshufb            m2,         [tab_Tm]
    pmaddubsw         m2,         m0
    pmaddubsw         m5,         m1
    paddw             m2,         m5
    pmaddwd           m2,         m7
    vbroadcasti128    m5,         [r0 + 24]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m2,         m5
    pmulhrsw          m2,         [pw_512]
    packuswb          m4,         m2
    vpermq            m4,         m4,     11011000b
    vextracti128      xm5,        m4,     1
    pshufd            xm4,        xm4,    11011000b
    pshufd            xm5,        xm5,    11011000b
    movu              [r2],       xm4
    movu              [r2 + 16],  xm5
    lea               r0,         [r0 + r1]
    lea               r2,         [r2 + r3]
    dec               r4d
    jnz               .loop
    RET
%endmacro

%macro IPFILTER_LUMA_64x_avx2 2
INIT_YMM avx2
cglobal interp_8tap_horiz_pp_%1x%2, 4,6,8
    sub               r0,    3
    mov               r4d,   r4m
%ifdef PIC
    lea               r5,        [tab_LumaCoeff]
    vpbroadcastd      m0,        [r5 + r4 * 8]
    vpbroadcastd      m1,        [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,        [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,        [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,        [tab_Tm + 16]
    vpbroadcastd      m7,        [pw_1]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,   %2
.loop:
    ; Row 0
    vbroadcasti128    m4,        [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,        m4,    m3
    pshufb            m4,        [tab_Tm]
    pmaddubsw         m4,        m0
    pmaddubsw         m5,        m1
    paddw             m4,        m5
    pmaddwd           m4,        m7
    vbroadcasti128    m5,        [r0 + 8]
    pshufb            m6,        m5,    m3
    pshufb            m5,        [tab_Tm]
    pmaddubsw         m5,        m0
    pmaddubsw         m6,        m1
    paddw             m5,        m6
    pmaddwd           m5,        m7
    packssdw          m4,        m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,        [pw_512]
    vbroadcasti128    m2,        [r0 + 16]
    pshufb            m5,        m2,    m3
    pshufb            m2,        [tab_Tm]
    pmaddubsw         m2,        m0
    pmaddubsw         m5,        m1
    paddw             m2,        m5
    pmaddwd           m2,        m7
    vbroadcasti128    m5,        [r0 + 24]
    pshufb            m6,        m5,    m3
    pshufb            m5,        [tab_Tm]
    pmaddubsw         m5,        m0
    pmaddubsw         m6,        m1
    paddw             m5,        m6
    pmaddwd           m5,        m7
    packssdw          m2,        m5
    pmulhrsw          m2,        [pw_512]
    packuswb          m4,        m2
    vpermq            m4,        m4,    11011000b
    vextracti128      xm5,       m4,    1
    pshufd            xm4,       xm4,   11011000b
    pshufd            xm5,       xm5,   11011000b
    movu              [r2],      xm4
    movu              [r2 + 16], xm5

    vbroadcasti128    m4,        [r0 + 32]
    pshufb            m5,        m4,    m3
    pshufb            m4,        [tab_Tm]
    pmaddubsw         m4,        m0
    pmaddubsw         m5,        m1
    paddw             m4,        m5
    pmaddwd           m4,        m7
    vbroadcasti128    m5,        [r0 + 40]
    pshufb            m6,        m5,    m3
    pshufb            m5,        [tab_Tm]
    pmaddubsw         m5,        m0
    pmaddubsw         m6,        m1
    paddw             m5,        m6
    pmaddwd           m5,        m7
    packssdw          m4,        m5
    pmulhrsw          m4,        [pw_512]
    vbroadcasti128    m2,        [r0 + 48]
    pshufb            m5,        m2,    m3
    pshufb            m2,        [tab_Tm]
    pmaddubsw         m2,        m0
    pmaddubsw         m5,        m1
    paddw             m2,        m5
    pmaddwd           m2,        m7
    vbroadcasti128    m5,        [r0 + 56]
    pshufb            m6,        m5,    m3
    pshufb            m5,        [tab_Tm]
    pmaddubsw         m5,        m0
    pmaddubsw         m6,        m1
    paddw             m5,        m6
    pmaddwd           m5,        m7
    packssdw          m2,        m5
    pmulhrsw          m2,        [pw_512]
    packuswb          m4,        m2
    vpermq            m4,        m4,    11011000b
    vextracti128      xm5,       m4,    1
    pshufd            xm4,       xm4,   11011000b
    pshufd            xm5,       xm5,   11011000b
    movu              [r2 +32],  xm4
    movu              [r2 + 48], xm5

    lea               r0,        [r0 + r1]
    lea               r2,        [r2 + r3]
    dec               r4d
    jnz               .loop
    RET
%endmacro

INIT_YMM avx2
cglobal interp_8tap_horiz_pp_48x64, 4,6,8
    sub               r0,         3
    mov               r4d,        r4m
%ifdef PIC
    lea               r5,         [tab_LumaCoeff]
    vpbroadcastd      m0,         [r5 + r4 * 8]
    vpbroadcastd      m1,         [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,         [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,         [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,         [tab_Tm + 16]
    vpbroadcastd      m7,         [pw_1]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,        64
.loop:
    ; Row 0
    vbroadcasti128    m4,         [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,         m4,     m3
    pshufb            m4,         [tab_Tm]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7
    vbroadcasti128    m5,         [r0 + 8]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m4,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,         [pw_512]

    vbroadcasti128    m2,         [r0 + 16]
    pshufb            m5,         m2,     m3
    pshufb            m2,         [tab_Tm]
    pmaddubsw         m2,         m0
    pmaddubsw         m5,         m1
    paddw             m2,         m5
    pmaddwd           m2,         m7
    vbroadcasti128    m5,         [r0 + 24]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m2,         m5
    pmulhrsw          m2,         [pw_512]
    packuswb          m4,         m2
    vpermq            m4,         m4,     11011000b
    vextracti128      xm5,        m4,     1
    pshufd            xm4,        xm4,    11011000b
    pshufd            xm5,        xm5,    11011000b
    movu              [r2],       xm4
    movu              [r2 + 16],  xm5

    vbroadcasti128    m4,         [r0 + 32]
    pshufb            m5,         m4,     m3
    pshufb            m4,         [tab_Tm]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7
    vbroadcasti128    m5,         [r0 + 40]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [tab_Tm]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m4,         m5
    pmulhrsw          m4,         [pw_512]
    packuswb          m4,         m4
    vpermq            m4,         m4,     11011000b
    pshufd            xm4,        xm4,    11011000b
    movu              [r2 + 32],  xm4

    lea               r0,         [r0 + r1]
    lea               r2,         [r2 + r3]
    dec               r4d
    jnz               .loop
    RET

INIT_YMM avx2 
cglobal interp_4tap_horiz_pp_4x4, 4,6,6
    mov             r4d, r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vpbroadcastd      m2,           [pw_1]
    vbroadcasti128    m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec                r0

    ; Row 0-1
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 2-3
    lea               r0,           [r0 + r1 * 2]
    vbroadcasti128    m4,           [r0]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    vinserti128       m4,           m4,      [r0 + r1],     1
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    pmulhrsw          m3,           [pw_512]
    vextracti128      xm4,          m3,     1
    packuswb          xm3,          xm4

    lea               r0,           [r3 * 3]
    movd              [r2],         xm3
    pextrd            [r2+r3],      xm3,     2
    pextrd            [r2+r3*2],    xm3,     1
    pextrd            [r2+r0],      xm3,     3
    RET

INIT_YMM avx2 
cglobal interp_4tap_horiz_pp_2x4, 4, 6, 3
    mov               r4d,           r4m

%ifdef PIC
    lea               r5,            [tab_ChromaCoeff]
    vpbroadcastd      m0,            [r5 + r4 * 4]
%else
    vpbroadcastd      m0,            [tab_ChromaCoeff + r4 * 4]
%endif

    dec               r0
    lea               r4,            [r1 * 3]
    movq              xm1,           [r0]
    movhps            xm1,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    movhps            xm2,           [r0 + r4]
    vinserti128       m1,            m1,          xm2,          1
    pshufb            m1,            [interp4_hpp_shuf]
    pmaddubsw         m1,            m0
    pmaddwd           m1,            [pw_1]
    vextracti128      xm2,           m1,          1
    packssdw          xm1,           xm2
    pmulhrsw          xm1,           [pw_512]
    packuswb          xm1,           xm1

    lea               r4,            [r3 * 3]
    pextrw            [r2],          xm1,         0
    pextrw            [r2 + r3],     xm1,         1
    pextrw            [r2 + r3 * 2], xm1,         2
    pextrw            [r2 + r4],     xm1,         3
    RET

INIT_YMM avx2 
cglobal interp_4tap_horiz_pp_2x8, 4, 6, 6
    mov               r4d,           r4m

%ifdef PIC
    lea               r5,            [tab_ChromaCoeff]
    vpbroadcastd      m0,            [r5 + r4 * 4]
%else
    vpbroadcastd      m0,            [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m4,            [interp4_hpp_shuf]
    mova              m5,            [pw_1]
    dec               r0
    lea               r4,            [r1 * 3]
    movq              xm1,           [r0]
    movhps            xm1,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    movhps            xm2,           [r0 + r4]
    vinserti128       m1,            m1,          xm2,          1
    lea               r0,            [r0 + r1 * 4]
    movq              xm3,           [r0]
    movhps            xm3,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    movhps            xm2,           [r0 + r4]
    vinserti128       m3,            m3,          xm2,          1

    pshufb            m1,            m4
    pshufb            m3,            m4
    pmaddubsw         m1,            m0
    pmaddubsw         m3,            m0
    pmaddwd           m1,            m5
    pmaddwd           m3,            m5
    packssdw          m1,            m3
    pmulhrsw          m1,            [pw_512]
    vextracti128      xm2,           m1,          1
    packuswb          xm1,           xm2

    lea               r4,            [r3 * 3]
    pextrw            [r2],          xm1,         0
    pextrw            [r2 + r3],     xm1,         1
    pextrw            [r2 + r3 * 2], xm1,         4
    pextrw            [r2 + r4],     xm1,         5
    lea               r2,            [r2 + r3 * 4]
    pextrw            [r2],          xm1,         2
    pextrw            [r2 + r3],     xm1,         3
    pextrw            [r2 + r3 * 2], xm1,         6
    pextrw            [r2 + r4],     xm1,         7
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_32x32, 4,6,7
    mov             r4d, r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]
    mova              m6,           [pw_512]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          32

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    vbroadcasti128    m4,           [r0 + 16]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + 20]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    movu              [r2],         m3
    lea               r2,           [r2 + r3]
    lea               r0,           [r0 + r1]
    dec               r4d
    jnz               .loop
    RET


INIT_YMM avx2
cglobal interp_4tap_horiz_pp_16x16, 4, 6, 7
    mov               r4d,          r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m6,           [pw_512]
    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          8

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + r1 + 4]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    vextracti128      xm4,          m3,       1
    movu              [r2],         xm3
    movu              [r2 + r3],    xm4
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]
    dec               r4d
    jnz               .loop
    RET
;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
    IPFILTER_LUMA 4, 4, pp
    IPFILTER_LUMA 4, 8, pp
    IPFILTER_LUMA 12, 16, pp
    IPFILTER_LUMA 4, 16, pp

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_8x8, 4,6,6
    mov               r4d,    r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    movu              m1,           [tab_Tm]
    vpbroadcastd      m2,           [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    sub               r0,           1
    mov               r4d,          2

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           [pw_512]
    lea               r0,           [r0 + r1 * 2]

    ; Row 2
    vbroadcasti128    m4,           [r0 ]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    ; Row 3
    vbroadcasti128    m5,           [r0 + r1]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           [pw_512]

    packuswb          m3,           m4
    mova              m5,           [interp_4tap_8x8_horiz_shuf]
    vpermd            m3,           m5,     m3
    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movhps            [r2 + r3],    xm3
    lea               r2,           [r2 + r3 * 2]
    movq              [r2],         xm4
    movhps            [r2 + r3],    xm4
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1*2]
    dec               r4d
    jnz               .loop
    RET

    IPFILTER_LUMA_AVX2 16, 4
    IPFILTER_LUMA_AVX2 16, 8
    IPFILTER_LUMA_AVX2 16, 12 
    IPFILTER_LUMA_AVX2 16, 16
    IPFILTER_LUMA_AVX2 16, 32
    IPFILTER_LUMA_AVX2 16, 64

    IPFILTER_LUMA_32x_avx2 32 , 8
    IPFILTER_LUMA_32x_avx2 32 , 16
    IPFILTER_LUMA_32x_avx2 32 , 24
    IPFILTER_LUMA_32x_avx2 32 , 32
    IPFILTER_LUMA_32x_avx2 32 , 64

    IPFILTER_LUMA_64x_avx2 64 , 64
    IPFILTER_LUMA_64x_avx2 64 , 48
    IPFILTER_LUMA_64x_avx2 64 , 32
    IPFILTER_LUMA_64x_avx2 64 , 16

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_8x2, 4, 6, 5
    mov               r4d,          r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,           [tab_Tm]
    mova              m2,           [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           [pw_512]
    vextracti128      xm4,          m3,          1
    packuswb          xm3,          xm4
    pshufd            xm3,          xm3,         11011000b
    movq              [r2],         xm3
    movhps            [r2 + r3],    xm3
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_8x6, 4, 6, 7
    mov               r4d,           r4m

%ifdef PIC
    lea               r5,            [tab_ChromaCoeff]
    vpbroadcastd      m0,            [r5 + r4 * 4]
%else
    vpbroadcastd      m0,            [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,            [tab_Tm]
    mova              m2,            [pw_1]
    mova              m6,            [pw_512]
    lea               r4,            [r1 * 3]
    lea               r5,            [r3 * 3]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    ; Row 0
    vbroadcasti128    m3,            [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,            m1
    pmaddubsw         m3,            m0
    pmaddwd           m3,            m2

    ; Row 1
    vbroadcasti128    m4,            [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,            m1
    pmaddubsw         m4,            m0
    pmaddwd           m4,            m2
    packssdw          m3,            m4
    pmulhrsw          m3,            m6

    ; Row 2
    vbroadcasti128    m4,            [r0 + r1 * 2]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,            m1
    pmaddubsw         m4,            m0
    pmaddwd           m4,            m2

    ; Row 3
    vbroadcasti128    m5,            [r0 + r4]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,            m1
    pmaddubsw         m5,            m0
    pmaddwd           m5,            m2
    packssdw          m4,            m5
    pmulhrsw          m4,            m6

    packuswb          m3,            m4
    mova              m5,            [interp8_hps_shuf]
    vpermd            m3,            m5,          m3
    vextracti128      xm4,           m3,          1
    movq              [r2],          xm3
    movhps            [r2 + r3],     xm3
    movq              [r2 + r3 * 2], xm4
    movhps            [r2 + r5],     xm4
    lea               r2,            [r2 + r3 * 4]
    lea               r0,            [r0 + r1 * 4]
    ; Row 4
    vbroadcasti128    m3,            [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,            m1
    pmaddubsw         m3,            m0
    pmaddwd           m3,            m2

    ; Row 5
    vbroadcasti128    m4,            [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,            m1
    pmaddubsw         m4,            m0
    pmaddwd           m4,            m2
    packssdw          m3,            m4
    pmulhrsw          m3,            m6
    vextracti128      xm4,           m3,          1
    packuswb          xm3,           xm4
    pshufd            xm3,           xm3,         11011000b
    movq              [r2],          xm3
    movhps            [r2 + r3],     xm3
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_6x8, 4, 6, 7
    mov               r4d,               r4m

%ifdef PIC
    lea               r5,                [tab_ChromaCoeff]
    vpbroadcastd      m0,                [r5 + r4 * 4]
%else
    vpbroadcastd      m0,                [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,                [tab_Tm]
    mova              m2,                [pw_1]
    mova              m6,                [pw_512]
    lea               r4,                [r1 * 3]
    lea               r5,                [r3 * 3]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
%rep 2
    ; Row 0
    vbroadcasti128    m3,                [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,                m1
    pmaddubsw         m3,                m0
    pmaddwd           m3,                m2

    ; Row 1
    vbroadcasti128    m4,                [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,                m1
    pmaddubsw         m4,                m0
    pmaddwd           m4,                m2
    packssdw          m3,                m4
    pmulhrsw          m3,                m6

    ; Row 2
    vbroadcasti128    m4,                [r0 + r1 * 2]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,                m1
    pmaddubsw         m4,                m0
    pmaddwd           m4,                m2

    ; Row 3
    vbroadcasti128    m5,                [r0 + r4]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,                m1
    pmaddubsw         m5,                m0
    pmaddwd           m5,                m2
    packssdw          m4,                m5
    pmulhrsw          m4,                m6

    packuswb          m3,                m4
    vextracti128      xm4,               m3,          1
    movd              [r2],              xm3
    pextrw            [r2 + 4],          xm4,         0
    pextrd            [r2 + r3],         xm3,         1
    pextrw            [r2 + r3 + 4],     xm4,         2
    pextrd            [r2 + r3 * 2],     xm3,         2
    pextrw            [r2 + r3 * 2 + 4], xm4,         4
    pextrd            [r2 + r5],         xm3,         3
    pextrw            [r2 + r5 + 4],     xm4,         6
    lea               r2,                [r2 + r3 * 4]
    lea               r0,                [r0 + r1 * 4]
%endrep
    RET

;-----------------------------------------------------------------------------------------------------------------------------
;void interp_horiz_ps_c(const pixel* src, intptr_t srcStride, int16_t* dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------

%macro IPFILTER_LUMA_PS_4xN_AVX2 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_horiz_ps_4x%1, 6,7,6
    mov                         r5d,               r5m
    mov                         r4d,               r4m
%ifdef PIC
    lea                         r6,                [tab_LumaCoeff]
    vpbroadcastq                m0,                [r6 + r4 * 8]
%else
    vpbroadcastq                m0,                [tab_LumaCoeff + r4 * 8]
%endif
    mova                        m1,                [tab_Lm]
    add                         r3d,               r3d
    vbroadcasti128              m2,                [pw_2000]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - pw_2000

    sub                         r0,                3
    test                        r5d,               r5d
    mov                         r5d,               %1                           ; loop count variable - height
    jz                         .preloop
    lea                         r6,                [r1 * 3]                     ; r8 = (N / 2 - 1) * srcStride
    sub                         r0,                r6                           ; r0(src) - 3 * srcStride
    add                         r5d,               7                            ; need extra 7 rows, just set a specially flag here, blkheight += N - 1  (7 - 3 = 4 ; since the last three rows not in loop)

.preloop:
    lea                         r6,                [r3 * 3]
.loop
    ; Row 0-1
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m3,                m1                           ; shuffled based on the col order tab_Lm
    pmaddubsw                   m3,                m0
    vbroadcasti128              m4,                [r0 + r1]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m1
    pmaddubsw                   m4,                m0
    phaddw                      m3,                m4                           ; DWORD [R1D R1C R0D R0C R1B R1A R0B R0A]

    ; Row 2-3
    lea                         r0,                [r0 + r1 * 2]                ;3rd row(i.e 2nd row)
    vbroadcasti128              m4,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m1
    pmaddubsw                   m4,                m0
    vbroadcasti128              m5,                [r0 + r1]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m5,                m1
    pmaddubsw                   m5,                m0
    phaddw                      m4,                m5                           ; DWORD [R3D R3C R2D R2C R3B R3A R2B R2A]
    phaddw                      m3,                m4                           ; all rows and col completed.

    mova                        m5,                [interp8_hps_shuf]
    vpermd                      m3,                m5,               m3
    psubw                       m3,                m2

    vextracti128                xm4,               m3,               1
    movq                        [r2],              xm3                          ;row 0
    movhps                      [r2 + r3],         xm3                          ;row 1
    movq                        [r2 + r3 * 2],     xm4                          ;row 2
    movhps                      [r2 + r6],         xm4                          ;row 3

    lea                         r0,                [r0 + r1 * 2]                ; first loop src ->5th row(i.e 4)
    lea                         r2,                [r2 + r3 * 4]                ; first loop dst ->5th row(i.e 4)
    sub                         r5d,               4
    jz                         .end
    cmp                         r5d,               4
    jge                        .loop

    ; Row 8-9
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m3,                m1
    pmaddubsw                   m3,                m0
    vbroadcasti128              m4,                [r0 + r1]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m1
    pmaddubsw                   m4,                m0
    phaddw                      m3,                m4                           ; DWORD [R1D R1C R0D R0C R1B R1A R0B R0A]

    ; Row 10
    vbroadcasti128              m4,                [r0 + r1 * 2]                ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m1
    pmaddubsw                   m4,                m0
    phaddw                      m4,                m4                           ; DWORD [R3D R3C R2D R2C R3B R3A R2B R2A]
    phaddw                      m3,                m4 

    vpermd                      m3,                m5,            m3            ; m5 don't broken in above
    psubw                       m3,                m2

    vextracti128                xm4,               m3,            1
    movq                        [r2],              xm3
    movhps                      [r2 + r3],         xm3
    movq                        [r2 + r3 * 2],     xm4
.end
    RET
%endif
%endmacro

    IPFILTER_LUMA_PS_4xN_AVX2 4
    IPFILTER_LUMA_PS_4xN_AVX2 8
    IPFILTER_LUMA_PS_4xN_AVX2 16

%macro IPFILTER_LUMA_PS_8xN_AVX2 1
; TODO: verify and enable on X86 mode
%if ARCH_X86_64 == 1
; void filter_hps(const pixel* src, intptr_t srcStride, int16_t* dst, intptr_t dstStride, int coeffIdx, int isRowExt)
INIT_YMM avx2
cglobal interp_8tap_horiz_ps_8x%1, 4,7,6
    mov                         r5d,        r5m
    mov                         r4d,        r4m
    shl                         r4d,        7
%ifdef PIC
    lea                         r6,         [pb_LumaCoeffVer]
    add                         r6,         r4
%else
    lea                         r6,         [pb_LumaCoeffVer + r4]
%endif
    add                         r3d,        r3d
    vpbroadcastd                m0,         [pw_2000]
    sub                         r0,         3
    lea                         r4,         [pb_8tap_hps_0]
    vbroadcasti128              m5,         [r4 + 0 * mmsize]

    ; check row count extend for interpolateHV
    test                        r5d,        r5d;
    mov                         r5d,        %1
    jz                         .enter_loop
    lea                         r4,         [r1 * 3]                        ; r8 = (N / 2 - 1) * srcStride
    sub                         r0,         r4                              ; r0(src)-r8
    add                         r5d,        8-1-2                           ; blkheight += N - 1  (7 - 3 = 4 ; since the last three rows not in loop)

.enter_loop:
    lea                         r4,         [pb_8tap_hps_0]

    ; ***** register map *****
    ; m0 - pw_2000
    ; r4 - base pointer of shuffle order table
    ; r5 - count of loop
    ; r6 - point to LumaCoeff
.loop:

    ; Row 0-1
    movu                        xm1,        [r0]
    movu                        xm2,        [r0 + r1]
    vinserti128                 m1,         m1,         xm2, 1
    pshufb                      m2,         m1,         m5                  ; [0 1 1 2 2 3 3 4 ...]
    pshufb                      m3,         m1,         [r4 + 1 * mmsize]   ; [2 3 3 4 4 5 5 6 ...]
    pshufb                      m4,         m1,         [r4 + 2 * mmsize]   ; [4 5 5 6 6 7 7 8 ...]
    pshufb                      m1,         m1,         [r4 + 3 * mmsize]   ; [6 7 7 8 8 9 9 A ...]
    pmaddubsw                   m2,         [r6 + 0 * mmsize]
    pmaddubsw                   m3,         [r6 + 1 * mmsize]
    pmaddubsw                   m4,         [r6 + 2 * mmsize]
    pmaddubsw                   m1,         [r6 + 3 * mmsize]
    paddw                       m2,         m3
    paddw                       m1,         m4
    paddw                       m1,         m2
    psubw                       m1,         m0

    vextracti128                xm2,        m1,         1
    movu                        [r2],       xm1                             ; row 0
    movu                        [r2 + r3],  xm2                             ; row 1

    lea                         r0,         [r0 + r1 * 2]                   ; first loop src ->5th row(i.e 4)
    lea                         r2,         [r2 + r3 * 2]                   ; first loop dst ->5th row(i.e 4)
    sub                         r5d,        2
    jg                         .loop
    jz                         .end             

    ; last row
    movu                        xm1,        [r0]
    pshufb                      xm2,        xm1,         xm5                ; [0 1 1 2 2 3 3 4 ...]
    pshufb                      xm3,        xm1,         [r4 + 1 * mmsize]  ; [2 3 3 4 4 5 5 6 ...]
    pshufb                      xm4,        xm1,         [r4 + 2 * mmsize]  ; [4 5 5 6 6 7 7 8 ...]
    pshufb                      xm1,        xm1,         [r4 + 3 * mmsize]  ; [6 7 7 8 8 9 9 A ...]
    pmaddubsw                   xm2,        [r6 + 0 * mmsize]
    pmaddubsw                   xm3,        [r6 + 1 * mmsize]
    pmaddubsw                   xm4,        [r6 + 2 * mmsize]
    pmaddubsw                   xm1,        [r6 + 3 * mmsize]
    paddw                       xm2,        xm3
    paddw                       xm1,        xm4
    paddw                       xm1,        xm2
    psubw                       xm1,        xm0
    movu                        [r2],       xm1                          ;row 0
.end
    RET
%endif
%endmacro ; IPFILTER_LUMA_PS_8xN_AVX2

IPFILTER_LUMA_PS_8xN_AVX2  4
IPFILTER_LUMA_PS_8xN_AVX2  8
IPFILTER_LUMA_PS_8xN_AVX2 16
IPFILTER_LUMA_PS_8xN_AVX2 32


%macro IPFILTER_LUMA_PS_16x_AVX2 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_horiz_ps_%1x%2, 6, 10, 7
    mov                         r5d,               r5m
    mov                         r4d,               r4m
%ifdef PIC
    lea                         r6,                [tab_LumaCoeff]
    vpbroadcastq                m0,                [r6 + r4 * 8]
%else
    vpbroadcastq                m0,                [tab_LumaCoeff + r4 * 8]
%endif
    mova                        m6,                [tab_Lm + 32]
    mova                        m1,                [tab_Lm]
    mov                         r9,                %2                           ;height
    add                         r3d,               r3d
    vbroadcasti128              m2,                [pw_2000]

    ; register map
    ; m0      - interpolate coeff
    ; m1 , m6 - shuffle order table
    ; m2      - pw_2000

    xor                         r7,                r7                          ; loop count variable
    sub                         r0,                3
    test                        r5d,               r5d
    jz                          .label
    lea                         r8,                [r1 * 3]                     ; r8 = (N / 2 - 1) * srcStride
    sub                         r0,                r8                           ; r0(src)-r8
    add                         r9,                7                            ; blkheight += N - 1  (7 - 1 = 6 ; since the last one row not in loop)

.label
    ; Row 0
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m3,             m6           ; row 0 (col 4 to 7)
    pshufb                      m3,                m1                           ; shuffled based on the col order tab_Lm row 0 (col 0 to 3)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    phaddw                      m3,                m4                           ; DWORD [R1D R1C R0D R0C R1B R1A R0B R0A]

    vbroadcasti128              m4,                [r0 + 8]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m5,                m4,            m6            ;row 1 (col 4 to 7)
    pshufb                      m4,                m1                           ;row 1 (col 0 to 3)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    phaddw                      m4,                m5                           ; DWORD [R3D R3C R2D R2C R3B R3A R2B R2A]
    phaddw                      m3,                m4                           ; all rows and col completed.

    mova                        m5,                [interp8_hps_shuf]
    vpermd                      m3,                m5,               m3
    psubw                       m3,                m2

    movu                        [r2],              m3                          ;row 0

    lea                         r0,                [r0 + r1]                ; first loop src ->5th row(i.e 4)
    lea                         r2,                [r2 + r3]                ; first loop dst ->5th row(i.e 4)
    dec                         r9d
    jnz                         .label

RET
%endif
%endmacro


IPFILTER_LUMA_PS_16x_AVX2 16 , 16
IPFILTER_LUMA_PS_16x_AVX2 16 , 8
IPFILTER_LUMA_PS_16x_AVX2 16 , 12
IPFILTER_LUMA_PS_16x_AVX2 16 , 4
IPFILTER_LUMA_PS_16x_AVX2 16 , 32
IPFILTER_LUMA_PS_16x_AVX2 16 , 64


;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro IPFILTER_LUMA_PP_W8 2
INIT_XMM sse4
cglobal interp_8tap_horiz_pp_%1x%2, 4,6,7
    mov         r4d, r4m

%ifdef PIC
    lea         r5, [tab_LumaCoeff]
    movh        m3, [r5 + r4 * 8]
%else
    movh        m3, [tab_LumaCoeff + r4 * 8]
%endif
    pshufd      m0, m3, 0                       ; m0 = coeff-L
    pshufd      m1, m3, 0x55                    ; m1 = coeff-H
    lea         r5, [tab_Tm]                    ; r5 = shuffle
    mova        m2, [pw_512]                    ; m2 = 512

    mov         r4d, %2
.loopH:
%assign x 0
%rep %1 / 8
    movu        m3, [r0 - 3 + x]                ; m3 = [F E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb      m4, m3, [r5 + 0*16]             ; m4 = [6 5 4 3 5 4 3 2 4 3 2 1 3 2 1 0]
    pshufb      m5, m3, [r5 + 1*16]             ; m5 = [A 9 8 7 9 8 7 6 8 7 6 5 7 6 5 4]
    pshufb          m3, [r5 + 2*16]             ; m3 = [E D C B D C B A C B A 9 B A 9 8]
    pmaddubsw   m4, m0
    pmaddubsw   m6, m5, m1
    pmaddubsw   m5, m0
    pmaddubsw   m3, m1
    paddw       m4, m6
    paddw       m5, m3
    phaddw      m4, m5
    pmulhrsw    m4, m2
    packuswb    m4, m4
    movh        [r2 + x], m4
%assign x x+8
%endrep

    add       r0, r1
    add       r2, r3

    dec       r4d
    jnz      .loopH
    RET
%endmacro

IPFILTER_LUMA_PP_W8      8,  4
IPFILTER_LUMA_PP_W8      8,  8
IPFILTER_LUMA_PP_W8      8, 16
IPFILTER_LUMA_PP_W8      8, 32
IPFILTER_LUMA_PP_W8     16,  4
IPFILTER_LUMA_PP_W8     16,  8
IPFILTER_LUMA_PP_W8     16, 12
IPFILTER_LUMA_PP_W8     16, 16
IPFILTER_LUMA_PP_W8     16, 32
IPFILTER_LUMA_PP_W8     16, 64
IPFILTER_LUMA_PP_W8     24, 32
IPFILTER_LUMA_PP_W8     32,  8
IPFILTER_LUMA_PP_W8     32, 16
IPFILTER_LUMA_PP_W8     32, 24
IPFILTER_LUMA_PP_W8     32, 32
IPFILTER_LUMA_PP_W8     32, 64
IPFILTER_LUMA_PP_W8     48, 64
IPFILTER_LUMA_PP_W8     64, 16
IPFILTER_LUMA_PP_W8     64, 32
IPFILTER_LUMA_PP_W8     64, 48
IPFILTER_LUMA_PP_W8     64, 64

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
    IPFILTER_LUMA 4, 4, ps
    IPFILTER_LUMA 8, 8, ps
    IPFILTER_LUMA 8, 4, ps
    IPFILTER_LUMA 4, 8, ps
    IPFILTER_LUMA 16, 16, ps
    IPFILTER_LUMA 16, 8, ps
    IPFILTER_LUMA 8, 16, ps
    IPFILTER_LUMA 16, 12, ps
    IPFILTER_LUMA 12, 16, ps
    IPFILTER_LUMA 16, 4, ps
    IPFILTER_LUMA 4, 16, ps
    IPFILTER_LUMA 32, 32, ps
    IPFILTER_LUMA 32, 16, ps
    IPFILTER_LUMA 16, 32, ps
    IPFILTER_LUMA 32, 24, ps
    IPFILTER_LUMA 24, 32, ps
    IPFILTER_LUMA 32, 8, ps
    IPFILTER_LUMA 8, 32, ps
    IPFILTER_LUMA 64, 64, ps
    IPFILTER_LUMA 64, 32, ps
    IPFILTER_LUMA 32, 64, ps
    IPFILTER_LUMA 64, 48, ps
    IPFILTER_LUMA 48, 64, ps
    IPFILTER_LUMA 64, 16, ps
    IPFILTER_LUMA 16, 64, ps

;-----------------------------------------------------------------------------
; Interpolate HV
;-----------------------------------------------------------------------------
%macro FILTER_HV8_START 7 ; (t0, t1, t2, t3, t4, off_src, off_coeff) -> (t3, t5), (t4, t1), [2]
    mova        %5, [r0 +  (%6 + 0) * 16]
    mova        %1, [r0 +  (%6 + 1) * 16]
    mova        %2, [r0 +  (%6 + 2) * 16]
    punpcklwd   %3, %5, %1
    punpckhwd   %5, %1
    pmaddwd     %3, [r5 + (%7) * 16]   ; R3 = L[0+1] -- Row 0
    pmaddwd     %5, [r5 + (%7) * 16]   ; R0 = H[0+1]
    punpcklwd   %4, %1, %2
    punpckhwd   %1, %2
    pmaddwd     %4, [r5 + (%7) * 16]   ; R4 = L[1+2] -- Row 1
    pmaddwd     %1, [r5 + (%7) * 16]   ; R1 = H[1+2]
%endmacro ; FILTER_HV8_START

%macro FILTER_HV8_MID 10 ; (Row3, prevRow, sum0L, sum1L, sum0H, sum1H, t6, t7, off_src, off_coeff) -> [6]
    mova        %8, [r0 +  (%9 + 0) * 16]
    mova        %1, [r0 +  (%9 + 1) * 16]
    punpcklwd   %7, %2, %8
    punpckhwd   %2, %8
    pmaddwd     %7, [r5 + %10 * 16]
    pmaddwd     %2, [r5 + %10 * 16]
    paddd       %3, %7              ; R3 = L[0+1+2+3] -- Row 0
    paddd       %5, %2              ; R0 = H[0+1+2+3]
    punpcklwd   %7, %8, %1
    punpckhwd   %8, %1
    pmaddwd     %7, [r5 + %10 * 16]
    pmaddwd     %8, [r5 + %10 * 16]
    paddd       %4, %7              ; R4 = L[1+2+3+4] -- Row 1
    paddd       %6, %8              ; R1 = H[1+2+3+4]
%endmacro ; FILTER_HV8_MID

; Round and Saturate
%macro FILTER_HV8_END 4 ; output in [1, 3]
    paddd       %1, [tab_c_526336]
    paddd       %2, [tab_c_526336]
    paddd       %3, [tab_c_526336]
    paddd       %4, [tab_c_526336]
    psrad       %1, 12
    psrad       %2, 12
    psrad       %3, 12
    psrad       %4, 12
    packssdw    %1, %2
    packssdw    %3, %4

    ; TODO: is merge better? I think this way is short dependency link
    packuswb    %1, %3
%endmacro ; FILTER_HV8_END

;-----------------------------------------------------------------------------
; void interp_8tap_hv_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int idxX, int idxY)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_8tap_hv_pp_8x8, 4, 7, 8, 0-15*16
%define coef        m7
%define stk_buf     rsp

    mov         r4d,        r4m
    mov         r5d,        r5m

%ifdef PIC
    lea         r6,         [tab_LumaCoeff]
    movh        coef,       [r6 + r4 * 8]
%else
    movh        coef,       [tab_LumaCoeff + r4 * 8]
%endif
    punpcklqdq  coef,       coef

    ; move to row -3
    lea         r6,         [r1 + r1 * 2]
    sub         r0,         r6

    xor         r6,         r6
    mov         r4,         rsp

.loopH:
    FILTER_H8_W8 m0, m1, m2, m3, coef, [pw_512], [r0 - 3]
    psubw       m1,         [pw_2000]
    mova        [r4],       m1

    add         r0,         r1
    add         r4,         16
    inc         r6
    cmp         r6,         8+7
    jnz         .loopH

    ; ready to phase V
    ; Here all of mN is free

    ; load coeff table
    shl         r5,         6
    lea         r6,         [tab_LumaCoeffV]
    lea         r5,         [r5 + r6]

    ; load intermedia buffer
    mov         r0,         stk_buf

    ; register mapping
    ; r0 - src
    ; r5 - coeff
    ; r6 - loop_i

    ; let's go
    xor         r6,         r6

    ; TODO: this loop have more than 70 instructions, I think it is more than Intel loop decode cache
.loopV:

    FILTER_HV8_START    m1, m2, m3, m4, m0,             0, 0
    FILTER_HV8_MID      m6, m2, m3, m4, m0, m1, m7, m5, 3, 1
    FILTER_HV8_MID      m5, m6, m3, m4, m0, m1, m7, m2, 5, 2
    FILTER_HV8_MID      m6, m5, m3, m4, m0, m1, m7, m2, 7, 3
    FILTER_HV8_END      m3, m0, m4, m1

    movh        [r2],       m3
    movhps      [r2 + r3],  m3

    lea         r0,         [r0 + 16 * 2]
    lea         r2,         [r2 + r3 * 2]

    inc         r6
    cmp         r6,         8/2
    jnz         .loopV

    RET

;-----------------------------------------------------------------------------
;void interp_4tap_vert_pp_2x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_pp_2x4, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif
lea         r4,        [r1 * 3]
lea         r5,        [r0 + 4 * r1]
pshufb      m0,        [tab_Cm]
mova        m1,        [pw_512]

movd        m2,        [r0]
movd        m3,        [r0 + r1]
movd        m4,        [r0 + 2 * r1]
movd        m5,        [r0 + r4]

punpcklbw   m2,        m3
punpcklbw   m6,        m4,        m5
punpcklbw   m2,        m6

pmaddubsw   m2,        m0

movd        m6,        [r5]

punpcklbw   m3,        m4
punpcklbw   m7,        m5,        m6
punpcklbw   m3,        m7

pmaddubsw   m3,        m0

phaddw      m2,        m3

pmulhrsw    m2,        m1

movd        m7,        [r5 + r1]

punpcklbw   m4,        m5
punpcklbw   m3,        m6,        m7
punpcklbw   m4,        m3

pmaddubsw   m4,        m0

movd        m3,        [r5 + 2 * r1]

punpcklbw   m5,        m6
punpcklbw   m7,        m3
punpcklbw   m5,        m7

pmaddubsw   m5,        m0

phaddw      m4,        m5

pmulhrsw    m4,        m1
packuswb    m2,        m4

pextrw      [r2],      m2, 0
pextrw      [r2 + r3], m2, 2
lea         r2,        [r2 + 2 * r3]
pextrw      [r2],      m2, 4
pextrw      [r2 + r3], m2, 6

RET

%macro FILTER_VER_CHROMA_AVX2_2x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_2x4, 4, 6, 2
    mov             r4d, r4m
    shl             r4d, 5
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeff_V]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeff_V + r4]
%endif

    lea             r4, [r1 * 3]

    pinsrw          xm1, [r0], 0
    pinsrw          xm1, [r0 + r1], 1
    pinsrw          xm1, [r0 + r1 * 2], 2
    pinsrw          xm1, [r0 + r4], 3
    lea             r0, [r0 + r1 * 4]
    pinsrw          xm1, [r0], 4
    pinsrw          xm1, [r0 + r1], 5
    pinsrw          xm1, [r0 + r1 * 2], 6

    pshufb          xm0, xm1, [interp_vert_shuf]
    pshufb          xm1, [interp_vert_shuf + 32]
    vinserti128     m0, m0, xm1, 1
    pmaddubsw       m0, [r5]
    vextracti128    xm1, m0, 1
    paddw           xm0, xm1
%ifidn %1,pp
    pmulhrsw        xm0, [pw_512]
    packuswb        xm0, xm0
    lea             r4, [r3 * 3]
    pextrw          [r2], xm0, 0
    pextrw          [r2 + r3], xm0, 1
    pextrw          [r2 + r3 * 2], xm0, 2
    pextrw          [r2 + r4], xm0, 3
%else
    add             r3d, r3d
    lea             r4, [r3 * 3]
    psubw           xm0, [pw_2000]
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    pextrd          [r2 + r3 * 2], xm0, 2
    pextrd          [r2 + r4], xm0, 3
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_2x4 pp
FILTER_VER_CHROMA_AVX2_2x4 ps

%macro FILTER_VER_CHROMA_AVX2_2x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_2x8, 4, 6, 2
    mov             r4d, r4m
    shl             r4d, 6
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]

    pinsrw          xm1, [r0], 0
    pinsrw          xm1, [r0 + r1], 1
    pinsrw          xm1, [r0 + r1 * 2], 2
    pinsrw          xm1, [r0 + r4], 3
    lea             r0, [r0 + r1 * 4]
    pinsrw          xm1, [r0], 4
    pinsrw          xm1, [r0 + r1], 5
    pinsrw          xm1, [r0 + r1 * 2], 6
    pinsrw          xm1, [r0 + r4], 7
    movhlps         xm0, xm1
    lea             r0, [r0 + r1 * 4]
    pinsrw          xm0, [r0], 4
    pinsrw          xm0, [r0 + r1], 5
    pinsrw          xm0, [r0 + r1 * 2], 6
    vinserti128     m1, m1, xm0, 1

    pshufb          m0, m1, [interp_vert_shuf]
    pshufb          m1, [interp_vert_shuf + 32]
    pmaddubsw       m0, [r5]
    pmaddubsw       m1, [r5 + 1 * mmsize]
    paddw           m0, m1
%ifidn %1,pp
    pmulhrsw        m0, [pw_512]
    vextracti128    xm1, m0, 1
    packuswb        xm0, xm1
    lea             r4, [r3 * 3]
    pextrw          [r2], xm0, 0
    pextrw          [r2 + r3], xm0, 1
    pextrw          [r2 + r3 * 2], xm0, 2
    pextrw          [r2 + r4], xm0, 3
    lea             r2, [r2 + r3 * 4]
    pextrw          [r2], xm0, 4
    pextrw          [r2 + r3], xm0, 5
    pextrw          [r2 + r3 * 2], xm0, 6
    pextrw          [r2 + r4], xm0, 7
%else
    add             r3d, r3d
    lea             r4, [r3 * 3]
    psubw           m0, [pw_2000]
    vextracti128    xm1, m0, 1
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    pextrd          [r2 + r3 * 2], xm0, 2
    pextrd          [r2 + r4], xm0, 3
    lea             r2, [r2 + r3 * 4]
    movd            [r2], xm1
    pextrd          [r2 + r3], xm1, 1
    pextrd          [r2 + r3 * 2], xm1, 2
    pextrd          [r2 + r4], xm1, 3
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_2x8 pp
FILTER_VER_CHROMA_AVX2_2x8 ps

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_2x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W2_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_2x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m0,        [tab_Cm]

mova        m1,        [pw_512]

mov         r4d,       %2
lea         r5,        [3 * r1]

.loop:
movd        m2,        [r0]
movd        m3,        [r0 + r1]
movd        m4,        [r0 + 2 * r1]
movd        m5,        [r0 + r5]

punpcklbw   m2,        m3
punpcklbw   m6,        m4,        m5
punpcklbw   m2,        m6

pmaddubsw   m2,        m0

lea         r0,        [r0 + 4 * r1]
movd        m6,        [r0]

punpcklbw   m3,        m4
punpcklbw   m7,        m5,        m6
punpcklbw   m3,        m7

pmaddubsw   m3,        m0

phaddw      m2,        m3

pmulhrsw    m2,        m1

movd        m7,        [r0 + r1]

punpcklbw   m4,        m5
punpcklbw   m3,        m6,        m7
punpcklbw   m4,        m3

pmaddubsw   m4,        m0

movd        m3,        [r0 + 2 * r1]

punpcklbw   m5,        m6
punpcklbw   m7,        m3
punpcklbw   m5,        m7

pmaddubsw   m5,        m0

phaddw      m4,        m5

pmulhrsw    m4,        m1
packuswb    m2,        m4

pextrw      [r2],      m2, 0
pextrw      [r2 + r3], m2, 2
lea         r2,        [r2 + 2 * r3]
pextrw      [r2],      m2, 4
pextrw      [r2 + r3], m2, 6

lea         r2,        [r2 + 2 * r3]

sub         r4,        4
jnz        .loop
RET
%endmacro

FILTER_V4_W2_H4 2, 8

FILTER_V4_W2_H4 2, 16

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_4x2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_pp_4x2, 4, 6, 6

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m0,        [tab_Cm]
lea         r5,        [r0 + 2 * r1]

movd        m2,        [r0]
movd        m3,        [r0 + r1]
movd        m4,        [r5]
movd        m5,        [r5 + r1]

punpcklbw   m2,        m3
punpcklbw   m1,        m4,        m5
punpcklbw   m2,        m1

pmaddubsw   m2,        m0

movd        m1,        [r0 + 4 * r1]

punpcklbw   m3,        m4
punpcklbw   m5,        m1
punpcklbw   m3,        m5

pmaddubsw   m3,        m0

phaddw      m2,        m3

pmulhrsw    m2,        [pw_512]
packuswb    m2,        m2
movd        [r2],      m2
pextrd      [r2 + r3], m2,  1

RET

%macro FILTER_VER_CHROMA_AVX2_4x2 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x2, 4, 6, 4
    mov             r4d, r4m
    shl             r4d, 5
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeff_V]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeff_V + r4]
%endif

    lea             r4, [r1 * 3]

    movd            xm1, [r0]
    movd            xm2, [r0 + r1]
    punpcklbw       xm1, xm2
    movd            xm3, [r0 + r1 * 2]
    punpcklbw       xm2, xm3
    movlhps         xm1, xm2
    movd            xm0, [r0 + r4]
    punpcklbw       xm3, xm0
    movd            xm2, [r0 + r1 * 4]
    punpcklbw       xm0, xm2
    movlhps         xm3, xm0
    vinserti128     m1, m1, xm3, 1                          ; m1 = row[x x x 4 3 2 1 0]

    pmaddubsw       m1, [r5]
    vextracti128    xm3, m1, 1
    paddw           xm1, xm3
%ifidn %1,pp
    pmulhrsw        xm1, [pw_512]
    packuswb        xm1, xm1
    movd            [r2], xm1
    pextrd          [r2 + r3], xm1, 1
%else
    add             r3d, r3d
    psubw           xm1, [pw_2000]
    movq            [r2], xm1
    movhps          [r2 + r3], xm1
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_4x2 pp
FILTER_VER_CHROMA_AVX2_4x2 ps

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_4x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_pp_4x4, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m0,        [tab_Cm]
mova        m1,        [pw_512]
lea         r5,        [r0 + 4 * r1]
lea         r4,        [r1 * 3]

movd        m2,        [r0]
movd        m3,        [r0 + r1]
movd        m4,        [r0 + 2 * r1]
movd        m5,        [r0 + r4]

punpcklbw   m2,        m3
punpcklbw   m6,        m4,        m5
punpcklbw   m2,        m6

pmaddubsw   m2,        m0

movd        m6,        [r5]

punpcklbw   m3,        m4
punpcklbw   m7,        m5,        m6
punpcklbw   m3,        m7

pmaddubsw   m3,        m0

phaddw      m2,        m3

pmulhrsw    m2,        m1

movd        m7,        [r5 + r1]

punpcklbw   m4,        m5
punpcklbw   m3,        m6,        m7
punpcklbw   m4,        m3

pmaddubsw   m4,        m0

movd        m3,        [r5 + 2 * r1]

punpcklbw   m5,        m6
punpcklbw   m7,        m3
punpcklbw   m5,        m7

pmaddubsw   m5,        m0

phaddw      m4,        m5

pmulhrsw    m4,        m1

packuswb    m2,        m4
movd        [r2],      m2
pextrd      [r2 + r3], m2, 1
lea         r2,        [r2 + 2 * r3]
pextrd      [r2],      m2, 2
pextrd      [r2 + r3], m2, 3
RET
%macro FILTER_VER_CHROMA_AVX2_4x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x4, 4, 6, 3
    mov             r4d, r4m
    shl             r4d, 6
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]

    movd            xm1, [r0]
    pinsrd          xm1, [r0 + r1], 1
    pinsrd          xm1, [r0 + r1 * 2], 2
    pinsrd          xm1, [r0 + r4], 3                       ; m1 = row[3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm2, [r0]
    pinsrd          xm2, [r0 + r1], 1
    pinsrd          xm2, [r0 + r1 * 2], 2                   ; m2 = row[x 6 5 4]
    vinserti128     m1, m1, xm2, 1                          ; m1 = row[x 6 5 4 3 2 1 0]
    mova            m2, [interp4_vpp_shuf1]
    vpermd          m0, m2, m1                              ; m0 = row[4 3 3 2 2 1 1 0]
    mova            m2, [interp4_vpp_shuf1 + mmsize]
    vpermd          m1, m2, m1                              ; m1 = row[6 5 5 4 4 3 3 2]

    mova            m2, [interp4_vpp_shuf]
    pshufb          m0, m0, m2
    pshufb          m1, m1, m2
    pmaddubsw       m0, [r5]
    pmaddubsw       m1, [r5 + mmsize]
    paddw           m0, m1                                  ; m0 = WORD ROW[3 2 1 0]
%ifidn %1,pp
    pmulhrsw        m0, [pw_512]
    vextracti128    xm1, m0, 1
    packuswb        xm0, xm1
    lea             r5, [r3 * 3]
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    pextrd          [r2 + r3 * 2], xm0, 2
    pextrd          [r2 + r5], xm0, 3
%else
    add             r3d, r3d
    psubw           m0, [pw_2000]
    vextracti128    xm1, m0, 1
    lea             r5, [r3 * 3]
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm1
    movhps          [r2 + r5], xm1
%endif
    RET
%endmacro
FILTER_VER_CHROMA_AVX2_4x4 pp
FILTER_VER_CHROMA_AVX2_4x4 ps

%macro FILTER_VER_CHROMA_AVX2_4x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x8, 4, 6, 5
    mov             r4d, r4m
    shl             r4d, 6
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]

    movd            xm1, [r0]
    pinsrd          xm1, [r0 + r1], 1
    pinsrd          xm1, [r0 + r1 * 2], 2
    pinsrd          xm1, [r0 + r4], 3                       ; m1 = row[3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm2, [r0]
    pinsrd          xm2, [r0 + r1], 1
    pinsrd          xm2, [r0 + r1 * 2], 2
    pinsrd          xm2, [r0 + r4], 3                       ; m2 = row[7 6 5 4]
    vinserti128     m1, m1, xm2, 1                          ; m1 = row[7 6 5 4 3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm3, [r0]
    pinsrd          xm3, [r0 + r1], 1
    pinsrd          xm3, [r0 + r1 * 2], 2                   ; m3 = row[x 10 9 8]
    vinserti128     m2, m2, xm3, 1                          ; m2 = row[x 10 9 8 7 6 5 4]
    mova            m3, [interp4_vpp_shuf1]
    vpermd          m0, m3, m1                              ; m0 = row[4 3 3 2 2 1 1 0]
    vpermd          m4, m3, m2                              ; m4 = row[8 7 7 6 6 5 5 4]
    mova            m3, [interp4_vpp_shuf1 + mmsize]
    vpermd          m1, m3, m1                              ; m1 = row[6 5 5 4 4 3 3 2]
    vpermd          m2, m3, m2                              ; m2 = row[10 9 9 8 8 7 7 6]

    mova            m3, [interp4_vpp_shuf]
    pshufb          m0, m0, m3
    pshufb          m1, m1, m3
    pshufb          m2, m2, m3
    pshufb          m4, m4, m3
    pmaddubsw       m0, [r5]
    pmaddubsw       m4, [r5]
    pmaddubsw       m1, [r5 + mmsize]
    pmaddubsw       m2, [r5 + mmsize]
    paddw           m0, m1                                  ; m0 = WORD ROW[3 2 1 0]
    paddw           m4, m2                                  ; m4 = WORD ROW[7 6 5 4]
%ifidn %1,pp
    pmulhrsw        m0, [pw_512]
    pmulhrsw        m4, [pw_512]
    packuswb        m0, m4
    vextracti128    xm1, m0, 1
    lea             r5, [r3 * 3]
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    movd            [r2 + r3 * 2], xm1
    pextrd          [r2 + r5], xm1, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm0, 3
    pextrd          [r2 + r3 * 2], xm1, 2
    pextrd          [r2 + r5], xm1, 3
%else
    add             r3d, r3d
    psubw           m0, [pw_2000]
    psubw           m4, [pw_2000]
    vextracti128    xm1, m0, 1
    vextracti128    xm2, m4, 1
    lea             r5, [r3 * 3]
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm1
    movhps          [r2 + r5], xm1
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r5], xm2
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_4x8 pp
FILTER_VER_CHROMA_AVX2_4x8 ps

%macro FILTER_VER_CHROMA_AVX2_4x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_4x16, 4, 6, 9
    mov             r4d, r4m
    shl             r4d, 6
    sub             r0, r1

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]

    movd            xm1, [r0]
    pinsrd          xm1, [r0 + r1], 1
    pinsrd          xm1, [r0 + r1 * 2], 2
    pinsrd          xm1, [r0 + r4], 3                       ; m1 = row[3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm2, [r0]
    pinsrd          xm2, [r0 + r1], 1
    pinsrd          xm2, [r0 + r1 * 2], 2
    pinsrd          xm2, [r0 + r4], 3                       ; m2 = row[7 6 5 4]
    vinserti128     m1, m1, xm2, 1                          ; m1 = row[7 6 5 4 3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm3, [r0]
    pinsrd          xm3, [r0 + r1], 1
    pinsrd          xm3, [r0 + r1 * 2], 2
    pinsrd          xm3, [r0 + r4], 3                       ; m3 = row[11 10 9 8]
    vinserti128     m2, m2, xm3, 1                          ; m2 = row[11 10 9 8 7 6 5 4]
    lea             r0, [r0 + r1 * 4]
    movd            xm4, [r0]
    pinsrd          xm4, [r0 + r1], 1
    pinsrd          xm4, [r0 + r1 * 2], 2
    pinsrd          xm4, [r0 + r4], 3                       ; m4 = row[15 14 13 12]
    vinserti128     m3, m3, xm4, 1                          ; m3 = row[15 14 13 12 11 10 9 8]
    lea             r0, [r0 + r1 * 4]
    movd            xm5, [r0]
    pinsrd          xm5, [r0 + r1], 1
    pinsrd          xm5, [r0 + r1 * 2], 2                   ; m5 = row[x 18 17 16]
    vinserti128     m4, m4, xm5, 1                          ; m4 = row[x 18 17 16 15 14 13 12]
    mova            m5, [interp4_vpp_shuf1]
    vpermd          m0, m5, m1                              ; m0 = row[4 3 3 2 2 1 1 0]
    vpermd          m6, m5, m2                              ; m6 = row[8 7 7 6 6 5 5 4]
    vpermd          m7, m5, m3                              ; m7 = row[12 11 11 10 10 9 9 8]
    vpermd          m8, m5, m4                              ; m8 = row[16 15 15 14 14 13 13 12]
    mova            m5, [interp4_vpp_shuf1 + mmsize]
    vpermd          m1, m5, m1                              ; m1 = row[6 5 5 4 4 3 3 2]
    vpermd          m2, m5, m2                              ; m2 = row[10 9 9 8 8 7 7 6]
    vpermd          m3, m5, m3                              ; m3 = row[14 13 13 12 12 11 11 10]
    vpermd          m4, m5, m4                              ; m4 = row[18 17 17 16 16 15 15 14]

    mova            m5, [interp4_vpp_shuf]
    pshufb          m0, m0, m5
    pshufb          m1, m1, m5
    pshufb          m2, m2, m5
    pshufb          m4, m4, m5
    pshufb          m3, m3, m5
    pshufb          m6, m6, m5
    pshufb          m7, m7, m5
    pshufb          m8, m8, m5
    pmaddubsw       m0, [r5]
    pmaddubsw       m6, [r5]
    pmaddubsw       m7, [r5]
    pmaddubsw       m8, [r5]
    pmaddubsw       m1, [r5 + mmsize]
    pmaddubsw       m2, [r5 + mmsize]
    pmaddubsw       m3, [r5 + mmsize]
    pmaddubsw       m4, [r5 + mmsize]
    paddw           m0, m1                                  ; m0 = WORD ROW[3 2 1 0]
    paddw           m6, m2                                  ; m6 = WORD ROW[7 6 5 4]
    paddw           m7, m3                                  ; m7 = WORD ROW[11 10 9 8]
    paddw           m8, m4                                  ; m8 = WORD ROW[15 14 13 12]
%ifidn %1,pp
    mova            m5, [pw_512]
    pmulhrsw        m0, m5
    pmulhrsw        m6, m5
    pmulhrsw        m7, m5
    pmulhrsw        m8, m5
    packuswb        m0, m6
    packuswb        m7, m8
    vextracti128    xm1, m0, 1
    vextracti128    xm2, m7, 1
    lea             r5, [r3 * 3]
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    movd            [r2 + r3 * 2], xm1
    pextrd          [r2 + r5], xm1, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm0, 3
    pextrd          [r2 + r3 * 2], xm1, 2
    pextrd          [r2 + r5], xm1, 3
    lea             r2, [r2 + r3 * 4]
    movd            [r2], xm7
    pextrd          [r2 + r3], xm7, 1
    movd            [r2 + r3 * 2], xm2
    pextrd          [r2 + r5], xm2, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm7, 2
    pextrd          [r2 + r3], xm7, 3
    pextrd          [r2 + r3 * 2], xm2, 2
    pextrd          [r2 + r5], xm2, 3
%else
    add             r3d, r3d
    mova            m5, [pw_2000]
    psubw           m0, m5
    psubw           m6, m5
    psubw           m7, m5
    psubw           m8, m5
    vextracti128    xm1, m0, 1
    vextracti128    xm2, m6, 1
    vextracti128    xm3, m7, 1
    vextracti128    xm4, m8, 1
    lea             r5, [r3 * 3]
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm1
    movhps          [r2 + r5], xm1
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm6
    movhps          [r2 + r3], xm6
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r5], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm7
    movhps          [r2 + r3], xm7
    movq            [r2 + r3 * 2], xm3
    movhps          [r2 + r5], xm3
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm8
    movhps          [r2 + r3], xm8
    movq            [r2 + r3 * 2], xm4
    movhps          [r2 + r5], xm4
%endif
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_4x16 pp
FILTER_VER_CHROMA_AVX2_4x16 ps

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W4_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_%1x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m0,        [tab_Cm]

mova        m1,        [pw_512]

mov         r4d,       %2

lea         r5,        [3 * r1]

.loop:
movd        m2,        [r0]
movd        m3,        [r0 + r1]
movd        m4,        [r0 + 2 * r1]
movd        m5,        [r0 + r5]

punpcklbw   m2,        m3
punpcklbw   m6,        m4,        m5
punpcklbw   m2,        m6

pmaddubsw   m2,        m0

lea         r0,        [r0 + 4 * r1]
movd        m6,        [r0]

punpcklbw   m3,        m4
punpcklbw   m7,        m5,        m6
punpcklbw   m3,        m7

pmaddubsw   m3,        m0

phaddw      m2,        m3

pmulhrsw    m2,        m1

movd        m7,        [r0 + r1]

punpcklbw   m4,        m5
punpcklbw   m3,        m6,        m7
punpcklbw   m4,        m3

pmaddubsw   m4,        m0

movd        m3,        [r0 + 2 * r1]

punpcklbw   m5,        m6
punpcklbw   m7,        m3
punpcklbw   m5,        m7

pmaddubsw   m5,        m0

phaddw      m4,        m5

pmulhrsw    m4,        m1
packuswb    m2,        m4
movd        [r2],      m2
pextrd      [r2 + r3], m2,  1
lea         r2,        [r2 + 2 * r3]
pextrd      [r2],      m2, 2
pextrd      [r2 + r3], m2, 3

lea         r2,        [r2 + 2 * r3]

sub         r4,        4
jnz        .loop
RET
%endmacro

FILTER_V4_W4_H4 4,  8
FILTER_V4_W4_H4 4, 16

FILTER_V4_W4_H4 4, 32

%macro FILTER_V4_W8_H2 0
punpcklbw   m1,        m2
punpcklbw   m7,        m3,        m0

pmaddubsw   m1,        m6
pmaddubsw   m7,        m5

paddw       m1,        m7

pmulhrsw    m1,        m4
packuswb    m1,        m1
%endmacro

%macro FILTER_V4_W8_H3 0
punpcklbw   m2,        m3
punpcklbw   m7,        m0,        m1

pmaddubsw   m2,        m6
pmaddubsw   m7,        m5

paddw       m2,        m7

pmulhrsw    m2,        m4
packuswb    m2,        m2
%endmacro

%macro FILTER_V4_W8_H4 0
punpcklbw   m3,        m0
punpcklbw   m7,        m1,        m2

pmaddubsw   m3,        m6
pmaddubsw   m7,        m5

paddw       m3,        m7

pmulhrsw    m3,        m4
packuswb    m3,        m3
%endmacro

%macro FILTER_V4_W8_H5 0
punpcklbw   m0,        m1
punpcklbw   m7,        m2,        m3

pmaddubsw   m0,        m6
pmaddubsw   m7,        m5

paddw       m0,        m7

pmulhrsw    m0,        m4
packuswb    m0,        m0
%endmacro

%macro FILTER_V4_W8_8x2 2
FILTER_V4_W8 %1, %2
movq        m0,        [r0 + 4 * r1]

FILTER_V4_W8_H2

movh        [r2 + r3], m1
%endmacro

%macro FILTER_V4_W8_8x4 2
FILTER_V4_W8_8x2 %1, %2
;8x3
lea         r6,        [r0 + 4 * r1]
movq        m1,        [r6 + r1]

FILTER_V4_W8_H3

movh        [r2 + 2 * r3], m2

;8x4
movq        m2,        [r6 + 2 * r1]

FILTER_V4_W8_H4

lea         r5,        [r2 + 2 * r3]
movh        [r5 + r3], m3
%endmacro

%macro FILTER_V4_W8_8x6 2
FILTER_V4_W8_8x4 %1, %2
;8x5
lea         r6,        [r6 + 2 * r1]
movq        m3,        [r6 + r1]

FILTER_V4_W8_H5

movh        [r2 + 4 * r3], m0

;8x6
movq        m0,        [r0 + 8 * r1]

FILTER_V4_W8_H2

lea         r5,        [r2 + 4 * r3]
movh        [r5 + r3], m1
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W8 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_%1x%2, 4, 7, 8

mov         r4d,       r4m

sub         r0,        r1
movq        m0,        [r0]
movq        m1,        [r0 + r1]
movq        m2,        [r0 + 2 * r1]
lea         r5,        [r0 + 2 * r1]
movq        m3,        [r5 + r1]

punpcklbw   m0,        m1
punpcklbw   m4,        m2,          m3

%ifdef PIC
lea         r6,        [tab_ChromaCoeff]
movd        m5,        [r6 + r4 * 4]
%else
movd        m5,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m6,        m5,       [tab_Vm]
pmaddubsw   m0,        m6

pshufb      m5,        [tab_Vm + 16]
pmaddubsw   m4,        m5

paddw       m0,        m4

mova        m4,        [pw_512]

pmulhrsw    m0,        m4
packuswb    m0,        m0
movh        [r2],      m0
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_8x2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
FILTER_V4_W8_8x2 8, 2

RET

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_8x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
FILTER_V4_W8_8x4 8, 4

RET

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_8x6(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
FILTER_V4_W8_8x6 8, 6

RET

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_4x2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_ps_4x2, 4, 6, 6

mov         r4d, r4m
sub         r0, r1
add         r3d, r3d

%ifdef PIC
lea         r5, [tab_ChromaCoeff]
movd        m0, [r5 + r4 * 4]
%else
movd        m0, [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m0, [tab_Cm]

movd        m2, [r0]
movd        m3, [r0 + r1]
lea         r5, [r0 + 2 * r1]
movd        m4, [r5]
movd        m5, [r5 + r1]

punpcklbw   m2, m3
punpcklbw   m1, m4, m5
punpcklbw   m2, m1

pmaddubsw   m2, m0

movd        m1, [r0 + 4 * r1]

punpcklbw   m3, m4
punpcklbw   m5, m1
punpcklbw   m3, m5

pmaddubsw   m3, m0

phaddw      m2, m3

psubw       m2, [pw_2000]
movh        [r2], m2
movhps      [r2 + r3], m2

RET

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_4x4(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_ps_4x4, 4, 6, 7

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m0, [tab_Cm]

    lea        r4, [r1 * 3]
    lea        r5, [r0 + 4 * r1]

    movd       m2, [r0]
    movd       m3, [r0 + r1]
    movd       m4, [r0 + 2 * r1]
    movd       m5, [r0 + r4]

    punpcklbw  m2, m3
    punpcklbw  m6, m4, m5
    punpcklbw  m2, m6

    pmaddubsw  m2, m0

    movd       m6, [r5]

    punpcklbw  m3, m4
    punpcklbw  m1, m5, m6
    punpcklbw  m3, m1

    pmaddubsw  m3, m0

    phaddw     m2, m3

    mova       m1, [pw_2000]

    psubw      m2, m1
    movh       [r2], m2
    movhps     [r2 + r3], m2

    movd       m2, [r5 + r1]

    punpcklbw  m4, m5
    punpcklbw  m3, m6, m2
    punpcklbw  m4, m3

    pmaddubsw  m4, m0

    movd       m3, [r5 + 2 * r1]

    punpcklbw  m5, m6
    punpcklbw  m2, m3
    punpcklbw  m5, m2

    pmaddubsw  m5, m0

    phaddw     m4, m5

    psubw      m4, m1
    lea        r2, [r2 + 2 * r3]
    movh       [r2], m4
    movhps     [r2 + r3], m4

    RET

;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_%1x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W4_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m0, [tab_Cm]

    mova       m1, [pw_2000]

    mov        r4d, %2/4
    lea        r5, [3 * r1]

.loop:
    movd       m2, [r0]
    movd       m3, [r0 + r1]
    movd       m4, [r0 + 2 * r1]
    movd       m5, [r0 + r5]

    punpcklbw  m2, m3
    punpcklbw  m6, m4, m5
    punpcklbw  m2, m6

    pmaddubsw  m2, m0

    lea        r0, [r0 + 4 * r1]
    movd       m6, [r0]

    punpcklbw  m3, m4
    punpcklbw  m7, m5, m6
    punpcklbw  m3, m7

    pmaddubsw  m3, m0

    phaddw     m2, m3

    psubw      m2, m1
    movh       [r2], m2
    movhps     [r2 + r3], m2

    movd       m2, [r0 + r1]

    punpcklbw  m4, m5
    punpcklbw  m3, m6, m2
    punpcklbw  m4, m3

    pmaddubsw  m4, m0

    movd       m3, [r0 + 2 * r1]

    punpcklbw  m5, m6
    punpcklbw  m2, m3
    punpcklbw  m5, m2

    pmaddubsw  m5, m0

    phaddw     m4, m5

    psubw      m4, m1
    lea        r2, [r2 + 2 * r3]
    movh       [r2], m4
    movhps     [r2 + r3], m4

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W4_H4 4, 8
FILTER_V_PS_W4_H4 4, 16

FILTER_V_PS_W4_H4 4, 32

;--------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_8x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W8_H8_H16_H2 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 6, 7

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m5, [r5 + r4 * 4]
%else
    movd       m5, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m6, m5, [tab_Vm]
    pshufb     m5, [tab_Vm + 16]
    mova       m4, [pw_2000]

    mov        r4d, %2/2
    lea        r5, [3 * r1]

.loopH:
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    movq       m2, [r0 + 2 * r1]
    movq       m3, [r0 + r5]

    punpcklbw  m0, m1
    punpcklbw  m1, m2
    punpcklbw  m2, m3

    pmaddubsw  m0, m6
    pmaddubsw  m2, m5

    paddw      m0, m2

    psubw      m0, m4
    movu       [r2], m0

    movq       m0, [r0 + 4 * r1]

    punpcklbw  m3, m0

    pmaddubsw  m1, m6
    pmaddubsw  m3, m5

    paddw      m1, m3
    psubw      m1, m4

    movu       [r2 + r3], m1

    lea        r0, [r0 + 2 * r1]
    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz       .loopH

    RET
%endmacro

FILTER_V_PS_W8_H8_H16_H2 8, 2
FILTER_V_PS_W8_H8_H16_H2 8, 4
FILTER_V_PS_W8_H8_H16_H2 8, 6

FILTER_V_PS_W8_H8_H16_H2 8, 12
FILTER_V_PS_W8_H8_H16_H2 8, 64

;--------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_8x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W8_H8_H16_H32 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m5, [r5 + r4 * 4]
%else
    movd       m5, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m6, m5, [tab_Vm]
    pshufb     m5, [tab_Vm + 16]
    mova       m4, [pw_2000]

    mov        r4d, %2/4
    lea        r5, [3 * r1]

.loop:
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    movq       m2, [r0 + 2 * r1]
    movq       m3, [r0 + r5]

    punpcklbw  m0, m1
    punpcklbw  m1, m2
    punpcklbw  m2, m3

    pmaddubsw  m0, m6
    pmaddubsw  m7, m2, m5

    paddw      m0, m7

    psubw       m0, m4
    movu       [r2], m0

    lea        r0, [r0 + 4 * r1]
    movq       m0, [r0]

    punpcklbw  m3, m0

    pmaddubsw  m1, m6
    pmaddubsw  m7, m3, m5

    paddw      m1, m7

    psubw      m1, m4
    movu       [r2 + r3], m1

    movq       m1, [r0 + r1]

    punpcklbw  m0, m1

    pmaddubsw  m2, m6
    pmaddubsw  m0, m5

    paddw      m2, m0

    psubw      m2, m4
    lea        r2, [r2 + 2 * r3]
    movu       [r2], m2

    movq       m2, [r0 + 2 * r1]

    punpcklbw  m1, m2

    pmaddubsw  m3, m6
    pmaddubsw  m1, m5

    paddw      m3, m1
    psubw      m3, m4

    movu       [r2 + r3], m3

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W8_H8_H16_H32 8,  8
FILTER_V_PS_W8_H8_H16_H32 8, 16
FILTER_V_PS_W8_H8_H16_H32 8, 32

;------------------------------------------------------------------------------------------------------------
;void interp_4tap_vert_ps_6x8(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W6 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_6x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m5, [r5 + r4 * 4]
%else
    movd       m5, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m6, m5, [tab_Vm]
    pshufb     m5, [tab_Vm + 16]
    mova       m4, [pw_2000]
    lea        r5, [3 * r1]
    mov        r4d, %2/4

.loop:
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    movq       m2, [r0 + 2 * r1]
    movq       m3, [r0 + r5]

    punpcklbw  m0, m1
    punpcklbw  m1, m2
    punpcklbw  m2, m3

    pmaddubsw  m0, m6
    pmaddubsw  m7, m2, m5

    paddw      m0, m7
    psubw      m0, m4

    movh       [r2], m0
    pshufd     m0, m0, 2
    movd       [r2 + 8], m0

    lea        r0, [r0 + 4 * r1]
    movq       m0, [r0]
    punpcklbw  m3, m0

    pmaddubsw  m1, m6
    pmaddubsw  m7, m3, m5

    paddw      m1, m7
    psubw      m1, m4

    movh       [r2 + r3], m1
    pshufd     m1, m1, 2
    movd       [r2 + r3 + 8], m1

    movq       m1, [r0 + r1]
    punpcklbw  m0, m1

    pmaddubsw  m2, m6
    pmaddubsw  m0, m5

    paddw      m2, m0
    psubw      m2, m4

    lea        r2,[r2 + 2 * r3]
    movh       [r2], m2
    pshufd     m2, m2, 2
    movd       [r2 + 8], m2

    movq       m2,[r0 + 2 * r1]
    punpcklbw  m1, m2

    pmaddubsw  m3, m6
    pmaddubsw  m1, m5

    paddw      m3, m1
    psubw      m3, m4

    movh       [r2 + r3], m3
    pshufd     m3, m3, 2
    movd       [r2 + r3 + 8], m3

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W6 6, 8
FILTER_V_PS_W6 6, 16

;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_12x16(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W12 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_12x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m1, m0, [tab_Vm]
    pshufb     m0, [tab_Vm + 16]

    mov        r4d, %2/2

.loop:
    movu       m2, [r0]
    movu       m3, [r0 + r1]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    lea        r0, [r0 + 2 * r1]
    movu       m5, [r0]
    movu       m7, [r0 + r1]

    punpcklbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m4, m6

    punpckhbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m2, m6

    mova       m6, [pw_2000]

    psubw      m4, m6
    psubw      m2, m6

    movu       [r2], m4
    movh       [r2 + 16], m2

    punpcklbw  m4, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m4, m1
    pmaddubsw  m3, m1

    movu       m2, [r0 + 2 * r1]

    punpcklbw  m5, m7, m2
    punpckhbw  m7, m2

    pmaddubsw  m5, m0
    pmaddubsw  m7, m0

    paddw      m4, m5
    paddw      m3, m7

    psubw      m4, m6
    psubw      m3, m6

    movu       [r2 + r3], m4
    movh       [r2 + r3 + 16], m3

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W12 12, 16
FILTER_V_PS_W12 12, 32

;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_16x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W16 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m1, m0, [tab_Vm]
    pshufb     m0, [tab_Vm + 16]
    mov        r4d, %2/2

.loop:
    movu       m2, [r0]
    movu       m3, [r0 + r1]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    lea        r0, [r0 + 2 * r1]
    movu       m5, [r0]
    movu       m7, [r0 + r1]

    punpcklbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m4, m6

    punpckhbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m2, m6

    mova       m6, [pw_2000]

    psubw      m4, m6
    psubw      m2, m6

    movu       [r2], m4
    movu       [r2 + 16], m2

    punpcklbw  m4, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m4, m1
    pmaddubsw  m3, m1

    movu       m5, [r0 + 2 * r1]

    punpcklbw  m2, m7, m5
    punpckhbw  m7, m5

    pmaddubsw  m2, m0
    pmaddubsw  m7, m0

    paddw      m4, m2
    paddw      m3, m7

    psubw      m4, m6
    psubw      m3, m6

    movu       [r2 + r3], m4
    movu       [r2 + r3 + 16], m3

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W16 16,  4
FILTER_V_PS_W16 16,  8
FILTER_V_PS_W16 16, 12
FILTER_V_PS_W16 16, 16
FILTER_V_PS_W16 16, 32

FILTER_V_PS_W16 16, 24
FILTER_V_PS_W16 16, 64

;--------------------------------------------------------------------------------------------------------------
;void interp_4tap_vert_ps_24x32(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_V4_PS_W24 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_24x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m1, m0, [tab_Vm]
    pshufb     m0, [tab_Vm + 16]

    mov        r4d, %2/2

.loop:
    movu       m2, [r0]
    movu       m3, [r0 + r1]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    lea        r5, [r0 + 2 * r1]

    movu       m5, [r5]
    movu       m7, [r5 + r1]

    punpcklbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m4, m6

    punpckhbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m2, m6

    mova       m6, [pw_2000]

    psubw      m4, m6
    psubw      m2, m6

    movu       [r2], m4
    movu       [r2 + 16], m2

    punpcklbw  m4, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m4, m1
    pmaddubsw  m3, m1

    movu       m2, [r5 + 2 * r1]

    punpcklbw  m5, m7, m2
    punpckhbw  m7, m2

    pmaddubsw  m5, m0
    pmaddubsw  m7, m0

    paddw      m4, m5
    paddw      m3, m7

    psubw      m4, m6
    psubw      m3, m6

    movu       [r2 + r3], m4
    movu       [r2 + r3 + 16], m3

    movq       m2, [r0 + 16]
    movq       m3, [r0 + r1 + 16]
    movq       m4, [r5 + 16]
    movq       m5, [r5 + r1 + 16]

    punpcklbw  m2, m3
    punpcklbw  m7, m4, m5

    pmaddubsw  m2, m1
    pmaddubsw  m7, m0

    paddw      m2, m7
    psubw      m2, m6

    movu       [r2 + 32], m2

    movq       m2, [r5 + 2 * r1 + 16]

    punpcklbw  m3, m4
    punpcklbw  m5, m2

    pmaddubsw  m3, m1
    pmaddubsw  m5, m0

    paddw      m3, m5
    psubw      m3,  m6

    movu       [r2 + r3 + 32], m3

    mov        r0, r5
    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V4_PS_W24 24, 32

FILTER_V4_PS_W24 24, 64

;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_32x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W32 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m1, m0, [tab_Vm]
    pshufb     m0, [tab_Vm + 16]

    mova       m7, [pw_2000]

    mov        r4d, %2

.loop:
    movu       m2, [r0]
    movu       m3, [r0 + r1]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    lea        r5, [r0 + 2 * r1]
    movu       m3, [r5]
    movu       m5, [r5 + r1]

    punpcklbw  m6, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m6, m0
    pmaddubsw  m3, m0

    paddw      m4, m6
    paddw      m2, m3

    psubw      m4, m7
    psubw      m2, m7

    movu       [r2], m4
    movu       [r2 + 16], m2

    movu       m2, [r0 + 16]
    movu       m3, [r0 + r1 + 16]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    movu       m3, [r5 + 16]
    movu       m5, [r5 + r1 + 16]

    punpcklbw  m6, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m6, m0
    pmaddubsw  m3, m0

    paddw      m4, m6
    paddw      m2, m3

    psubw      m4, m7
    psubw      m2, m7

    movu       [r2 + 32], m4
    movu       [r2 + 48], m2

    lea        r0, [r0 + r1]
    lea        r2, [r2 + r3]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W32 32,  8
FILTER_V_PS_W32 32, 16
FILTER_V_PS_W32 32, 24
FILTER_V_PS_W32 32, 32

FILTER_V_PS_W32 32, 48
FILTER_V_PS_W32 32, 64

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W8_H8_H16_H32 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_%1x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m5,        [r5 + r4 * 4]
%else
movd        m5,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m6,        m5,       [tab_Vm]
pshufb      m5,        [tab_Vm + 16]
mova        m4,        [pw_512]
lea         r5,        [r1 * 3]

mov         r4d,       %2

.loop:
movq        m0,        [r0]
movq        m1,        [r0 + r1]
movq        m2,        [r0 + 2 * r1]
movq        m3,        [r0 + r5]

punpcklbw   m0,        m1
punpcklbw   m1,        m2
punpcklbw   m2,        m3

pmaddubsw   m0,        m6
pmaddubsw   m7,        m2, m5

paddw       m0,        m7

pmulhrsw    m0,        m4
packuswb    m0,        m0
movh        [r2],      m0

lea         r0,        [r0 + 4 * r1]
movq        m0,        [r0]

punpcklbw   m3,        m0

pmaddubsw   m1,        m6
pmaddubsw   m7,        m3, m5

paddw       m1,        m7

pmulhrsw    m1,        m4
packuswb    m1,        m1
movh        [r2 + r3], m1

movq        m1,        [r0 + r1]

punpcklbw   m0,        m1

pmaddubsw   m2,        m6
pmaddubsw   m0,        m5

paddw       m2,        m0

pmulhrsw    m2,        m4

movq        m7,        [r0 + 2 * r1]
punpcklbw   m1,        m7

pmaddubsw   m3,        m6
pmaddubsw   m1,        m5

paddw       m3,        m1

pmulhrsw    m3,        m4
packuswb    m2,        m3

lea         r2,        [r2 + 2 * r3]
movh        [r2],      m2
movhps      [r2 + r3], m2

lea         r2,        [r2 + 2 * r3]

sub         r4,         4
jnz        .loop
RET
%endmacro

FILTER_V4_W8_H8_H16_H32 8,  8
FILTER_V4_W8_H8_H16_H32 8, 16
FILTER_V4_W8_H8_H16_H32 8, 32

FILTER_V4_W8_H8_H16_H32 8, 12
FILTER_V4_W8_H8_H16_H32 8, 64

%macro PROCESS_CHROMA_AVX2_W8_8R 0
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m5, m1, xm2, 1                  ; m5 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1                        ; m4 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m2, m3, xm4, 1                  ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3                        ; m1 = [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4                        ; m3 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50]
    vinserti128     m1, m1, xm3, 1                  ; m1 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50] - [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    pmaddubsw       m1, [r5]
    movq            xm3, [r0 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3                        ; m4 = [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 8
    punpcklbw       xm3, xm0                        ; m3 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70]
    vinserti128     m4, m4, xm3, 1                  ; m4 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70] - [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    pmaddubsw       m3, m4, [r5 + 1 * mmsize]
    paddw           m1, m3
    pmaddubsw       m4, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3                        ; m0 = [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    movq            xm6, [r0 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6                        ; m3 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90]
    vinserti128     m0, m0, xm3, 1                  ; m0 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90] - [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    pmaddubsw       m0, [r5 + 1 * mmsize]
    paddw           m4, m0
%endmacro

%macro FILTER_VER_CHROMA_AVX2_8x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x8, 4, 6, 7
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
    PROCESS_CHROMA_AVX2_W8_8R
%ifidn %1,pp
    lea             r4, [r3 * 3]
    mova            m3, [pw_512]
    pmulhrsw        m5, m3                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m3                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m3                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r4], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm1
    movq            [r2 + r3], xm4
    movhps          [r2 + r3 * 2], xm1
    movhps          [r2 + r4], xm4
%else
    add             r3d, r3d
    vbroadcasti128  m3, [pw_2000]
    lea             r4, [r3 * 3]
    psubw           m5, m3                          ; m5 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    psubw           m1, m3                          ; m1 = word: row 4, row 5
    psubw           m4, m3                          ; m4 = word: row 6, row 7
    vextracti128    xm6, m5, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm0, m1, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm6
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm4
    vextracti128    xm4, m4, 1
    movu            [r2 + r4], xm4
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x8 pp
FILTER_VER_CHROMA_AVX2_8x8 ps

%macro FILTER_VER_CHROMA_AVX2_8x6 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x6, 4, 6, 6
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1

    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m5, m1, xm2, 1                  ; m5 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1                        ; m4 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m2, m3, xm4, 1                  ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3                        ; m1 = [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4                        ; m3 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50]
    vinserti128     m1, m1, xm3, 1                  ; m1 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50] - [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    pmaddubsw       m1, [r5]
    movq            xm3, [r0 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3                        ; m4 = [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 8
    punpcklbw       xm3, xm0                        ; m3 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70]
    vinserti128     m4, m4, xm3, 1                  ; m4 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70] - [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    pmaddubsw       m4, [r5 + 1 * mmsize]
    paddw           m1, m4
%ifidn %1,pp
    lea             r4, [r3 * 3]
    mova            m3, [pw_512]
    pmulhrsw        m5, m3                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m3                          ; m1 = word: row 4, row 5
    packuswb        m5, m2
    packuswb        m1, m1
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r4], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm1
    movq            [r2 + r3], xm4
%else
    add             r3d, r3d
    mova            m3, [pw_2000]
    lea             r4, [r3 * 3]
    psubw           m5, m3                          ; m5 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    psubw           m1, m3                          ; m1 = word: row 4, row 5
    vextracti128    xm4, m5, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm0, m1, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm4
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm1
    movu            [r2 + r3], xm0
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x6 pp
FILTER_VER_CHROMA_AVX2_8x6 ps

%macro PROCESS_CHROMA_AVX2_W8_16R 1
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3
    vinserti128     m5, m1, xm2, 1
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1
    vinserti128     m2, m3, xm4, 1
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    pmaddubsw       m1, [r5]
    movq            xm3, [r0 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 8
    punpcklbw       xm3, xm0
    vinserti128     m4, m4, xm3, 1
    pmaddubsw       m3, m4, [r5 + 1 * mmsize]
    paddw           m1, m3
    pmaddubsw       m4, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3
    movq            xm6, [r0 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6
    vinserti128     m0, m0, xm3, 1
    pmaddubsw       m3, m0, [r5 + 1 * mmsize]
    paddw           m4, m3
    pmaddubsw       m0, [r5]
%ifidn %1,pp
    pmulhrsw        m5, m7                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m7                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m7                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m7                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm1
    movq            [r2 + r3], xm4
    movhps          [r2 + r3 * 2], xm1
    movhps          [r2 + r6], xm4
%else
    psubw           m5, m7                          ; m5 = word: row 0, row 1
    psubw           m2, m7                          ; m2 = word: row 2, row 3
    psubw           m1, m7                          ; m1 = word: row 4, row 5
    psubw           m4, m7                          ; m4 = word: row 6, row 7
    vextracti128    xm3, m5, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm3
    vextracti128    xm3, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r2, [r2 + r3 * 4]
    vextracti128    xm5, m1, 1
    vextracti128    xm3, m4, 1
    movu            [r2], xm1
    movu            [r2 + r3], xm5
    movu            [r2 + r3 * 2], xm4
    movu            [r2 + r6], xm3
%endif
    movq            xm3, [r0 + r4]                  ; m3 = row 11
    punpcklbw       xm6, xm3
    lea             r0, [r0 + r1 * 4]
    movq            xm5, [r0]                       ; m5 = row 12
    punpcklbw       xm3, xm5
    vinserti128     m6, m6, xm3, 1
    pmaddubsw       m3, m6, [r5 + 1 * mmsize]
    paddw           m0, m3
    pmaddubsw       m6, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 13
    punpcklbw       xm5, xm3
    movq            xm2, [r0 + r1 * 2]              ; m2 = row 14
    punpcklbw       xm3, xm2
    vinserti128     m5, m5, xm3, 1
    pmaddubsw       m3, m5, [r5 + 1 * mmsize]
    paddw           m6, m3
    pmaddubsw       m5, [r5]
    movq            xm3, [r0 + r4]                  ; m3 = row 15
    punpcklbw       xm2, xm3
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 16
    punpcklbw       xm3, xm1
    vinserti128     m2, m2, xm3, 1
    pmaddubsw       m3, m2, [r5 + 1 * mmsize]
    paddw           m5, m3
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 17
    punpcklbw       xm1, xm3
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpcklbw       xm3, xm4
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5 + 1 * mmsize]
    paddw           m2, m1
    lea             r2, [r2 + r3 * 4]
%ifidn %1,pp
    pmulhrsw        m0, m7                          ; m0 = word: row 8, row 9
    pmulhrsw        m6, m7                          ; m6 = word: row 10, row 11
    pmulhrsw        m5, m7                          ; m5 = word: row 12, row 13
    pmulhrsw        m2, m7                          ; m2 = word: row 14, row 15
    packuswb        m0, m6
    packuswb        m5, m2
    vextracti128    xm6, m0, 1
    vextracti128    xm2, m5, 1
    movq            [r2], xm0
    movq            [r2 + r3], xm6
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm6
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm2
%else
    psubw           m0, m7                          ; m0 = word: row 8, row 9
    psubw           m6, m7                          ; m6 = word: row 10, row 11
    psubw           m5, m7                          ; m5 = word: row 12, row 13
    psubw           m2, m7                          ; m2 = word: row 14, row 15
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m6, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm3
    lea             r2, [r2 + r3 * 4]
    vextracti128    xm1, m5, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif
%endmacro

%macro FILTER_VER_CHROMA_AVX2_8x16 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x16, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    mova            m7, [pw_2000]
%endif
    lea             r6, [r3 * 3]
    PROCESS_CHROMA_AVX2_W8_16R %1
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x16 pp
FILTER_VER_CHROMA_AVX2_8x16 ps

%macro FILTER_VER_CHROMA_AVX2_8x32 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x32, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    mova            m7, [pw_2000]
%endif
    lea             r6, [r3 * 3]
%rep 2
    PROCESS_CHROMA_AVX2_W8_16R %1
    lea             r2, [r2 + r3 * 4]
%endrep
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x32 pp
FILTER_VER_CHROMA_AVX2_8x32 ps

%macro PROCESS_CHROMA_AVX2_W8_4R 0
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m0, m1, xm2, 1                  ; m0 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m0, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1                        ; m4 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m2, m3, xm4, 1                  ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3                        ; m1 = [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4                        ; m3 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50]
    vinserti128     m1, m1, xm3, 1                  ; m1 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50] - [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    pmaddubsw       m1, [r5 + 1 * mmsize]
    paddw           m2, m1
%endmacro

%macro FILTER_VER_CHROMA_AVX2_8x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x4, 4, 6, 5
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
    PROCESS_CHROMA_AVX2_W8_4R
%ifidn %1,pp
    lea             r4, [r3 * 3]
    mova            m3, [pw_512]
    pmulhrsw        m0, m3                          ; m0 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    packuswb        m0, m2
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r4], xm2
%else
    add             r3d, r3d
    vbroadcasti128  m3, [pw_2000]
    lea             r4, [r3 * 3]
    psubw           m0, m3                          ; m0 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    vextracti128    xm1, m0, 1
    vextracti128    xm4, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm4
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x4 pp
FILTER_VER_CHROMA_AVX2_8x4 ps

%macro FILTER_VER_CHROMA_AVX2_8x2 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x2, 4, 6, 4
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1

    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m1, m1, xm2, 1                  ; m1 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m1, [r5]
    movq            xm2, [r0 + r4]                  ; m2 = row 3
    punpcklbw       xm3, xm2                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    movq            xm0, [r0 + r1 * 4]              ; m0 = row 4
    punpcklbw       xm2, xm0                        ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m3, m3, xm2, 1                  ; m3 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m1, m3
%ifidn %1,pp
    pmulhrsw        m1, [pw_512]                    ; m1 = word: row 0, row 1
    packuswb        m1, m1
    vextracti128    xm0, m1, 1
    movq            [r2], xm1
    movq            [r2 + r3], xm0
%else
    add             r3d, r3d
    psubw           m1, [pw_2000]                   ; m1 = word: row 0, row 1
    vextracti128    xm0, m1, 1
    movu            [r2], xm1
    movu            [r2 + r3], xm0
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_8x2 pp
FILTER_VER_CHROMA_AVX2_8x2 ps

%macro FILTER_VER_CHROMA_AVX2_6x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_6x8, 4, 6, 7
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
    PROCESS_CHROMA_AVX2_W8_8R
%ifidn %1,pp
    lea             r4, [r3 * 3]
    mova            m3, [pw_512]
    pmulhrsw        m5, m3                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m3                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m3                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movd            [r2], xm5
    pextrw          [r2 + 4], xm5, 2
    movd            [r2 + r3], xm2
    pextrw          [r2 + r3 + 4], xm2, 2
    pextrd          [r2 + r3 * 2], xm5, 2
    pextrw          [r2 + r3 * 2 + 4], xm5, 6
    pextrd          [r2 + r4], xm2, 2
    pextrw          [r2 + r4 + 4], xm2, 6
    lea             r2, [r2 + r3 * 4]
    movd            [r2], xm1
    pextrw          [r2 + 4], xm1, 2
    movd            [r2 + r3], xm4
    pextrw          [r2 + r3 + 4], xm4, 2
    pextrd          [r2 + r3 * 2], xm1, 2
    pextrw          [r2 + r3 * 2 + 4], xm1, 6
    pextrd          [r2 + r4], xm4, 2
    pextrw          [r2 + r4 + 4], xm4, 6
%else
    add             r3d, r3d
    vbroadcasti128  m3, [pw_2000]
    lea             r4, [r3 * 3]
    psubw           m5, m3                          ; m5 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    psubw           m1, m3                          ; m1 = word: row 4, row 5
    psubw           m4, m3                          ; m4 = word: row 6, row 7
    vextracti128    xm6, m5, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm0, m1, 1
    movq            [r2], xm5
    pextrd          [r2 + 8], xm5, 2
    movq            [r2 + r3], xm6
    pextrd          [r2 + r3 + 8], xm6, 2
    movq            [r2 + r3 * 2], xm2
    pextrd          [r2 + r3 * 2 + 8], xm2, 2
    movq            [r2 + r4], xm3
    pextrd          [r2 + r4 + 8], xm3, 2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm1
    pextrd          [r2 + 8], xm1, 2
    movq            [r2 + r3], xm0
    pextrd          [r2 + r3 + 8], xm0, 2
    movq            [r2 + r3 * 2], xm4
    pextrd          [r2 + r3 * 2 + 8], xm4, 2
    vextracti128    xm4, m4, 1
    movq            [r2 + r4], xm4
    pextrd          [r2 + r4 + 8], xm4, 2
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_6x8 pp
FILTER_VER_CHROMA_AVX2_6x8 ps

;-----------------------------------------------------------------------------
;void interp_4tap_vert_pp_6x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W6_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_6x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m5,        [r5 + r4 * 4]
%else
movd        m5,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m6,        m5,       [tab_Vm]
pshufb      m5,        [tab_Vm + 16]
mova        m4,        [pw_512]

mov         r4d,       %2
lea         r5,        [3 * r1]

.loop:
movq        m0,        [r0]
movq        m1,        [r0 + r1]
movq        m2,        [r0 + 2 * r1]
movq        m3,        [r0 + r5]

punpcklbw   m0,        m1
punpcklbw   m1,        m2
punpcklbw   m2,        m3

pmaddubsw   m0,        m6
pmaddubsw   m7,        m2, m5

paddw       m0,        m7

pmulhrsw    m0,        m4
packuswb    m0,        m0
movd        [r2],      m0
pextrw      [r2 + 4],  m0,    2

lea         r0,        [r0 + 4 * r1]

movq        m0,        [r0]
punpcklbw   m3,        m0

pmaddubsw   m1,        m6
pmaddubsw   m7,        m3, m5

paddw       m1,        m7

pmulhrsw    m1,        m4
packuswb    m1,        m1
movd        [r2 + r3],      m1
pextrw      [r2 + r3 + 4],  m1,    2

movq        m1,        [r0 + r1]
punpcklbw   m7,        m0,        m1

pmaddubsw   m2,        m6
pmaddubsw   m7,        m5

paddw       m2,        m7

pmulhrsw    m2,        m4
packuswb    m2,        m2
lea         r2,        [r2 + 2 * r3]
movd        [r2],      m2
pextrw      [r2 + 4],  m2,    2

movq        m2,        [r0 + 2 * r1]
punpcklbw   m1,        m2

pmaddubsw   m3,        m6
pmaddubsw   m1,        m5

paddw       m3,        m1

pmulhrsw    m3,        m4
packuswb    m3,        m3

movd        [r2 + r3],        m3
pextrw      [r2 + r3 + 4],    m3,    2

lea         r2,        [r2 + 2 * r3]

sub         r4,         4
jnz        .loop
RET
%endmacro

FILTER_V4_W6_H4 6, 8

FILTER_V4_W6_H4 6, 16

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_12x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W12_H2 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_12x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m1,        m0,       [tab_Vm]
pshufb      m0,        [tab_Vm + 16]

mov         r4d,       %2

.loop:
movu        m2,        [r0]
movu        m3,        [r0 + r1]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

lea         r0,        [r0 + 2 * r1]
movu        m5,        [r0]
movu        m7,        [r0 + r1]

punpcklbw   m6,        m5,        m7
pmaddubsw   m6,        m0
paddw       m4,        m6

punpckhbw   m6,        m5,        m7
pmaddubsw   m6,        m0
paddw       m2,        m6

mova        m6,        [pw_512]

pmulhrsw    m4,        m6
pmulhrsw    m2,        m6

packuswb    m4,        m2

movh         [r2],     m4
pextrd       [r2 + 8], m4,  2

punpcklbw   m4,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m4,        m1
pmaddubsw   m3,        m1

movu        m5,        [r0 + 2 * r1]

punpcklbw   m2,        m7,        m5
punpckhbw   m7,        m5

pmaddubsw   m2,        m0
pmaddubsw   m7,        m0

paddw       m4,        m2
paddw       m3,        m7

pmulhrsw    m4,        m6
pmulhrsw    m3,        m6

packuswb    m4,        m3

movh        [r2 + r3],      m4
pextrd      [r2 + r3 + 8],  m4,  2

lea         r2,        [r2 + 2 * r3]

sub         r4,        2
jnz        .loop
RET
%endmacro

FILTER_V4_W12_H2 12, 16

FILTER_V4_W12_H2 12, 32

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_16x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W16_H2 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_16x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m1,        m0,       [tab_Vm]
pshufb      m0,        [tab_Vm + 16]

mov         r4d,       %2/2

.loop:
movu        m2,        [r0]
movu        m3,        [r0 + r1]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

lea         r0,        [r0 + 2 * r1]
movu        m5,        [r0]
movu        m6,        [r0 + r1]

punpckhbw   m7,        m5,        m6
pmaddubsw   m7,        m0
paddw       m2,        m7

punpcklbw   m7,        m5,        m6
pmaddubsw   m7,        m0
paddw       m4,        m7

mova        m7,        [pw_512]

pmulhrsw    m4,        m7
pmulhrsw    m2,        m7

packuswb    m4,        m2

movu        [r2],      m4

punpcklbw   m4,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m4,        m1
pmaddubsw   m3,        m1

movu        m5,        [r0 + 2 * r1]

punpcklbw   m2,        m6,        m5
punpckhbw   m6,        m5

pmaddubsw   m2,        m0
pmaddubsw   m6,        m0

paddw       m4,        m2
paddw       m3,        m6

pmulhrsw    m4,        m7
pmulhrsw    m3,        m7

packuswb    m4,        m3

movu        [r2 + r3],      m4

lea         r2,        [r2 + 2 * r3]

dec         r4d
jnz        .loop
RET
%endmacro

FILTER_V4_W16_H2 16,  4
FILTER_V4_W16_H2 16,  8
FILTER_V4_W16_H2 16, 12
FILTER_V4_W16_H2 16, 16
FILTER_V4_W16_H2 16, 32

FILTER_V4_W16_H2 16, 24
FILTER_V4_W16_H2 16, 64

%macro FILTER_VER_CHROMA_AVX2_16x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_16x16, 4, 6, 15
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    mova            m12, [r5]
    mova            m13, [r5 + mmsize]
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%endif
    lea             r5, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, m12
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, m12
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, m13
    paddw           m0, m4
    pmaddubsw       m2, m12
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, m13
    paddw           m1, m5
    pmaddubsw       m3, m12
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, m13
    paddw           m2, m6
    pmaddubsw       m4, m12
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, m13
    paddw           m3, m7
    pmaddubsw       m5, m12
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, m13
    paddw           m4, m8
    pmaddubsw       m6, m12
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, m13
    paddw           m5, m9
    pmaddubsw       m7, m12
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, m13
    paddw           m6, m10
    pmaddubsw       m8, m12
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, m13
    paddw           m7, m11
    pmaddubsw       m9, m12

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    packuswb        m6, m7
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r5], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r5], xm7
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r5], m3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m4
    movu            [r2 + r3], m5
    movu            [r2 + r3 * 2], m6
    movu            [r2 + r5], m7
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm6, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm6, 1
    pmaddubsw       m6, m10, m13
    paddw           m8, m6
    pmaddubsw       m10, m12
    lea             r0, [r0 + r1 * 4]
    movu            xm6, [r0]                       ; m6 = row 12
    punpckhbw       xm7, xm11, xm6
    punpcklbw       xm11, xm6
    vinserti128     m11, m11, xm7, 1
    pmaddubsw       m7, m11, m13
    paddw           m9, m7
    pmaddubsw       m11, m12

    movu            xm7, [r0 + r1]                  ; m7 = row 13
    punpckhbw       xm0, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm0, 1
    pmaddubsw       m0, m6, m13
    paddw           m10, m0
    pmaddubsw       m6, m12
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm7, xm0
    punpcklbw       xm7, xm0
    vinserti128     m7, m7, xm1, 1
    pmaddubsw       m1, m7, m13
    paddw           m11, m1
    pmaddubsw       m7, m12
    movu            xm1, [r0 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, m13
    paddw           m6, m2
    pmaddubsw       m0, m12
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, m13
    paddw           m7, m3
    pmaddubsw       m1, m12
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m2, m13
    paddw           m0, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m3, m13
    paddw           m1, m3

%ifidn %1,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    pmulhrsw        m6, m14                         ; m6 = word: row 12
    pmulhrsw        m7, m14                         ; m7 = word: row 13
    pmulhrsw        m0, m14                         ; m0 = word: row 14
    pmulhrsw        m1, m14                         ; m1 = word: row 15
    packuswb        m8, m9
    packuswb        m10, m11
    packuswb        m6, m7
    packuswb        m0, m1
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m6, m6, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm7, m6, 1
    vextracti128    xm1, m0, 1
    movu            [r2], xm8
    movu            [r2 + r3], xm9
    movu            [r2 + r3 * 2], xm10
    movu            [r2 + r5], xm11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm6
    movu            [r2 + r3], xm7
    movu            [r2 + r3 * 2], xm0
    movu            [r2 + r5], xm1
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    psubw           m6, m14                         ; m6 = word: row 12
    psubw           m7, m14                         ; m7 = word: row 13
    psubw           m0, m14                         ; m0 = word: row 14
    psubw           m1, m14                         ; m1 = word: row 15
    movu            [r2], m8
    movu            [r2 + r3], m9
    movu            [r2 + r3 * 2], m10
    movu            [r2 + r5], m11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m6
    movu            [r2 + r3], m7
    movu            [r2 + r3 * 2], m0
    movu            [r2 + r5], m1
%endif
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_16x16 pp
FILTER_VER_CHROMA_AVX2_16x16 ps
%macro FILTER_VER_CHROMA_AVX2_16x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_16x8, 4, 7, 7
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m6, [pw_512]
%else
    add             r3d, r3d
    mova            m6, [pw_2000]
%endif
    lea             r6, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
%ifidn %1,pp
    pmulhrsw        m0, m6                          ; m0 = word: row 0
    pmulhrsw        m1, m6                          ; m1 = word: row 1
    packuswb        m0, m1
    vpermq          m0, m0, 11011000b
    vextracti128    xm1, m0, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
%else
    psubw           m0, m6                          ; m0 = word: row 0
    psubw           m1, m6                          ; m1 = word: row 1
    movu            [r2], m0
    movu            [r2 + r3], m1
%endif

    movu            xm0, [r0 + r1]                  ; m0 = row 5
    punpckhbw       xm1, xm4, xm0
    punpcklbw       xm4, xm0
    vinserti128     m4, m4, xm1, 1
    pmaddubsw       m1, m4, [r5 + mmsize]
    paddw           m2, m1
    pmaddubsw       m4, [r5]
    movu            xm1, [r0 + r1 * 2]              ; m1 = row 6
    punpckhbw       xm5, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm5, 1
    pmaddubsw       m5, m0, [r5 + mmsize]
    paddw           m3, m5
    pmaddubsw       m0, [r5]
%ifidn %1,pp
    pmulhrsw        m2, m6                          ; m2 = word: row 2
    pmulhrsw        m3, m6                          ; m3 = word: row 3
    packuswb        m2, m3
    vpermq          m2, m2, 11011000b
    vextracti128    xm3, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%else
    psubw           m2, m6                          ; m2 = word: row 2
    psubw           m3, m6                          ; m3 = word: row 3
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
%endif

    movu            xm2, [r0 + r4]                  ; m2 = row 7
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + mmsize]
    paddw           m4, m3
    pmaddubsw       m1, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm3, [r0]                       ; m3 = row 8
    punpckhbw       xm5, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm5, 1
    pmaddubsw       m5, m2, [r5 + mmsize]
    paddw           m0, m5
    pmaddubsw       m2, [r5]
    lea             r2, [r2 + r3 * 4]
%ifidn %1,pp
    pmulhrsw        m4, m6                          ; m4 = word: row 4
    pmulhrsw        m0, m6                          ; m0 = word: row 5
    packuswb        m4, m0
    vpermq          m4, m4, 11011000b
    vextracti128    xm0, m4, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm0
%else
    psubw           m4, m6                          ; m4 = word: row 4
    psubw           m0, m6                          ; m0 = word: row 5
    movu            [r2], m4
    movu            [r2 + r3], m0
%endif

    movu            xm5, [r0 + r1]                  ; m5 = row 9
    punpckhbw       xm4, xm3, xm5
    punpcklbw       xm3, xm5
    vinserti128     m3, m3, xm4, 1
    pmaddubsw       m3, [r5 + mmsize]
    paddw           m1, m3
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 10
    punpckhbw       xm0, xm5, xm4
    punpcklbw       xm5, xm4
    vinserti128     m5, m5, xm0, 1
    pmaddubsw       m5, [r5 + mmsize]
    paddw           m2, m5
%ifidn %1,pp
    pmulhrsw        m1, m6                          ; m1 = word: row 6
    pmulhrsw        m2, m6                          ; m2 = word: row 7
    packuswb        m1, m2
    vpermq          m1, m1, 11011000b
    vextracti128    xm2, m1, 1
    movu            [r2 + r3 * 2], xm1
    movu            [r2 + r6], xm2
%else
    psubw           m1, m6                          ; m1 = word: row 6
    psubw           m2, m6                          ; m2 = word: row 7
    movu            [r2 + r3 * 2], m1
    movu            [r2 + r6], m2
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_16x8 pp
FILTER_VER_CHROMA_AVX2_16x8 ps

%macro FILTER_VER_CHROMA_AVX2_16x12 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_16x12, 4, 6, 10
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    mova            m8, [r5]
    mova            m9, [r5 + mmsize]
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m7, [pw_2000]
%endif
    lea             r5, [r3 * 3]

    movu            xm0, [r0]
    vinserti128     m0, m0, [r0 + r1 * 2], 1
    movu            xm1, [r0 + r1]
    vinserti128     m1, m1, [r0 + r4], 1

    punpcklbw       m2, m0, m1
    punpckhbw       m3, m0, m1
    vperm2i128      m4, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    pmaddubsw       m4, m8
    pmaddubsw       m3, m2, m9
    paddw           m4, m3
    pmaddubsw       m2, m8

    vextracti128    xm0, m0, 1
    lea             r0, [r0 + r1 * 4]
    vinserti128     m0, m0, [r0], 1

    punpcklbw       m5, m1, m0
    punpckhbw       m3, m1, m0
    vperm2i128      m6, m5, m3, 0x20
    vperm2i128      m5, m5, m3, 0x31
    pmaddubsw       m6, m8
    pmaddubsw       m3, m5, m9
    paddw           m6, m3
    pmaddubsw       m5, m8
%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 0
    pmulhrsw        m6, m7                         ; m6 = word: row 1
    packuswb        m4, m6
    vpermq          m4, m4, 11011000b
    vextracti128    xm6, m4, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm6
%else
    psubw           m4, m7                         ; m4 = word: row 0
    psubw           m6, m7                         ; m6 = word: row 1
    movu            [r2], m4
    movu            [r2 + r3], m6
%endif

    movu            xm4, [r0 + r1 * 2]
    vinserti128     m4, m4, [r0 + r1], 1
    vextracti128    xm1, m4, 1
    vinserti128     m0, m0, xm1, 0

    punpcklbw       m6, m0, m4
    punpckhbw       m1, m0, m4
    vperm2i128      m0, m6, m1, 0x20
    vperm2i128      m6, m6, m1, 0x31
    pmaddubsw       m1, m0, m9
    paddw           m5, m1
    pmaddubsw       m0, m8
    pmaddubsw       m1, m6, m9
    paddw           m2, m1
    pmaddubsw       m6, m8

%ifidn %1,pp
    pmulhrsw        m2, m7                         ; m2 = word: row 2
    pmulhrsw        m5, m7                         ; m5 = word: row 3
    packuswb        m2, m5
    vpermq          m2, m2, 11011000b
    vextracti128    xm5, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r5], xm5
%else
    psubw           m2, m7                         ; m2 = word: row 2
    psubw           m5, m7                         ; m5 = word: row 3
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r5], m5
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]
    lea             r0, [r0 + r1 * 4]
    vinserti128     m1, m1, [r0], 1
    vinserti128     m4, m4, xm1, 1

    punpcklbw       m2, m4, m1
    punpckhbw       m5, m4, m1
    vperm2i128      m3, m2, m5, 0x20
    vperm2i128      m2, m2, m5, 0x31
    pmaddubsw       m5, m3, m9
    paddw           m6, m5
    pmaddubsw       m3, m8
    pmaddubsw       m5, m2, m9
    paddw           m0, m5
    pmaddubsw       m2, m8

%ifidn %1,pp
    pmulhrsw        m6, m7                         ; m6 = word: row 4
    pmulhrsw        m0, m7                         ; m0 = word: row 5
    packuswb        m6, m0
    vpermq          m6, m6, 11011000b
    vextracti128    xm0, m6, 1
    movu            [r2], xm6
    movu            [r2 + r3], xm0
%else
    psubw           m6, m7                         ; m6 = word: row 4
    psubw           m0, m7                         ; m0 = word: row 5
    movu            [r2], m6
    movu            [r2 + r3], m0
%endif

    movu            xm6, [r0 + r1 * 2]
    vinserti128     m6, m6, [r0 + r1], 1
    vextracti128    xm0, m6, 1
    vinserti128     m1, m1, xm0, 0

    punpcklbw       m4, m1, m6
    punpckhbw       m5, m1, m6
    vperm2i128      m0, m4, m5, 0x20
    vperm2i128      m5, m4, m5, 0x31
    pmaddubsw       m4, m0, m9
    paddw           m2, m4
    pmaddubsw       m0, m8
    pmaddubsw       m4, m5, m9
    paddw           m3, m4
    pmaddubsw       m5, m8

%ifidn %1,pp
    pmulhrsw        m3, m7                         ; m3 = word: row 6
    pmulhrsw        m2, m7                         ; m2 = word: row 7
    packuswb        m3, m2
    vpermq          m3, m3, 11011000b
    vextracti128    xm2, m3, 1
    movu            [r2 + r3 * 2], xm3
    movu            [r2 + r5], xm2
%else
    psubw           m3, m7                         ; m3 = word: row 6
    psubw           m2, m7                         ; m2 = word: row 7
    movu            [r2 + r3 * 2], m3
    movu            [r2 + r5], m2
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm3, [r0 + r4]
    lea             r0, [r0 + r1 * 4]
    vinserti128     m3, m3, [r0], 1
    vinserti128     m6, m6, xm3, 1

    punpcklbw       m2, m6, m3
    punpckhbw       m1, m6, m3
    vperm2i128      m4, m2, m1, 0x20
    vperm2i128      m2, m2, m1, 0x31
    pmaddubsw       m1, m4, m9
    paddw           m5, m1
    pmaddubsw       m4, m8
    pmaddubsw       m1, m2, m9
    paddw           m0, m1
    pmaddubsw       m2, m8

%ifidn %1,pp
    pmulhrsw        m5, m7                         ; m5 = word: row 8
    pmulhrsw        m0, m7                         ; m0 = word: row 9
    packuswb        m5, m0
    vpermq          m5, m5, 11011000b
    vextracti128    xm0, m5, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm0
%else
    psubw           m5, m7                         ; m5 = word: row 8
    psubw           m0, m7                         ; m0 = word: row 9
    movu            [r2], m5
    movu            [r2 + r3], m0
%endif

    movu            xm5, [r0 + r1 * 2]
    vinserti128     m5, m5, [r0 + r1], 1
    vextracti128    xm0, m5, 1
    vinserti128     m3, m3, xm0, 0

    punpcklbw       m1, m3, m5
    punpckhbw       m0, m3, m5
    vperm2i128      m6, m1, m0, 0x20
    vperm2i128      m0, m1, m0, 0x31
    pmaddubsw       m1, m6, m9
    paddw           m2, m1
    pmaddubsw       m1, m0, m9
    paddw           m4, m1

%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 10
    pmulhrsw        m2, m7                         ; m2 = word: row 11
    packuswb        m4, m2
    vpermq          m4, m4, 11011000b
    vextracti128    xm2, m4, 1
    movu            [r2 + r3 * 2], xm4
    movu            [r2 + r5], xm2
%else
    psubw           m4, m7                         ; m4 = word: row 10
    psubw           m2, m7                         ; m2 = word: row 11
    movu            [r2 + r3 * 2], m4
    movu            [r2 + r5], m2
%endif
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_16x12 pp
FILTER_VER_CHROMA_AVX2_16x12 ps

%macro FILTER_VER_CHROMA_AVX2_16x32 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_16x32, 4, 8, 8
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    mova            m7, [pw_2000]
%endif
    lea             r6, [r3 * 3]
    mov             r7d, 2
.loopH:
    movu            xm0, [r0]
    vinserti128     m0, m0, [r0 + r1 * 2], 1
    movu            xm1, [r0 + r1]
    vinserti128     m1, m1, [r0 + r4], 1

    punpcklbw       m2, m0, m1
    punpckhbw       m3, m0, m1
    vperm2i128      m4, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    pmaddubsw       m4, [r5]
    pmaddubsw       m3, m2, [r5 + mmsize]
    paddw           m4, m3
    pmaddubsw       m2, [r5]

    vextracti128    xm0, m0, 1
    lea             r0, [r0 + r1 * 4]
    vinserti128     m0, m0, [r0], 1

    punpcklbw       m5, m1, m0
    punpckhbw       m3, m1, m0
    vperm2i128      m6, m5, m3, 0x20
    vperm2i128      m5, m5, m3, 0x31
    pmaddubsw       m6, [r5]
    pmaddubsw       m3, m5, [r5 + mmsize]
    paddw           m6, m3
    pmaddubsw       m5, [r5]
%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 0
    pmulhrsw        m6, m7                         ; m6 = word: row 1
    packuswb        m4, m6
    vpermq          m4, m4, 11011000b
    vextracti128    xm6, m4, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm6
%else
    psubw           m4, m7                         ; m4 = word: row 0
    psubw           m6, m7                         ; m6 = word: row 1
    movu            [r2], m4
    movu            [r2 + r3], m6
%endif

    movu            xm4, [r0 + r1 * 2]
    vinserti128     m4, m4, [r0 + r1], 1
    vextracti128    xm1, m4, 1
    vinserti128     m0, m0, xm1, 0

    punpcklbw       m6, m0, m4
    punpckhbw       m1, m0, m4
    vperm2i128      m0, m6, m1, 0x20
    vperm2i128      m6, m6, m1, 0x31
    pmaddubsw       m1, m0, [r5 + mmsize]
    paddw           m5, m1
    pmaddubsw       m0, [r5]
    pmaddubsw       m1, m6, [r5 + mmsize]
    paddw           m2, m1
    pmaddubsw       m6, [r5]

%ifidn %1,pp
    pmulhrsw        m2, m7                         ; m2 = word: row 2
    pmulhrsw        m5, m7                         ; m5 = word: row 3
    packuswb        m2, m5
    vpermq          m2, m2, 11011000b
    vextracti128    xm5, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm5
%else
    psubw           m2, m7                         ; m2 = word: row 2
    psubw           m5, m7                         ; m5 = word: row 3
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m5
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]
    lea             r0, [r0 + r1 * 4]
    vinserti128     m1, m1, [r0], 1
    vinserti128     m4, m4, xm1, 1

    punpcklbw       m2, m4, m1
    punpckhbw       m5, m4, m1
    vperm2i128      m3, m2, m5, 0x20
    vperm2i128      m2, m2, m5, 0x31
    pmaddubsw       m5, m3, [r5 + mmsize]
    paddw           m6, m5
    pmaddubsw       m3, [r5]
    pmaddubsw       m5, m2, [r5 + mmsize]
    paddw           m0, m5
    pmaddubsw       m2, [r5]

%ifidn %1,pp
    pmulhrsw        m6, m7                         ; m6 = word: row 4
    pmulhrsw        m0, m7                         ; m0 = word: row 5
    packuswb        m6, m0
    vpermq          m6, m6, 11011000b
    vextracti128    xm0, m6, 1
    movu            [r2], xm6
    movu            [r2 + r3], xm0
%else
    psubw           m6, m7                         ; m6 = word: row 4
    psubw           m0, m7                         ; m0 = word: row 5
    movu            [r2], m6
    movu            [r2 + r3], m0
%endif

    movu            xm6, [r0 + r1 * 2]
    vinserti128     m6, m6, [r0 + r1], 1
    vextracti128    xm0, m6, 1
    vinserti128     m1, m1, xm0, 0

    punpcklbw       m4, m1, m6
    punpckhbw       m5, m1, m6
    vperm2i128      m0, m4, m5, 0x20
    vperm2i128      m5, m4, m5, 0x31
    pmaddubsw       m4, m0, [r5 + mmsize]
    paddw           m2, m4
    pmaddubsw       m0, [r5]
    pmaddubsw       m4, m5, [r5 + mmsize]
    paddw           m3, m4
    pmaddubsw       m5, [r5]

%ifidn %1,pp
    pmulhrsw        m3, m7                         ; m3 = word: row 6
    pmulhrsw        m2, m7                         ; m2 = word: row 7
    packuswb        m3, m2
    vpermq          m3, m3, 11011000b
    vextracti128    xm2, m3, 1
    movu            [r2 + r3 * 2], xm3
    movu            [r2 + r6], xm2
%else
    psubw           m3, m7                         ; m3 = word: row 6
    psubw           m2, m7                         ; m2 = word: row 7
    movu            [r2 + r3 * 2], m3
    movu            [r2 + r6], m2
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm3, [r0 + r4]
    lea             r0, [r0 + r1 * 4]
    vinserti128     m3, m3, [r0], 1
    vinserti128     m6, m6, xm3, 1

    punpcklbw       m2, m6, m3
    punpckhbw       m1, m6, m3
    vperm2i128      m4, m2, m1, 0x20
    vperm2i128      m2, m2, m1, 0x31
    pmaddubsw       m1, m4, [r5 + mmsize]
    paddw           m5, m1
    pmaddubsw       m4, [r5]
    pmaddubsw       m1, m2, [r5 + mmsize]
    paddw           m0, m1
    pmaddubsw       m2, [r5]

%ifidn %1,pp
    pmulhrsw        m5, m7                         ; m5 = word: row 8
    pmulhrsw        m0, m7                         ; m0 = word: row 9
    packuswb        m5, m0
    vpermq          m5, m5, 11011000b
    vextracti128    xm0, m5, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm0
%else
    psubw           m5, m7                         ; m5 = word: row 8
    psubw           m0, m7                         ; m0 = word: row 9
    movu            [r2], m5
    movu            [r2 + r3], m0
%endif

    movu            xm5, [r0 + r1 * 2]
    vinserti128     m5, m5, [r0 + r1], 1
    vextracti128    xm0, m5, 1
    vinserti128     m3, m3, xm0, 0

    punpcklbw       m1, m3, m5
    punpckhbw       m0, m3, m5
    vperm2i128      m6, m1, m0, 0x20
    vperm2i128      m0, m1, m0, 0x31
    pmaddubsw       m1, m6, [r5 + mmsize]
    paddw           m2, m1
    pmaddubsw       m6, [r5]
    pmaddubsw       m1, m0, [r5 + mmsize]
    paddw           m4, m1
    pmaddubsw       m0, [r5]

%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 10
    pmulhrsw        m2, m7                         ; m2 = word: row 11
    packuswb        m4, m2
    vpermq          m4, m4, 11011000b
    vextracti128    xm2, m4, 1
    movu            [r2 + r3 * 2], xm4
    movu            [r2 + r6], xm2
%else
    psubw           m4, m7                         ; m4 = word: row 10
    psubw           m2, m7                         ; m2 = word: row 11
    movu            [r2 + r3 * 2], m4
    movu            [r2 + r6], m2
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm3, [r0 + r4]
    lea             r0, [r0 + r1 * 4]
    vinserti128     m3, m3, [r0], 1
    vinserti128     m5, m5, xm3, 1

    punpcklbw       m2, m5, m3
    punpckhbw       m1, m5, m3
    vperm2i128      m4, m2, m1, 0x20
    vperm2i128      m2, m2, m1, 0x31
    pmaddubsw       m1, m4, [r5 + mmsize]
    paddw           m0, m1
    pmaddubsw       m4, [r5]
    pmaddubsw       m1, m2, [r5 + mmsize]
    paddw           m6, m1
    pmaddubsw       m2, [r5]

%ifidn %1,pp
    pmulhrsw        m0, m7                         ; m0 = word: row 12
    pmulhrsw        m6, m7                         ; m6 = word: row 13
    packuswb        m0, m6
    vpermq          m0, m0, 11011000b
    vextracti128    xm6, m0, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm6
%else
    psubw           m0, m7                         ; m0 = word: row 12
    psubw           m6, m7                         ; m6 = word: row 13
    movu            [r2], m0
    movu            [r2 + r3], m6
%endif

    movu            xm5, [r0 + r1 * 2]
    vinserti128     m5, m5, [r0 + r1], 1
    vextracti128    xm0, m5, 1
    vinserti128     m3, m3, xm0, 0

    punpcklbw       m1, m3, m5
    punpckhbw       m0, m3, m5
    vperm2i128      m6, m1, m0, 0x20
    vperm2i128      m0, m1, m0, 0x31
    pmaddubsw       m6, [r5 + mmsize]
    paddw           m2, m6
    pmaddubsw       m0, [r5 + mmsize]
    paddw           m4, m0

%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 14
    pmulhrsw        m2, m7                         ; m2 = word: row 15
    packuswb        m4, m2
    vpermq          m4, m4, 11011000b
    vextracti128    xm2, m4, 1
    movu            [r2 + r3 * 2], xm4
    movu            [r2 + r6], xm2
%else
    psubw           m4, m7                         ; m4 = word: row 14
    psubw           m2, m7                         ; m2 = word: row 15
    movu            [r2 + r3 * 2], m4
    movu            [r2 + r6], m2
%endif
    lea             r2, [r2 + r3 * 4]
    dec             r7d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_16x32 pp
FILTER_VER_CHROMA_AVX2_16x32 ps

%macro FILTER_VER_CHROMA_AVX2_24x32 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_24x32, 4, 9, 10
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    mova            m8, [r5]
    mova            m9, [r5 + mmsize]
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m7, [pw_2000]
%endif
    lea             r6, [r3 * 3]
    mov             r5d, 2
.loopH:
    movu            xm0, [r0]
    vinserti128     m0, m0, [r0 + r1 * 2], 1
    movu            xm1, [r0 + r1]
    vinserti128     m1, m1, [r0 + r4], 1

    punpcklbw       m2, m0, m1
    punpckhbw       m3, m0, m1
    vperm2i128      m4, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    pmaddubsw       m4, m8
    pmaddubsw       m3, m2, m9
    paddw           m4, m3
    pmaddubsw       m2, m8

    vextracti128    xm0, m0, 1
    lea             r7, [r0 + r1 * 4]
    vinserti128     m0, m0, [r7], 1

    punpcklbw       m5, m1, m0
    punpckhbw       m3, m1, m0
    vperm2i128      m6, m5, m3, 0x20
    vperm2i128      m5, m5, m3, 0x31
    pmaddubsw       m6, m8
    pmaddubsw       m3, m5, m9
    paddw           m6, m3
    pmaddubsw       m5, m8
%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 0
    pmulhrsw        m6, m7                         ; m6 = word: row 1
    packuswb        m4, m6
    vpermq          m4, m4, 11011000b
    vextracti128    xm6, m4, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm6
%else
    psubw           m4, m7                         ; m4 = word: row 0
    psubw           m6, m7                         ; m6 = word: row 1
    movu            [r2], m4
    movu            [r2 + r3], m6
%endif

    movu            xm4, [r7 + r1 * 2]
    vinserti128     m4, m4, [r7 + r1], 1
    vextracti128    xm1, m4, 1
    vinserti128     m0, m0, xm1, 0

    punpcklbw       m6, m0, m4
    punpckhbw       m1, m0, m4
    vperm2i128      m0, m6, m1, 0x20
    vperm2i128      m6, m6, m1, 0x31
    pmaddubsw       m1, m0, m9
    paddw           m5, m1
    pmaddubsw       m0, m8
    pmaddubsw       m1, m6, m9
    paddw           m2, m1
    pmaddubsw       m6, m8

%ifidn %1,pp
    pmulhrsw        m2, m7                         ; m2 = word: row 2
    pmulhrsw        m5, m7                         ; m5 = word: row 3
    packuswb        m2, m5
    vpermq          m2, m2, 11011000b
    vextracti128    xm5, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm5
%else
    psubw           m2, m7                         ; m2 = word: row 2
    psubw           m5, m7                         ; m5 = word: row 3
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m5
%endif
    lea             r8, [r2 + r3 * 4]

    movu            xm1, [r7 + r4]
    lea             r7, [r7 + r1 * 4]
    vinserti128     m1, m1, [r7], 1
    vinserti128     m4, m4, xm1, 1

    punpcklbw       m2, m4, m1
    punpckhbw       m5, m4, m1
    vperm2i128      m3, m2, m5, 0x20
    vperm2i128      m2, m2, m5, 0x31
    pmaddubsw       m5, m3, m9
    paddw           m6, m5
    pmaddubsw       m3, m8
    pmaddubsw       m5, m2, m9
    paddw           m0, m5
    pmaddubsw       m2, m8

%ifidn %1,pp
    pmulhrsw        m6, m7                         ; m6 = word: row 4
    pmulhrsw        m0, m7                         ; m0 = word: row 5
    packuswb        m6, m0
    vpermq          m6, m6, 11011000b
    vextracti128    xm0, m6, 1
    movu            [r8], xm6
    movu            [r8 + r3], xm0
%else
    psubw           m6, m7                         ; m6 = word: row 4
    psubw           m0, m7                         ; m0 = word: row 5
    movu            [r8], m6
    movu            [r8 + r3], m0
%endif

    movu            xm6, [r7 + r1 * 2]
    vinserti128     m6, m6, [r7 + r1], 1
    vextracti128    xm0, m6, 1
    vinserti128     m1, m1, xm0, 0

    punpcklbw       m4, m1, m6
    punpckhbw       m5, m1, m6
    vperm2i128      m0, m4, m5, 0x20
    vperm2i128      m5, m4, m5, 0x31
    pmaddubsw       m4, m0, m9
    paddw           m2, m4
    pmaddubsw       m0, m8
    pmaddubsw       m4, m5, m9
    paddw           m3, m4
    pmaddubsw       m5, m8

%ifidn %1,pp
    pmulhrsw        m3, m7                         ; m3 = word: row 6
    pmulhrsw        m2, m7                         ; m2 = word: row 7
    packuswb        m3, m2
    vpermq          m3, m3, 11011000b
    vextracti128    xm2, m3, 1
    movu            [r8 + r3 * 2], xm3
    movu            [r8 + r6], xm2
%else
    psubw           m3, m7                         ; m3 = word: row 6
    psubw           m2, m7                         ; m2 = word: row 7
    movu            [r8 + r3 * 2], m3
    movu            [r8 + r6], m2
%endif
    lea             r8, [r8 + r3 * 4]

    movu            xm3, [r7 + r4]
    lea             r7, [r7 + r1 * 4]
    vinserti128     m3, m3, [r7], 1
    vinserti128     m6, m6, xm3, 1

    punpcklbw       m2, m6, m3
    punpckhbw       m1, m6, m3
    vperm2i128      m4, m2, m1, 0x20
    vperm2i128      m2, m2, m1, 0x31
    pmaddubsw       m1, m4, m9
    paddw           m5, m1
    pmaddubsw       m4, m8
    pmaddubsw       m1, m2, m9
    paddw           m0, m1
    pmaddubsw       m2, m8

%ifidn %1,pp
    pmulhrsw        m5, m7                         ; m5 = word: row 8
    pmulhrsw        m0, m7                         ; m0 = word: row 9
    packuswb        m5, m0
    vpermq          m5, m5, 11011000b
    vextracti128    xm0, m5, 1
    movu            [r8], xm5
    movu            [r8 + r3], xm0
%else
    psubw           m5, m7                         ; m5 = word: row 8
    psubw           m0, m7                         ; m0 = word: row 9
    movu            [r8], m5
    movu            [r8 + r3], m0
%endif

    movu            xm5, [r7 + r1 * 2]
    vinserti128     m5, m5, [r7 + r1], 1
    vextracti128    xm0, m5, 1
    vinserti128     m3, m3, xm0, 0

    punpcklbw       m1, m3, m5
    punpckhbw       m0, m3, m5
    vperm2i128      m6, m1, m0, 0x20
    vperm2i128      m0, m1, m0, 0x31
    pmaddubsw       m1, m6, m9
    paddw           m2, m1
    pmaddubsw       m6, m8
    pmaddubsw       m1, m0, m9
    paddw           m4, m1
    pmaddubsw       m0, m8

%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 10
    pmulhrsw        m2, m7                         ; m2 = word: row 11
    packuswb        m4, m2
    vpermq          m4, m4, 11011000b
    vextracti128    xm2, m4, 1
    movu            [r8 + r3 * 2], xm4
    movu            [r8 + r6], xm2
%else
    psubw           m4, m7                         ; m4 = word: row 10
    psubw           m2, m7                         ; m2 = word: row 11
    movu            [r8 + r3 * 2], m4
    movu            [r8 + r6], m2
%endif
    lea             r8, [r8 + r3 * 4]

    movu            xm3, [r7 + r4]
    lea             r7, [r7 + r1 * 4]
    vinserti128     m3, m3, [r7], 1
    vinserti128     m5, m5, xm3, 1

    punpcklbw       m2, m5, m3
    punpckhbw       m1, m5, m3
    vperm2i128      m4, m2, m1, 0x20
    vperm2i128      m2, m2, m1, 0x31
    pmaddubsw       m1, m4, m9
    paddw           m0, m1
    pmaddubsw       m4, m8
    pmaddubsw       m1, m2, m9
    paddw           m6, m1
    pmaddubsw       m2, m8

%ifidn %1,pp
    pmulhrsw        m0, m7                         ; m0 = word: row 12
    pmulhrsw        m6, m7                         ; m6 = word: row 13
    packuswb        m0, m6
    vpermq          m0, m0, 11011000b
    vextracti128    xm6, m0, 1
    movu            [r8], xm0
    movu            [r8 + r3], xm6
%else
    psubw           m0, m7                         ; m0 = word: row 12
    psubw           m6, m7                         ; m6 = word: row 13
    movu            [r8], m0
    movu            [r8 + r3], m6
%endif

    movu            xm5, [r7 + r1 * 2]
    vinserti128     m5, m5, [r7 + r1], 1
    vextracti128    xm0, m5, 1
    vinserti128     m3, m3, xm0, 0

    punpcklbw       m1, m3, m5
    punpckhbw       m0, m3, m5
    vperm2i128      m6, m1, m0, 0x20
    vperm2i128      m0, m1, m0, 0x31
    pmaddubsw       m6, m9
    paddw           m2, m6
    pmaddubsw       m0, m9
    paddw           m4, m0

%ifidn %1,pp
    pmulhrsw        m4, m7                         ; m4 = word: row 14
    pmulhrsw        m2, m7                         ; m2 = word: row 15
    packuswb        m4, m2
    vpermq          m4, m4, 11011000b
    vextracti128    xm2, m4, 1
    movu            [r8 + r3 * 2], xm4
    movu            [r8 + r6], xm2
    add             r2, 16
%else
    psubw           m4, m7                         ; m4 = word: row 14
    psubw           m2, m7                         ; m2 = word: row 15
    movu            [r8 + r3 * 2], m4
    movu            [r8 + r6], m2
    add             r2, 32
%endif
    add             r0, 16
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3
    vinserti128     m5, m1, xm2, 1
    pmaddubsw       m5, m8
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4
    lea             r7, [r0 + r1 * 4]
    movq            xm1, [r7]                       ; m1 = row 4
    punpcklbw       xm4, xm1
    vinserti128     m2, m3, xm4, 1
    pmaddubsw       m0, m2, m9
    paddw           m5, m0
    pmaddubsw       m2, m8
    movq            xm3, [r7 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3
    movq            xm4, [r7 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m0, m1, m9
    paddw           m2, m0
    pmaddubsw       m1, m8
    movq            xm3, [r7 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm0, [r7]                       ; m0 = row 8
    punpcklbw       xm3, xm0
    vinserti128     m4, m4, xm3, 1
    pmaddubsw       m3, m4, m9
    paddw           m1, m3
    pmaddubsw       m4, m8
    movq            xm3, [r7 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3
    movq            xm6, [r7 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6
    vinserti128     m0, m0, xm3, 1
    pmaddubsw       m3, m0, m9
    paddw           m4, m3
    pmaddubsw       m0, m8

%ifidn %1,pp
    pmulhrsw        m5, m7                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m7                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m7                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m7                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm2
    lea             r8, [r2 + r3 * 4]
    movq            [r8], xm1
    movq            [r8 + r3], xm4
    movhps          [r8 + r3 * 2], xm1
    movhps          [r8 + r6], xm4
%else
    psubw           m5, m7                          ; m5 = word: row 0, row 1
    psubw           m2, m7                          ; m2 = word: row 2, row 3
    psubw           m1, m7                          ; m1 = word: row 4, row 5
    psubw           m4, m7                          ; m4 = word: row 6, row 7
    vextracti128    xm3, m5, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm3
    vextracti128    xm3, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    vextracti128    xm3, m1, 1
    lea             r8, [r2 + r3 * 4]
    movu            [r8], xm1
    movu            [r8 + r3], xm3
    vextracti128    xm3, m4, 1
    movu            [r8 + r3 * 2], xm4
    movu            [r8 + r6], xm3
%endif
    lea             r8, [r8 + r3 * 4]

    movq            xm3, [r7 + r4]                  ; m3 = row 11
    punpcklbw       xm6, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm5, [r7]                       ; m5 = row 12
    punpcklbw       xm3, xm5
    vinserti128     m6, m6, xm3, 1
    pmaddubsw       m3, m6, m9
    paddw           m0, m3
    pmaddubsw       m6, m8
    movq            xm3, [r7 + r1]                  ; m3 = row 13
    punpcklbw       xm5, xm3
    movq            xm2, [r7 + r1 * 2]              ; m2 = row 14
    punpcklbw       xm3, xm2
    vinserti128     m5, m5, xm3, 1
    pmaddubsw       m3, m5, m9
    paddw           m6, m3
    pmaddubsw       m5, m8
    movq            xm3, [r7 + r4]                  ; m3 = row 15
    punpcklbw       xm2, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm1, [r7]                       ; m1 = row 16
    punpcklbw       xm3, xm1
    vinserti128     m2, m2, xm3, 1
    pmaddubsw       m3, m2, m9
    paddw           m5, m3
    pmaddubsw       m2, m8
    movq            xm3, [r7 + r1]                  ; m3 = row 17
    punpcklbw       xm1, xm3
    movq            xm4, [r7 + r1 * 2]              ; m4 = row 18
    punpcklbw       xm3, xm4
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, m9
    paddw           m2, m3
%ifidn %1,pp
    pmulhrsw        m0, m7                          ; m0 = word: row 8, row 9
    pmulhrsw        m6, m7                          ; m6 = word: row 10, row 11
    pmulhrsw        m5, m7                          ; m5 = word: row 12, row 13
    pmulhrsw        m2, m7                          ; m2 = word: row 14, row 15
    packuswb        m0, m6
    packuswb        m5, m2
    vextracti128    xm6, m0, 1
    vextracti128    xm2, m5, 1
    movq            [r8], xm0
    movq            [r8 + r3], xm6
    movhps          [r8 + r3 * 2], xm0
    movhps          [r8 + r6], xm6
    lea             r8, [r8 + r3 * 4]
    movq            [r8], xm5
    movq            [r8 + r3], xm2
    movhps          [r8 + r3 * 2], xm5
    movhps          [r8 + r6], xm2
    lea             r2, [r8 + r3 * 4 - 16]
%else
    psubw           m0, m7                          ; m0 = word: row 8, row 9
    psubw           m6, m7                          ; m6 = word: row 10, row 11
    psubw           m5, m7                          ; m5 = word: row 12, row 13
    psubw           m2, m7                          ; m2 = word: row 14, row 15
    vextracti128    xm3, m0, 1
    movu            [r8], xm0
    movu            [r8 + r3], xm3
    vextracti128    xm3, m6, 1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm3
    vextracti128    xm3, m5, 1
    lea             r8, [r8 + r3 * 4]
    movu            [r8], xm5
    movu            [r8 + r3], xm3
    vextracti128    xm3, m2, 1
    movu            [r8 + r3 * 2], xm2
    movu            [r8 + r6], xm3
    lea             r2, [r8 + r3 * 4 - 32]
%endif
    lea             r0, [r7 - 16]
    dec             r5d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_24x32 pp
FILTER_VER_CHROMA_AVX2_24x32 ps

%macro FILTER_VER_CHROMA_AVX2_16x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_16x4, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    mova            m7, [pw_2000]
%endif

    movu            xm0, [r0]
    vinserti128     m0, m0, [r0 + r1 * 2], 1
    movu            xm1, [r0 + r1]
    vinserti128     m1, m1, [r0 + r4], 1

    punpcklbw       m2, m0, m1
    punpckhbw       m3, m0, m1
    vperm2i128      m4, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    pmaddubsw       m4, [r5]
    pmaddubsw       m3, m2, [r5 + mmsize]
    paddw           m4, m3
    pmaddubsw       m2, [r5]

    vextracti128    xm0, m0, 1
    lea             r0, [r0 + r1 * 4]
    vinserti128     m0, m0, [r0], 1

    punpcklbw       m5, m1, m0
    punpckhbw       m3, m1, m0
    vperm2i128      m6, m5, m3, 0x20
    vperm2i128      m5, m5, m3, 0x31
    pmaddubsw       m6, [r5]
    pmaddubsw       m3, m5, [r5 + mmsize]
    paddw           m6, m3
    pmaddubsw       m5, [r5]
%ifidn %1,pp
    pmulhrsw        m4, m7                          ; m4 = word: row 0
    pmulhrsw        m6, m7                          ; m6 = word: row 1
    packuswb        m4, m6
    vpermq          m4, m4, 11011000b
    vextracti128    xm6, m4, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm6
%else
    psubw           m4, m7                          ; m4 = word: row 0
    psubw           m6, m7                          ; m6 = word: row 1
    movu            [r2], m4
    movu            [r2 + r3], m6
%endif
    lea             r2, [r2 + r3 * 2]

    movu            xm4, [r0 + r1 * 2]
    vinserti128     m4, m4, [r0 + r1], 1
    vextracti128    xm1, m4, 1
    vinserti128     m0, m0, xm1, 0

    punpcklbw       m6, m0, m4
    punpckhbw       m1, m0, m4
    vperm2i128      m0, m6, m1, 0x20
    vperm2i128      m6, m6, m1, 0x31
    pmaddubsw       m0, [r5 + mmsize]
    paddw           m5, m0
    pmaddubsw       m6, [r5 + mmsize]
    paddw           m2, m6

%ifidn %1,pp
    pmulhrsw        m2, m7                          ; m2 = word: row 2
    pmulhrsw        m5, m7                          ; m5 = word: row 3
    packuswb        m2, m5
    vpermq          m2, m2, 11011000b
    vextracti128    xm5, m2, 1
    movu            [r2], xm2
    movu            [r2 + r3], xm5
%else
    psubw           m2, m7                          ; m2 = word: row 2
    psubw           m5, m7                          ; m5 = word: row 3
    movu            [r2], m2
    movu            [r2 + r3], m5
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_16x4 pp
FILTER_VER_CHROMA_AVX2_16x4 ps

%macro FILTER_VER_CHROMA_AVX2_12x16 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_12x16, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m7, [pw_2000]
%endif
    lea             r6, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
%ifidn %1,pp
    pmulhrsw        m0, m7                          ; m0 = word: row 0
    pmulhrsw        m1, m7                          ; m1 = word: row 1
    packuswb        m0, m1
    vextracti128    xm1, m0, 1
    movq            [r2], xm0
    movd            [r2 + 8], xm1
    movhps          [r2 + r3], xm0
    pextrd          [r2 + r3 + 8], xm1, 2
%else
    psubw           m0, m7                          ; m0 = word: row 0
    psubw           m1, m7                          ; m1 = word: row 1
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    movq            [r2 + 16], xm0
    movu            [r2 + r3], xm1
    vextracti128    xm1, m1, 1
    movq            [r2 + r3 + 16], xm1
%endif

    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm0, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm0, 1
    pmaddubsw       m0, m5, [r5 + 1 * mmsize]
    paddw           m3, m0
    pmaddubsw       m5, [r5]
%ifidn %1,pp
    pmulhrsw        m2, m7                          ; m2 = word: row 2
    pmulhrsw        m3, m7                          ; m3 = word: row 3
    packuswb        m2, m3
    vextracti128    xm3, m2, 1
    movq            [r2 + r3 * 2], xm2
    movd            [r2 + r3 * 2 + 8], xm3
    movhps          [r2 + r6], xm2
    pextrd          [r2 + r6 + 8], xm3, 2
%else
    psubw           m2, m7                          ; m2 = word: row 2
    psubw           m3, m7                          ; m3 = word: row 3
    movu            [r2 + r3 * 2], xm2
    vextracti128    xm2, m2, 1
    movq            [r2 + r3 * 2 + 16], xm2
    movu            [r2 + r6], xm3
    vextracti128    xm3, m3, 1
    movq            [r2 + r6 + 16], xm3
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm0, [r0 + r4]                  ; m0 = row 7
    punpckhbw       xm3, xm6, xm0
    punpcklbw       xm6, xm0
    vinserti128     m6, m6, xm3, 1
    pmaddubsw       m3, m6, [r5 + 1 * mmsize]
    paddw           m4, m3
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm3, [r0]                       ; m3 = row 8
    punpckhbw       xm1, xm0, xm3
    punpcklbw       xm0, xm3
    vinserti128     m0, m0, xm1, 1
    pmaddubsw       m1, m0, [r5 + 1 * mmsize]
    paddw           m5, m1
    pmaddubsw       m0, [r5]
%ifidn %1,pp
    pmulhrsw        m4, m7                          ; m4 = word: row 4
    pmulhrsw        m5, m7                          ; m5 = word: row 5
    packuswb        m4, m5
    vextracti128    xm5, m4, 1
    movq            [r2], xm4
    movd            [r2 + 8], xm5
    movhps          [r2 + r3], xm4
    pextrd          [r2 + r3 + 8], xm5, 2
%else
    psubw           m4, m7                          ; m4 = word: row 4
    psubw           m5, m7                          ; m5 = word: row 5
    movu            [r2], xm4
    vextracti128    xm4, m4, 1
    movq            [r2 + 16], xm4
    movu            [r2 + r3], xm5
    vextracti128    xm5, m5, 1
    movq            [r2 + r3 + 16], xm5
%endif

    movu            xm1, [r0 + r1]                  ; m1 = row 9
    punpckhbw       xm2, xm3, xm1
    punpcklbw       xm3, xm1
    vinserti128     m3, m3, xm2, 1
    pmaddubsw       m2, m3, [r5 + 1 * mmsize]
    paddw           m6, m2
    pmaddubsw       m3, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 10
    punpckhbw       xm4, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm4, 1
    pmaddubsw       m4, m1, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m1, [r5]

%ifidn %1,pp
    pmulhrsw        m6, m7                          ; m6 = word: row 6
    pmulhrsw        m0, m7                          ; m0 = word: row 7
    packuswb        m6, m0
    vextracti128    xm0, m6, 1
    movq            [r2 + r3 * 2], xm6
    movd            [r2 + r3 * 2 + 8], xm0
    movhps          [r2 + r6], xm6
    pextrd          [r2 + r6 + 8], xm0, 2
%else
    psubw           m6, m7                          ; m6 = word: row 6
    psubw           m0, m7                          ; m0 = word: row 7
    movu            [r2 + r3 * 2], xm6
    vextracti128    xm6, m6, 1
    movq            [r2 + r3 * 2 + 16], xm6
    movu            [r2 + r6], xm0
    vextracti128    xm0, m0, 1
    movq            [r2 + r6 + 16], xm0
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm4, [r0 + r4]                  ; m4 = row 11
    punpckhbw       xm6, xm2, xm4
    punpcklbw       xm2, xm4
    vinserti128     m2, m2, xm6, 1
    pmaddubsw       m6, m2, [r5 + 1 * mmsize]
    paddw           m3, m6
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm6, [r0]                       ; m6 = row 12
    punpckhbw       xm0, xm4, xm6
    punpcklbw       xm4, xm6
    vinserti128     m4, m4, xm0, 1
    pmaddubsw       m0, m4, [r5 + 1 * mmsize]
    paddw           m1, m0
    pmaddubsw       m4, [r5]
%ifidn %1,pp
    pmulhrsw        m3, m7                          ; m3 = word: row 8
    pmulhrsw        m1, m7                          ; m1 = word: row 9
    packuswb        m3, m1
    vextracti128    xm1, m3, 1
    movq            [r2], xm3
    movd            [r2 + 8], xm1
    movhps          [r2 + r3], xm3
    pextrd          [r2 + r3 + 8], xm1, 2
%else
    psubw           m3, m7                          ; m3 = word: row 8
    psubw           m1, m7                          ; m1 = word: row 9
    movu            [r2], xm3
    vextracti128    xm3, m3, 1
    movq            [r2 + 16], xm3
    movu            [r2 + r3], xm1
    vextracti128    xm1, m1, 1
    movq            [r2 + r3 + 16], xm1
%endif

    movu            xm0, [r0 + r1]                  ; m0 = row 13
    punpckhbw       xm1, xm6, xm0
    punpcklbw       xm6, xm0
    vinserti128     m6, m6, xm1, 1
    pmaddubsw       m1, m6, [r5 + 1 * mmsize]
    paddw           m2, m1
    pmaddubsw       m6, [r5]
    movu            xm1, [r0 + r1 * 2]              ; m1 = row 14
    punpckhbw       xm5, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm5, 1
    pmaddubsw       m5, m0, [r5 + 1 * mmsize]
    paddw           m4, m5
    pmaddubsw       m0, [r5]
%ifidn %1,pp
    pmulhrsw        m2, m7                          ; m2 = word: row 10
    pmulhrsw        m4, m7                          ; m4 = word: row 11
    packuswb        m2, m4
    vextracti128    xm4, m2, 1
    movq            [r2 + r3 * 2], xm2
    movd            [r2 + r3 * 2 + 8], xm4
    movhps          [r2 + r6], xm2
    pextrd          [r2 + r6 + 8], xm4, 2
%else
    psubw           m2, m7                          ; m2 = word: row 10
    psubw           m4, m7                          ; m4 = word: row 11
    movu            [r2 + r3 * 2], xm2
    vextracti128    xm2, m2, 1
    movq            [r2 + r3 * 2 + 16], xm2
    movu            [r2 + r6], xm4
    vextracti128    xm4, m4, 1
    movq            [r2 + r6 + 16], xm4
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm5, [r0 + r4]                  ; m5 = row 15
    punpckhbw       xm2, xm1, xm5
    punpcklbw       xm1, xm5
    vinserti128     m1, m1, xm2, 1
    pmaddubsw       m2, m1, [r5 + 1 * mmsize]
    paddw           m6, m2
    pmaddubsw       m1, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm5, xm2
    punpcklbw       xm5, xm2
    vinserti128     m5, m5, xm3, 1
    pmaddubsw       m3, m5, [r5 + 1 * mmsize]
    paddw           m0, m3
    pmaddubsw       m5, [r5]
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m2, [r5 + 1 * mmsize]
    paddw           m1, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm2, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm2, 1
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m5, m3

%ifidn %1,pp
    pmulhrsw        m6, m7                          ; m6 = word: row 12
    pmulhrsw        m0, m7                          ; m0 = word: row 13
    pmulhrsw        m1, m7                          ; m1 = word: row 14
    pmulhrsw        m5, m7                          ; m5 = word: row 15
    packuswb        m6, m0
    packuswb        m1, m5
    vextracti128    xm0, m6, 1
    vextracti128    xm5, m1, 1
    movq            [r2], xm6
    movd            [r2 + 8], xm0
    movhps          [r2 + r3], xm6
    pextrd          [r2 + r3 + 8], xm0, 2
    movq            [r2 + r3 * 2], xm1
    movd            [r2 + r3 * 2 + 8], xm5
    movhps          [r2 + r6], xm1
    pextrd          [r2 + r6 + 8], xm5, 2
%else
    psubw           m6, m7                          ; m6 = word: row 12
    psubw           m0, m7                          ; m0 = word: row 13
    psubw           m1, m7                          ; m1 = word: row 14
    psubw           m5, m7                          ; m5 = word: row 15
    movu            [r2], xm6
    vextracti128    xm6, m6, 1
    movq            [r2 + 16], xm6
    movu            [r2 + r3], xm0
    vextracti128    xm0, m0, 1
    movq            [r2 + r3 + 16], xm0
    movu            [r2 + r3 * 2], xm1
    vextracti128    xm1, m1, 1
    movq            [r2 + r3 * 2 + 16], xm1
    movu            [r2 + r6], xm5
    vextracti128    xm5, m5, 1
    movq            [r2 + r6 + 16], xm5
%endif
    RET
%endmacro

FILTER_VER_CHROMA_AVX2_12x16 pp
FILTER_VER_CHROMA_AVX2_12x16 ps

;-----------------------------------------------------------------------------
;void interp_4tap_vert_pp_24x32(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W24 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_24x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m1,        m0,       [tab_Vm]
pshufb      m0,        [tab_Vm + 16]

mov         r4d,       %2

.loop:
movu        m2,        [r0]
movu        m3,        [r0 + r1]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

lea         r5,        [r0 + 2 * r1]
movu        m5,        [r5]
movu        m7,        [r5 + r1]

punpcklbw   m6,        m5,        m7
pmaddubsw   m6,        m0
paddw       m4,        m6

punpckhbw   m6,        m5,        m7
pmaddubsw   m6,        m0
paddw       m2,        m6

mova        m6,        [pw_512]

pmulhrsw    m4,        m6
pmulhrsw    m2,        m6

packuswb    m4,        m2

movu        [r2],      m4

punpcklbw   m4,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m4,        m1
pmaddubsw   m3,        m1

movu        m2,        [r5 + 2 * r1]

punpcklbw   m5,        m7,        m2
punpckhbw   m7,        m2

pmaddubsw   m5,        m0
pmaddubsw   m7,        m0

paddw       m4,        m5
paddw       m3,        m7

pmulhrsw    m4,        m6
pmulhrsw    m3,        m6

packuswb    m4,        m3

movu        [r2 + r3],      m4

movq        m2,        [r0 + 16]
movq        m3,        [r0 + r1 + 16]
movq        m4,        [r5 + 16]
movq        m5,        [r5 + r1 + 16]

punpcklbw   m2,        m3
punpcklbw   m4,        m5

pmaddubsw   m2,        m1
pmaddubsw   m4,        m0

paddw       m2,        m4

pmulhrsw    m2,        m6

movq        m3,        [r0 + r1 + 16]
movq        m4,        [r5 + 16]
movq        m5,        [r5 + r1 + 16]
movq        m7,        [r5 + 2 * r1 + 16]

punpcklbw   m3,        m4
punpcklbw   m5,        m7

pmaddubsw   m3,        m1
pmaddubsw   m5,        m0

paddw       m3,        m5

pmulhrsw    m3,        m6
packuswb    m2,        m3

movh        [r2 + 16], m2
movhps      [r2 + r3 + 16], m2

mov         r0,        r5
lea         r2,        [r2 + 2 * r3]

sub         r4,        2
jnz        .loop
RET
%endmacro

FILTER_V4_W24 24, 32

FILTER_V4_W24 24, 64

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_32x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W32 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_%1x%2, 4, 6, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m1,        m0,       [tab_Vm]
pshufb      m0,        [tab_Vm + 16]

mova        m7,        [pw_512]

mov         r4d,       %2

.loop:
movu        m2,        [r0]
movu        m3,        [r0 + r1]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

lea         r5,        [r0 + 2 * r1]
movu        m3,        [r5]
movu        m5,        [r5 + r1]

punpcklbw   m6,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m6,        m0
pmaddubsw   m3,        m0

paddw       m4,        m6
paddw       m2,        m3

pmulhrsw    m4,        m7
pmulhrsw    m2,        m7

packuswb    m4,        m2

movu        [r2],      m4

movu        m2,        [r0 + 16]
movu        m3,        [r0 + r1 + 16]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

movu        m3,        [r5 + 16]
movu        m5,        [r5 + r1 + 16]

punpcklbw   m6,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m6,        m0
pmaddubsw   m3,        m0

paddw       m4,        m6
paddw       m2,        m3

pmulhrsw    m4,        m7
pmulhrsw    m2,        m7

packuswb    m4,        m2

movu        [r2 + 16], m4

lea         r0,        [r0 + r1]
lea         r2,        [r2 + r3]

dec         r4
jnz        .loop
RET
%endmacro

FILTER_V4_W32 32,  8
FILTER_V4_W32 32, 16
FILTER_V4_W32 32, 24
FILTER_V4_W32 32, 32

FILTER_V4_W32 32, 48
FILTER_V4_W32 32, 64

%macro FILTER_VER_CHROMA_AVX2_32xN 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_32x%2, 4, 7, 13
    mov             r4d, r4m
    shl             r4d, 6

%ifdef PIC
    lea             r5, [tab_ChromaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_ChromaCoeffVer_32 + r4]
%endif

    mova            m10, [r5]
    mova            m11, [r5 + mmsize]
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,pp
    mova            m12, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m12, [pw_2000]
%endif
    lea             r5, [r3 * 3]
    mov             r6d, %2 / 4
.loopW:
    movu            m0, [r0]                        ; m0 = row 0
    movu            m1, [r0 + r1]                   ; m1 = row 1
    punpcklbw       m2, m0, m1
    punpckhbw       m3, m0, m1
    pmaddubsw       m2, m10
    pmaddubsw       m3, m10
    movu            m0, [r0 + r1 * 2]               ; m0 = row 2
    punpcklbw       m4, m1, m0
    punpckhbw       m5, m1, m0
    pmaddubsw       m4, m10
    pmaddubsw       m5, m10
    movu            m1, [r0 + r4]                   ; m1 = row 3
    punpcklbw       m6, m0, m1
    punpckhbw       m7, m0, m1
    pmaddubsw       m8, m6, m11
    pmaddubsw       m9, m7, m11
    pmaddubsw       m6, m10
    pmaddubsw       m7, m10
    paddw           m2, m8
    paddw           m3, m9
%ifidn %1,pp
    pmulhrsw        m2, m12
    pmulhrsw        m3, m12
    packuswb        m2, m3
    movu            [r2], m2
%else
    psubw           m2, m12
    psubw           m3, m12
    vperm2i128      m0, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    movu            [r2], m0
    movu            [r2 + mmsize], m2
%endif
    lea             r0, [r0 + r1 * 4]
    movu            m0, [r0]                        ; m0 = row 4
    punpcklbw       m2, m1, m0
    punpckhbw       m3, m1, m0
    pmaddubsw       m8, m2, m11
    pmaddubsw       m9, m3, m11
    pmaddubsw       m2, m10
    pmaddubsw       m3, m10
    paddw           m4, m8
    paddw           m5, m9
%ifidn %1,pp
    pmulhrsw        m4, m12
    pmulhrsw        m5, m12
    packuswb        m4, m5
    movu            [r2 + r3], m4
%else
    psubw           m4, m12
    psubw           m5, m12
    vperm2i128      m1, m4, m5, 0x20
    vperm2i128      m4, m4, m5, 0x31
    movu            [r2 + r3], m1
    movu            [r2 + r3 + mmsize], m4
%endif

    movu            m1, [r0 + r1]                   ; m1 = row 5
    punpcklbw       m4, m0, m1
    punpckhbw       m5, m0, m1
    pmaddubsw       m4, m11
    pmaddubsw       m5, m11
    paddw           m6, m4
    paddw           m7, m5
%ifidn %1,pp
    pmulhrsw        m6, m12
    pmulhrsw        m7, m12
    packuswb        m6, m7
    movu            [r2 + r3 * 2], m6
%else
    psubw           m6, m12
    psubw           m7, m12
    vperm2i128      m0, m6, m7, 0x20
    vperm2i128      m6, m6, m7, 0x31
    movu            [r2 + r3 * 2], m0
    movu            [r2 + r3 * 2 + mmsize], m6
%endif

    movu            m0, [r0 + r1 * 2]               ; m0 = row 6
    punpcklbw       m6, m1, m0
    punpckhbw       m7, m1, m0
    pmaddubsw       m6, m11
    pmaddubsw       m7, m11
    paddw           m2, m6
    paddw           m3, m7
%ifidn %1,pp
    pmulhrsw        m2, m12
    pmulhrsw        m3, m12
    packuswb        m2, m3
    movu            [r2 + r5], m2
%else
    psubw           m2, m12
    psubw           m3, m12
    vperm2i128      m0, m2, m3, 0x20
    vperm2i128      m2, m2, m3, 0x31
    movu            [r2 + r5], m0
    movu            [r2 + r5 + mmsize], m2
%endif
    lea             r2, [r2 + r3 * 4]
    dec             r6d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_CHROMA_AVX2_32xN pp, 32
FILTER_VER_CHROMA_AVX2_32xN pp, 24
FILTER_VER_CHROMA_AVX2_32xN pp, 16
FILTER_VER_CHROMA_AVX2_32xN pp, 8
FILTER_VER_CHROMA_AVX2_32xN ps, 32
FILTER_VER_CHROMA_AVX2_32xN ps, 24
FILTER_VER_CHROMA_AVX2_32xN ps, 16
FILTER_VER_CHROMA_AVX2_32xN ps, 8

;-----------------------------------------------------------------------------
; void interp_4tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_V4_W16n_H2 2
INIT_XMM sse4
cglobal interp_4tap_vert_pp_%1x%2, 4, 7, 8

mov         r4d,       r4m
sub         r0,        r1

%ifdef PIC
lea         r5,        [tab_ChromaCoeff]
movd        m0,        [r5 + r4 * 4]
%else
movd        m0,        [tab_ChromaCoeff + r4 * 4]
%endif

pshufb      m1,        m0,       [tab_Vm]
pshufb      m0,        [tab_Vm + 16]

mov         r4d,       %2/2

.loop:

mov         r6d,       %1/16

.loopW:

movu        m2,        [r0]
movu        m3,        [r0 + r1]

punpcklbw   m4,        m2,        m3
punpckhbw   m2,        m3

pmaddubsw   m4,        m1
pmaddubsw   m2,        m1

lea         r5,        [r0 + 2 * r1]
movu        m5,        [r5]
movu        m6,        [r5 + r1]

punpckhbw   m7,        m5,        m6
pmaddubsw   m7,        m0
paddw       m2,        m7

punpcklbw   m7,        m5,        m6
pmaddubsw   m7,        m0
paddw       m4,        m7

mova        m7,        [pw_512]

pmulhrsw    m4,        m7
pmulhrsw    m2,        m7

packuswb    m4,        m2

movu        [r2],      m4

punpcklbw   m4,        m3,        m5
punpckhbw   m3,        m5

pmaddubsw   m4,        m1
pmaddubsw   m3,        m1

movu        m5,        [r5 + 2 * r1]

punpcklbw   m2,        m6,        m5
punpckhbw   m6,        m5

pmaddubsw   m2,        m0
pmaddubsw   m6,        m0

paddw       m4,        m2
paddw       m3,        m6

pmulhrsw    m4,        m7
pmulhrsw    m3,        m7

packuswb    m4,        m3

movu        [r2 + r3],      m4

add         r0,        16
add         r2,        16
dec         r6d
jnz         .loopW

lea         r0,        [r0 + r1 * 2 - %1]
lea         r2,        [r2 + r3 * 2 - %1]

dec         r4d
jnz        .loop
RET
%endmacro

FILTER_V4_W16n_H2 64, 64
FILTER_V4_W16n_H2 64, 32
FILTER_V4_W16n_H2 64, 48
FILTER_V4_W16n_H2 48, 64
FILTER_V4_W16n_H2 64, 16
;-----------------------------------------------------------------------------
; void pixelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;-----------------------------------------------------------------------------
%macro PIXEL_WH_4xN 2
INIT_XMM ssse3
cglobal pixelToShort_%1x%2, 3, 7, 6

    ; load width and height
    mov         r3d, %1
    mov         r4d, %2
    ; load constant
    mova        m4, [pb_128]
    mova        m5, [tab_c_64_n64]
.loopH:
    xor         r5d, r5d

.loopW:
    mov         r6, r0
    movh        m0, [r6]
    punpcklbw   m0, m4
    pmaddubsw   m0, m5

    movh        m1, [r6 + r1]
    punpcklbw   m1, m4
    pmaddubsw   m1, m5

    movh        m2, [r6 + r1 * 2]
    punpcklbw   m2, m4
    pmaddubsw   m2, m5

    lea         r6, [r6 + r1 * 2]
    movh        m3, [r6 + r1]
    punpcklbw   m3, m4
    pmaddubsw   m3, m5

    add         r5, 8
    cmp         r5, r3
    jg          .width4
    movu        [r2 + r5 * 2 + FENC_STRIDE * 0 - 16], m0
    movu        [r2 + r5 * 2 + FENC_STRIDE * 2 - 16], m1
    movu        [r2 + r5 * 2 + FENC_STRIDE * 4 - 16], m2
    movu        [r2 + r5 * 2 + FENC_STRIDE * 6 - 16], m3
    je          .nextH
    jmp         .loopW

.width4:
    movh        [r2 + r5 * 2 + FENC_STRIDE * 0 - 16], m0
    movh        [r2 + r5 * 2 + FENC_STRIDE * 2 - 16], m1
    movh        [r2 + r5 * 2 + FENC_STRIDE * 4 - 16], m2
    movh        [r2 + r5 * 2 + FENC_STRIDE * 6 - 16], m3

.nextH:
    lea         r0, [r0 + r1 * 4]
    add         r2, FENC_STRIDE * 8

    sub         r4d, 4
    jnz         .loopH
    RET
%endmacro
PIXEL_WH_4xN 4, 4
PIXEL_WH_4xN 4, 8
PIXEL_WH_4xN 4, 16

;-----------------------------------------------------------------------------
; void pixelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;-----------------------------------------------------------------------------
%macro PIXEL_WH_8xN 2
INIT_XMM ssse3
cglobal pixelToShort_%1x%2, 3, 7, 6

    ; load width and height
    mov         r3d, %1
    mov         r4d, %2

    ; load constant
    mova        m4, [pb_128]
    mova        m5, [tab_c_64_n64]

.loopH
    xor         r5d, r5d
.loopW
    lea         r6, [r0 + r5]

    movh        m0, [r6]
    punpcklbw   m0, m4
    pmaddubsw   m0, m5

    movh        m1, [r6 + r1]
    punpcklbw   m1, m4
    pmaddubsw   m1, m5

    movh        m2, [r6 + r1 * 2]
    punpcklbw   m2, m4
    pmaddubsw   m2, m5

    lea         r6, [r6 + r1 * 2]
    movh        m3, [r6 + r1]
    punpcklbw   m3, m4
    pmaddubsw   m3, m5

    add         r5, 8
    cmp         r5, r3

    movu        [r2 + FENC_STRIDE * 0], m0
    movu        [r2 + FENC_STRIDE * 2], m1
    movu        [r2 + FENC_STRIDE * 4], m2
    movu        [r2 + FENC_STRIDE * 6], m3

    je          .nextH
    jmp         .loopW


.nextH:
    lea         r0, [r0 + r1 * 4]
    add         r2, FENC_STRIDE * 8

    sub         r4d, 4
    jnz         .loopH
    RET
%endmacro
PIXEL_WH_8xN 8, 8
PIXEL_WH_8xN 8, 4
PIXEL_WH_8xN 8, 16
PIXEL_WH_8xN 8, 32


;-----------------------------------------------------------------------------
; void pixelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;-----------------------------------------------------------------------------
%macro PIXEL_WH_16xN 2
INIT_XMM ssse3
cglobal pixelToShort_%1x%2, 3, 7, 6

    ; load width and height
    mov         r3d, %1
    mov         r4d, %2

    ; load constant
    mova        m4, [pb_128]
    mova        m5, [tab_c_64_n64]

.loopH:
    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5]

    movh        m0, [r6]
    punpcklbw   m0, m4
    pmaddubsw   m0, m5

    movh        m1, [r6 + r1]
    punpcklbw   m1, m4
    pmaddubsw   m1, m5

    movh        m2, [r6 + r1 * 2]
    punpcklbw   m2, m4
    pmaddubsw   m2, m5

    lea         r6, [r6 + r1 * 2]
    movh        m3, [r6 + r1]
    punpcklbw   m3, m4
    pmaddubsw   m3, m5

    add         r5, 8
    cmp         r5, r3

    movu        [r2 + r5 * 2 + FENC_STRIDE * 0 - 16], m0
    movu        [r2 + r5 * 2 + FENC_STRIDE * 2 - 16], m1
    movu        [r2 + r5 * 2 + FENC_STRIDE * 4 - 16], m2
    movu        [r2 + r5 * 2 + FENC_STRIDE * 6 - 16], m3
    je          .nextH
    jmp         .loopW


.nextH:
    lea         r0, [r0 + r1 * 4]
    add         r2, FENC_STRIDE * 8

    sub         r4d, 4
    jnz         .loopH

    RET
%endmacro
PIXEL_WH_16xN 16, 16
PIXEL_WH_16xN 16, 8
PIXEL_WH_16xN 16, 4
PIXEL_WH_16xN 16, 12
PIXEL_WH_16xN 16, 32
PIXEL_WH_16xN 16, 64

;-----------------------------------------------------------------------------
; void pixelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;-----------------------------------------------------------------------------
%macro PIXEL_WH_32xN 2
INIT_XMM ssse3
cglobal pixelToShort_%1x%2, 3, 7, 6

    ; load width and height
    mov         r3d, %1
    mov         r4d, %2

    ; load constant
    mova        m4, [pb_128]
    mova        m5, [tab_c_64_n64]

.loopH:
    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5]

    movh        m0, [r6]
    punpcklbw   m0, m4
    pmaddubsw   m0, m5

    movh        m1, [r6 + r1]
    punpcklbw   m1, m4
    pmaddubsw   m1, m5

    movh        m2, [r6 + r1 * 2]
    punpcklbw   m2, m4
    pmaddubsw   m2, m5

    lea         r6, [r6 + r1 * 2]
    movh        m3, [r6 + r1]
    punpcklbw   m3, m4
    pmaddubsw   m3, m5

    add         r5, 8
    cmp         r5, r3

    movu        [r2 + r5 * 2 + FENC_STRIDE * 0 - 16], m0
    movu        [r2 + r5 * 2 + FENC_STRIDE * 2 - 16], m1
    movu        [r2 + r5 * 2 + FENC_STRIDE * 4 - 16], m2
    movu        [r2 + r5 * 2 + FENC_STRIDE * 6 - 16], m3
    je          .nextH
    jmp         .loopW


.nextH:
    lea         r0, [r0 + r1 * 4]
    add         r2, FENC_STRIDE * 8

    sub         r4d, 4
    jnz         .loopH

    RET
%endmacro
PIXEL_WH_32xN 32, 32
PIXEL_WH_32xN 32, 8
PIXEL_WH_32xN 32, 16
PIXEL_WH_32xN 32, 24
PIXEL_WH_32xN 32, 64

;-----------------------------------------------------------------------------
; void pixelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;-----------------------------------------------------------------------------
%macro PIXEL_WH_64xN 2
INIT_XMM ssse3
cglobal pixelToShort_%1x%2, 3, 7, 6

    ; load width and height
    mov         r3d, %1
    mov         r4d, %2

    ; load constant
    mova        m4, [pb_128]
    mova        m5, [tab_c_64_n64]

.loopH:
    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5]

    movh        m0, [r6]
    punpcklbw   m0, m4
    pmaddubsw   m0, m5

    movh        m1, [r6 + r1]
    punpcklbw   m1, m4
    pmaddubsw   m1, m5

    movh        m2, [r6 + r1 * 2]
    punpcklbw   m2, m4
    pmaddubsw   m2, m5

    lea         r6, [r6 + r1 * 2]
    movh        m3, [r6 + r1]
    punpcklbw   m3, m4
    pmaddubsw   m3, m5

    add         r5, 8
    cmp         r5, r3

    movu        [r2 + r5 * 2 + FENC_STRIDE * 0 - 16], m0
    movu        [r2 + r5 * 2 + FENC_STRIDE * 2 - 16], m1
    movu        [r2 + r5 * 2 + FENC_STRIDE * 4 - 16], m2
    movu        [r2 + r5 * 2 + FENC_STRIDE * 6 - 16], m3
    je          .nextH
    jmp         .loopW


.nextH:
    lea         r0, [r0 + r1 * 4]
    add         r2, FENC_STRIDE * 8

    sub         r4d, 4
    jnz         .loopH

    RET
%endmacro
PIXEL_WH_64xN 64, 64
PIXEL_WH_64xN 64, 16
PIXEL_WH_64xN 64, 32
PIXEL_WH_64xN 64, 48

%macro PROCESS_LUMA_W4_4R 0
    movd        m0, [r0]
    movd        m1, [r0 + r1]
    punpcklbw   m2, m0, m1                     ; m2=[0 1]

    lea         r0, [r0 + 2 * r1]
    movd        m0, [r0]
    punpcklbw   m1, m0                         ; m1=[1 2]
    punpcklqdq  m2, m1                         ; m2=[0 1 1 2]
    pmaddubsw   m4, m2, [r6 + 0 * 16]          ; m4=[0+1 1+2]

    movd        m1, [r0 + r1]
    punpcklbw   m5, m0, m1                     ; m2=[2 3]
    lea         r0, [r0 + 2 * r1]
    movd        m0, [r0]
    punpcklbw   m1, m0                         ; m1=[3 4]
    punpcklqdq  m5, m1                         ; m5=[2 3 3 4]
    pmaddubsw   m2, m5, [r6 + 1 * 16]          ; m2=[2+3 3+4]
    paddw       m4, m2                         ; m4=[0+1+2+3 1+2+3+4]                   Row1-2
    pmaddubsw   m5, [r6 + 0 * 16]              ; m5=[2+3 3+4]                           Row3-4

    movd        m1, [r0 + r1]
    punpcklbw   m2, m0, m1                     ; m2=[4 5]
    lea         r0, [r0 + 2 * r1]
    movd        m0, [r0]
    punpcklbw   m1, m0                         ; m1=[5 6]
    punpcklqdq  m2, m1                         ; m2=[4 5 5 6]
    pmaddubsw   m1, m2, [r6 + 2 * 16]          ; m1=[4+5 5+6]
    paddw       m4, m1                         ; m4=[0+1+2+3+4+5 1+2+3+4+5+6]           Row1-2
    pmaddubsw   m2, [r6 + 1 * 16]              ; m2=[4+5 5+6]
    paddw       m5, m2                         ; m5=[2+3+4+5 3+4+5+6]                   Row3-4

    movd        m1, [r0 + r1]
    punpcklbw   m2, m0, m1                     ; m2=[6 7]
    lea         r0, [r0 + 2 * r1]
    movd        m0, [r0]
    punpcklbw   m1, m0                         ; m1=[7 8]
    punpcklqdq  m2, m1                         ; m2=[6 7 7 8]
    pmaddubsw   m1, m2, [r6 + 3 * 16]          ; m1=[6+7 7+8]
    paddw       m4, m1                         ; m4=[0+1+2+3+4+5+6+7 1+2+3+4+5+6+7+8]   Row1-2 end
    pmaddubsw   m2, [r6 + 2 * 16]              ; m2=[6+7 7+8]
    paddw       m5, m2                         ; m5=[2+3+4+5+6+7 3+4+5+6+7+8]           Row3-4

    movd        m1, [r0 + r1]
    punpcklbw   m2, m0, m1                     ; m2=[8 9]
    movd        m0, [r0 + 2 * r1]
    punpcklbw   m1, m0                         ; m1=[9 10]
    punpcklqdq  m2, m1                         ; m2=[8 9 9 10]
    pmaddubsw   m2, [r6 + 3 * 16]              ; m2=[8+9 9+10]
    paddw       m5, m2                         ; m5=[2+3+4+5+6+7+8+9 3+4+5+6+7+8+9+10]  Row3-4 end
%endmacro

%macro PROCESS_LUMA_W8_4R 0
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklbw  m0, m1
    pmaddubsw  m7, m0, [r6 + 0 *16]            ;m7=[0+1]               Row1

    lea        r0, [r0 + 2 * r1]
    movq       m0, [r0]
    punpcklbw  m1, m0
    pmaddubsw  m6, m1, [r6 + 0 *16]            ;m6=[1+2]               Row2

    movq       m1, [r0 + r1]
    punpcklbw  m0, m1
    pmaddubsw  m5, m0, [r6 + 0 *16]            ;m5=[2+3]               Row3
    pmaddubsw  m0, [r6 + 1 * 16]
    paddw      m7, m0                          ;m7=[0+1+2+3]           Row1

    lea        r0, [r0 + 2 * r1]
    movq       m0, [r0]
    punpcklbw  m1, m0
    pmaddubsw  m4, m1, [r6 + 0 *16]            ;m4=[3+4]               Row4
    pmaddubsw  m1, [r6 + 1 * 16]
    paddw      m6, m1                          ;m6 = [1+2+3+4]         Row2

    movq       m1, [r0 + r1]
    punpcklbw  m0, m1
    pmaddubsw  m2, m0, [r6 + 1 * 16]
    pmaddubsw  m0, [r6 + 2 * 16]
    paddw      m7, m0                          ;m7=[0+1+2+3+4+5]       Row1
    paddw      m5, m2                          ;m5=[2+3+4+5]           Row3

    lea        r0, [r0 + 2 * r1]
    movq       m0, [r0]
    punpcklbw  m1, m0
    pmaddubsw  m2, m1, [r6 + 1 * 16]
    pmaddubsw  m1, [r6 + 2 * 16]
    paddw      m6, m1                          ;m6=[1+2+3+4+5+6]       Row2
    paddw      m4, m2                          ;m4=[3+4+5+6]           Row4

    movq       m1, [r0 + r1]
    punpcklbw  m0, m1
    pmaddubsw  m2, m0, [r6 + 2 * 16]
    pmaddubsw  m0, [r6 + 3 * 16]
    paddw      m7, m0                          ;m7=[0+1+2+3+4+5+6+7]   Row1 end
    paddw      m5, m2                          ;m5=[2+3+4+5+6+7]       Row3

    lea        r0, [r0 + 2 * r1]
    movq       m0, [r0]
    punpcklbw  m1, m0
    pmaddubsw  m2, m1, [r6 + 2 * 16]
    pmaddubsw  m1, [r6 + 3 * 16]
    paddw      m6, m1                          ;m6=[1+2+3+4+5+6+7+8]   Row2 end
    paddw      m4, m2                          ;m4=[3+4+5+6+7+8]       Row4

    movq       m1, [r0 + r1]
    punpcklbw  m0, m1
    pmaddubsw  m0, [r6 + 3 * 16]
    paddw      m5, m0                          ;m5=[2+3+4+5+6+7+8+9]   Row3 end

    movq       m0, [r0 + 2 * r1]
    punpcklbw  m1, m0
    pmaddubsw  m1, [r6 + 3 * 16]
    paddw      m4, m1                          ;m4=[3+4+5+6+7+8+9+10]  Row4 end
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_%3_4x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_4xN 3
INIT_XMM sse4
cglobal interp_8tap_vert_%3_%1x%2, 5, 7, 6
    lea       r5, [3 * r1]
    sub       r0, r5
    shl       r4d, 6
%ifidn %3,ps
    add       r3d, r3d
%endif

%ifdef PIC
    lea       r5, [tab_LumaCoeffVer]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffVer + r4]
%endif

%ifidn %3,pp
    mova      m3, [pw_512]
%else
    mova      m3, [pw_2000]
%endif

    mov       r4d, %2/4
    lea       r5, [4 * r1]

.loopH:
    PROCESS_LUMA_W4_4R

%ifidn %3,pp
    pmulhrsw  m4, m3
    pmulhrsw  m5, m3

    packuswb  m4, m5

    movd      [r2], m4
    pextrd    [r2 + r3], m4, 1
    lea       r2, [r2 + 2 * r3]
    pextrd    [r2], m4, 2
    pextrd    [r2 + r3], m4, 3
%else
    psubw     m4, m3
    psubw     m5, m3

    movlps    [r2], m4
    movhps    [r2 + r3], m4
    lea       r2, [r2 + 2 * r3]
    movlps    [r2], m5
    movhps    [r2 + r3], m5
%endif

    sub       r0, r5
    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro


INIT_YMM avx2
cglobal interp_8tap_vert_pp_4x4, 4,6,8
    mov             r4d, r4m
    lea             r5, [r1 * 3]
    sub             r0, r5

    ; TODO: VPGATHERDD
    movd            xm1, [r0]                       ; m1 = row0
    movd            xm2, [r0 + r1]                  ; m2 = row1
    punpcklbw       xm1, xm2                        ; m1 = [13 03 12 02 11 01 10 00]

    movd            xm3, [r0 + r1 * 2]              ; m3 = row2
    punpcklbw       xm2, xm3                        ; m2 = [23 13 22 12 21 11 20 10]
    movd            xm4, [r0 + r5]
    punpcklbw       xm3, xm4                        ; m3 = [33 23 32 22 31 21 30 20]
    punpcklwd       xm1, xm3                        ; m1 = [33 23 13 03 32 22 12 02 31 21 11 01 30 20 10 00]

    lea             r0, [r0 + r1 * 4]
    movd            xm5, [r0]                       ; m5 = row4
    punpcklbw       xm4, xm5                        ; m4 = [43 33 42 32 41 31 40 30]
    punpcklwd       xm2, xm4                        ; m2 = [43 33 21 13 42 32 22 12 41 31 21 11 40 30 20 10]
    vinserti128     m1, m1, xm2, 1                  ; m1 = [43 33 21 13 42 32 22 12 41 31 21 11 40 30 20 10] - [33 23 13 03 32 22 12 02 31 21 11 01 30 20 10 00]
    movd            xm2, [r0 + r1]                  ; m2 = row5
    punpcklbw       xm5, xm2                        ; m5 = [53 43 52 42 51 41 50 40]
    punpcklwd       xm3, xm5                        ; m3 = [53 43 44 23 52 42 32 22 51 41 31 21 50 40 30 20]
    movd            xm6, [r0 + r1 * 2]              ; m6 = row6
    punpcklbw       xm2, xm6                        ; m2 = [63 53 62 52 61 51 60 50]
    punpcklwd       xm4, xm2                        ; m4 = [63 53 43 33 62 52 42 32 61 51 41 31 60 50 40 30]
    vinserti128     m3, m3, xm4, 1                  ; m3 = [63 53 43 33 62 52 42 32 61 51 41 31 60 50 40 30] - [53 43 44 23 52 42 32 22 51 41 31 21 50 40 30 20]
    movd            xm4, [r0 + r5]                  ; m4 = row7
    punpcklbw       xm6, xm4                        ; m6 = [73 63 72 62 71 61 70 60]
    punpcklwd       xm5, xm6                        ; m5 = [73 63 53 43 72 62 52 42 71 61 51 41 70 60 50 40]

    lea             r0, [r0 + r1 * 4]
    movd            xm7, [r0]                       ; m7 = row8
    punpcklbw       xm4, xm7                        ; m4 = [83 73 82 72 81 71 80 70]
    punpcklwd       xm2, xm4                        ; m2 = [83 73 63 53 82 72 62 52 81 71 61 51 80 70 60 50]
    vinserti128     m5, m5, xm2, 1                  ; m5 = [83 73 63 53 82 72 62 52 81 71 61 51 80 70 60 50] - [73 63 53 43 72 62 52 42 71 61 51 41 70 60 50 40]
    movd            xm2, [r0 + r1]                  ; m2 = row9
    punpcklbw       xm7, xm2                        ; m7 = [93 83 92 82 91 81 90 80]
    punpcklwd       xm6, xm7                        ; m6 = [93 83 73 63 92 82 72 62 91 81 71 61 90 80 70 60]
    movd            xm7, [r0 + r1 * 2]              ; m7 = rowA
    punpcklbw       xm2, xm7                        ; m2 = [A3 93 A2 92 A1 91 A0 90]
    punpcklwd       xm4, xm2                        ; m4 = [A3 93 83 73 A2 92 82 72 A1 91 81 71 A0 90 80 70]
    vinserti128     m6, m6, xm4, 1                  ; m6 = [A3 93 83 73 A2 92 82 72 A1 91 81 71 A0 90 80 70] - [93 83 73 63 92 82 72 62 91 81 71 61 90 80 70 60]

    ; load filter coeff
%ifdef PIC
    lea             r5, [tab_LumaCoeff]
    vpbroadcastd    m0, [r5 + r4 * 8 + 0]
    vpbroadcastd    m2, [r5 + r4 * 8 + 4]
%else
    vpbroadcastd    m0, [tab_LumaCoeff + r4 * 8 + 0]
    vpbroadcastd    m2, [tab_LumaCoeff + r4 * 8 + 4]
%endif

    pmaddubsw       m1, m0
    pmaddubsw       m3, m0
    pmaddubsw       m5, m2
    pmaddubsw       m6, m2
    vbroadcasti128  m0, [pw_1]
    pmaddwd         m1, m0
    pmaddwd         m3, m0
    pmaddwd         m5, m0
    pmaddwd         m6, m0
    paddd           m1, m5                          ; m1 = DQWORD ROW[1 0]
    paddd           m3, m6                          ; m3 = DQWORD ROW[3 2]
    packssdw        m1, m3                          ; m1 =  QWORD ROW[3 1 2 0]

    ; TODO: does it overflow?
    pmulhrsw        m1, [pw_512]
    vextracti128    xm2, m1, 1
    packuswb        xm1, xm2                        ; m1 =  DWORD ROW[3 1 2 0]
    movd            [r2], xm1
    pextrd          [r2 + r3], xm1, 2
    pextrd          [r2 + r3 * 2], xm1, 1
    lea             r4, [r3 * 3]
    pextrd          [r2 + r4], xm1, 3
    RET

INIT_YMM avx2
cglobal interp_8tap_vert_ps_4x4, 4, 6, 5
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

    add             r3d, r3d

    movd            xm1, [r0]
    pinsrd          xm1, [r0 + r1], 1
    pinsrd          xm1, [r0 + r1 * 2], 2
    pinsrd          xm1, [r0 + r4], 3                       ; m1 = row[3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm2, [r0]
    pinsrd          xm2, [r0 + r1], 1
    pinsrd          xm2, [r0 + r1 * 2], 2
    pinsrd          xm2, [r0 + r4], 3                       ; m2 = row[7 6 5 4]
    vinserti128     m1, m1, xm2, 1                          ; m1 = row[7 6 5 4 3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm3, [r0]
    pinsrd          xm3, [r0 + r1], 1
    pinsrd          xm3, [r0 + r1 * 2], 2                   ; m3 = row[x 10 9 8]
    vinserti128     m2, m2, xm3, 1                          ; m2 = row[x 10 9 8 7 6 5 4]
    mova            m3, [interp4_vpp_shuf1]
    vpermd          m0, m3, m1                              ; m0 = row[4 3 3 2 2 1 1 0]
    vpermd          m4, m3, m2                              ; m4 = row[8 7 7 6 6 5 5 4]
    mova            m3, [interp4_vpp_shuf1 + mmsize]
    vpermd          m1, m3, m1                              ; m1 = row[6 5 5 4 4 3 3 2]
    vpermd          m2, m3, m2                              ; m2 = row[10 9 9 8 8 7 7 6]

    mova            m3, [interp4_vpp_shuf]
    pshufb          m0, m0, m3
    pshufb          m1, m1, m3
    pshufb          m4, m4, m3
    pshufb          m2, m2, m3
    pmaddubsw       m0, [r5]
    pmaddubsw       m1, [r5 + mmsize]
    pmaddubsw       m4, [r5 + 2 * mmsize]
    pmaddubsw       m2, [r5 + 3 * mmsize]
    paddw           m0, m1
    paddw           m0, m4
    paddw           m0, m2                                  ; m0 = WORD ROW[3 2 1 0]

    psubw           m0, [pw_2000]
    vextracti128    xm2, m0, 1
    lea             r5, [r3 * 3]
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r5], xm2
    RET

%macro FILTER_VER_LUMA_AVX2_4xN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%3_%1x%2, 4, 9, 10
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
    lea             r6, [r1 * 4]
%ifidn %3,pp
    mova            m6, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m6, [pw_2000]
%endif
    lea             r8, [r3 * 3]
    mova            m5, [interp4_vpp_shuf]
    mova            m0, [interp4_vpp_shuf1]
    mova            m7, [interp4_vpp_shuf1 + mmsize]
    mov             r7d, %2 / 8
.loop:
    movd            xm1, [r0]
    pinsrd          xm1, [r0 + r1], 1
    pinsrd          xm1, [r0 + r1 * 2], 2
    pinsrd          xm1, [r0 + r4], 3                       ; m1 = row[3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm2, [r0]
    pinsrd          xm2, [r0 + r1], 1
    pinsrd          xm2, [r0 + r1 * 2], 2
    pinsrd          xm2, [r0 + r4], 3                       ; m2 = row[7 6 5 4]
    vinserti128     m1, m1, xm2, 1                          ; m1 = row[7 6 5 4 3 2 1 0]
    lea             r0, [r0 + r1 * 4]
    movd            xm3, [r0]
    pinsrd          xm3, [r0 + r1], 1
    pinsrd          xm3, [r0 + r1 * 2], 2
    pinsrd          xm3, [r0 + r4], 3                       ; m3 = row[11 10 9 8]
    vinserti128     m2, m2, xm3, 1                          ; m2 = row[11 10 9 8 7 6 5 4]
    lea             r0, [r0 + r1 * 4]
    movd            xm4, [r0]
    pinsrd          xm4, [r0 + r1], 1
    pinsrd          xm4, [r0 + r1 * 2], 2                   ; m4 = row[x 14 13 12]
    vinserti128     m3, m3, xm4, 1                          ; m3 = row[x 14 13 12 11 10 9 8]
    vpermd          m8, m0, m1                              ; m8 = row[4 3 3 2 2 1 1 0]
    vpermd          m4, m0, m2                              ; m4 = row[8 7 7 6 6 5 5 4]
    vpermd          m1, m7, m1                              ; m1 = row[6 5 5 4 4 3 3 2]
    vpermd          m2, m7, m2                              ; m2 = row[10 9 9 8 8 7 7 6]
    vpermd          m9, m0, m3                              ; m9 = row[12 11 11 10 10 9 9 8]
    vpermd          m3, m7, m3                              ; m3 = row[14 13 13 12 12 11 11 10]

    pshufb          m8, m8, m5
    pshufb          m1, m1, m5
    pshufb          m4, m4, m5
    pshufb          m9, m9, m5
    pshufb          m2, m2, m5
    pshufb          m3, m3, m5
    pmaddubsw       m8, [r5]
    pmaddubsw       m1, [r5 + mmsize]
    pmaddubsw       m9, [r5 + 2 * mmsize]
    pmaddubsw       m3, [r5 + 3 * mmsize]
    paddw           m8, m1
    paddw           m9, m3
    pmaddubsw       m1, m4, [r5 + 2 * mmsize]
    pmaddubsw       m3, m2, [r5 + 3 * mmsize]
    pmaddubsw       m4, [r5]
    pmaddubsw       m2, [r5 + mmsize]
    paddw           m3, m1
    paddw           m2, m4
    paddw           m8, m3                                  ; m8 = WORD ROW[3 2 1 0]
    paddw           m9, m2                                  ; m9 = WORD ROW[7 6 5 4]

%ifidn %3,pp
    pmulhrsw        m8, m6
    pmulhrsw        m9, m6
    packuswb        m8, m9
    vextracti128    xm1, m8, 1
    movd            [r2], xm8
    pextrd          [r2 + r3], xm8, 1
    movd            [r2 + r3 * 2], xm1
    pextrd          [r2 + r8], xm1, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm8, 2
    pextrd          [r2 + r3], xm8, 3
    pextrd          [r2 + r3 * 2], xm1, 2
    pextrd          [r2 + r8], xm1, 3
%else
    psubw           m8, m6
    psubw           m9, m6
    vextracti128    xm1, m8, 1
    vextracti128    xm2, m9, 1
    movq            [r2], xm8
    movhps          [r2 + r3], xm8
    movq            [r2 + r3 * 2], xm1
    movhps          [r2 + r8], xm1
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm9
    movhps          [r2 + r3], xm9
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r8], xm2
%endif
    lea             r2, [r2 + r3 * 4]
    sub             r0, r6
    dec             r7d
    jnz             .loop
    RET
%endif
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_4x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 4, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_4x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 8, pp
FILTER_VER_LUMA_AVX2_4xN 4, 8, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_4x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 16, pp
FILTER_VER_LUMA_AVX2_4xN 4, 16, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_4x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 4, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_4x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 8, ps
FILTER_VER_LUMA_AVX2_4xN 4, 8, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_4x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_4xN 4, 16, ps
FILTER_VER_LUMA_AVX2_4xN 4, 16, ps

%macro PROCESS_LUMA_AVX2_W8_8R 0
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m5, m1, xm2, 1                  ; m5 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1                        ; m4 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m2, m3, xm4, 1                  ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3                        ; m1 = [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4                        ; m3 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50]
    vinserti128     m1, m1, xm3, 1                  ; m1 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50] - [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m5, m3
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    pmaddubsw       m1, [r5]
    movq            xm3, [r0 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3                        ; m4 = [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 8
    punpcklbw       xm3, xm0                        ; m3 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70]
    vinserti128     m4, m4, xm3, 1                  ; m4 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70] - [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    pmaddubsw       m3, m4, [r5 + 3 * mmsize]
    paddw           m5, m3
    pmaddubsw       m3, m4, [r5 + 2 * mmsize]
    paddw           m2, m3
    pmaddubsw       m3, m4, [r5 + 1 * mmsize]
    paddw           m1, m3
    pmaddubsw       m4, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3                        ; m0 = [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    movq            xm6, [r0 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6                        ; m3 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90]
    vinserti128     m0, m0, xm3, 1                  ; m0 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90] - [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    pmaddubsw       m3, m0, [r5 + 3 * mmsize]
    paddw           m2, m3
    pmaddubsw       m3, m0, [r5 + 2 * mmsize]
    paddw           m1, m3
    pmaddubsw       m0, [r5 + 1 * mmsize]
    paddw           m4, m0

    movq            xm3, [r0 + r4]                  ; m3 = row 11
    punpcklbw       xm6, xm3                        ; m6 = [B7 A7 B6 A6 B5 A5 B4 A4 B3 A3 B2 A2 B1 A1 B0 A0]
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 12
    punpcklbw       xm3, xm0                        ; m3 = [C7 B7 C6 B6 C5 B5 C4 B4 C3 B3 C2 B2 C1 B1 C0 B0]
    vinserti128     m6, m6, xm3, 1                  ; m6 = [C7 B7 C6 B6 C5 B5 C4 B4 C3 B3 C2 B2 C1 B1 C0 B0] - [B7 A7 B6 A6 B5 A5 B4 A4 B3 A3 B2 A2 B1 A1 B0 A0]
    pmaddubsw       m3, m6, [r5 + 3 * mmsize]
    paddw           m1, m3
    pmaddubsw       m6, [r5 + 2 * mmsize]
    paddw           m4, m6
    movq            xm3, [r0 + r1]                  ; m3 = row 13
    punpcklbw       xm0, xm3                        ; m0 = [D7 C7 D6 C6 D5 C5 D4 C4 D3 C3 D2 C2 D1 C1 D0 C0]
    movq            xm6, [r0 + r1 * 2]              ; m6 = row 14
    punpcklbw       xm3, xm6                        ; m3 = [E7 D7 E6 D6 E5 D5 E4 D4 E3 D3 E2 D2 E1 D1 E0 D0]
    vinserti128     m0, m0, xm3, 1                  ; m0 = [E7 D7 E6 D6 E5 D5 E4 D4 E3 D3 E2 D2 E1 D1 E0 D0] - [D7 C7 D6 C6 D5 C5 D4 C4 D3 C3 D2 C2 D1 C1 D0 C0]
    pmaddubsw       m0, [r5 + 3 * mmsize]
    paddw           m4, m0
%endmacro

%macro PROCESS_LUMA_AVX2_W8_4R 0
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2                        ; m1 = [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3                        ; m2 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10]
    vinserti128     m5, m1, xm2, 1                  ; m5 = [27 17 26 16 25 15 24 14 23 13 22 12 21 11 20 10] - [17 07 16 06 15 05 14 04 13 03 12 02 11 01 10 00]
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4                        ; m3 = [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    lea             r0, [r0 + r1 * 4]
    movq            xm1, [r0]                       ; m1 = row 4
    punpcklbw       xm4, xm1                        ; m4 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30]
    vinserti128     m2, m3, xm4, 1                  ; m2 = [47 37 46 36 45 35 44 34 43 33 42 32 41 31 40 30] - [37 27 36 26 35 25 34 24 33 23 32 22 31 21 30 20]
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r0 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3                        ; m1 = [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    movq            xm4, [r0 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4                        ; m3 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50]
    vinserti128     m1, m1, xm3, 1                  ; m1 = [67 57 66 56 65 55 64 54 63 53 62 52 61 51 60 50] - [57 47 56 46 55 45 54 44 53 43 52 42 51 41 50 40]
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m5, m3
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    movq            xm3, [r0 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3                        ; m4 = [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    lea             r0, [r0 + r1 * 4]
    movq            xm0, [r0]                       ; m0 = row 8
    punpcklbw       xm3, xm0                        ; m3 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70]
    vinserti128     m4, m4, xm3, 1                  ; m4 = [87 77 86 76 85 75 84 74 83 73 82 72 81 71 80 70] - [77 67 76 66 75 65 74 64 73 63 72 62 71 61 70 60]
    pmaddubsw       m3, m4, [r5 + 3 * mmsize]
    paddw           m5, m3
    pmaddubsw       m3, m4, [r5 + 2 * mmsize]
    paddw           m2, m3
    movq            xm3, [r0 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3                        ; m0 = [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    movq            xm6, [r0 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6                        ; m3 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90]
    vinserti128     m0, m0, xm3, 1                  ; m0 = [A7 97 A6 96 A5 95 A4 94 A3 93 A2 92 A1 91 A0 90] - [97 87 96 86 95 85 94 84 93 83 92 82 91 81 90 80]
    pmaddubsw       m3, m0, [r5 + 3 * mmsize]
    paddw           m2, m3
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_%3_8x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_8xN 3
INIT_XMM sse4
cglobal interp_8tap_vert_%3_%1x%2, 5, 7, 8
    lea       r5, [3 * r1]
    sub       r0, r5
    shl       r4d, 6

%ifidn %3,ps
    add       r3d, r3d
%endif

%ifdef PIC
    lea       r5, [tab_LumaCoeffVer]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffVer + r4]
%endif

 %ifidn %3,pp
    mova      m3, [pw_512]
%else
    mova      m3, [pw_2000]
%endif

    mov       r4d, %2/4
    lea       r5, [4 * r1]

.loopH:
    PROCESS_LUMA_W8_4R

%ifidn %3,pp
    pmulhrsw  m7, m3
    pmulhrsw  m6, m3
    pmulhrsw  m5, m3
    pmulhrsw  m4, m3

    packuswb  m7, m6
    packuswb  m5, m4

    movlps    [r2], m7
    movhps    [r2 + r3], m7
    lea       r2, [r2 + 2 * r3]
    movlps    [r2], m5
    movhps    [r2 + r3], m5
%else
    psubw     m7, m3
    psubw     m6, m3
    psubw     m5, m3
    psubw     m4, m3

    movu      [r2], m7
    movu      [r2 + r3], m6
    lea       r2, [r2 + 2 * r3]
    movu      [r2], m5
    movu      [r2 + r3], m4
%endif

    sub       r0, r5
    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

%macro FILTER_VER_LUMA_AVX2_8xN 3
INIT_YMM avx2
cglobal interp_8tap_vert_%3_%1x%2, 4, 7, 8, 0-gprsize
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
    lea             r6, [r1 * 4]
%ifidn %3,pp
    mova            m7, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m7, [pw_2000]
%endif
    mov             word [rsp], %2 / 8

.loop:
    PROCESS_LUMA_AVX2_W8_8R
%ifidn %3,pp
    pmulhrsw        m5, m7                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m7                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m7                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m7                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    lea             r2, [r2 + r3 * 2]
    movhps          [r2], xm5
    movhps          [r2 + r3], xm2
    lea             r2, [r2 + r3 * 2]
    movq            [r2], xm1
    movq            [r2 + r3], xm4
    lea             r2, [r2 + r3 * 2]
    movhps          [r2], xm1
    movhps          [r2 + r3], xm4
%else
    psubw           m5, m7                          ; m5 = word: row 0, row 1
    psubw           m2, m7                          ; m2 = word: row 2, row 3
    psubw           m1, m7                          ; m1 = word: row 4, row 5
    psubw           m4, m7                          ; m4 = word: row 6, row 7
    vextracti128    xm6, m5, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm0, m1, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm6
    lea             r2, [r2 + r3 * 2]
    movu            [r2], xm2
    movu            [r2 + r3], xm3
    lea             r2, [r2 + r3 * 2]
    movu            [r2], xm1
    movu            [r2 + r3], xm0
    lea             r2, [r2 + r3 * 2]
    movu            [r2], xm4
    vextracti128    xm4, m4, 1
    movu            [r2 + r3], xm4
%endif
    lea             r2, [r2 + r3 * 2]
    sub             r0, r6
    dec             word [rsp]
    jnz             .loop
    RET
%endmacro

%macro FILTER_VER_LUMA_AVX2_8x8 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_8x8, 4, 6, 7
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
    PROCESS_LUMA_AVX2_W8_8R
%ifidn %1,pp
    mova            m3, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m3, [pw_2000]
%endif
    lea             r4, [r3 * 3]
%ifidn %1,pp
    pmulhrsw        m5, m3                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    pmulhrsw        m1, m3                          ; m1 = word: row 4, row 5
    pmulhrsw        m4, m3                          ; m4 = word: row 6, row 7
    packuswb        m5, m2
    packuswb        m1, m4
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r4], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm1
    movq            [r2 + r3], xm4
    movhps          [r2 + r3 * 2], xm1
    movhps          [r2 + r4], xm4
%else
    psubw           m5, m3                          ; m5 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    psubw           m1, m3                          ; m1 = word: row 4, row 5
    psubw           m4, m3                          ; m4 = word: row 6, row 7
    vextracti128    xm6, m5, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm0, m1, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm6
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm4
    vextracti128    xm4, m4, 1
    movu            [r2 + r4], xm4
%endif
    RET
%endmacro

%macro FILTER_VER_LUMA_AVX2_8x4 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_8x4, 4, 6, 7
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
    PROCESS_LUMA_AVX2_W8_4R
%ifidn %1,pp
    mova            m3, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m3, [pw_2000]
%endif
    lea             r4, [r3 * 3]
%ifidn %1,pp
    pmulhrsw        m5, m3                          ; m5 = word: row 0, row 1
    pmulhrsw        m2, m3                          ; m2 = word: row 2, row 3
    packuswb        m5, m2
    vextracti128    xm2, m5, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r4], xm2
%else
    psubw           m5, m3                          ; m5 = word: row 0, row 1
    psubw           m2, m3                          ; m2 = word: row 2, row 3
    movu            [r2], xm5
    vextracti128    xm5, m5, 1
    movu            [r2 + r3], xm5
    movu            [r2 + r3 * 2], xm2
    vextracti128    xm2, m2, 1
    movu            [r2 + r4], xm2
%endif
    RET
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_8x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 4, pp
FILTER_VER_LUMA_AVX2_8x4 pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_8x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 8, pp
FILTER_VER_LUMA_AVX2_8x8 pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_8x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 16, pp
FILTER_VER_LUMA_AVX2_8xN 8, 16, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_8x32(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 32, pp
FILTER_VER_LUMA_AVX2_8xN 8, 32, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_8x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 4, ps
FILTER_VER_LUMA_AVX2_8x4 ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_8x8(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 8, ps
FILTER_VER_LUMA_AVX2_8x8 ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_8x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 16, ps
FILTER_VER_LUMA_AVX2_8xN 8, 16, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_8x32(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_8xN 8, 32, ps
FILTER_VER_LUMA_AVX2_8xN 8, 32, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_%3_12x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_12xN 3
INIT_XMM sse4
cglobal interp_8tap_vert_%3_%1x%2, 5, 7, 8
    lea       r5, [3 * r1]
    sub       r0, r5
    shl       r4d, 6
%ifidn %3,ps
    add       r3d, r3d
%endif

%ifdef PIC
    lea       r5, [tab_LumaCoeffVer]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffVer + r4]
%endif

 %ifidn %3,pp
    mova      m3, [pw_512]
%else
    mova      m3, [pw_2000]
%endif

    mov       r4d, %2/4

.loopH:
    PROCESS_LUMA_W8_4R

%ifidn %3,pp
    pmulhrsw  m7, m3
    pmulhrsw  m6, m3
    pmulhrsw  m5, m3
    pmulhrsw  m4, m3

    packuswb  m7, m6
    packuswb  m5, m4

    movlps    [r2], m7
    movhps    [r2 + r3], m7
    lea       r5, [r2 + 2 * r3]
    movlps    [r5], m5
    movhps    [r5 + r3], m5
%else
    psubw     m7, m3
    psubw     m6, m3
    psubw     m5, m3
    psubw     m4, m3

    movu      [r2], m7
    movu      [r2 + r3], m6
    lea       r5, [r2 + 2 * r3]
    movu      [r5], m5
    movu      [r5 + r3], m4
%endif

    lea       r5, [8 * r1 - 8]
    sub       r0, r5
%ifidn %3,pp 
    add       r2, 8
%else
    add       r2, 16
%endif

    PROCESS_LUMA_W4_4R

%ifidn %3,pp
    pmulhrsw  m4, m3
    pmulhrsw  m5, m3

    packuswb  m4, m5

    movd      [r2], m4
    pextrd    [r2 + r3], m4, 1
    lea       r5, [r2 + 2 * r3]
    pextrd    [r5], m4, 2
    pextrd    [r5 + r3], m4, 3
%else
    psubw     m4, m3
    psubw     m5, m3

    movlps    [r2], m4
    movhps    [r2 + r3], m4
    lea       r5, [r2 + 2 * r3]
    movlps    [r5], m5
    movhps    [r5 + r3], m5
%endif

    lea       r5, [4 * r1 + 8]
    sub       r0, r5
%ifidn %3,pp
    lea       r2, [r2 + 4 * r3 - 8]
%else
    lea       r2, [r2 + 4 * r3 - 16]
%endif

    dec       r4d
    jnz       .loopH

    RET
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_12x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_12xN 12, 16, pp

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_12x16(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
FILTER_VER_LUMA_12xN 12, 16, ps

%macro FILTER_VER_LUMA_AVX2_12x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_12x16, 4, 7, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%endif
    lea             r6, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    pmaddubsw       m8, [r5]
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    pmaddubsw       m9, [r5]
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    pmaddubsw       m12, m10, [r5 + 1 * mmsize]
    paddw           m8, m12
    pmaddubsw       m10, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    pmaddubsw       m13, m11, [r5 + 1 * mmsize]
    paddw           m9, m13
    pmaddubsw       m11, [r5]

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movq            [r2], xm0
    pextrd          [r2 + 8], xm0, 2
    movq            [r2 + r3], xm1
    pextrd          [r2 + r3 + 8], xm1, 2
    movq            [r2 + r3 * 2], xm2
    pextrd          [r2 + r3 * 2 + 8], xm2, 2
    movq            [r2 + r6], xm3
    pextrd          [r2 + r6 + 8], xm3, 2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    pextrd          [r2 + 8], xm4, 2
    movq            [r2 + r3], xm5
    pextrd          [r2 + r3 + 8], xm5, 2
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    movq            [r2 + 16], xm0
    movu            [r2 + r3], xm1
    vextracti128    xm1, m1, 1
    movq            [r2 + r3 + 16], xm1
    movu            [r2 + r3 * 2], xm2
    vextracti128    xm2, m2, 1
    movq            [r2 + r3 * 2 + 16], xm2
    movu            [r2 + r6], xm3
    vextracti128    xm3, m3, 1
    movq            [r2 + r6 + 16], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    vextracti128    xm4, m4, 1
    movq            [r2 + 16], xm4
    movu            [r2 + r3], xm5
    vextracti128    xm5, m5, 1
    movq            [r2 + r3 + 16], xm5
%endif

    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    pmaddubsw       m0, m12, [r5 + 2 * mmsize]
    paddw           m8, m0
    pmaddubsw       m0, m12, [r5 + 1 * mmsize]
    paddw           m10, m0
    pmaddubsw       m12, [r5]
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
    pmaddubsw       m1, m13, [r5 + 2 * mmsize]
    paddw           m9, m1
    pmaddubsw       m1, m13, [r5 + 1 * mmsize]
    paddw           m11, m1
    pmaddubsw       m13, [r5]

%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movq            [r2 + r3 * 2], xm6
    pextrd          [r2 + r3 * 2 + 8], xm6, 2
    movq            [r2 + r6], xm7
    pextrd          [r2 + r6 + 8], xm7, 2
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2 + r3 * 2], xm6
    vextracti128    xm6, m6, 1
    movq            [r2 + r3 * 2 + 16], xm6
    movu            [r2 + r6], xm7
    vextracti128    xm7, m7, 1
    movq            [r2 + r6 + 16], xm7
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, [r5 + 3 * mmsize]
    paddw           m8, m2
    pmaddubsw       m2, m0, [r5 + 2 * mmsize]
    paddw           m10, m2
    pmaddubsw       m2, m0, [r5 + 1 * mmsize]
    paddw           m12, m2
    pmaddubsw       m0, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 3 * mmsize]
    paddw           m9, m3
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m11, m3
    pmaddubsw       m3, m1, [r5 + 1 * mmsize]
    paddw           m13, m3
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 3 * mmsize]
    paddw           m10, m4
    pmaddubsw       m4, m2, [r5 + 2 * mmsize]
    paddw           m12, m4
    pmaddubsw       m2, [r5 + 1 * mmsize]
    paddw           m0, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 3 * mmsize]
    paddw           m11, m5
    pmaddubsw       m5, m3, [r5 + 2 * mmsize]
    paddw           m13, m5
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m1, m3
    movu            xm5, [r0 + r4]                  ; m5 = row 19
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 3 * mmsize]
    paddw           m12, m6
    pmaddubsw       m4, [r5 + 2 * mmsize]
    paddw           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm6, [r0]                       ; m6 = row 20
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 3 * mmsize]
    paddw           m13, m7
    pmaddubsw       m5, [r5 + 2 * mmsize]
    paddw           m1, m5
    movu            xm7, [r0 + r1]                  ; m7 = row 21
    punpckhbw       xm2, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddubsw       m6, [r5 + 3 * mmsize]
    paddw           m0, m6
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 22
    punpckhbw       xm3, xm7, xm2
    punpcklbw       xm7, xm2
    vinserti128     m7, m7, xm3, 1
    pmaddubsw       m7, [r5 + 3 * mmsize]
    paddw           m1, m7

%ifidn %1,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    pmulhrsw        m12, m14                        ; m12 = word: row 12
    pmulhrsw        m13, m14                        ; m13 = word: row 13
    pmulhrsw        m0, m14                         ; m0 = word: row 14
    pmulhrsw        m1, m14                         ; m1 = word: row 15
    packuswb        m8, m9
    packuswb        m10, m11
    packuswb        m12, m13
    packuswb        m0, m1
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm1, m0, 1
    movq            [r2], xm8
    pextrd          [r2 + 8], xm8, 2
    movq            [r2 + r3], xm9
    pextrd          [r2 + r3 + 8], xm9, 2
    movq            [r2 + r3 * 2], xm10
    pextrd          [r2 + r3 * 2 + 8], xm10, 2
    movq            [r2 + r6], xm11
    pextrd          [r2 + r6 + 8], xm11, 2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm12
    pextrd          [r2 + 8], xm12, 2
    movq            [r2 + r3], xm13
    pextrd          [r2 + r3 + 8], xm13, 2
    movq            [r2 + r3 * 2], xm0
    pextrd          [r2 + r3 * 2 + 8], xm0, 2
    movq            [r2 + r6], xm1
    pextrd          [r2 + r6 + 8], xm1, 2
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    psubw           m12, m14                        ; m12 = word: row 12
    psubw           m13, m14                        ; m13 = word: row 13
    psubw           m0, m14                         ; m0 = word: row 14
    psubw           m1, m14                         ; m1 = word: row 15
    movu            [r2], xm8
    vextracti128    xm8, m8, 1
    movq            [r2 + 16], xm8
    movu            [r2 + r3], xm9
    vextracti128    xm9, m9, 1
    movq            [r2 + r3 + 16], xm9
    movu            [r2 + r3 * 2], xm10
    vextracti128    xm10, m10, 1
    movq            [r2 + r3 * 2 + 16], xm10
    movu            [r2 + r6], xm11
    vextracti128    xm11, m11, 1
    movq            [r2 + r6 + 16], xm11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm12
    vextracti128    xm12, m12, 1
    movq            [r2 + 16], xm12
    movu            [r2 + r3], xm13
    vextracti128    xm13, m13, 1
    movq            [r2 + r3 + 16], xm13
    movu            [r2 + r3 * 2], xm0
    vextracti128    xm0, m0, 1
    movq            [r2 + r3 * 2 + 16], xm0
    movu            [r2 + r6], xm1
    vextracti128    xm1, m1, 1
    movq            [r2 + r6 + 16], xm1
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_12x16 pp
FILTER_VER_LUMA_AVX2_12x16 ps

%macro FILTER_VER_LUMA_AVX2_16x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_16x16, 4, 7, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%endif
    lea             r6, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    pmaddubsw       m8, [r5]
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    pmaddubsw       m9, [r5]
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    pmaddubsw       m12, m10, [r5 + 1 * mmsize]
    paddw           m8, m12
    pmaddubsw       m10, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    pmaddubsw       m13, m11, [r5 + 1 * mmsize]
    paddw           m9, m13
    pmaddubsw       m11, [r5]

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m4
    movu            [r2 + r3], m5
%endif

    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    pmaddubsw       m0, m12, [r5 + 2 * mmsize]
    paddw           m8, m0
    pmaddubsw       m0, m12, [r5 + 1 * mmsize]
    paddw           m10, m0
    pmaddubsw       m12, [r5]
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
    pmaddubsw       m1, m13, [r5 + 2 * mmsize]
    paddw           m9, m1
    pmaddubsw       m1, m13, [r5 + 1 * mmsize]
    paddw           m11, m1
    pmaddubsw       m13, [r5]

%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2 + r3 * 2], m6
    movu            [r2 + r6], m7
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, [r5 + 3 * mmsize]
    paddw           m8, m2
    pmaddubsw       m2, m0, [r5 + 2 * mmsize]
    paddw           m10, m2
    pmaddubsw       m2, m0, [r5 + 1 * mmsize]
    paddw           m12, m2
    pmaddubsw       m0, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 3 * mmsize]
    paddw           m9, m3
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m11, m3
    pmaddubsw       m3, m1, [r5 + 1 * mmsize]
    paddw           m13, m3
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 3 * mmsize]
    paddw           m10, m4
    pmaddubsw       m4, m2, [r5 + 2 * mmsize]
    paddw           m12, m4
    pmaddubsw       m2, [r5 + 1 * mmsize]
    paddw           m0, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 3 * mmsize]
    paddw           m11, m5
    pmaddubsw       m5, m3, [r5 + 2 * mmsize]
    paddw           m13, m5
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m1, m3
    movu            xm5, [r0 + r4]                  ; m5 = row 19
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 3 * mmsize]
    paddw           m12, m6
    pmaddubsw       m4, [r5 + 2 * mmsize]
    paddw           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm6, [r0]                       ; m6 = row 20
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 3 * mmsize]
    paddw           m13, m7
    pmaddubsw       m5, [r5 + 2 * mmsize]
    paddw           m1, m5
    movu            xm7, [r0 + r1]                  ; m7 = row 21
    punpckhbw       xm2, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddubsw       m6, [r5 + 3 * mmsize]
    paddw           m0, m6
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 22
    punpckhbw       xm3, xm7, xm2
    punpcklbw       xm7, xm2
    vinserti128     m7, m7, xm3, 1
    pmaddubsw       m7, [r5 + 3 * mmsize]
    paddw           m1, m7

%ifidn %1,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    pmulhrsw        m12, m14                        ; m12 = word: row 12
    pmulhrsw        m13, m14                        ; m13 = word: row 13
    pmulhrsw        m0, m14                         ; m0 = word: row 14
    pmulhrsw        m1, m14                         ; m1 = word: row 15
    packuswb        m8, m9
    packuswb        m10, m11
    packuswb        m12, m13
    packuswb        m0, m1
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm1, m0, 1
    movu            [r2], xm8
    movu            [r2 + r3], xm9
    movu            [r2 + r3 * 2], xm10
    movu            [r2 + r6], xm11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm12
    movu            [r2 + r3], xm13
    movu            [r2 + r3 * 2], xm0
    movu            [r2 + r6], xm1
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    psubw           m12, m14                        ; m12 = word: row 12
    psubw           m13, m14                        ; m13 = word: row 13
    psubw           m0, m14                         ; m0 = word: row 14
    psubw           m1, m14                         ; m1 = word: row 15
    movu            [r2], m8
    movu            [r2 + r3], m9
    movu            [r2 + r3 * 2], m10
    movu            [r2 + r6], m11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m12
    movu            [r2 + r3], m13
    movu            [r2 + r3 * 2], m0
    movu            [r2 + r6], m1
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_16x16 pp
FILTER_VER_LUMA_AVX2_16x16 ps

%macro FILTER_VER_LUMA_AVX2_16x12 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_16x12, 4, 7, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%endif
    lea             r6, [r3 * 3]

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    pmaddubsw       m8, [r5]
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    pmaddubsw       m9, [r5]
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    pmaddubsw       m12, m10, [r5 + 1 * mmsize]
    paddw           m8, m12
    pmaddubsw       m10, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    pmaddubsw       m13, m11, [r5 + 1 * mmsize]
    paddw           m9, m13
    pmaddubsw       m11, [r5]

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m4
    movu            [r2 + r3], m5
%endif

    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    pmaddubsw       m0, m12, [r5 + 2 * mmsize]
    paddw           m8, m0
    pmaddubsw       m0, m12, [r5 + 1 * mmsize]
    paddw           m10, m0
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
    pmaddubsw       m1, m13, [r5 + 2 * mmsize]
    paddw           m9, m1
    pmaddubsw       m1, m13, [r5 + 1 * mmsize]
    paddw           m11, m1

%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2 + r3 * 2], m6
    movu            [r2 + r6], m7
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, [r5 + 3 * mmsize]
    paddw           m8, m2
    pmaddubsw       m2, m0, [r5 + 2 * mmsize]
    paddw           m10, m2
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 3 * mmsize]
    paddw           m9, m3
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m11, m3
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 3 * mmsize]
    paddw           m10, m4
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 3 * mmsize]
    paddw           m11, m5

%ifidn %1,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    packuswb        m8, m9
    packuswb        m10, m11
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    movu            [r2], xm8
    movu            [r2 + r3], xm9
    movu            [r2 + r3 * 2], xm10
    movu            [r2 + r6], xm11
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    movu            [r2], m8
    movu            [r2 + r3], m9
    movu            [r2 + r3 * 2], m10
    movu            [r2 + r6], m11
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_16x12 pp
FILTER_VER_LUMA_AVX2_16x12 ps

%macro FILTER_VER_LUMA_AVX2_16x8 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_16x8, 4, 6, 15
    mov             r4d, r4m
    shl             r4d, 7
%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%endif
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    lea             r4, [r3 * 3]
%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r4], m3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m4
    movu            [r2 + r3], m5
%endif
    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r4], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2 + r3 * 2], m6
    movu            [r2 + r4], m7
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_16x8 pp
FILTER_VER_LUMA_AVX2_16x8 ps

%macro FILTER_VER_LUMA_AVX2_16x4 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_16x4, 4, 6, 13
    mov             r4d, r4m
    shl             r4d, 7
%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,pp
    mova            m12, [pw_512]
%else
    add             r3d, r3d
    vbroadcasti128  m12, [pw_2000]
%endif
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
%ifidn %1,pp
    pmulhrsw        m0, m12                         ; m0 = word: row 0
    pmulhrsw        m1, m12                         ; m1 = word: row 1
    pmulhrsw        m2, m12                         ; m2 = word: row 2
    pmulhrsw        m3, m12                         ; m3 = word: row 3
    packuswb        m0, m1
    packuswb        m2, m3
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    lea             r4, [r3 * 3]
    movu            [r2 + r4], xm3
%else
    psubw           m0, m12                         ; m0 = word: row 0
    psubw           m1, m12                         ; m1 = word: row 1
    psubw           m2, m12                         ; m2 = word: row 2
    psubw           m3, m12                         ; m3 = word: row 3
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    lea             r4, [r3 * 3]
    movu            [r2 + r4], m3
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_16x4 pp
FILTER_VER_LUMA_AVX2_16x4 ps
%macro FILTER_VER_LUMA_AVX2_16xN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%3_%1x%2, 4, 9, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %3,ps
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%else
    mova            m14, [pw_512]
%endif
    lea             r6, [r3 * 3]
    lea             r7, [r1 * 4]
    mov             r8d, %2 / 16

.loop:
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    pmaddubsw       m8, [r5]
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    pmaddubsw       m9, [r5]
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    pmaddubsw       m12, m10, [r5 + 1 * mmsize]
    paddw           m8, m12
    pmaddubsw       m10, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    pmaddubsw       m13, m11, [r5 + 1 * mmsize]
    paddw           m9, m13
    pmaddubsw       m11, [r5]

%ifidn %3,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m4
    movu            [r2 + r3], m5
%endif

    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    pmaddubsw       m0, m12, [r5 + 2 * mmsize]
    paddw           m8, m0
    pmaddubsw       m0, m12, [r5 + 1 * mmsize]
    paddw           m10, m0
    pmaddubsw       m12, [r5]
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
    pmaddubsw       m1, m13, [r5 + 2 * mmsize]
    paddw           m9, m1
    pmaddubsw       m1, m13, [r5 + 1 * mmsize]
    paddw           m11, m1
    pmaddubsw       m13, [r5]

%ifidn %3,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r2 + r3 * 2], m6
    movu            [r2 + r6], m7
%endif

    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, [r5 + 3 * mmsize]
    paddw           m8, m2
    pmaddubsw       m2, m0, [r5 + 2 * mmsize]
    paddw           m10, m2
    pmaddubsw       m2, m0, [r5 + 1 * mmsize]
    paddw           m12, m2
    pmaddubsw       m0, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 3 * mmsize]
    paddw           m9, m3
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m11, m3
    pmaddubsw       m3, m1, [r5 + 1 * mmsize]
    paddw           m13, m3
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 3 * mmsize]
    paddw           m10, m4
    pmaddubsw       m4, m2, [r5 + 2 * mmsize]
    paddw           m12, m4
    pmaddubsw       m2, [r5 + 1 * mmsize]
    paddw           m0, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 3 * mmsize]
    paddw           m11, m5
    pmaddubsw       m5, m3, [r5 + 2 * mmsize]
    paddw           m13, m5
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m1, m3
    movu            xm5, [r0 + r4]                  ; m5 = row 19
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 3 * mmsize]
    paddw           m12, m6
    pmaddubsw       m4, [r5 + 2 * mmsize]
    paddw           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm6, [r0]                       ; m6 = row 20
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 3 * mmsize]
    paddw           m13, m7
    pmaddubsw       m5, [r5 + 2 * mmsize]
    paddw           m1, m5
    movu            xm7, [r0 + r1]                  ; m7 = row 21
    punpckhbw       xm2, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddubsw       m6, [r5 + 3 * mmsize]
    paddw           m0, m6
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 22
    punpckhbw       xm3, xm7, xm2
    punpcklbw       xm7, xm2
    vinserti128     m7, m7, xm3, 1
    pmaddubsw       m7, [r5 + 3 * mmsize]
    paddw           m1, m7

%ifidn %3,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    pmulhrsw        m12, m14                        ; m12 = word: row 12
    pmulhrsw        m13, m14                        ; m13 = word: row 13
    pmulhrsw        m0, m14                         ; m0 = word: row 14
    pmulhrsw        m1, m14                         ; m1 = word: row 15
    packuswb        m8, m9
    packuswb        m10, m11
    packuswb        m12, m13
    packuswb        m0, m1
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm1, m0, 1
    movu            [r2], xm8
    movu            [r2 + r3], xm9
    movu            [r2 + r3 * 2], xm10
    movu            [r2 + r6], xm11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm12
    movu            [r2 + r3], xm13
    movu            [r2 + r3 * 2], xm0
    movu            [r2 + r6], xm1
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    psubw           m12, m14                        ; m12 = word: row 12
    psubw           m13, m14                        ; m13 = word: row 13
    psubw           m0, m14                         ; m0 = word: row 14
    psubw           m1, m14                         ; m1 = word: row 15
    movu            [r2], m8
    movu            [r2 + r3], m9
    movu            [r2 + r3 * 2], m10
    movu            [r2 + r6], m11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], m12
    movu            [r2 + r3], m13
    movu            [r2 + r3 * 2], m0
    movu            [r2 + r6], m1
%endif

    lea             r2, [r2 + r3 * 4]
    sub             r0, r7
    dec             r8d
    jnz             .loop
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_16xN 16, 32, pp
FILTER_VER_LUMA_AVX2_16xN 16, 64, pp
FILTER_VER_LUMA_AVX2_16xN 16, 32, ps
FILTER_VER_LUMA_AVX2_16xN 16, 64, ps

%macro PROCESS_LUMA_AVX2_W16_16R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r7 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    pmaddubsw       m8, [r5]
    movu            xm10, [r7 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    pmaddubsw       m9, [r5]
    movu            xm11, [r7 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    pmaddubsw       m12, m10, [r5 + 1 * mmsize]
    paddw           m8, m12
    pmaddubsw       m10, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm12, [r7]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13
    pmaddubsw       m13, m11, [r5 + 1 * mmsize]
    paddw           m9, m13
    pmaddubsw       m11, [r5]

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r8, [r2 + r3 * 4]
    movu            [r8], xm4
    movu            [r8 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
    lea             r8, [r2 + r3 * 4]
    movu            [r8], m4
    movu            [r8 + r3], m5
%endif

    movu            xm13, [r7 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    pmaddubsw       m0, m12, [r5 + 2 * mmsize]
    paddw           m8, m0
    pmaddubsw       m0, m12, [r5 + 1 * mmsize]
    paddw           m10, m0
    pmaddubsw       m12, [r5]
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1
    pmaddubsw       m1, m13, [r5 + 2 * mmsize]
    paddw           m9, m1
    pmaddubsw       m1, m13, [r5 + 1 * mmsize]
    paddw           m11, m1
    pmaddubsw       m13, [r5]

%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r8 + r3 * 2], m6
    movu            [r8 + r6], m7
%endif

    lea             r8, [r8 + r3 * 4]

    movu            xm1, [r7 + r4]                  ; m1 = row 15
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m2, m0, [r5 + 3 * mmsize]
    paddw           m8, m2
    pmaddubsw       m2, m0, [r5 + 2 * mmsize]
    paddw           m10, m2
    pmaddubsw       m2, m0, [r5 + 1 * mmsize]
    paddw           m12, m2
    pmaddubsw       m0, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm2, [r7]                       ; m2 = row 16
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 3 * mmsize]
    paddw           m9, m3
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m11, m3
    pmaddubsw       m3, m1, [r5 + 1 * mmsize]
    paddw           m13, m3
    pmaddubsw       m1, [r5]
    movu            xm3, [r7 + r1]                  ; m3 = row 17
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 3 * mmsize]
    paddw           m10, m4
    pmaddubsw       m4, m2, [r5 + 2 * mmsize]
    paddw           m12, m4
    pmaddubsw       m2, [r5 + 1 * mmsize]
    paddw           m0, m2
    movu            xm4, [r7 + r1 * 2]              ; m4 = row 18
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 3 * mmsize]
    paddw           m11, m5
    pmaddubsw       m5, m3, [r5 + 2 * mmsize]
    paddw           m13, m5
    pmaddubsw       m3, [r5 + 1 * mmsize]
    paddw           m1, m3
    movu            xm5, [r7 + r4]                  ; m5 = row 19
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 3 * mmsize]
    paddw           m12, m6
    pmaddubsw       m4, [r5 + 2 * mmsize]
    paddw           m0, m4
    lea             r7, [r7 + r1 * 4]
    movu            xm6, [r7]                       ; m6 = row 20
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 3 * mmsize]
    paddw           m13, m7
    pmaddubsw       m5, [r5 + 2 * mmsize]
    paddw           m1, m5
    movu            xm7, [r7 + r1]                  ; m7 = row 21
    punpckhbw       xm2, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddubsw       m6, [r5 + 3 * mmsize]
    paddw           m0, m6
    movu            xm2, [r7 + r1 * 2]              ; m2 = row 22
    punpckhbw       xm3, xm7, xm2
    punpcklbw       xm7, xm2
    vinserti128     m7, m7, xm3, 1
    pmaddubsw       m7, [r5 + 3 * mmsize]
    paddw           m1, m7

%ifidn %1,pp
    pmulhrsw        m8, m14                         ; m8 = word: row 8
    pmulhrsw        m9, m14                         ; m9 = word: row 9
    pmulhrsw        m10, m14                        ; m10 = word: row 10
    pmulhrsw        m11, m14                        ; m11 = word: row 11
    pmulhrsw        m12, m14                        ; m12 = word: row 12
    pmulhrsw        m13, m14                        ; m13 = word: row 13
    pmulhrsw        m0, m14                         ; m0 = word: row 14
    pmulhrsw        m1, m14                         ; m1 = word: row 15
    packuswb        m8, m9
    packuswb        m10, m11
    packuswb        m12, m13
    packuswb        m0, m1
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm1, m0, 1
    movu            [r8], xm8
    movu            [r8 + r3], xm9
    movu            [r8 + r3 * 2], xm10
    movu            [r8 + r6], xm11
    lea             r8, [r8 + r3 * 4]
    movu            [r8], xm12
    movu            [r8 + r3], xm13
    movu            [r8 + r3 * 2], xm0
    movu            [r8 + r6], xm1
%else
    psubw           m8, m14                         ; m8 = word: row 8
    psubw           m9, m14                         ; m9 = word: row 9
    psubw           m10, m14                        ; m10 = word: row 10
    psubw           m11, m14                        ; m11 = word: row 11
    psubw           m12, m14                        ; m12 = word: row 12
    psubw           m13, m14                        ; m13 = word: row 13
    psubw           m0, m14                         ; m0 = word: row 14
    psubw           m1, m14                         ; m1 = word: row 15
    movu            [r8], m8
    movu            [r8 + r3], m9
    movu            [r8 + r3 * 2], m10
    movu            [r8 + r6], m11
    lea             r8, [r8 + r3 * 4]
    movu            [r8], m12
    movu            [r8 + r3], m13
    movu            [r8 + r3 * 2], m0
    movu            [r8 + r6], m1
%endif
%endmacro

%macro PROCESS_LUMA_AVX2_W16_8R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhbw       xm2, xm0, xm1
    punpcklbw       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddubsw       m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhbw       xm3, xm1, xm2
    punpcklbw       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhbw       xm4, xm2, xm3
    punpcklbw       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddubsw       m4, m2, [r5 + 1 * mmsize]
    paddw           m0, m4
    pmaddubsw       m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhbw       xm5, xm3, xm4
    punpcklbw       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddubsw       m5, m3, [r5 + 1 * mmsize]
    paddw           m1, m5
    pmaddubsw       m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhbw       xm6, xm4, xm5
    punpcklbw       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddubsw       m6, m4, [r5 + 2 * mmsize]
    paddw           m0, m6
    pmaddubsw       m6, m4, [r5 + 1 * mmsize]
    paddw           m2, m6
    pmaddubsw       m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhbw       xm7, xm5, xm6
    punpcklbw       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddubsw       m7, m5, [r5 + 2 * mmsize]
    paddw           m1, m7
    pmaddubsw       m7, m5, [r5 + 1 * mmsize]
    paddw           m3, m7
    pmaddubsw       m5, [r5]
    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhbw       xm8, xm6, xm7
    punpcklbw       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddubsw       m8, m6, [r5 + 3 * mmsize]
    paddw           m0, m8
    pmaddubsw       m8, m6, [r5 + 2 * mmsize]
    paddw           m2, m8
    pmaddubsw       m8, m6, [r5 + 1 * mmsize]
    paddw           m4, m8
    pmaddubsw       m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhbw       xm9, xm7, xm8
    punpcklbw       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddubsw       m9, m7, [r5 + 3 * mmsize]
    paddw           m1, m9
    pmaddubsw       m9, m7, [r5 + 2 * mmsize]
    paddw           m3, m9
    pmaddubsw       m9, m7, [r5 + 1 * mmsize]
    paddw           m5, m9
    pmaddubsw       m7, [r5]
    movu            xm9, [r7 + r1]                  ; m9 = row 9
    punpckhbw       xm10, xm8, xm9
    punpcklbw       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddubsw       m10, m8, [r5 + 3 * mmsize]
    paddw           m2, m10
    pmaddubsw       m10, m8, [r5 + 2 * mmsize]
    paddw           m4, m10
    pmaddubsw       m10, m8, [r5 + 1 * mmsize]
    paddw           m6, m10
    movu            xm10, [r7 + r1 * 2]             ; m10 = row 10
    punpckhbw       xm11, xm9, xm10
    punpcklbw       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddubsw       m11, m9, [r5 + 3 * mmsize]
    paddw           m3, m11
    pmaddubsw       m11, m9, [r5 + 2 * mmsize]
    paddw           m5, m11
    pmaddubsw       m11, m9, [r5 + 1 * mmsize]
    paddw           m7, m11
    movu            xm11, [r7 + r4]                 ; m11 = row 11
    punpckhbw       xm12, xm10, xm11
    punpcklbw       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddubsw       m12, m10, [r5 + 3 * mmsize]
    paddw           m4, m12
    pmaddubsw       m12, m10, [r5 + 2 * mmsize]
    paddw           m6, m12
    lea             r7, [r7 + r1 * 4]
    movu            xm12, [r7]                      ; m12 = row 12
    punpckhbw       xm13, xm11, xm12
    punpcklbw       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddubsw       m13, m11, [r5 + 3 * mmsize]
    paddw           m5, m13
    pmaddubsw       m13, m11, [r5 + 2 * mmsize]
    paddw           m7, m13

%ifidn %1,pp
    pmulhrsw        m0, m14                         ; m0 = word: row 0
    pmulhrsw        m1, m14                         ; m1 = word: row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2
    pmulhrsw        m3, m14                         ; m3 = word: row 3
    pmulhrsw        m4, m14                         ; m4 = word: row 4
    pmulhrsw        m5, m14                         ; m5 = word: row 5
    packuswb        m0, m1
    packuswb        m2, m3
    packuswb        m4, m5
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    lea             r8, [r2 + r3 * 4]
    movu            [r8], xm4
    movu            [r8 + r3], xm5
%else
    psubw           m0, m14                         ; m0 = word: row 0
    psubw           m1, m14                         ; m1 = word: row 1
    psubw           m2, m14                         ; m2 = word: row 2
    psubw           m3, m14                         ; m3 = word: row 3
    psubw           m4, m14                         ; m4 = word: row 4
    psubw           m5, m14                         ; m5 = word: row 5
    movu            [r2], m0
    movu            [r2 + r3], m1
    movu            [r2 + r3 * 2], m2
    movu            [r2 + r6], m3
    lea             r8, [r2 + r3 * 4]
    movu            [r8], m4
    movu            [r8 + r3], m5
%endif

    movu            xm13, [r7 + r1]                 ; m13 = row 13
    punpckhbw       xm0, xm12, xm13
    punpcklbw       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddubsw       m0, m12, [r5 + 3 * mmsize]
    paddw           m6, m0
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 14
    punpckhbw       xm1, xm13, xm0
    punpcklbw       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddubsw       m1, m13, [r5 + 3 * mmsize]
    paddw           m7, m1

%ifidn %1,pp
    pmulhrsw        m6, m14                         ; m6 = word: row 6
    pmulhrsw        m7, m14                         ; m7 = word: row 7
    packuswb        m6, m7
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m6, 1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%else
    psubw           m6, m14                         ; m6 = word: row 6
    psubw           m7, m14                         ; m7 = word: row 7
    movu            [r8 + r3 * 2], m6
    movu            [r8 + r6], m7
%endif
%endmacro

%macro FILTER_VER_LUMA_AVX2_24x32 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_24x32, 4, 11, 15
    mov             r4d, r4m
    shl             r4d, 7
%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,ps
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%else
    mova            m14, [pw_512]
%endif
    lea             r6, [r3 * 3]
    lea             r10, [r1 * 4]
    mov             r9d, 2
.loopH:
    PROCESS_LUMA_AVX2_W16_16R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    movq            xm1, [r0]                       ; m1 = row 0
    movq            xm2, [r0 + r1]                  ; m2 = row 1
    punpcklbw       xm1, xm2
    movq            xm3, [r0 + r1 * 2]              ; m3 = row 2
    punpcklbw       xm2, xm3
    vinserti128     m5, m1, xm2, 1
    pmaddubsw       m5, [r5]
    movq            xm4, [r0 + r4]                  ; m4 = row 3
    punpcklbw       xm3, xm4
    lea             r7, [r0 + r1 * 4]
    movq            xm1, [r7]                       ; m1 = row 4
    punpcklbw       xm4, xm1
    vinserti128     m2, m3, xm4, 1
    pmaddubsw       m0, m2, [r5 + 1 * mmsize]
    paddw           m5, m0
    pmaddubsw       m2, [r5]
    movq            xm3, [r7 + r1]                  ; m3 = row 5
    punpcklbw       xm1, xm3
    movq            xm4, [r7 + r1 * 2]              ; m4 = row 6
    punpcklbw       xm3, xm4
    vinserti128     m1, m1, xm3, 1
    pmaddubsw       m3, m1, [r5 + 2 * mmsize]
    paddw           m5, m3
    pmaddubsw       m0, m1, [r5 + 1 * mmsize]
    paddw           m2, m0
    pmaddubsw       m1, [r5]
    movq            xm3, [r7 + r4]                  ; m3 = row 7
    punpcklbw       xm4, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm0, [r7]                       ; m0 = row 8
    punpcklbw       xm3, xm0
    vinserti128     m4, m4, xm3, 1
    pmaddubsw       m3, m4, [r5 + 3 * mmsize]
    paddw           m5, m3
    pmaddubsw       m3, m4, [r5 + 2 * mmsize]
    paddw           m2, m3
    pmaddubsw       m3, m4, [r5 + 1 * mmsize]
    paddw           m1, m3
    pmaddubsw       m4, [r5]
    movq            xm3, [r7 + r1]                  ; m3 = row 9
    punpcklbw       xm0, xm3
    movq            xm6, [r7 + r1 * 2]              ; m6 = row 10
    punpcklbw       xm3, xm6
    vinserti128     m0, m0, xm3, 1
    pmaddubsw       m3, m0, [r5 + 3 * mmsize]
    paddw           m2, m3
    pmaddubsw       m3, m0, [r5 + 2 * mmsize]
    paddw           m1, m3
    pmaddubsw       m3, m0, [r5 + 1 * mmsize]
    paddw           m4, m3
    pmaddubsw       m0, [r5]

    movq            xm3, [r7 + r4]                  ; m3 = row 11
    punpcklbw       xm6, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm7, [r7]                       ; m7 = row 12
    punpcklbw       xm3, xm7
    vinserti128     m6, m6, xm3, 1
    pmaddubsw       m3, m6, [r5 + 3 * mmsize]
    paddw           m1, m3
    pmaddubsw       m3, m6, [r5 + 2 * mmsize]
    paddw           m4, m3
    pmaddubsw       m3, m6, [r5 + 1 * mmsize]
    paddw           m0, m3
    pmaddubsw       m6, [r5]
    movq            xm3, [r7 + r1]                  ; m3 = row 13
    punpcklbw       xm7, xm3
    movq            xm8, [r7 + r1 * 2]              ; m8 = row 14
    punpcklbw       xm3, xm8
    vinserti128     m7, m7, xm3, 1
    pmaddubsw       m3, m7, [r5 + 3 * mmsize]
    paddw           m4, m3
    pmaddubsw       m3, m7, [r5 + 2 * mmsize]
    paddw           m0, m3
    pmaddubsw       m3, m7, [r5 + 1 * mmsize]
    paddw           m6, m3
    pmaddubsw       m7, [r5]
    movq            xm3, [r7 + r4]                  ; m3 = row 15
    punpcklbw       xm8, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm9, [r7]                       ; m9 = row 16
    punpcklbw       xm3, xm9
    vinserti128     m8, m8, xm3, 1
    pmaddubsw       m3, m8, [r5 + 3 * mmsize]
    paddw           m0, m3
    pmaddubsw       m3, m8, [r5 + 2 * mmsize]
    paddw           m6, m3
    pmaddubsw       m3, m8, [r5 + 1 * mmsize]
    paddw           m7, m3
    pmaddubsw       m8, [r5]
    movq            xm3, [r7 + r1]                  ; m3 = row 17
    punpcklbw       xm9, xm3
    movq            xm10, [r7 + r1 * 2]             ; m10 = row 18
    punpcklbw       xm3, xm10
    vinserti128     m9, m9, xm3, 1
    pmaddubsw       m3, m9, [r5 + 3 * mmsize]
    paddw           m6, m3
    pmaddubsw       m3, m9, [r5 + 2 * mmsize]
    paddw           m7, m3
    pmaddubsw       m3, m9, [r5 + 1 * mmsize]
    paddw           m8, m3
    movq            xm3, [r7 + r4]                  ; m3 = row 19
    punpcklbw       xm10, xm3
    lea             r7, [r7 + r1 * 4]
    movq            xm9, [r7]                       ; m9 = row 20
    punpcklbw       xm3, xm9
    vinserti128     m10, m10, xm3, 1
    pmaddubsw       m3, m10, [r5 + 3 * mmsize]
    paddw           m7, m3
    pmaddubsw       m3, m10, [r5 + 2 * mmsize]
    paddw           m8, m3
    movq            xm3, [r7 + r1]                  ; m3 = row 21
    punpcklbw       xm9, xm3
    movq            xm10, [r7 + r1 * 2]             ; m10 = row 22
    punpcklbw       xm3, xm10
    vinserti128     m9, m9, xm3, 1
    pmaddubsw       m3, m9, [r5 + 3 * mmsize]
    paddw           m8, m3
%ifidn %1,pp
    pmulhrsw        m5, m14                         ; m5 = word: row 0, row 1
    pmulhrsw        m2, m14                         ; m2 = word: row 2, row 3
    pmulhrsw        m1, m14                         ; m1 = word: row 4, row 5
    pmulhrsw        m4, m14                         ; m4 = word: row 6, row 7
    pmulhrsw        m0, m14                         ; m0 = word: row 8, row 9
    pmulhrsw        m6, m14                         ; m6 = word: row 10, row 11
    pmulhrsw        m7, m14                         ; m7 = word: row 12, row 13
    pmulhrsw        m8, m14                         ; m8 = word: row 14, row 15
    packuswb        m5, m2
    packuswb        m1, m4
    packuswb        m0, m6
    packuswb        m7, m8
    vextracti128    xm2, m5, 1
    vextracti128    xm4, m1, 1
    vextracti128    xm6, m0, 1
    vextracti128    xm8, m7, 1
    movq            [r2], xm5
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm2
    lea             r8, [r2 + r3 * 4]
    movq            [r8], xm1
    movq            [r8 + r3], xm4
    movhps          [r8 + r3 * 2], xm1
    movhps          [r8 + r6], xm4
    lea             r8, [r8 + r3 * 4]
    movq            [r8], xm0
    movq            [r8 + r3], xm6
    movhps          [r8 + r3 * 2], xm0
    movhps          [r8 + r6], xm6
    lea             r8, [r8 + r3 * 4]
    movq            [r8], xm7
    movq            [r8 + r3], xm8
    movhps          [r8 + r3 * 2], xm7
    movhps          [r8 + r6], xm8
%else
    psubw           m5, m14                         ; m5 = word: row 0, row 1
    psubw           m2, m14                         ; m2 = word: row 2, row 3
    psubw           m1, m14                         ; m1 = word: row 4, row 5
    psubw           m4, m14                         ; m4 = word: row 6, row 7
    psubw           m0, m14                         ; m0 = word: row 8, row 9
    psubw           m6, m14                         ; m6 = word: row 10, row 11
    psubw           m7, m14                         ; m7 = word: row 12, row 13
    psubw           m8, m14                         ; m8 = word: row 14, row 15
    vextracti128    xm3, m5, 1
    movu            [r2], xm5
    movu            [r2 + r3], xm3
    vextracti128    xm3, m2, 1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    vextracti128    xm3, m1, 1
    lea             r8, [r2 + r3 * 4]
    movu            [r8], xm1
    movu            [r8 + r3], xm3
    vextracti128    xm3, m4, 1
    movu            [r8 + r3 * 2], xm4
    movu            [r8 + r6], xm3
    vextracti128    xm3, m0, 1
    lea             r8, [r8 + r3 * 4]
    movu            [r8], xm0
    movu            [r8 + r3], xm3
    vextracti128    xm3, m6, 1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm3
    vextracti128    xm3, m7, 1
    lea             r8, [r8 + r3 * 4]
    movu            [r8], xm7
    movu            [r8 + r3], xm3
    vextracti128    xm3, m8, 1
    movu            [r8 + r3 * 2], xm8
    movu            [r8 + r6], xm3
%endif
    sub             r7, r10
    lea             r0, [r7 - 16]
%ifidn %1,pp
    lea             r2, [r8 + r3 * 4 - 16]
%else
    lea             r2, [r8 + r3 * 4 - 32]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_24x32 pp
FILTER_VER_LUMA_AVX2_24x32 ps

%macro FILTER_VER_LUMA_AVX2_32xN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%3_%1x%2, 4, 12, 15
    mov             r4d, r4m
    shl             r4d, 7
%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %3,ps
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%else
    mova            m14, [pw_512]
%endif
    lea             r6, [r3 * 3]
    lea             r11, [r1 * 4]
    mov             r9d, %2 / 16
.loopH:
    mov             r10d, %1 / 16
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %3
%ifidn %3,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r10d
    jnz             .loopW
    sub             r7, r11
    lea             r0, [r7 - 16]
%ifidn %3,pp
    lea             r2, [r8 + r3 * 4 - 16]
%else
    lea             r2, [r8 + r3 * 4 - 32]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_32xN 32, 32, pp
FILTER_VER_LUMA_AVX2_32xN 32, 64, pp
FILTER_VER_LUMA_AVX2_32xN 32, 32, ps
FILTER_VER_LUMA_AVX2_32xN 32, 64, ps

%macro FILTER_VER_LUMA_AVX2_32x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_32x16, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,ps
    add             r3d, r3d
    vbroadcasti128  m14, [pw_2000]
%else
    mova            m14, [pw_512]
%endif
    lea             r6, [r3 * 3]
    mov             r9d, 2
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_32x16 pp
FILTER_VER_LUMA_AVX2_32x16 ps
 
%macro FILTER_VER_LUMA_AVX2_32x24 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_32x24, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7
%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,ps
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    vbroadcasti128  m14, [pw_2000]
%endif
    mov             r9d, 2
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    lea             r9, [r1 * 4]
    sub             r7, r9
    lea             r0, [r7 - 16]
%ifidn %1,pp
    lea             r2, [r8 + r3 * 4 - 16]
%else
    lea             r2, [r8 + r3 * 4 - 32]
%endif
    mov             r9d, 2
.loop:
    PROCESS_LUMA_AVX2_W16_8R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r9d
    jnz             .loop
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_32x24 pp
FILTER_VER_LUMA_AVX2_32x24 ps

%macro FILTER_VER_LUMA_AVX2_32x8 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_32x8, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,ps
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
%ifidn %1,pp
    mova            m14, [pw_512]
%else
    vbroadcasti128  m14, [pw_2000]
%endif
    mov             r9d, 2
.loopW:
    PROCESS_LUMA_AVX2_W16_8R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_32x8 pp
FILTER_VER_LUMA_AVX2_32x8 ps

%macro FILTER_VER_LUMA_AVX2_48x64 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_48x64, 4, 12, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %1,ps
    add             r3d, r3d
%endif

    lea             r6, [r3 * 3]
    lea             r11, [r1 * 4]

%ifidn %1,pp
    mova            m14, [pw_512]
%else
    vbroadcasti128  m14, [pw_2000]
%endif

    mov             r9d, 4
.loopH:
    mov             r10d, 3
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r10d
    jnz             .loopW
    sub             r7, r11
    lea             r0, [r7 - 32]
%ifidn %1,pp
    lea             r2, [r8 + r3 * 4 - 32]
%else
    lea             r2, [r8 + r3 * 4 - 64]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_48x64 pp
FILTER_VER_LUMA_AVX2_48x64 ps

%macro FILTER_VER_LUMA_AVX2_64xN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%3_%1x%2, 4, 12, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %3,ps
    add             r3d, r3d
%endif

    lea             r6, [r3 * 3]
    lea             r11, [r1 * 4]

%ifidn %3,pp
    mova            m14, [pw_512]
%else
    vbroadcasti128  m14, [pw_2000]
%endif

    mov             r9d, %2 / 16
.loopH:
    mov             r10d, %1 / 16
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %3
%ifidn %3,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r10d
    jnz             .loopW
    sub             r7, r11
    lea             r0, [r7 - 48]
%ifidn %3,pp
    lea             r2, [r8 + r3 * 4 - 48]
%else
    lea             r2, [r8 + r3 * 4 - 96]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_64xN 64, 32, pp
FILTER_VER_LUMA_AVX2_64xN 64, 48, pp
FILTER_VER_LUMA_AVX2_64xN 64, 64, pp
FILTER_VER_LUMA_AVX2_64xN 64, 32, ps
FILTER_VER_LUMA_AVX2_64xN 64, 48, ps
FILTER_VER_LUMA_AVX2_64xN 64, 64, ps

%macro FILTER_VER_LUMA_AVX2_64x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_64x16, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7

%ifdef PIC
    lea             r5, [tab_LumaCoeffVer_32]
    add             r5, r4
%else
    lea             r5, [tab_LumaCoeffVer_32 + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %1,ps
    add             r3d, r3d
%endif

    lea             r6, [r3 * 3]

%ifidn %1,pp
    mova            m14, [pw_512]
%else
    vbroadcasti128  m14, [pw_2000]
%endif

    mov             r9d, 4
.loopW:
    PROCESS_LUMA_AVX2_W16_16R %1
%ifidn %1,pp
    add             r2, 16
%else
    add             r2, 32
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_64x16 pp
FILTER_VER_LUMA_AVX2_64x16 ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA 3
INIT_XMM sse4
cglobal interp_8tap_vert_%3_%1x%2, 5, 7, 8 ,0-gprsize
    lea       r5, [3 * r1]
    sub       r0, r5
    shl       r4d, 6
%ifidn %3,ps
    add       r3d, r3d
%endif

%ifdef PIC
    lea       r5, [tab_LumaCoeffVer]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffVer + r4]
%endif

%ifidn %3,pp
    mova      m3, [pw_512]
%else
    mova      m3, [pw_2000]
%endif
    mov       dword [rsp], %2/4

.loopH:
    mov       r4d, (%1/8)
.loopW:
    PROCESS_LUMA_W8_4R
%ifidn %3,pp
    pmulhrsw  m7, m3
    pmulhrsw  m6, m3
    pmulhrsw  m5, m3
    pmulhrsw  m4, m3

    packuswb  m7, m6
    packuswb  m5, m4

    movlps    [r2], m7
    movhps    [r2 + r3], m7
    lea       r5, [r2 + 2 * r3]
    movlps    [r5], m5
    movhps    [r5 + r3], m5
%else
    psubw     m7, m3
    psubw     m6, m3
    psubw     m5, m3
    psubw     m4, m3

    movu      [r2], m7
    movu      [r2 + r3], m6
    lea       r5, [r2 + 2 * r3]
    movu      [r5], m5
    movu      [r5 + r3], m4
%endif

    lea       r5, [8 * r1 - 8]
    sub       r0, r5
%ifidn %3,pp
    add       r2, 8
%else
    add       r2, 16
%endif
    dec       r4d
    jnz       .loopW

    lea       r0, [r0 + 4 * r1 - %1]
%ifidn %3,pp
    lea       r2, [r2 + 4 * r3 - %1]
%else
    lea       r2, [r2 + 4 * r3 - 2 * %1]
%endif

    dec       dword [rsp]
    jnz       .loopH

    RET
%endmacro

FILTER_VER_LUMA 16, 4, pp
FILTER_VER_LUMA 16, 8, pp
FILTER_VER_LUMA 16, 12, pp
FILTER_VER_LUMA 16, 16, pp
FILTER_VER_LUMA 16, 32, pp
FILTER_VER_LUMA 16, 64, pp
FILTER_VER_LUMA 24, 32, pp
FILTER_VER_LUMA 32, 8, pp
FILTER_VER_LUMA 32, 16, pp
FILTER_VER_LUMA 32, 24, pp
FILTER_VER_LUMA 32, 32, pp
FILTER_VER_LUMA 32, 64, pp
FILTER_VER_LUMA 48, 64, pp
FILTER_VER_LUMA 64, 16, pp
FILTER_VER_LUMA 64, 32, pp
FILTER_VER_LUMA 64, 48, pp
FILTER_VER_LUMA 64, 64, pp

FILTER_VER_LUMA 16, 4, ps
FILTER_VER_LUMA 16, 8, ps
FILTER_VER_LUMA 16, 12, ps
FILTER_VER_LUMA 16, 16, ps
FILTER_VER_LUMA 16, 32, ps
FILTER_VER_LUMA 16, 64, ps
FILTER_VER_LUMA 24, 32, ps
FILTER_VER_LUMA 32, 8, ps
FILTER_VER_LUMA 32, 16, ps
FILTER_VER_LUMA 32, 24, ps
FILTER_VER_LUMA 32, 32, ps
FILTER_VER_LUMA 32, 64, ps
FILTER_VER_LUMA 48, 64, ps
FILTER_VER_LUMA 64, 16, ps
FILTER_VER_LUMA 64, 32, ps
FILTER_VER_LUMA 64, 48, ps
FILTER_VER_LUMA 64, 64, ps

%macro PROCESS_LUMA_SP_W4_4R 0
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r6 + 0 *16]                ;m0=[0+1]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m1, m4                          ;m1=[1 2]
    pmaddwd    m1, [r6 + 0 *16]                ;m1=[1+2]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[2 3]
    pmaddwd    m2, m4, [r6 + 0 *16]            ;m2=[2+3]  Row3
    pmaddwd    m4, [r6 + 1 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[3 4]
    pmaddwd    m3, m5, [r6 + 0 *16]            ;m3=[3+4]  Row4
    pmaddwd    m5, [r6 + 1 * 16]
    paddd      m1, m5                          ;m1 = [1+2+3+4]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[4 5]
    pmaddwd    m6, m4, [r6 + 1 * 16]
    paddd      m2, m6                          ;m2=[2+3+4+5]  Row3
    pmaddwd    m4, [r6 + 2 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3+4+5]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[5 6]
    pmaddwd    m6, m5, [r6 + 1 * 16]
    paddd      m3, m6                          ;m3=[3+4+5+6]  Row4
    pmaddwd    m5, [r6 + 2 * 16]
    paddd      m1, m5                          ;m1=[1+2+3+4+5+6]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[6 7]
    pmaddwd    m6, m4, [r6 + 2 * 16]
    paddd      m2, m6                          ;m2=[2+3+4+5+6+7]  Row3
    pmaddwd    m4, [r6 + 3 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3+4+5+6+7]  Row1 end

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[7 8]
    pmaddwd    m6, m5, [r6 + 2 * 16]
    paddd      m3, m6                          ;m3=[3+4+5+6+7+8]  Row4
    pmaddwd    m5, [r6 + 3 * 16]
    paddd      m1, m5                          ;m1=[1+2+3+4+5+6+7+8]  Row2 end

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[8 9]
    pmaddwd    m4, [r6 + 3 * 16]
    paddd      m2, m4                          ;m2=[2+3+4+5+6+7+8+9]  Row3 end

    movq       m4, [r0 + 2 * r1]
    punpcklwd  m5, m4                          ;m5=[9 10]
    pmaddwd    m5, [r6 + 3 * 16]
    paddd      m3, m5                          ;m3=[3+4+5+6+7+8+9+10]  Row4 end
%endmacro

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_sp_%1x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_SP 2
INIT_XMM sse4
cglobal interp_8tap_vert_sp_%1x%2, 5, 7, 8 ,0-gprsize

    add       r1d, r1d
    lea       r5, [r1 + 2 * r1]
    sub       r0, r5
    shl       r4d, 6

%ifdef PIC
    lea       r5, [tab_LumaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffV + r4]
%endif

    mova      m7, [tab_c_526336]

    mov       dword [rsp], %2/4
.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_LUMA_SP_W4_4R

    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7

    psrad     m0, 12
    psrad     m1, 12
    psrad     m2, 12
    psrad     m3, 12

    packssdw  m0, m1
    packssdw  m2, m3

    packuswb  m0, m2

    movd      [r2], m0
    pextrd    [r2 + r3], m0, 1
    lea       r5, [r2 + 2 * r3]
    pextrd    [r5], m0, 2
    pextrd    [r5 + r3], m0, 3

    lea       r5, [8 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 4

    dec       r4d
    jnz       .loopW

    lea       r0, [r0 + 4 * r1 - 2 * %1]
    lea       r2, [r2 + 4 * r3 - %1]

    dec       dword [rsp]
    jnz       .loopH

    RET
%endmacro

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_sp_%1x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
    FILTER_VER_LUMA_SP 4, 4
    FILTER_VER_LUMA_SP 8, 8
    FILTER_VER_LUMA_SP 8, 4
    FILTER_VER_LUMA_SP 4, 8
    FILTER_VER_LUMA_SP 16, 16
    FILTER_VER_LUMA_SP 16, 8
    FILTER_VER_LUMA_SP 8, 16
    FILTER_VER_LUMA_SP 16, 12
    FILTER_VER_LUMA_SP 12, 16
    FILTER_VER_LUMA_SP 16, 4
    FILTER_VER_LUMA_SP 4, 16
    FILTER_VER_LUMA_SP 32, 32
    FILTER_VER_LUMA_SP 32, 16
    FILTER_VER_LUMA_SP 16, 32
    FILTER_VER_LUMA_SP 32, 24
    FILTER_VER_LUMA_SP 24, 32
    FILTER_VER_LUMA_SP 32, 8
    FILTER_VER_LUMA_SP 8, 32
    FILTER_VER_LUMA_SP 64, 64
    FILTER_VER_LUMA_SP 64, 32
    FILTER_VER_LUMA_SP 32, 64
    FILTER_VER_LUMA_SP 64, 48
    FILTER_VER_LUMA_SP 48, 64
    FILTER_VER_LUMA_SP 64, 16
    FILTER_VER_LUMA_SP 16, 64

; TODO: combin of U and V is more performance, but need more register
; TODO: use two path for height alignment to 4 and otherwise may improvement 10% performance, but code is more complex, so I disable it
INIT_XMM ssse3
cglobal chroma_p2s, 3, 7, 4

    ; load width and height
    mov         r3d, r3m
    mov         r4d, r4m

    ; load constant
    mova        m2, [pb_128]
    mova        m3, [tab_c_64_n64]

.loopH:

    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5]

    movh        m0, [r6]
    punpcklbw   m0, m2
    pmaddubsw   m0, m3

    movh        m1, [r6 + r1]
    punpcklbw   m1, m2
    pmaddubsw   m1, m3

    add         r5d, 8
    cmp         r5d, r3d
    lea         r6, [r2 + r5 * 2]
    jg          .width4
    movu        [r6 + FENC_STRIDE / 2 * 0 - 16], m0
    movu        [r6 + FENC_STRIDE / 2 * 2 - 16], m1
    je          .nextH
    jmp         .loopW

.width4:
    test        r3d, 4
    jz          .width2
    test        r3d, 2
    movh        [r6 + FENC_STRIDE / 2 * 0 - 16], m0
    movh        [r6 + FENC_STRIDE / 2 * 2 - 16], m1
    lea         r6, [r6 + 8]
    pshufd      m0, m0, 2
    pshufd      m1, m1, 2
    jz          .nextH

.width2:
    movd        [r6 + FENC_STRIDE / 2 * 0 - 16], m0
    movd        [r6 + FENC_STRIDE / 2 * 2 - 16], m1

.nextH:
    lea         r0, [r0 + r1 * 2]
    add         r2, FENC_STRIDE / 2 * 4

    sub         r4d, 2
    jnz         .loopH

    RET

%macro PROCESS_CHROMA_SP_W4_4R 0
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r6 + 0 *16]                ;m0=[0+1]         Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m1, m4                          ;m1=[1 2]
    pmaddwd    m1, [r6 + 0 *16]                ;m1=[1+2]         Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[2 3]
    pmaddwd    m2, m4, [r6 + 0 *16]            ;m2=[2+3]         Row3
    pmaddwd    m4, [r6 + 1 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3]     Row1 done

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[3 4]
    pmaddwd    m3, m5, [r6 + 0 *16]            ;m3=[3+4]         Row4
    pmaddwd    m5, [r6 + 1 * 16]
    paddd      m1, m5                          ;m1 = [1+2+3+4]   Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[4 5]
    pmaddwd    m4, [r6 + 1 * 16]
    paddd      m2, m4                          ;m2=[2+3+4+5]     Row3

    movq       m4, [r0 + 2 * r1]
    punpcklwd  m5, m4                          ;m5=[5 6]
    pmaddwd    m5, [r6 + 1 * 16]
    paddd      m3, m5                          ;m3=[3+4+5+6]     Row4
%endmacro

;--------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_sp_%1x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SP 2
INIT_XMM sse4
cglobal interp_4tap_vert_sp_%1x%2, 5, 7, 7 ,0-gprsize

    add       r1d, r1d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_ChromaCoeffV + r4]
%endif

    mova      m6, [tab_c_526336]

    mov       dword [rsp], %2/4

.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_CHROMA_SP_W4_4R

    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6

    psrad     m0, 12
    psrad     m1, 12
    psrad     m2, 12
    psrad     m3, 12

    packssdw  m0, m1
    packssdw  m2, m3

    packuswb  m0, m2

    movd      [r2], m0
    pextrd    [r2 + r3], m0, 1
    lea       r5, [r2 + 2 * r3]
    pextrd    [r5], m0, 2
    pextrd    [r5 + r3], m0, 3

    lea       r5, [4 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 4

    dec       r4d
    jnz       .loopW

    lea       r0, [r0 + 4 * r1 - 2 * %1]
    lea       r2, [r2 + 4 * r3 - %1]

    dec       dword [rsp]
    jnz       .loopH

    RET
%endmacro

    FILTER_VER_CHROMA_SP 4, 4
    FILTER_VER_CHROMA_SP 4, 8
    FILTER_VER_CHROMA_SP 16, 16
    FILTER_VER_CHROMA_SP 16, 8
    FILTER_VER_CHROMA_SP 16, 12
    FILTER_VER_CHROMA_SP 12, 16
    FILTER_VER_CHROMA_SP 16, 4
    FILTER_VER_CHROMA_SP 4, 16
    FILTER_VER_CHROMA_SP 32, 32
    FILTER_VER_CHROMA_SP 32, 16
    FILTER_VER_CHROMA_SP 16, 32
    FILTER_VER_CHROMA_SP 32, 24
    FILTER_VER_CHROMA_SP 24, 32
    FILTER_VER_CHROMA_SP 32, 8

    FILTER_VER_CHROMA_SP 16, 24
    FILTER_VER_CHROMA_SP 16, 64
    FILTER_VER_CHROMA_SP 12, 32
    FILTER_VER_CHROMA_SP 4, 32
    FILTER_VER_CHROMA_SP 32, 64
    FILTER_VER_CHROMA_SP 32, 48
    FILTER_VER_CHROMA_SP 24, 64

    FILTER_VER_CHROMA_SP 64, 64
    FILTER_VER_CHROMA_SP 64, 32
    FILTER_VER_CHROMA_SP 64, 48
    FILTER_VER_CHROMA_SP 48, 64
    FILTER_VER_CHROMA_SP 64, 16


%macro PROCESS_CHROMA_SP_W2_4R 1
    movd       m0, [r0]
    movd       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]

    lea        r0, [r0 + 2 * r1]
    movd       m2, [r0]
    punpcklwd  m1, m2                          ;m1=[1 2]
    punpcklqdq m0, m1                          ;m0=[0 1 1 2]
    pmaddwd    m0, [%1 + 0 *16]                ;m0=[0+1 1+2] Row 1-2

    movd       m1, [r0 + r1]
    punpcklwd  m2, m1                          ;m2=[2 3]

    lea        r0, [r0 + 2 * r1]
    movd       m3, [r0]
    punpcklwd  m1, m3                          ;m2=[3 4]
    punpcklqdq m2, m1                          ;m2=[2 3 3 4]

    pmaddwd    m4, m2, [%1 + 1 * 16]           ;m4=[2+3 3+4] Row 1-2
    pmaddwd    m2, [%1 + 0 * 16]               ;m2=[2+3 3+4] Row 3-4
    paddd      m0, m4                          ;m0=[0+1+2+3 1+2+3+4] Row 1-2

    movd       m1, [r0 + r1]
    punpcklwd  m3, m1                          ;m3=[4 5]

    movd       m4, [r0 + 2 * r1]
    punpcklwd  m1, m4                          ;m1=[5 6]
    punpcklqdq m3, m1                          ;m2=[4 5 5 6]
    pmaddwd    m3, [%1 + 1 * 16]               ;m3=[4+5 5+6] Row 3-4
    paddd      m2, m3                          ;m2=[2+3+4+5 3+4+5+6] Row 3-4
%endmacro

;-------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_sp_%1x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SP_W2_4R 2
INIT_XMM sse4
cglobal interp_4tap_vert_sp_%1x%2, 5, 6, 6

    add       r1d, r1d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r5, [r5 + r4]
%else
    lea       r5, [tab_ChromaCoeffV + r4]
%endif

    mova      m5, [tab_c_526336]

    mov       r4d, (%2/4)

.loopH:
    PROCESS_CHROMA_SP_W2_4R r5

    paddd     m0, m5
    paddd     m2, m5

    psrad     m0, 12
    psrad     m2, 12

    packssdw  m0, m2
    packuswb  m0, m0

    pextrw    [r2], m0, 0
    pextrw    [r2 + r3], m0, 1
    lea       r2, [r2 + 2 * r3]
    pextrw    [r2], m0, 2
    pextrw    [r2 + r3], m0, 3

    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SP_W2_4R 2, 4
FILTER_VER_CHROMA_SP_W2_4R 2, 8

FILTER_VER_CHROMA_SP_W2_4R 2, 16

;--------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_sp_4x2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_sp_4x2, 5, 6, 5

    add        r1d, r1d
    sub        r0, r1
    shl        r4d, 5

%ifdef PIC
    lea        r5, [tab_ChromaCoeffV]
    lea        r5, [r5 + r4]
%else
    lea        r5, [tab_ChromaCoeffV + r4]
%endif

    mova       m4, [tab_c_526336]

    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r5 + 0 *16]                ;m0=[0+1]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m2, [r0]
    punpcklwd  m1, m2                          ;m1=[1 2]
    pmaddwd    m1, [r5 + 0 *16]                ;m1=[1+2]  Row2

    movq       m3, [r0 + r1]
    punpcklwd  m2, m3                          ;m4=[2 3]
    pmaddwd    m2, [r5 + 1 * 16]
    paddd      m0, m2                          ;m0=[0+1+2+3]  Row1 done
    paddd      m0, m4
    psrad      m0, 12

    movq       m2, [r0 + 2 * r1]
    punpcklwd  m3, m2                          ;m5=[3 4]
    pmaddwd    m3, [r5 + 1 * 16]
    paddd      m1, m3                          ;m1 = [1+2+3+4]  Row2 done
    paddd      m1, m4
    psrad      m1, 12

    packssdw   m0, m1
    packuswb   m0, m0

    movd       [r2], m0
    pextrd     [r2 + r3], m0, 1

    RET

;-------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_sp_6x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SP_W6_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_sp_6x%2, 5, 7, 7

    add       r1d, r1d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_ChromaCoeffV + r4]
%endif

    mova      m6, [tab_c_526336]

    mov       r4d, %2/4

.loopH:
    PROCESS_CHROMA_SP_W4_4R

    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6

    psrad     m0, 12
    psrad     m1, 12
    psrad     m2, 12
    psrad     m3, 12

    packssdw  m0, m1
    packssdw  m2, m3

    packuswb  m0, m2

    movd      [r2], m0
    pextrd    [r2 + r3], m0, 1
    lea       r5, [r2 + 2 * r3]
    pextrd    [r5], m0, 2
    pextrd    [r5 + r3], m0, 3

    lea       r5, [4 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 4

    PROCESS_CHROMA_SP_W2_4R r6

    paddd     m0, m6
    paddd     m2, m6

    psrad     m0, 12
    psrad     m2, 12

    packssdw  m0, m2
    packuswb  m0, m0

    pextrw    [r2], m0, 0
    pextrw    [r2 + r3], m0, 1
    lea       r2, [r2 + 2 * r3]
    pextrw    [r2], m0, 2
    pextrw    [r2 + r3], m0, 3

    sub       r0, 2 * 4
    lea       r2, [r2 + 2 * r3 - 4]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SP_W6_H4 6, 8

FILTER_VER_CHROMA_SP_W6_H4 6, 16

%macro PROCESS_CHROMA_SP_W8_2R 0
    movu       m1, [r0]
    movu       m3, [r0 + r1]
    punpcklwd  m0, m1, m3
    pmaddwd    m0, [r5 + 0 * 16]                ;m0 = [0l+1l]  Row1l
    punpckhwd  m1, m3
    pmaddwd    m1, [r5 + 0 * 16]                ;m1 = [0h+1h]  Row1h

    movu       m4, [r0 + 2 * r1]
    punpcklwd  m2, m3, m4
    pmaddwd    m2, [r5 + 0 * 16]                ;m2 = [1l+2l]  Row2l
    punpckhwd  m3, m4
    pmaddwd    m3, [r5 + 0 * 16]                ;m3 = [1h+2h]  Row2h

    lea        r0, [r0 + 2 * r1]
    movu       m5, [r0 + r1]
    punpcklwd  m6, m4, m5
    pmaddwd    m6, [r5 + 1 * 16]                ;m6 = [2l+3l]  Row1l
    paddd      m0, m6                           ;m0 = [0l+1l+2l+3l]  Row1l sum
    punpckhwd  m4, m5
    pmaddwd    m4, [r5 + 1 * 16]                ;m6 = [2h+3h]  Row1h
    paddd      m1, m4                           ;m1 = [0h+1h+2h+3h]  Row1h sum

    movu       m4, [r0 + 2 * r1]
    punpcklwd  m6, m5, m4
    pmaddwd    m6, [r5 + 1 * 16]                ;m6 = [3l+4l]  Row2l
    paddd      m2, m6                           ;m2 = [1l+2l+3l+4l]  Row2l sum
    punpckhwd  m5, m4
    pmaddwd    m5, [r5 + 1 * 16]                ;m1 = [3h+4h]  Row2h
    paddd      m3, m5                           ;m3 = [1h+2h+3h+4h]  Row2h sum
%endmacro

;--------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_sp_8x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SP_W8_H2 2
INIT_XMM sse2
cglobal interp_4tap_vert_sp_%1x%2, 5, 6, 8

    add       r1d, r1d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r5, [r5 + r4]
%else
    lea       r5, [tab_ChromaCoeffV + r4]
%endif

    mova      m7, [tab_c_526336]

    mov       r4d, %2/2
.loopH:
    PROCESS_CHROMA_SP_W8_2R

    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7

    psrad     m0, 12
    psrad     m1, 12
    psrad     m2, 12
    psrad     m3, 12

    packssdw  m0, m1
    packssdw  m2, m3

    packuswb  m0, m2

    movlps    [r2], m0
    movhps    [r2 + r3], m0

    lea       r2, [r2 + 2 * r3]

    dec r4d
    jnz .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SP_W8_H2 8, 2
FILTER_VER_CHROMA_SP_W8_H2 8, 4
FILTER_VER_CHROMA_SP_W8_H2 8, 6
FILTER_VER_CHROMA_SP_W8_H2 8, 8
FILTER_VER_CHROMA_SP_W8_H2 8, 16
FILTER_VER_CHROMA_SP_W8_H2 8, 32

FILTER_VER_CHROMA_SP_W8_H2 8, 12
FILTER_VER_CHROMA_SP_W8_H2 8, 64


;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_2x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro FILTER_HORIZ_CHROMA_2xN 2
INIT_XMM sse4
cglobal interp_4tap_horiz_ps_%1x%2, 4, 7, 4, src, srcstride, dst, dststride
%define coef2  m3
%define Tm0    m2
%define t1     m1
%define t0     m0

    dec        srcq
    mov        r4d, r4m
    add        dststrided, dststrided

%ifdef PIC
    lea        r6, [tab_ChromaCoeff]
    movd       coef2, [r6 + r4 * 4]
%else
    movd       coef2, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufd     coef2, coef2, 0
    mova       t1, [pw_2000]
    mova       Tm0, [tab_Tm]

    mov        r4d, %2
    cmp        r5m, byte 0
    je         .loopH
    sub        srcq, srcstrideq
    add        r4d, 3

.loopH:
    movh       t0, [srcq]
    pshufb     t0, t0, Tm0
    pmaddubsw  t0, coef2
    phaddw     t0, t0
    psubw      t0, t1
    movd       [dstq], t0

    lea        srcq, [srcq + srcstrideq]
    lea        dstq, [dstq + dststrideq]

    dec        r4d
    jnz        .loopH

    RET
%endmacro

FILTER_HORIZ_CHROMA_2xN 2, 4
FILTER_HORIZ_CHROMA_2xN 2, 8

FILTER_HORIZ_CHROMA_2xN 2, 16

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_4x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro FILTER_HORIZ_CHROMA_4xN 2
INIT_XMM sse4
cglobal interp_4tap_horiz_ps_%1x%2, 4, 7, 4, src, srcstride, dst, dststride
%define coef2  m3
%define Tm0    m2
%define t1     m1
%define t0     m0

    dec        srcq
    mov        r4d, r4m
    add        dststrided, dststrided

%ifdef PIC
    lea        r6, [tab_ChromaCoeff]
    movd       coef2, [r6 + r4 * 4]
%else
    movd       coef2, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufd     coef2, coef2, 0
    mova       t1, [pw_2000]
    mova       Tm0, [tab_Tm]

    mov        r4d, %2
    cmp        r5m, byte 0
    je         .loopH
    sub        srcq, srcstrideq
    add        r4d, 3

.loopH:
    movh       t0, [srcq]
    pshufb     t0, t0, Tm0
    pmaddubsw  t0, coef2
    phaddw     t0, t0
    psubw      t0, t1
    movlps     [dstq], t0

    lea        srcq, [srcq + srcstrideq]
    lea        dstq, [dstq + dststrideq]

    dec        r4d
    jnz        .loopH
    RET
%endmacro

FILTER_HORIZ_CHROMA_4xN 4, 2
FILTER_HORIZ_CHROMA_4xN 4, 4
FILTER_HORIZ_CHROMA_4xN 4, 8
FILTER_HORIZ_CHROMA_4xN 4, 16

FILTER_HORIZ_CHROMA_4xN 4, 32

%macro PROCESS_CHROMA_W6 3
    movu       %1, [srcq]
    pshufb     %2, %1, Tm0
    pmaddubsw  %2, coef2
    pshufb     %1, %1, Tm1
    pmaddubsw  %1, coef2
    phaddw     %2, %1
    psubw      %2, %3
    movh       [dstq], %2
    pshufd     %2, %2, 2
    movd       [dstq + 8], %2
%endmacro

%macro PROCESS_CHROMA_W12 3
    movu       %1, [srcq]
    pshufb     %2, %1, Tm0
    pmaddubsw  %2, coef2
    pshufb     %1, %1, Tm1
    pmaddubsw  %1, coef2
    phaddw     %2, %1
    psubw      %2, %3
    movu       [dstq], %2
    movu       %1, [srcq + 8]
    pshufb     %1, %1, Tm0
    pmaddubsw  %1, coef2
    phaddw     %1, %1
    psubw      %1, %3
    movh       [dstq + 16], %1
%endmacro

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_6x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro FILTER_HORIZ_CHROMA 2
INIT_XMM sse4
cglobal interp_4tap_horiz_ps_%1x%2, 4, 7, 6, src, srcstride, dst, dststride
%define coef2    m5
%define Tm0      m4
%define Tm1      m3
%define t2       m2
%define t1       m1
%define t0       m0

    dec     srcq
    mov     r4d, r4m
    add     dststrided, dststrided

%ifdef PIC
    lea     r6, [tab_ChromaCoeff]
    movd    coef2, [r6 + r4 * 4]
%else
    movd    coef2, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufd  coef2, coef2, 0
    mova    t2, [pw_2000]
    mova    Tm0, [tab_Tm]
    mova    Tm1, [tab_Tm + 16]

    mov     r4d, %2
    cmp     r5m, byte 0
    je      .loopH
    sub     srcq, srcstrideq
    add     r4d, 3

.loopH:
    PROCESS_CHROMA_W%1  t0, t1, t2
    add     srcq, srcstrideq
    add     dstq, dststrideq

    dec     r4d
    jnz     .loopH

    RET
%endmacro

FILTER_HORIZ_CHROMA 6, 8
FILTER_HORIZ_CHROMA 12, 16

FILTER_HORIZ_CHROMA 6, 16
FILTER_HORIZ_CHROMA 12, 32

%macro PROCESS_CHROMA_W8 3
    movu        %1, [srcq]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    psubw       %2, %3
    movu        [dstq], %2
%endmacro

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_8x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro FILTER_HORIZ_CHROMA_8xN 2
INIT_XMM sse4
cglobal interp_4tap_horiz_ps_%1x%2, 4, 7, 6, src, srcstride, dst, dststride
%define coef2    m5
%define Tm0      m4
%define Tm1      m3
%define t2       m2
%define t1       m1
%define t0       m0

    dec     srcq
    mov     r4d, r4m
    add     dststrided, dststrided

%ifdef PIC
    lea     r6, [tab_ChromaCoeff]
    movd    coef2, [r6 + r4 * 4]
%else
    movd    coef2, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufd  coef2, coef2, 0
    mova    t2, [pw_2000]
    mova    Tm0, [tab_Tm]
    mova    Tm1, [tab_Tm + 16]

    mov     r4d, %2
    cmp     r5m, byte 0
    je      .loopH
    sub     srcq, srcstrideq
    add     r4d, 3

.loopH:
    PROCESS_CHROMA_W8  t0, t1, t2
    add     srcq, srcstrideq
    add     dstq, dststrideq

    dec     r4d
    jnz     .loopH

    RET
%endmacro

FILTER_HORIZ_CHROMA_8xN 8, 2
FILTER_HORIZ_CHROMA_8xN 8, 4
FILTER_HORIZ_CHROMA_8xN 8, 6
FILTER_HORIZ_CHROMA_8xN 8, 8
FILTER_HORIZ_CHROMA_8xN 8, 16
FILTER_HORIZ_CHROMA_8xN 8, 32

FILTER_HORIZ_CHROMA_8xN 8, 12
FILTER_HORIZ_CHROMA_8xN 8, 64

%macro PROCESS_CHROMA_W16 4
    movu        %1, [srcq]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    psubw       %2, %3
    psubw       %4, %3
    movu        [dstq], %2
    movu        [dstq + 16], %4
%endmacro

%macro PROCESS_CHROMA_W24 4
    movu        %1, [srcq]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    psubw       %2, %3
    psubw       %4, %3
    movu        [dstq], %2
    movu        [dstq + 16], %4
    movu        %1, [srcq + 16]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    psubw       %2, %3
    movu        [dstq + 32], %2
%endmacro

%macro PROCESS_CHROMA_W32 4
    movu        %1, [srcq]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    psubw       %2, %3
    psubw       %4, %3
    movu        [dstq], %2
    movu        [dstq + 16], %4
    movu        %1, [srcq + 16]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + 24]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    psubw       %2, %3
    psubw       %4, %3
    movu        [dstq + 32], %2
    movu        [dstq + 48], %4
%endmacro

%macro PROCESS_CHROMA_W16o 5
    movu        %1, [srcq + %5]
    pshufb      %2, %1, Tm0
    pmaddubsw   %2, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %2, %1
    movu        %1, [srcq + %5 + 8]
    pshufb      %4, %1, Tm0
    pmaddubsw   %4, coef2
    pshufb      %1, %1, Tm1
    pmaddubsw   %1, coef2
    phaddw      %4, %1
    psubw       %2, %3
    psubw       %4, %3
    movu        [dstq + %5 * 2], %2
    movu        [dstq + %5 * 2 + 16], %4
%endmacro

%macro PROCESS_CHROMA_W48 4
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 0
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 16
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 32
%endmacro

%macro PROCESS_CHROMA_W64 4
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 0
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 16
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 32
    PROCESS_CHROMA_W16o %1, %2, %3, %4, 48
%endmacro

;------------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_%1x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;------------------------------------------------------------------------------------------------------------------------------
%macro FILTER_HORIZ_CHROMA_WxN 2
INIT_XMM sse4
cglobal interp_4tap_horiz_ps_%1x%2, 4, 7, 7, src, srcstride, dst, dststride
%define coef2    m6
%define Tm0      m5
%define Tm1      m4
%define t3       m3
%define t2       m2
%define t1       m1
%define t0       m0

    dec     srcq
    mov     r4d, r4m
    add     dststrided, dststrided

%ifdef PIC
    lea     r6, [tab_ChromaCoeff]
    movd    coef2, [r6 + r4 * 4]
%else
    movd    coef2, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufd  coef2, coef2, 0
    mova    t2, [pw_2000]
    mova    Tm0, [tab_Tm]
    mova    Tm1, [tab_Tm + 16]

    mov     r4d, %2
    cmp     r5m, byte 0
    je      .loopH
    sub     srcq, srcstrideq
    add     r4d, 3

.loopH:
    PROCESS_CHROMA_W%1   t0, t1, t2, t3
    add     srcq, srcstrideq
    add     dstq, dststrideq

    dec     r4d
    jnz     .loopH

    RET
%endmacro

FILTER_HORIZ_CHROMA_WxN 16, 4
FILTER_HORIZ_CHROMA_WxN 16, 8
FILTER_HORIZ_CHROMA_WxN 16, 12
FILTER_HORIZ_CHROMA_WxN 16, 16
FILTER_HORIZ_CHROMA_WxN 16, 32
FILTER_HORIZ_CHROMA_WxN 24, 32
FILTER_HORIZ_CHROMA_WxN 32,  8
FILTER_HORIZ_CHROMA_WxN 32, 16
FILTER_HORIZ_CHROMA_WxN 32, 24
FILTER_HORIZ_CHROMA_WxN 32, 32

FILTER_HORIZ_CHROMA_WxN 16, 24
FILTER_HORIZ_CHROMA_WxN 16, 64
FILTER_HORIZ_CHROMA_WxN 24, 64
FILTER_HORIZ_CHROMA_WxN 32, 48
FILTER_HORIZ_CHROMA_WxN 32, 64

FILTER_HORIZ_CHROMA_WxN 64, 64
FILTER_HORIZ_CHROMA_WxN 64, 32
FILTER_HORIZ_CHROMA_WxN 64, 48
FILTER_HORIZ_CHROMA_WxN 48, 64
FILTER_HORIZ_CHROMA_WxN 64, 16


;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_%1x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W16n 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_%1x%2, 4, 7, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m1, m0, [tab_Vm]
    pshufb     m0, [tab_Vm + 16]
    mov        r4d, %2/2

.loop:

    mov         r6d,       %1/16

.loopW:

    movu       m2, [r0]
    movu       m3, [r0 + r1]

    punpcklbw  m4, m2, m3
    punpckhbw  m2, m3

    pmaddubsw  m4, m1
    pmaddubsw  m2, m1

    lea        r5, [r0 + 2 * r1]
    movu       m5, [r5]
    movu       m7, [r5 + r1]

    punpcklbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m4, m6

    punpckhbw  m6, m5, m7
    pmaddubsw  m6, m0
    paddw      m2, m6

    mova       m6, [pw_2000]

    psubw      m4, m6
    psubw      m2, m6

    movu       [r2], m4
    movu       [r2 + 16], m2

    punpcklbw  m4, m3, m5
    punpckhbw  m3, m5

    pmaddubsw  m4, m1
    pmaddubsw  m3, m1

    movu       m5, [r5 + 2 * r1]

    punpcklbw  m2, m7, m5
    punpckhbw  m7, m5

    pmaddubsw  m2, m0
    pmaddubsw  m7, m0

    paddw      m4, m2
    paddw      m3, m7

    psubw      m4, m6
    psubw      m3, m6

    movu       [r2 + r3], m4
    movu       [r2 + r3 + 16], m3

    add         r0,        16
    add         r2,        32
    dec         r6d
    jnz         .loopW

    lea         r0,        [r0 + r1 * 2 - %1]
    lea         r2,        [r2 + r3 * 2 - %1 * 2]

    dec        r4d
    jnz        .loop
    RET
%endmacro

FILTER_V_PS_W16n 64, 64
FILTER_V_PS_W16n 64, 32
FILTER_V_PS_W16n 64, 48
FILTER_V_PS_W16n 48, 64
FILTER_V_PS_W16n 64, 16


;------------------------------------------------------------------------------------------------------------
;void interp_4tap_vert_ps_2x4(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;------------------------------------------------------------------------------------------------------------
INIT_XMM sse4
cglobal interp_4tap_vert_ps_2x4, 4, 6, 7

    mov         r4d, r4m
    sub         r0, r1
    add         r3d, r3d

%ifdef PIC
    lea         r5, [tab_ChromaCoeff]
    movd        m0, [r5 + r4 * 4]
%else
    movd        m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb      m0, [tab_Cm]

    lea         r5, [3 * r1]

    movd        m2, [r0]
    movd        m3, [r0 + r1]
    movd        m4, [r0 + 2 * r1]
    movd        m5, [r0 + r5]

    punpcklbw   m2, m3
    punpcklbw   m6, m4, m5
    punpcklbw   m2, m6

    pmaddubsw   m2, m0

    lea         r0, [r0 + 4 * r1]
    movd        m6, [r0]

    punpcklbw   m3, m4
    punpcklbw   m1, m5, m6
    punpcklbw   m3, m1

    pmaddubsw   m3, m0
    phaddw      m2, m3

    mova        m1, [pw_2000]

    psubw       m2, m1

    movd        [r2], m2
    pextrd      [r2 + r3], m2, 2

    movd        m2, [r0 + r1]

    punpcklbw   m4, m5
    punpcklbw   m3, m6, m2
    punpcklbw   m4, m3

    pmaddubsw   m4, m0

    movd        m3, [r0 + 2 * r1]

    punpcklbw   m5, m6
    punpcklbw   m2, m3
    punpcklbw   m5, m2

    pmaddubsw   m5, m0
    phaddw      m4, m5
    psubw       m4, m1

    lea         r2, [r2 + 2 * r3]
    movd        [r2], m4
    pextrd      [r2 + r3], m4, 2

    RET

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ps_2x8(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
%macro FILTER_V_PS_W2 2
INIT_XMM sse4
cglobal interp_4tap_vert_ps_2x%2, 4, 6, 8

    mov        r4d, r4m
    sub        r0, r1
    add        r3d, r3d

%ifdef PIC
    lea        r5, [tab_ChromaCoeff]
    movd       m0, [r5 + r4 * 4]
%else
    movd       m0, [tab_ChromaCoeff + r4 * 4]
%endif

    pshufb     m0, [tab_Cm]

    mova       m1, [pw_2000]
    lea        r5, [3 * r1]
    mov        r4d, %2/4
.loop:
    movd       m2, [r0]
    movd       m3, [r0 + r1]
    movd       m4, [r0 + 2 * r1]
    movd       m5, [r0 + r5]

    punpcklbw  m2, m3
    punpcklbw  m6, m4, m5
    punpcklbw  m2, m6

    pmaddubsw  m2, m0

    lea        r0, [r0 + 4 * r1]
    movd       m6, [r0]

    punpcklbw  m3, m4
    punpcklbw  m7, m5, m6
    punpcklbw  m3, m7

    pmaddubsw  m3, m0

    phaddw     m2, m3
    psubw      m2, m1


    movd       [r2], m2
    pshufd     m2, m2, 2
    movd       [r2 + r3], m2

    movd       m2, [r0 + r1]

    punpcklbw  m4, m5
    punpcklbw  m3, m6, m2
    punpcklbw  m4, m3

    pmaddubsw  m4, m0

    movd       m3, [r0 + 2 * r1]

    punpcklbw  m5, m6
    punpcklbw  m2, m3
    punpcklbw  m5, m2

    pmaddubsw  m5, m0

    phaddw     m4, m5

    psubw      m4, m1

    lea        r2, [r2 + 2 * r3]
    movd       [r2], m4
    pshufd     m4 , m4 ,2
    movd       [r2 + r3], m4

    lea        r2, [r2 + 2 * r3]

    dec        r4d
    jnz        .loop

RET
%endmacro

FILTER_V_PS_W2 2, 8

FILTER_V_PS_W2 2, 16

;-----------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ss_%1x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SS 2
INIT_XMM sse2
cglobal interp_4tap_vert_ss_%1x%2, 5, 7, 6 ,0-gprsize

    add       r1d, r1d
    add       r3d, r3d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_ChromaCoeffV + r4]
%endif

    mov       dword [rsp], %2/4

.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_CHROMA_SP_W4_4R

    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3

    movlps    [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movlps    [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [4 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 2 * 4

    dec       r4d
    jnz       .loopW

    lea       r0, [r0 + 4 * r1 - 2 * %1]
    lea       r2, [r2 + 4 * r3 - 2 * %1]

    dec       dword [rsp]
    jnz       .loopH

    RET
%endmacro

    FILTER_VER_CHROMA_SS 4, 4
    FILTER_VER_CHROMA_SS 4, 8
    FILTER_VER_CHROMA_SS 16, 16
    FILTER_VER_CHROMA_SS 16, 8
    FILTER_VER_CHROMA_SS 16, 12
    FILTER_VER_CHROMA_SS 12, 16
    FILTER_VER_CHROMA_SS 16, 4
    FILTER_VER_CHROMA_SS 4, 16
    FILTER_VER_CHROMA_SS 32, 32
    FILTER_VER_CHROMA_SS 32, 16
    FILTER_VER_CHROMA_SS 16, 32
    FILTER_VER_CHROMA_SS 32, 24
    FILTER_VER_CHROMA_SS 24, 32
    FILTER_VER_CHROMA_SS 32, 8

    FILTER_VER_CHROMA_SS 16, 24
    FILTER_VER_CHROMA_SS 12, 32
    FILTER_VER_CHROMA_SS 4, 32
    FILTER_VER_CHROMA_SS 32, 64
    FILTER_VER_CHROMA_SS 16, 64
    FILTER_VER_CHROMA_SS 32, 48
    FILTER_VER_CHROMA_SS 24, 64

    FILTER_VER_CHROMA_SS 64, 64
    FILTER_VER_CHROMA_SS 64, 32
    FILTER_VER_CHROMA_SS 64, 48
    FILTER_VER_CHROMA_SS 48, 64
    FILTER_VER_CHROMA_SS 64, 16

%macro FILTER_VER_CHROMA_S_AVX2_4x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x4, 4, 6, 7
    mov             r4d, r4m
    add             r1d, r1d
    shl             r4d, 6
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m6, [pd_526336]
%else
    add             r3d, r3d
%endif

    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m4, [r5 + 1 * mmsize]
    paddd           m2, m4

%ifidn %1,sp
    paddd           m0, m6
    paddd           m2, m6
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2
    vextracti128    xm2, m0, 1
    lea             r4, [r3 * 3]

%ifidn %1,sp
    packuswb        xm0, xm2
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 2
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r4], xm0, 3
%else
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r4], xm2
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_4x4 sp
FILTER_VER_CHROMA_S_AVX2_4x4 ss

%macro FILTER_VER_CHROMA_S_AVX2_4x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x8, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m5, m4, [r5 + 1 * mmsize]
    paddd           m2, m5
    pmaddwd         m4, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm1, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm6, [r0]
    punpcklwd       xm3, xm6
    vinserti128     m1, m1, xm3, 1                  ; m1 = [8 7 7 6]
    pmaddwd         m5, m1, [r5 + 1 * mmsize]
    paddd           m4, m5
    pmaddwd         m1, [r5]
    movq            xm3, [r0 + r1]
    punpcklwd       xm6, xm3
    movq            xm5, [r0 + 2 * r1]
    punpcklwd       xm3, xm5
    vinserti128     m6, m6, xm3, 1                  ; m6 = [A 9 9 8]
    pmaddwd         m6, [r5 + 1 * mmsize]
    paddd           m1, m6
    lea             r4, [r3 * 3]

%ifidn %1,sp
    paddd           m0, m7
    paddd           m2, m7
    paddd           m4, m7
    paddd           m1, m7
    psrad           m0, 12
    psrad           m2, 12
    psrad           m4, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m2, 6
    psrad           m4, 6
    psrad           m1, 6
%endif
    packssdw        m0, m2
    packssdw        m4, m1
%ifidn %1,sp
    packuswb        m0, m4
    vextracti128    xm2, m0, 1
    movd            [r2], xm0
    movd            [r2 + r3], xm2
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r4], xm2, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm2, 2
    pextrd          [r2 + r3 * 2], xm0, 3
    pextrd          [r2 + r4], xm2, 3
%else
    vextracti128    xm2, m0, 1
    vextracti128    xm1, m4, 1
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r4], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    movq            [r2 + r3], xm1
    movhps          [r2 + r3 * 2], xm4
    movhps          [r2 + r4], xm1
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_4x8 sp
FILTER_VER_CHROMA_S_AVX2_4x8 ss

%macro PROCESS_CHROMA_AVX2_W4_16R 1
    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m5, m4, [r5 + 1 * mmsize]
    paddd           m2, m5
    pmaddwd         m4, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm1, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm6, [r0]
    punpcklwd       xm3, xm6
    vinserti128     m1, m1, xm3, 1                  ; m1 = [8 7 7 6]
    pmaddwd         m5, m1, [r5 + 1 * mmsize]
    paddd           m4, m5
    pmaddwd         m1, [r5]
    movq            xm3, [r0 + r1]
    punpcklwd       xm6, xm3
    movq            xm5, [r0 + 2 * r1]
    punpcklwd       xm3, xm5
    vinserti128     m6, m6, xm3, 1                  ; m6 = [10 9 9 8]
    pmaddwd         m3, m6, [r5 + 1 * mmsize]
    paddd           m1, m3
    pmaddwd         m6, [r5]

%ifidn %1,sp
    paddd           m0, m7
    paddd           m2, m7
    paddd           m4, m7
    paddd           m1, m7
    psrad           m4, 12
    psrad           m1, 12
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
    psrad           m4, 6
    psrad           m1, 6
%endif
    packssdw        m0, m2
    packssdw        m4, m1
%ifidn %1,sp
    packuswb        m0, m4
    vextracti128    xm4, m0, 1
    movd            [r2], xm0
    movd            [r2 + r3], xm4
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r6], xm4, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm4, 2
    pextrd          [r2 + r3 * 2], xm0, 3
    pextrd          [r2 + r6], xm4, 3
%else
    vextracti128    xm2, m0, 1
    vextracti128    xm1, m4, 1
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    movq            [r2 + r3], xm1
    movhps          [r2 + r3 * 2], xm4
    movhps          [r2 + r6], xm1
%endif

    movq            xm2, [r0 + r4]
    punpcklwd       xm5, xm2
    lea             r0, [r0 + 4 * r1]
    movq            xm0, [r0]
    punpcklwd       xm2, xm0
    vinserti128     m5, m5, xm2, 1                  ; m5 = [12 11 11 10]
    pmaddwd         m2, m5, [r5 + 1 * mmsize]
    paddd           m6, m2
    pmaddwd         m5, [r5]
    movq            xm2, [r0 + r1]
    punpcklwd       xm0, xm2
    movq            xm3, [r0 + 2 * r1]
    punpcklwd       xm2, xm3
    vinserti128     m0, m0, xm2, 1                  ; m0 = [14 13 13 12]
    pmaddwd         m2, m0, [r5 + 1 * mmsize]
    paddd           m5, m2
    pmaddwd         m0, [r5]
    movq            xm4, [r0 + r4]
    punpcklwd       xm3, xm4
    lea             r0, [r0 + 4 * r1]
    movq            xm1, [r0]
    punpcklwd       xm4, xm1
    vinserti128     m3, m3, xm4, 1                  ; m3 = [16 15 15 14]
    pmaddwd         m4, m3, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m3, [r5]
    movq            xm4, [r0 + r1]
    punpcklwd       xm1, xm4
    movq            xm2, [r0 + 2 * r1]
    punpcklwd       xm4, xm2
    vinserti128     m1, m1, xm4, 1                  ; m1 = [18 17 17 16]
    pmaddwd         m1, [r5 + 1 * mmsize]
    paddd           m3, m1

%ifidn %1,sp
    paddd           m6, m7
    paddd           m5, m7
    paddd           m0, m7
    paddd           m3, m7
    psrad           m6, 12
    psrad           m5, 12
    psrad           m0, 12
    psrad           m3, 12
%else
    psrad           m6, 6
    psrad           m5, 6
    psrad           m0, 6
    psrad           m3, 6
%endif
    packssdw        m6, m5
    packssdw        m0, m3
    lea             r2, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m6, m0
    vextracti128    xm0, m6, 1
    movd            [r2], xm6
    movd            [r2 + r3], xm0
    pextrd          [r2 + r3 * 2], xm6, 1
    pextrd          [r2 + r6], xm0, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm6, 2
    pextrd          [r2 + r3], xm0, 2
    pextrd          [r2 + r3 * 2], xm6, 3
    pextrd          [r2 + r6], xm0, 3
%else
    vextracti128    xm5, m6, 1
    vextracti128    xm3, m0, 1
    movq            [r2], xm6
    movq            [r2 + r3], xm5
    movhps          [r2 + r3 * 2], xm6
    movhps          [r2 + r6], xm5
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm0
    movq            [r2 + r3], xm3
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm3
%endif
%endmacro

%macro FILTER_VER_CHROMA_S_AVX2_4x16 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x16, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    PROCESS_CHROMA_AVX2_W4_16R %1
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_4x16 sp
FILTER_VER_CHROMA_S_AVX2_4x16 ss

%macro FILTER_VER_CHROMA_S_AVX2_4x2 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_4x2, 4, 6, 6
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m5, [pd_526336]
%else
    add             r3d, r3d
%endif
    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    movq            xm4, [r0 + 4 * r1]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
%ifidn %1,sp
    paddd           m0, m5
    psrad           m0, 12
%else
    psrad           m0, 6
%endif
    vextracti128    xm1, m0, 1
    packssdw        xm0, xm1
%ifidn %1,sp
    packuswb        xm0, xm0
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
%else
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_4x2 sp
FILTER_VER_CHROMA_S_AVX2_4x2 ss

%macro FILTER_VER_CHROMA_S_AVX2_2x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_2x4, 4, 6, 6
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m5, [pd_526336]
%else
    add             r3d, r3d
%endif
    movd            xm0, [r0]
    movd            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movd            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    punpcklqdq      xm0, xm1                        ; m0 = [2 1 1 0]
    movd            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movd            xm4, [r0]
    punpcklwd       xm3, xm4
    punpcklqdq      xm2, xm3                        ; m2 = [4 3 3 2]
    vinserti128     m0, m0, xm2, 1                  ; m0 = [4 3 3 2 2 1 1 0]
    movd            xm1, [r0 + r1]
    punpcklwd       xm4, xm1
    movd            xm3, [r0 + r1 * 2]
    punpcklwd       xm1, xm3
    punpcklqdq      xm4, xm1                        ; m4 = [6 5 5 4]
    vinserti128     m2, m2, xm4, 1                  ; m2 = [6 5 5 4 4 3 3 2]
    pmaddwd         m0, [r5]
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
%ifidn %1,sp
    paddd           m0, m5
    psrad           m0, 12
%else
    psrad           m0, 6
%endif
    vextracti128    xm1, m0, 1
    packssdw        xm0, xm1
    lea             r4, [r3 * 3]
%ifidn %1,sp
    packuswb        xm0, xm0
    pextrw          [r2], xm0, 0
    pextrw          [r2 + r3], xm0, 1
    pextrw          [r2 + 2 * r3], xm0, 2
    pextrw          [r2 + r4], xm0, 3
%else
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    pextrd          [r2 + 2 * r3], xm0, 2
    pextrd          [r2 + r4], xm0, 3
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_2x4 sp
FILTER_VER_CHROMA_S_AVX2_2x4 ss

%macro FILTER_VER_CHROMA_S_AVX2_8x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x8, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    pmaddwd         m3, [r5]
    paddd           m1, m5
%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m7
    paddd           m3, m7
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3

    movu            xm1, [r0 + r4]                  ; m1 = row 7
    punpckhwd       xm3, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm3, 1
    pmaddwd         m3, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m3

    lea             r4, [r3 * 3]
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
%endif
    lea             r2, [r2 + r3 * 4]
    lea             r0, [r0 + r1 * 4]
    movu            xm0, [r0]                       ; m0 = row 8
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    pmaddwd         m1, [r5]
    paddd           m5, m2
%ifidn %1,sp
    paddd           m4, m7
    paddd           m5, m7
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5

    movu            xm2, [r0 + r1]                  ; m2 = row 9
    punpckhwd       xm5, xm0, xm2
    punpcklwd       xm0, xm2
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m0, [r5 + 1 * mmsize]
    paddd           m6, m0
    movu            xm5, [r0 + r1 * 2]              ; m5 = row 10
    punpckhwd       xm0, xm2, xm5
    punpcklwd       xm2, xm5
    vinserti128     m2, m2, xm0, 1
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m1, m2

%ifidn %1,sp
    paddd           m6, m7
    paddd           m1, m7
    psrad           m6, 12
    psrad           m1, 12
%else
    psrad           m6, 6
    psrad           m1, 6
%endif
    packssdw        m6, m1
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm6
    movhps          [r2 + r4], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm5, m4, 1
    vextracti128    xm1, m6, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm5
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r4], xm1
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_8x8 sp
FILTER_VER_CHROMA_S_AVX2_8x8 ss

%macro PROCESS_CHROMA_S_AVX2_W8_16R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m3, m7
    pmaddwd         m5, [r5]
%ifidn %1,sp
    paddd           m0, m9
    paddd           m1, m9
    paddd           m2, m9
    paddd           m3, m9
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif

    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    paddd           m4, m8
    pmaddwd         m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhwd       xm0, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm0, 1
    pmaddwd         m0, m7, [r5 + 1 * mmsize]
    paddd           m5, m0
    pmaddwd         m7, [r5]
    movu            xm0, [r7 + r1]                  ; m0 = row 9
    punpckhwd       xm1, xm8, xm0
    punpcklwd       xm8, xm0
    vinserti128     m8, m8, xm1, 1
    pmaddwd         m1, m8, [r5 + 1 * mmsize]
    paddd           m6, m1
    pmaddwd         m8, [r5]
    movu            xm1, [r7 + r1 * 2]              ; m1 = row 10
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m2, m0, [r5 + 1 * mmsize]
    paddd           m7, m2
    pmaddwd         m0, [r5]
%ifidn %1,sp
    paddd           m4, m9
    paddd           m5, m9
    psrad           m4, 12
    psrad           m5, 12
    paddd           m6, m9
    paddd           m7, m9
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m4, 6
    psrad           m5, 6
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m4, m5
    packssdw        m6, m7
    lea             r8, [r2 + r3 * 4]
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm5, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm5
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%endif

    movu            xm2, [r7 + r4]                  ; m2 = row 11
    punpckhwd       xm4, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm4, 1
    pmaddwd         m4, m1, [r5 + 1 * mmsize]
    paddd           m8, m4
    pmaddwd         m1, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 12
    punpckhwd       xm5, xm2, xm4
    punpcklwd       xm2, xm4
    vinserti128     m2, m2, xm5, 1
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    paddd           m0, m5
    pmaddwd         m2, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 13
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m1, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 14
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m2, m7
    pmaddwd         m5, [r5]
%ifidn %1,sp
    paddd           m8, m9
    paddd           m0, m9
    paddd           m1, m9
    paddd           m2, m9
    psrad           m8, 12
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
%else
    psrad           m8, 6
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
%endif
    packssdw        m8, m0
    packssdw        m1, m2
    lea             r8, [r8 + r3 * 4]
%ifidn %1,sp
    packuswb        m8, m1
    vpermd          m8, m3, m8
    vextracti128    xm1, m8, 1
    movq            [r8], xm8
    movhps          [r8 + r3], xm8
    movq            [r8 + r3 * 2], xm1
    movhps          [r8 + r6], xm1
%else
    vpermq          m8, m8, 11011000b
    vpermq          m1, m1, 11011000b
    vextracti128    xm0, m8, 1
    vextracti128    xm2, m1, 1
    movu            [r8], xm8
    movu            [r8 + r3], xm0
    movu            [r8 + r3 * 2], xm1
    movu            [r8 + r6], xm2
%endif
    lea             r8, [r8 + r3 * 4]

    movu            xm7, [r7 + r4]                  ; m7 = row 15
    punpckhwd       xm2, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddwd         m2, m6, [r5 + 1 * mmsize]
    paddd           m4, m2
    pmaddwd         m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm2, [r7]                       ; m2 = row 16
    punpckhwd       xm1, xm7, xm2
    punpcklwd       xm7, xm2
    vinserti128     m7, m7, xm1, 1
    pmaddwd         m1, m7, [r5 + 1 * mmsize]
    paddd           m5, m1
    pmaddwd         m7, [r5]
    movu            xm1, [r7 + r1]                  ; m1 = row 17
    punpckhwd       xm0, xm2, xm1
    punpcklwd       xm2, xm1
    vinserti128     m2, m2, xm0, 1
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m6, m2
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 18
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m1, [r5 + 1 * mmsize]
    paddd           m7, m1

%ifidn %1,sp
    paddd           m4, m9
    paddd           m5, m9
    paddd           m6, m9
    paddd           m7, m9
    psrad           m4, 12
    psrad           m5, 12
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m4, 6
    psrad           m5, 6
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m4, m5
    packssdw        m6, m7
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm5, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm5
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%endif
%endmacro

%macro FILTER_VER_CHROMA_S_AVX2_Nx16 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_%2x16, 4, 10, 10
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m9, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, %2 / 8
.loopW:
    PROCESS_CHROMA_S_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_Nx16 sp, 16
FILTER_VER_CHROMA_S_AVX2_Nx16 sp, 32
FILTER_VER_CHROMA_S_AVX2_Nx16 ss, 16
FILTER_VER_CHROMA_S_AVX2_Nx16 ss, 32

%macro FILTER_VER_CHROMA_S_AVX2_NxN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%3_%1x%2, 4, 11, 10
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif
    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %3,sp
    mova            m9, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, %2 / 16
.loopH:
    mov             r10d, %1 / 8
.loopW:
    PROCESS_CHROMA_S_AVX2_W8_16R %3
%ifidn %3,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r10d
    jnz             .loopW
    lea             r0, [r7 - 2 * %1 + 16]
%ifidn %3,sp
    lea             r2, [r8 + r3 * 4 - %1 + 8]
%else
    lea             r2, [r8 + r3 * 4 - 2 * %1 + 16]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_NxN 16, 32, sp
FILTER_VER_CHROMA_S_AVX2_NxN 24, 32, sp
FILTER_VER_CHROMA_S_AVX2_NxN 32, 32, sp
FILTER_VER_CHROMA_S_AVX2_NxN 16, 32, ss
FILTER_VER_CHROMA_S_AVX2_NxN 24, 32, ss
FILTER_VER_CHROMA_S_AVX2_NxN 32, 32, ss

%macro PROCESS_CHROMA_S_AVX2_W8_4R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m4, [r5 + 1 * mmsize]
    paddd           m2, m4
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm4, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm4, 1
    pmaddwd         m5, [r5 + 1 * mmsize]
    paddd           m3, m5
%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    paddd           m2, m7
    paddd           m3, m7
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
%endif
%endmacro

%macro FILTER_VER_CHROMA_S_AVX2_8x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x4, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    PROCESS_CHROMA_S_AVX2_W8_4R %1
    lea             r4, [r3 * 3]
%ifidn %1,sp
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
%else
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_8x4 sp
FILTER_VER_CHROMA_S_AVX2_8x4 ss

%macro FILTER_VER_CHROMA_S_AVX2_12x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_12x16, 4, 9, 10
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m9, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    PROCESS_CHROMA_S_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    mova            m7, m9
    PROCESS_CHROMA_AVX2_W4_16R %1
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_12x16 sp
FILTER_VER_CHROMA_S_AVX2_12x16 ss

%macro FILTER_VER_CHROMA_S_AVX2_16x12 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_16x12, 4, 9, 9
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m8, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
%rep 2
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
%ifidn %1,sp
    paddd           m0, m8
    paddd           m1, m8
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m8
    paddd           m3, m8
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif
    lea             r8, [r2 + r3 * 4]

    movu            xm1, [r7 + r4]                  ; m1 = row 7
    punpckhwd       xm0, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm0, 1
    pmaddwd         m0, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m0
    lea             r7, [r7 + r1 * 4]
    movu            xm0, [r7]                       ; m0 = row 8
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    pmaddwd         m1, [r5]
    paddd           m5, m2
%ifidn %1,sp
    paddd           m4, m8
    paddd           m5, m8
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5

    movu            xm2, [r7 + r1]                  ; m2 = row 9
    punpckhwd       xm5, xm0, xm2
    punpcklwd       xm0, xm2
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m5, m0, [r5 + 1 * mmsize]
    paddd           m6, m5
    pmaddwd         m0, [r5]
    movu            xm5, [r7 + r1 * 2]              ; m5 = row 10
    punpckhwd       xm7, xm2, xm5
    punpcklwd       xm2, xm5
    vinserti128     m2, m2, xm7, 1
    pmaddwd         m7, m2, [r5 + 1 * mmsize]
    paddd           m1, m7
    pmaddwd         m2, [r5]

%ifidn %1,sp
    paddd           m6, m8
    paddd           m1, m8
    psrad           m6, 12
    psrad           m1, 12
%else
    psrad           m6, 6
    psrad           m1, 6
%endif
    packssdw        m6, m1
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m4, 1
    vextracti128    xm1, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm7
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm1
%endif
    lea             r8, [r8 + r3 * 4]

    movu            xm7, [r7 + r4]                  ; m7 = row 11
    punpckhwd       xm1, xm5, xm7
    punpcklwd       xm5, xm7
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    paddd           m0, m1
    pmaddwd         m5, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm1, [r7]                       ; m1 = row 12
    punpckhwd       xm4, xm7, xm1
    punpcklwd       xm7, xm1
    vinserti128     m7, m7, xm4, 1
    pmaddwd         m4, m7, [r5 + 1 * mmsize]
    paddd           m2, m4
    pmaddwd         m7, [r5]
%ifidn %1,sp
    paddd           m0, m8
    paddd           m2, m8
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2

    movu            xm4, [r7 + r1]                  ; m4 = row 13
    punpckhwd       xm2, xm1, xm4
    punpcklwd       xm1, xm4
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m1, [r5 + 1 * mmsize]
    paddd           m5, m1
    movu            xm2, [r7 + r1 * 2]              ; m2 = row 14
    punpckhwd       xm6, xm4, xm2
    punpcklwd       xm4, xm2
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m4, [r5 + 1 * mmsize]
    paddd           m7, m4
%ifidn %1,sp
    paddd           m5, m8
    paddd           m7, m8
    psrad           m5, 12
    psrad           m7, 12
%else
    psrad           m5, 6
    psrad           m7, 6
%endif
    packssdw        m5, m7
%ifidn %1,sp
    packuswb        m0, m5
    vpermd          m0, m3, m0
    vextracti128    xm5, m0, 1
    movq            [r8], xm0
    movhps          [r8 + r3], xm0
    movq            [r8 + r3 * 2], xm5
    movhps          [r8 + r6], xm5
    add             r2, 8
%else
    vpermq          m0, m0, 11011000b
    vpermq          m5, m5, 11011000b
    vextracti128    xm7, m0, 1
    vextracti128    xm6, m5, 1
    movu            [r8], xm0
    movu            [r8 + r3], xm7
    movu            [r8 + r3 * 2], xm5
    movu            [r8 + r6], xm6
    add             r2, 16
%endif
    add             r0, 16
%endrep
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_16x12 sp
FILTER_VER_CHROMA_S_AVX2_16x12 ss

%macro FILTER_VER_CHROMA_S_AVX2_16x4 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_16x4, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
%rep 2
    PROCESS_CHROMA_S_AVX2_W8_4R %1
    lea             r6, [r3 * 3]
%ifidn %1,sp
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
    add             r2, 8
%else
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    add             r2, 16
%endif
    lea             r6, [4 * r1 - 16]
    sub             r0, r6
%endrep
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_16x4 sp
FILTER_VER_CHROMA_S_AVX2_16x4 ss

%macro PROCESS_CHROMA_S_AVX2_W8_8R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m7
    paddd           m3, m7
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif
    lea             r8, [r2 + r3 * 4]

    movu            xm1, [r7 + r4]                  ; m1 = row 7
    punpckhwd       xm0, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm0, 1
    pmaddwd         m0, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m0
    lea             r7, [r7 + r1 * 4]
    movu            xm0, [r7]                       ; m0 = row 8
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    pmaddwd         m1, [r5]
    paddd           m5, m2
%ifidn %1,sp
    paddd           m4, m7
    paddd           m5, m7
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5

    movu            xm2, [r7 + r1]                  ; m2 = row 9
    punpckhwd       xm5, xm0, xm2
    punpcklwd       xm0, xm2
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m0, [r5 + 1 * mmsize]
    paddd           m6, m0
    movu            xm5, [r7 + r1 * 2]              ; m5 = row 10
    punpckhwd       xm0, xm2, xm5
    punpcklwd       xm2, xm5
    vinserti128     m2, m2, xm0, 1
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m1, m2

%ifidn %1,sp
    paddd           m6, m7
    paddd           m1, m7
    psrad           m6, 12
    psrad           m1, 12
%else
    psrad           m6, 6
    psrad           m1, 6
%endif
    packssdw        m6, m1
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m4, 1
    vextracti128    xm1, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm7
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm1
%endif
%endmacro

%macro FILTER_VER_CHROMA_S_AVX2_Nx8 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_%2x8, 4, 9, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
%rep %2 / 8
    PROCESS_CHROMA_S_AVX2_W8_8R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
%endrep
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_Nx8 sp, 32
FILTER_VER_CHROMA_S_AVX2_Nx8 sp, 16
FILTER_VER_CHROMA_S_AVX2_Nx8 ss, 32
FILTER_VER_CHROMA_S_AVX2_Nx8 ss, 16

%macro FILTER_VER_CHROMA_S_AVX2_8x2 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x2, 4, 6, 6
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m5, [pd_526336]
%else
    add             r3d, r3d
%endif

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
    movu            xm4, [r0 + r1 * 4]              ; m4 = row 4
    punpckhwd       xm2, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm2, 1
    pmaddwd         m3, [r5 + 1 * mmsize]
    paddd           m1, m3
%ifidn %1,sp
    paddd           m0, m5
    paddd           m1, m5
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1
%ifidn %1,sp
    vextracti128    xm1, m0, 1
    packuswb        xm0, xm1
    pshufd          xm0, xm0, 11011000b
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
%else
    vpermq          m0, m0, 11011000b
    vextracti128    xm1, m0, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_8x2 sp
FILTER_VER_CHROMA_S_AVX2_8x2 ss

%macro FILTER_VER_CHROMA_S_AVX2_8x6 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_8x6, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    pmaddwd         m3, [r5]
    paddd           m1, m5
%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m7
    paddd           m3, m7
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3

    movu            xm1, [r0 + r4]                  ; m1 = row 7
    punpckhwd       xm3, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm3, 1
    pmaddwd         m6, [r5 + 1 * mmsize]
    paddd           m4, m6
    movu            xm6, [r0 + r1 * 4]              ; m6 = row 8
    punpckhwd       xm3, xm1, xm6
    punpcklwd       xm1, xm6
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5 + 1 * mmsize]
    paddd           m5, m1
%ifidn %1,sp
    paddd           m4, m7
    paddd           m5, m7
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5
    lea             r4, [r3 * 3]
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    vextracti128    xm5, m4, 1
    packuswb        xm4, xm5
    pshufd          xm4, xm4, 11011000b
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vpermq          m4, m4, 11011000b
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    vextracti128    xm5, m4, 1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm4
    movu            [r2 + r3], xm5
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_8x6 sp
FILTER_VER_CHROMA_S_AVX2_8x6 ss

%macro FILTER_VER_CHROMA_S_AVX2_8xN 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_8x%2, 4, 7, 9
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m8, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
%rep %2 / 16
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
%ifidn %1,sp
    paddd           m0, m8
    paddd           m1, m8
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m8
    paddd           m3, m8
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m3, [interp8_hps_shuf]
    vpermd          m0, m3, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    movu            [r2], xm0
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2 + r3], xm0
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm1, [r0 + r4]                  ; m1 = row 7
    punpckhwd       xm0, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm0, 1
    pmaddwd         m0, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m0
    lea             r0, [r0 + r1 * 4]
    movu            xm0, [r0]                       ; m0 = row 8
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    pmaddwd         m1, [r5]
    paddd           m5, m2
%ifidn %1,sp
    paddd           m4, m8
    paddd           m5, m8
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5

    movu            xm2, [r0 + r1]                  ; m2 = row 9
    punpckhwd       xm5, xm0, xm2
    punpcklwd       xm0, xm2
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m5, m0, [r5 + 1 * mmsize]
    paddd           m6, m5
    pmaddwd         m0, [r5]
    movu            xm5, [r0 + r1 * 2]              ; m5 = row 10
    punpckhwd       xm7, xm2, xm5
    punpcklwd       xm2, xm5
    vinserti128     m2, m2, xm7, 1
    pmaddwd         m7, m2, [r5 + 1 * mmsize]
    paddd           m1, m7
    pmaddwd         m2, [r5]

%ifidn %1,sp
    paddd           m6, m8
    paddd           m1, m8
    psrad           m6, 12
    psrad           m1, 12
%else
    psrad           m6, 6
    psrad           m1, 6
%endif
    packssdw        m6, m1
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m3, m4
    vextracti128    xm6, m4, 1
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm6
    movhps          [r2 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm7, m4, 1
    vextracti128    xm1, m6, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm7
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm1
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm7, [r0 + r4]                  ; m7 = row 11
    punpckhwd       xm1, xm5, xm7
    punpcklwd       xm5, xm7
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    paddd           m0, m1
    pmaddwd         m5, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm1, [r0]                       ; m1 = row 12
    punpckhwd       xm4, xm7, xm1
    punpcklwd       xm7, xm1
    vinserti128     m7, m7, xm4, 1
    pmaddwd         m4, m7, [r5 + 1 * mmsize]
    paddd           m2, m4
    pmaddwd         m7, [r5]
%ifidn %1,sp
    paddd           m0, m8
    paddd           m2, m8
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2

    movu            xm4, [r0 + r1]                  ; m4 = row 13
    punpckhwd       xm2, xm1, xm4
    punpcklwd       xm1, xm4
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    paddd           m5, m2
    pmaddwd         m1, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 14
    punpckhwd       xm6, xm4, xm2
    punpcklwd       xm4, xm2
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m7, m6
    pmaddwd         m4, [r5]
%ifidn %1,sp
    paddd           m5, m8
    paddd           m7, m8
    psrad           m5, 12
    psrad           m7, 12
%else
    psrad           m5, 6
    psrad           m7, 6
%endif
    packssdw        m5, m7
%ifidn %1,sp
    packuswb        m0, m5
    vpermd          m0, m3, m0
    vextracti128    xm5, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm5
    movhps          [r2 + r6], xm5
%else
    vpermq          m0, m0, 11011000b
    vpermq          m5, m5, 11011000b
    vextracti128    xm7, m0, 1
    vextracti128    xm6, m5, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm7
    movu            [r2 + r3 * 2], xm5
    movu            [r2 + r6], xm6
%endif
    lea             r2, [r2 + r3 * 4]

    movu            xm6, [r0 + r4]                  ; m6 = row 15
    punpckhwd       xm5, xm2, xm6
    punpcklwd       xm2, xm6
    vinserti128     m2, m2, xm5, 1
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm0, [r0]                       ; m0 = row 16
    punpckhwd       xm5, xm6, xm0
    punpcklwd       xm6, xm0
    vinserti128     m6, m6, xm5, 1
    pmaddwd         m5, m6, [r5 + 1 * mmsize]
    paddd           m4, m5
    pmaddwd         m6, [r5]
%ifidn %1,sp
    paddd           m1, m8
    paddd           m4, m8
    psrad           m1, 12
    psrad           m4, 12
%else
    psrad           m1, 6
    psrad           m4, 6
%endif
    packssdw        m1, m4

    movu            xm5, [r0 + r1]                  ; m5 = row 17
    punpckhwd       xm4, xm0, xm5
    punpcklwd       xm0, xm5
    vinserti128     m0, m0, xm4, 1
    pmaddwd         m0, [r5 + 1 * mmsize]
    paddd           m2, m0
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhwd       xm0, xm5, xm4
    punpcklwd       xm5, xm4
    vinserti128     m5, m5, xm0, 1
    pmaddwd         m5, [r5 + 1 * mmsize]
    paddd           m6, m5
%ifidn %1,sp
    paddd           m2, m8
    paddd           m6, m8
    psrad           m2, 12
    psrad           m6, 12
%else
    psrad           m2, 6
    psrad           m6, 6
%endif
    packssdw        m2, m6
%ifidn %1,sp
    packuswb        m1, m2
    vpermd          m1, m3, m1
    vextracti128    xm2, m1, 1
    movq            [r2], xm1
    movhps          [r2 + r3], xm1
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m1, m1, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm6, m1, 1
    vextracti128    xm4, m2, 1
    movu            [r2], xm1
    movu            [r2 + r3], xm6
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm4
%endif
    lea             r2, [r2 + r3 * 4]
%endrep
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_8xN sp, 16
FILTER_VER_CHROMA_S_AVX2_8xN sp, 32
FILTER_VER_CHROMA_S_AVX2_8xN ss, 16
FILTER_VER_CHROMA_S_AVX2_8xN ss, 32

%macro FILTER_VER_CHROMA_S_AVX2_32x24 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_4tap_vert_%1_32x24, 4, 10, 10
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m9, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, 4
.loopW:
    PROCESS_CHROMA_S_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
%ifidn %1,sp
    lea             r2, [r8 + r3 * 4 - 24]
%else
    lea             r2, [r8 + r3 * 4 - 48]
%endif
    lea             r0, [r7 - 48]
    mova            m7, m9
    mov             r9d, 4
.loop:
    PROCESS_CHROMA_S_AVX2_W8_8R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loop
    RET
%endif
%endmacro

FILTER_VER_CHROMA_S_AVX2_32x24 sp
FILTER_VER_CHROMA_S_AVX2_32x24 ss

%macro FILTER_VER_CHROMA_S_AVX2_2x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_2x8, 4, 6, 7
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d
    sub             r0, r1

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
%ifidn %1,sp
    mova            m6, [pd_526336]
%else
    add             r3d, r3d
%endif
    movd            xm0, [r0]
    movd            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movd            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    punpcklqdq      xm0, xm1                        ; m0 = [2 1 1 0]
    movd            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movd            xm4, [r0]
    punpcklwd       xm3, xm4
    punpcklqdq      xm2, xm3                        ; m2 = [4 3 3 2]
    vinserti128     m0, m0, xm2, 1                  ; m0 = [4 3 3 2 2 1 1 0]
    movd            xm1, [r0 + r1]
    punpcklwd       xm4, xm1
    movd            xm3, [r0 + r1 * 2]
    punpcklwd       xm1, xm3
    punpcklqdq      xm4, xm1                        ; m4 = [6 5 5 4]
    vinserti128     m2, m2, xm4, 1                  ; m2 = [6 5 5 4 4 3 3 2]
    pmaddwd         m0, [r5]
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
    movd            xm1, [r0 + r4]
    punpcklwd       xm3, xm1
    lea             r0, [r0 + 4 * r1]
    movd            xm2, [r0]
    punpcklwd       xm1, xm2
    punpcklqdq      xm3, xm1                        ; m3 = [8 7 7 6]
    vinserti128     m4, m4, xm3, 1                  ; m4 = [8 7 7 6 6 5 5 4]
    movd            xm1, [r0 + r1]
    punpcklwd       xm2, xm1
    movd            xm5, [r0 + r1 * 2]
    punpcklwd       xm1, xm5
    punpcklqdq      xm2, xm1                        ; m2 = [10 9 9 8]
    vinserti128     m3, m3, xm2, 1                  ; m3 = [10 9 9 8 8 7 7 6]
    pmaddwd         m4, [r5]
    pmaddwd         m3, [r5 + 1 * mmsize]
    paddd           m4, m3
%ifidn %1,sp
    paddd           m0, m6
    paddd           m4, m6
    psrad           m0, 12
    psrad           m4, 12
%else
    psrad           m0, 6
    psrad           m4, 6
%endif
    packssdw        m0, m4
    vextracti128    xm4, m0, 1
    lea             r4, [r3 * 3]
%ifidn %1,sp
    packuswb        xm0, xm4
    pextrw          [r2], xm0, 0
    pextrw          [r2 + r3], xm0, 1
    pextrw          [r2 + 2 * r3], xm0, 4
    pextrw          [r2 + r4], xm0, 5
    lea             r2, [r2 + r3 * 4]
    pextrw          [r2], xm0, 2
    pextrw          [r2 + r3], xm0, 3
    pextrw          [r2 + 2 * r3], xm0, 6
    pextrw          [r2 + r4], xm0, 7
%else
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 1
    movd            [r2 + 2 * r3], xm4
    pextrd          [r2 + r4], xm4, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm0, 3
    pextrd          [r2 + 2 * r3], xm4, 2
    pextrd          [r2 + r4], xm4, 3
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_2x8 sp
FILTER_VER_CHROMA_S_AVX2_2x8 ss

%macro FILTER_VER_CHROMA_S_AVX2_6x8 1
INIT_YMM avx2
cglobal interp_4tap_vert_%1_6x8, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 6
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_ChromaCoeffV]
    add             r5, r4
%else
    lea             r5, [pw_ChromaCoeffV + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r1
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    pmaddwd         m3, [r5]
    paddd           m1, m5
%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m0, m1

    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm1, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm1, 1
    pmaddwd         m1, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m1
%ifidn %1,sp
    paddd           m2, m7
    paddd           m3, m7
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m2, m3

    movu            xm1, [r0 + r4]                  ; m1 = row 7
    punpckhwd       xm3, xm6, xm1
    punpcklwd       xm6, xm1
    vinserti128     m6, m6, xm3, 1
    pmaddwd         m3, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m3

    lea             r4, [r3 * 3]
%ifidn %1,sp
    packuswb        m0, m2
    vextracti128    xm2, m0, 1
    movd            [r2], xm0
    pextrw          [r2 + 4], xm2, 0
    pextrd          [r2 + r3], xm0, 1
    pextrw          [r2 + r3 + 4], xm2, 2
    pextrd          [r2 + r3 * 2], xm0, 2
    pextrw          [r2 + r3 * 2 + 4], xm2, 4
    pextrd          [r2 + r4], xm0, 3
    pextrw          [r2 + r4 + 4], xm2, 6
%else
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
    vextracti128    xm0, m0, 1
    vextracti128    xm3, m2, 1
    movd            [r2 + 8], xm0
    pextrd          [r2 + r3 + 8], xm0, 2
    movd            [r2 + r3 * 2 + 8], xm3
    pextrd          [r2 + r4 + 8], xm3, 2
%endif
    lea             r2, [r2 + r3 * 4]
    lea             r0, [r0 + r1 * 4]
    movu            xm0, [r0]                       ; m0 = row 8
    punpckhwd       xm2, xm1, xm0
    punpcklwd       xm1, xm0
    vinserti128     m1, m1, xm2, 1
    pmaddwd         m2, m1, [r5 + 1 * mmsize]
    pmaddwd         m1, [r5]
    paddd           m5, m2
%ifidn %1,sp
    paddd           m4, m7
    paddd           m5, m7
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m4, m5

    movu            xm2, [r0 + r1]                  ; m2 = row 9
    punpckhwd       xm5, xm0, xm2
    punpcklwd       xm0, xm2
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m0, [r5 + 1 * mmsize]
    paddd           m6, m0
    movu            xm5, [r0 + r1 * 2]              ; m5 = row 10
    punpckhwd       xm0, xm2, xm5
    punpcklwd       xm2, xm5
    vinserti128     m2, m2, xm0, 1
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m1, m2

%ifidn %1,sp
    paddd           m6, m7
    paddd           m1, m7
    psrad           m6, 12
    psrad           m1, 12
%else
    psrad           m6, 6
    psrad           m1, 6
%endif
    packssdw        m6, m1
%ifidn %1,sp
    packuswb        m4, m6
    vextracti128    xm6, m4, 1
    movd            [r2], xm4
    pextrw          [r2 + 4], xm6, 0
    pextrd          [r2 + r3], xm4, 1
    pextrw          [r2 + r3 + 4], xm6, 2
    pextrd          [r2 + r3 * 2], xm4, 2
    pextrw          [r2 + r3 * 2 + 4], xm6, 4
    pextrd          [r2 + r4], xm4, 3
    pextrw          [r2 + r4 + 4], xm6, 6
%else
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm6
    movhps          [r2 + r4], xm6
    vextracti128    xm5, m4, 1
    vextracti128    xm1, m6, 1
    movd            [r2 + 8], xm5
    pextrd          [r2 + r3 + 8], xm5, 2
    movd            [r2 + r3 * 2 + 8], xm1
    pextrd          [r2 + r4 + 8], xm1, 2
%endif
    RET
%endmacro

FILTER_VER_CHROMA_S_AVX2_6x8 sp
FILTER_VER_CHROMA_S_AVX2_6x8 ss

;---------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_ss_%1x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SS_W2_4R 2
INIT_XMM sse4
cglobal interp_4tap_vert_ss_%1x%2, 5, 6, 5

    add       r1d, r1d
    add       r3d, r3d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r5, [r5 + r4]
%else
    lea       r5, [tab_ChromaCoeffV + r4]
%endif

    mov       r4d, (%2/4)

.loopH:
    PROCESS_CHROMA_SP_W2_4R r5

    psrad     m0, 6
    psrad     m2, 6

    packssdw  m0, m2

    movd      [r2], m0
    pextrd    [r2 + r3], m0, 1
    lea       r2, [r2 + 2 * r3]
    pextrd    [r2], m0, 2
    pextrd    [r2 + r3], m0, 3

    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SS_W2_4R 2, 4
FILTER_VER_CHROMA_SS_W2_4R 2, 8

FILTER_VER_CHROMA_SS_W2_4R 2, 16

;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ss_4x2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
INIT_XMM sse2
cglobal interp_4tap_vert_ss_4x2, 5, 6, 4

    add        r1d, r1d
    add        r3d, r3d
    sub        r0, r1
    shl        r4d, 5

%ifdef PIC
    lea        r5, [tab_ChromaCoeffV]
    lea        r5, [r5 + r4]
%else
    lea        r5, [tab_ChromaCoeffV + r4]
%endif

    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r5 + 0 *16]                ;m0=[0+1]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m2, [r0]
    punpcklwd  m1, m2                          ;m1=[1 2]
    pmaddwd    m1, [r5 + 0 *16]                ;m1=[1+2]  Row2

    movq       m3, [r0 + r1]
    punpcklwd  m2, m3                          ;m4=[2 3]
    pmaddwd    m2, [r5 + 1 * 16]
    paddd      m0, m2                          ;m0=[0+1+2+3]  Row1 done
    psrad      m0, 6

    movq       m2, [r0 + 2 * r1]
    punpcklwd  m3, m2                          ;m5=[3 4]
    pmaddwd    m3, [r5 + 1 * 16]
    paddd      m1, m3                          ;m1=[1+2+3+4]  Row2 done
    psrad      m1, 6

    packssdw   m0, m1

    movlps     [r2], m0
    movhps     [r2 + r3], m0

    RET

;-------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_ss_6x8(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SS_W6_H4 2
INIT_XMM sse4
cglobal interp_4tap_vert_ss_6x%2, 5, 7, 6

    add       r1d, r1d
    add       r3d, r3d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_ChromaCoeffV + r4]
%endif

    mov       r4d, %2/4

.loopH:
    PROCESS_CHROMA_SP_W4_4R

    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3

    movlps    [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movlps    [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [4 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 2 * 4

    PROCESS_CHROMA_SP_W2_4R r6

    psrad     m0, 6
    psrad     m2, 6

    packssdw  m0, m2

    movd      [r2], m0
    pextrd    [r2 + r3], m0, 1
    lea       r2, [r2 + 2 * r3]
    pextrd    [r2], m0, 2
    pextrd    [r2 + r3], m0, 3

    sub       r0, 2 * 4
    lea       r2, [r2 + 2 * r3 - 2 * 4]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SS_W6_H4 6, 8

FILTER_VER_CHROMA_SS_W6_H4 6, 16


;----------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_ss_8x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;----------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SS_W8_H2 2
INIT_XMM sse2
cglobal interp_4tap_vert_ss_%1x%2, 5, 6, 7

    add       r1d, r1d
    add       r3d, r3d
    sub       r0, r1
    shl       r4d, 5

%ifdef PIC
    lea       r5, [tab_ChromaCoeffV]
    lea       r5, [r5 + r4]
%else
    lea       r5, [tab_ChromaCoeffV + r4]
%endif

    mov       r4d, %2/2
.loopH:
    PROCESS_CHROMA_SP_W8_2R

    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3

    movu      [r2], m0
    movu      [r2 + r3], m2

    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_SS_W8_H2 8, 2
FILTER_VER_CHROMA_SS_W8_H2 8, 4
FILTER_VER_CHROMA_SS_W8_H2 8, 6
FILTER_VER_CHROMA_SS_W8_H2 8, 8
FILTER_VER_CHROMA_SS_W8_H2 8, 16
FILTER_VER_CHROMA_SS_W8_H2 8, 32

FILTER_VER_CHROMA_SS_W8_H2 8, 12
FILTER_VER_CHROMA_SS_W8_H2 8, 64

;-----------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ss_%1x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_SS 2
INIT_XMM sse2
cglobal interp_8tap_vert_ss_%1x%2, 5, 7, 7 ,0-gprsize

    add        r1d, r1d
    add        r3d, r3d
    lea        r5, [3 * r1]
    sub        r0, r5
    shl        r4d, 6

%ifdef PIC
    lea        r5, [tab_LumaCoeffV]
    lea        r6, [r5 + r4]
%else
    lea        r6, [tab_LumaCoeffV + r4]
%endif

    mov        dword [rsp], %2/4
.loopH:
    mov        r4d, (%1/4)
.loopW:
    movq       m0, [r0]
    movq       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r6 + 0 *16]                ;m0=[0+1]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m1, m4                          ;m1=[1 2]
    pmaddwd    m1, [r6 + 0 *16]                ;m1=[1+2]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[2 3]
    pmaddwd    m2, m4, [r6 + 0 *16]            ;m2=[2+3]  Row3
    pmaddwd    m4, [r6 + 1 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[3 4]
    pmaddwd    m3, m5, [r6 + 0 *16]            ;m3=[3+4]  Row4
    pmaddwd    m5, [r6 + 1 * 16]
    paddd      m1, m5                          ;m1 = [1+2+3+4]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[4 5]
    pmaddwd    m6, m4, [r6 + 1 * 16]
    paddd      m2, m6                          ;m2=[2+3+4+5]  Row3
    pmaddwd    m4, [r6 + 2 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3+4+5]  Row1

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[5 6]
    pmaddwd    m6, m5, [r6 + 1 * 16]
    paddd      m3, m6                          ;m3=[3+4+5+6]  Row4
    pmaddwd    m5, [r6 + 2 * 16]
    paddd      m1, m5                          ;m1=[1+2+3+4+5+6]  Row2

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[6 7]
    pmaddwd    m6, m4, [r6 + 2 * 16]
    paddd      m2, m6                          ;m2=[2+3+4+5+6+7]  Row3
    pmaddwd    m4, [r6 + 3 * 16]
    paddd      m0, m4                          ;m0=[0+1+2+3+4+5+6+7]  Row1 end
    psrad      m0, 6

    lea        r0, [r0 + 2 * r1]
    movq       m4, [r0]
    punpcklwd  m5, m4                          ;m5=[7 8]
    pmaddwd    m6, m5, [r6 + 2 * 16]
    paddd      m3, m6                          ;m3=[3+4+5+6+7+8]  Row4
    pmaddwd    m5, [r6 + 3 * 16]
    paddd      m1, m5                          ;m1=[1+2+3+4+5+6+7+8]  Row2 end
    psrad      m1, 6

    packssdw   m0, m1

    movlps     [r2], m0
    movhps     [r2 + r3], m0

    movq       m5, [r0 + r1]
    punpcklwd  m4, m5                          ;m4=[8 9]
    pmaddwd    m4, [r6 + 3 * 16]
    paddd      m2, m4                          ;m2=[2+3+4+5+6+7+8+9]  Row3 end
    psrad      m2, 6

    movq       m4, [r0 + 2 * r1]
    punpcklwd  m5, m4                          ;m5=[9 10]
    pmaddwd    m5, [r6 + 3 * 16]
    paddd      m3, m5                          ;m3=[3+4+5+6+7+8+9+10]  Row4 end
    psrad      m3, 6

    packssdw   m2, m3

    movlps     [r2 + 2 * r3], m2
    lea        r5, [3 * r3]
    movhps     [r2 + r5], m2

    lea        r5, [8 * r1 - 2 * 4]
    sub        r0, r5
    add        r2, 2 * 4

    dec        r4d
    jnz        .loopW

    lea        r0, [r0 + 4 * r1 - 2 * %1]
    lea        r2, [r2 + 4 * r3 - 2 * %1]

    dec        dword [rsp]
    jnz        .loopH

    RET
%endmacro

    FILTER_VER_LUMA_SS 4, 4
    FILTER_VER_LUMA_SS 8, 8
    FILTER_VER_LUMA_SS 8, 4
    FILTER_VER_LUMA_SS 4, 8
    FILTER_VER_LUMA_SS 16, 16
    FILTER_VER_LUMA_SS 16, 8
    FILTER_VER_LUMA_SS 8, 16
    FILTER_VER_LUMA_SS 16, 12
    FILTER_VER_LUMA_SS 12, 16
    FILTER_VER_LUMA_SS 16, 4
    FILTER_VER_LUMA_SS 4, 16
    FILTER_VER_LUMA_SS 32, 32
    FILTER_VER_LUMA_SS 32, 16
    FILTER_VER_LUMA_SS 16, 32
    FILTER_VER_LUMA_SS 32, 24
    FILTER_VER_LUMA_SS 24, 32
    FILTER_VER_LUMA_SS 32, 8
    FILTER_VER_LUMA_SS 8, 32
    FILTER_VER_LUMA_SS 64, 64
    FILTER_VER_LUMA_SS 64, 32
    FILTER_VER_LUMA_SS 32, 64
    FILTER_VER_LUMA_SS 64, 48
    FILTER_VER_LUMA_SS 48, 64
    FILTER_VER_LUMA_SS 64, 16
    FILTER_VER_LUMA_SS 16, 64

%macro FILTER_VER_LUMA_AVX2_4x4 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_4x4, 4, 6, 7
    mov             r4d, r4m
    add             r1d, r1d
    shl             r4d, 7

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %1,sp
    mova            m6, [pd_526336]
%else
    add             r3d, r3d
%endif

    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m5, m4, [r5 + 2 * mmsize]
    pmaddwd         m4, [r5 + 1 * mmsize]
    paddd           m0, m5
    paddd           m2, m4
    movq            xm3, [r0 + r4]
    punpcklwd       xm1, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m1, m1, xm3, 1                  ; m1 = [8 7 7 6]
    pmaddwd         m5, m1, [r5 + 3 * mmsize]
    pmaddwd         m1, [r5 + 2 * mmsize]
    paddd           m0, m5
    paddd           m2, m1
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + 2 * r1]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [A 9 9 8]
    pmaddwd         m4, [r5 + 3 * mmsize]
    paddd           m2, m4

%ifidn %1,sp
    paddd           m0, m6
    paddd           m2, m6
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2
    vextracti128    xm2, m0, 1
    lea             r4, [r3 * 3]

%ifidn %1,sp
    packuswb        xm0, xm2
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 2
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r4], xm0, 3
%else
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r4], xm2
%endif
    RET
%endmacro

FILTER_VER_LUMA_AVX2_4x4 sp
FILTER_VER_LUMA_AVX2_4x4 ss

%macro FILTER_VER_LUMA_AVX2_4x8 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_4x8, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]

    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m5, m4, [r5 + 2 * mmsize]
    paddd           m0, m5
    pmaddwd         m5, m4, [r5 + 1 * mmsize]
    paddd           m2, m5
    pmaddwd         m4, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm1, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm6, [r0]
    punpcklwd       xm3, xm6
    vinserti128     m1, m1, xm3, 1                  ; m1 = [8 7 7 6]
    pmaddwd         m5, m1, [r5 + 3 * mmsize]
    paddd           m0, m5
    pmaddwd         m5, m1, [r5 + 2 * mmsize]
    paddd           m2, m5
    pmaddwd         m5, m1, [r5 + 1 * mmsize]
    paddd           m4, m5
    pmaddwd         m1, [r5]
    movq            xm3, [r0 + r1]
    punpcklwd       xm6, xm3
    movq            xm5, [r0 + 2 * r1]
    punpcklwd       xm3, xm5
    vinserti128     m6, m6, xm3, 1                  ; m6 = [A 9 9 8]
    pmaddwd         m3, m6, [r5 + 3 * mmsize]
    paddd           m2, m3
    pmaddwd         m3, m6, [r5 + 2 * mmsize]
    paddd           m4, m3
    pmaddwd         m6, [r5 + 1 * mmsize]
    paddd           m1, m6

%ifidn %1,sp
    paddd           m0, m7
    paddd           m2, m7
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2

    movq            xm3, [r0 + r4]
    punpcklwd       xm5, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm2, [r0]
    punpcklwd       xm3, xm2
    vinserti128     m5, m5, xm3, 1                  ; m5 = [C B B A]
    pmaddwd         m3, m5, [r5 + 3 * mmsize]
    paddd           m4, m3
    pmaddwd         m5, [r5 + 2 * mmsize]
    paddd           m1, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm2, xm3
    movq            xm5, [r0 + 2 * r1]
    punpcklwd       xm3, xm5
    vinserti128     m2, m2, xm3, 1                  ; m2 = [E D D C]
    pmaddwd         m2, [r5 + 3 * mmsize]
    paddd           m1, m2

%ifidn %1,sp
    paddd           m4, m7
    paddd           m1, m7
    psrad           m4, 12
    psrad           m1, 12
%else
    psrad           m4, 6
    psrad           m1, 6
%endif
    packssdw        m4, m1

%ifidn %1,sp
    packuswb        m0, m4
    vextracti128    xm2, m0, 1
    movd            [r2], xm0
    movd            [r2 + r3], xm2
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r6], xm2, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm0, 2
    pextrd          [r2 + r3], xm2, 2
    pextrd          [r2 + r3 * 2], xm0, 3
    pextrd          [r2 + r6], xm2, 3
%else
    vextracti128    xm2, m0, 1
    vextracti128    xm1, m4, 1
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm2
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm4
    movq            [r2 + r3], xm1
    movhps          [r2 + r3 * 2], xm4
    movhps          [r2 + r6], xm1
%endif
    RET
%endmacro

FILTER_VER_LUMA_AVX2_4x8 sp
FILTER_VER_LUMA_AVX2_4x8 ss

%macro PROCESS_LUMA_AVX2_W4_16R 1
    movq            xm0, [r0]
    movq            xm1, [r0 + r1]
    punpcklwd       xm0, xm1
    movq            xm2, [r0 + r1 * 2]
    punpcklwd       xm1, xm2
    vinserti128     m0, m0, xm1, 1                  ; m0 = [2 1 1 0]
    pmaddwd         m0, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm2, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm4, [r0]
    punpcklwd       xm3, xm4
    vinserti128     m2, m2, xm3, 1                  ; m2 = [4 3 3 2]
    pmaddwd         m5, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m5
    movq            xm3, [r0 + r1]
    punpcklwd       xm4, xm3
    movq            xm1, [r0 + r1 * 2]
    punpcklwd       xm3, xm1
    vinserti128     m4, m4, xm3, 1                  ; m4 = [6 5 5 4]
    pmaddwd         m5, m4, [r5 + 2 * mmsize]
    paddd           m0, m5
    pmaddwd         m5, m4, [r5 + 1 * mmsize]
    paddd           m2, m5
    pmaddwd         m4, [r5]
    movq            xm3, [r0 + r4]
    punpcklwd       xm1, xm3
    lea             r0, [r0 + 4 * r1]
    movq            xm6, [r0]
    punpcklwd       xm3, xm6
    vinserti128     m1, m1, xm3, 1                  ; m1 = [8 7 7 6]
    pmaddwd         m5, m1, [r5 + 3 * mmsize]
    paddd           m0, m5
    pmaddwd         m5, m1, [r5 + 2 * mmsize]
    paddd           m2, m5
    pmaddwd         m5, m1, [r5 + 1 * mmsize]
    paddd           m4, m5
    pmaddwd         m1, [r5]
    movq            xm3, [r0 + r1]
    punpcklwd       xm6, xm3
    movq            xm5, [r0 + 2 * r1]
    punpcklwd       xm3, xm5
    vinserti128     m6, m6, xm3, 1                  ; m6 = [10 9 9 8]
    pmaddwd         m3, m6, [r5 + 3 * mmsize]
    paddd           m2, m3
    pmaddwd         m3, m6, [r5 + 2 * mmsize]
    paddd           m4, m3
    pmaddwd         m3, m6, [r5 + 1 * mmsize]
    paddd           m1, m3
    pmaddwd         m6, [r5]

%ifidn %1,sp
    paddd           m0, m7
    paddd           m2, m7
    psrad           m0, 12
    psrad           m2, 12
%else
    psrad           m0, 6
    psrad           m2, 6
%endif
    packssdw        m0, m2
    vextracti128    xm2, m0, 1
%ifidn %1,sp
    packuswb        xm0, xm2
    movd            [r2], xm0
    pextrd          [r2 + r3], xm0, 2
    pextrd          [r2 + r3 * 2], xm0, 1
    pextrd          [r2 + r6], xm0, 3
%else
    movq            [r2], xm0
    movq            [r2 + r3], xm2
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm2
%endif

    movq            xm2, [r0 + r4]
    punpcklwd       xm5, xm2
    lea             r0, [r0 + 4 * r1]
    movq            xm0, [r0]
    punpcklwd       xm2, xm0
    vinserti128     m5, m5, xm2, 1                  ; m5 = [12 11 11 10]
    pmaddwd         m2, m5, [r5 + 3 * mmsize]
    paddd           m4, m2
    pmaddwd         m2, m5, [r5 + 2 * mmsize]
    paddd           m1, m2
    pmaddwd         m2, m5, [r5 + 1 * mmsize]
    paddd           m6, m2
    pmaddwd         m5, [r5]
    movq            xm2, [r0 + r1]
    punpcklwd       xm0, xm2
    movq            xm3, [r0 + 2 * r1]
    punpcklwd       xm2, xm3
    vinserti128     m0, m0, xm2, 1                  ; m0 = [14 13 13 12]
    pmaddwd         m2, m0, [r5 + 3 * mmsize]
    paddd           m1, m2
    pmaddwd         m2, m0, [r5 + 2 * mmsize]
    paddd           m6, m2
    pmaddwd         m2, m0, [r5 + 1 * mmsize]
    paddd           m5, m2
    pmaddwd         m0, [r5]

%ifidn %1,sp
    paddd           m4, m7
    paddd           m1, m7
    psrad           m4, 12
    psrad           m1, 12
%else
    psrad           m4, 6
    psrad           m1, 6
%endif
    packssdw        m4, m1
    vextracti128    xm1, m4, 1
    lea             r2, [r2 + r3 * 4]
%ifidn %1,sp
    packuswb        xm4, xm1
    movd            [r2], xm4
    pextrd          [r2 + r3], xm4, 2
    pextrd          [r2 + r3 * 2], xm4, 1
    pextrd          [r2 + r6], xm4, 3
%else
    movq            [r2], xm4
    movq            [r2 + r3], xm1
    movhps          [r2 + r3 * 2], xm4
    movhps          [r2 + r6], xm1
%endif

    movq            xm4, [r0 + r4]
    punpcklwd       xm3, xm4
    lea             r0, [r0 + 4 * r1]
    movq            xm1, [r0]
    punpcklwd       xm4, xm1
    vinserti128     m3, m3, xm4, 1                  ; m3 = [16 15 15 14]
    pmaddwd         m4, m3, [r5 + 3 * mmsize]
    paddd           m6, m4
    pmaddwd         m4, m3, [r5 + 2 * mmsize]
    paddd           m5, m4
    pmaddwd         m4, m3, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m3, [r5]
    movq            xm4, [r0 + r1]
    punpcklwd       xm1, xm4
    movq            xm2, [r0 + 2 * r1]
    punpcklwd       xm4, xm2
    vinserti128     m1, m1, xm4, 1                  ; m1 = [18 17 17 16]
    pmaddwd         m4, m1, [r5 + 3 * mmsize]
    paddd           m5, m4
    pmaddwd         m4, m1, [r5 + 2 * mmsize]
    paddd           m0, m4
    pmaddwd         m1, [r5 + 1 * mmsize]
    paddd           m3, m1
    movq            xm4, [r0 + r4]
    punpcklwd       xm2, xm4
    lea             r0, [r0 + 4 * r1]
    movq            xm1, [r0]
    punpcklwd       xm4, xm1
    vinserti128     m2, m2, xm4, 1                  ; m2 = [20 19 19 18]
    pmaddwd         m4, m2, [r5 + 3 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5 + 2 * mmsize]
    paddd           m3, m2
    movq            xm4, [r0 + r1]
    punpcklwd       xm1, xm4
    movq            xm2, [r0 + 2 * r1]
    punpcklwd       xm4, xm2
    vinserti128     m1, m1, xm4, 1                  ; m1 = [22 21 21 20]
    pmaddwd         m1, [r5 + 3 * mmsize]
    paddd           m3, m1

%ifidn %1,sp
    paddd           m6, m7
    paddd           m5, m7
    paddd           m0, m7
    paddd           m3, m7
    psrad           m6, 12
    psrad           m5, 12
    psrad           m0, 12
    psrad           m3, 12
%else
    psrad           m6, 6
    psrad           m5, 6
    psrad           m0, 6
    psrad           m3, 6
%endif
    packssdw        m6, m5
    packssdw        m0, m3
    lea             r2, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m6, m0
    vextracti128    xm0, m6, 1
    movd            [r2], xm6
    movd            [r2 + r3], xm0
    pextrd          [r2 + r3 * 2], xm6, 1
    pextrd          [r2 + r6], xm0, 1
    lea             r2, [r2 + r3 * 4]
    pextrd          [r2], xm6, 2
    pextrd          [r2 + r3], xm0, 2
    pextrd          [r2 + r3 * 2], xm6, 3
    pextrd          [r2 + r6], xm0, 3
%else
    vextracti128    xm5, m6, 1
    vextracti128    xm3, m0, 1
    movq            [r2], xm6
    movq            [r2 + r3], xm5
    movhps          [r2 + r3 * 2], xm6
    movhps          [r2 + r6], xm5
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm0
    movq            [r2 + r3], xm3
    movhps          [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm3
%endif
%endmacro

%macro FILTER_VER_LUMA_AVX2_4x16 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_4x16, 4, 7, 8
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    PROCESS_LUMA_AVX2_W4_16R %1
    RET
%endmacro

FILTER_VER_LUMA_AVX2_4x16 sp
FILTER_VER_LUMA_AVX2_4x16 ss

%macro FILTER_VER_LUMA_S_AVX2_8x8 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_8x8, 4, 6, 12
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %1,sp
    mova            m11, [pd_526336]
%else
    add             r3d, r3d
%endif

    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    pmaddwd         m2, [r5]
    paddd           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    pmaddwd         m3, [r5]
    paddd           m1, m5
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 2 * mmsize]
    paddd           m1, m7
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    pmaddwd         m5, [r5]
    paddd           m3, m7
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 3 * mmsize]
    paddd           m0, m8
    pmaddwd         m8, m6, [r5 + 2 * mmsize]
    paddd           m2, m8
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    pmaddwd         m6, [r5]
    paddd           m4, m8
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhwd       xm9, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddwd         m9, m7, [r5 + 3 * mmsize]
    paddd           m1, m9
    pmaddwd         m9, m7, [r5 + 2 * mmsize]
    paddd           m3, m9
    pmaddwd         m9, m7, [r5 + 1 * mmsize]
    pmaddwd         m7, [r5]
    paddd           m5, m9
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m2, m10
    pmaddwd         m10, m8, [r5 + 2 * mmsize]
    pmaddwd         m8, [r5 + 1 * mmsize]
    paddd           m4, m10
    paddd           m6, m8
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhwd       xm8, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm8, 1
    pmaddwd         m8, m9, [r5 + 3 * mmsize]
    paddd           m3, m8
    pmaddwd         m8, m9, [r5 + 2 * mmsize]
    pmaddwd         m9, [r5 + 1 * mmsize]
    paddd           m5, m8
    paddd           m7, m9
    movu            xm8, [r0 + r4]                  ; m8 = row 11
    punpckhwd       xm9, xm10, xm8
    punpcklwd       xm10, xm8
    vinserti128     m10, m10, xm9, 1
    pmaddwd         m9, m10, [r5 + 3 * mmsize]
    pmaddwd         m10, [r5 + 2 * mmsize]
    paddd           m4, m9
    paddd           m6, m10

    lea             r4, [r3 * 3]
%ifidn %1,sp
    paddd           m0, m11
    paddd           m1, m11
    paddd           m2, m11
    paddd           m3, m11
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m1, [interp8_hps_shuf]
    vpermd          m0, m1, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
%endif

    lea             r0, [r0 + r1 * 4]
    movu            xm9, [r0]                       ; m9 = row 12
    punpckhwd       xm3, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm3, 1
    pmaddwd         m3, m8, [r5 + 3 * mmsize]
    pmaddwd         m8, [r5 + 2 * mmsize]
    paddd           m5, m3
    paddd           m7, m8
    movu            xm3, [r0 + r1]                  ; m3 = row 13
    punpckhwd       xm0, xm9, xm3
    punpcklwd       xm9, xm3
    vinserti128     m9, m9, xm0, 1
    pmaddwd         m9, [r5 + 3 * mmsize]
    paddd           m6, m9
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhwd       xm9, xm3, xm0
    punpcklwd       xm3, xm0
    vinserti128     m3, m3, xm9, 1
    pmaddwd         m3, [r5 + 3 * mmsize]
    paddd           m7, m3

%ifidn %1,sp
    paddd           m4, m11
    paddd           m5, m11
    paddd           m6, m11
    paddd           m7, m11
    psrad           m4, 12
    psrad           m5, 12
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m4, 6
    psrad           m5, 6
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m4, m5
    packssdw        m6, m7
    lea             r2, [r2 + r3 * 4]
%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m1, m4
    vextracti128    xm6, m4, 1
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm6
    movhps          [r2 + r4], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm5, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm5
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r4], xm7
%endif
    RET
%endif
%endmacro

FILTER_VER_LUMA_S_AVX2_8x8 sp
FILTER_VER_LUMA_S_AVX2_8x8 ss

%macro FILTER_VER_LUMA_S_AVX2_8xN 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_8x%2, 4, 9, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    lea             r7, [r1 * 4]
    mov             r8d, %2 / 16
.loopH:
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 2 * mmsize]
    paddd           m1, m7
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m3, m7
    pmaddwd         m5, [r5]
    movu            xm7, [r0 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 3 * mmsize]
    paddd           m0, m8
    pmaddwd         m8, m6, [r5 + 2 * mmsize]
    paddd           m2, m8
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    paddd           m4, m8
    pmaddwd         m6, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm8, [r0]                       ; m8 = row 8
    punpckhwd       xm9, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddwd         m9, m7, [r5 + 3 * mmsize]
    paddd           m1, m9
    pmaddwd         m9, m7, [r5 + 2 * mmsize]
    paddd           m3, m9
    pmaddwd         m9, m7, [r5 + 1 * mmsize]
    paddd           m5, m9
    pmaddwd         m7, [r5]
    movu            xm9, [r0 + r1]                  ; m9 = row 9
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m2, m10
    pmaddwd         m10, m8, [r5 + 2 * mmsize]
    paddd           m4, m10
    pmaddwd         m10, m8, [r5 + 1 * mmsize]
    paddd           m6, m10
    pmaddwd         m8, [r5]
    movu            xm10, [r0 + r1 * 2]             ; m10 = row 10
    punpckhwd       xm11, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddwd         m11, m9, [r5 + 3 * mmsize]
    paddd           m3, m11
    pmaddwd         m11, m9, [r5 + 2 * mmsize]
    paddd           m5, m11
    pmaddwd         m11, m9, [r5 + 1 * mmsize]
    paddd           m7, m11
    pmaddwd         m9, [r5]
    movu            xm11, [r0 + r4]                 ; m11 = row 11
    punpckhwd       xm12, xm10, xm11
    punpcklwd       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddwd         m12, m10, [r5 + 3 * mmsize]
    paddd           m4, m12
    pmaddwd         m12, m10, [r5 + 2 * mmsize]
    paddd           m6, m12
    pmaddwd         m12, m10, [r5 + 1 * mmsize]
    paddd           m8, m12
    pmaddwd         m10, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm12, [r0]                      ; m12 = row 12
    punpckhwd       xm13, xm11, xm12
    punpcklwd       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddwd         m13, m11, [r5 + 3 * mmsize]
    paddd           m5, m13
    pmaddwd         m13, m11, [r5 + 2 * mmsize]
    paddd           m7, m13
    pmaddwd         m13, m11, [r5 + 1 * mmsize]
    paddd           m9, m13
    pmaddwd         m11, [r5]

%ifidn %1,sp
    paddd           m0, m14
    paddd           m1, m14
    paddd           m2, m14
    paddd           m3, m14
    paddd           m4, m14
    paddd           m5, m14
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
    packssdw        m4, m5
%ifidn %1,sp
    packuswb        m0, m2
    mova            m1, [interp8_hps_shuf]
    vpermd          m0, m1, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif

    movu            xm13, [r0 + r1]                 ; m13 = row 13
    punpckhwd       xm0, xm12, xm13
    punpcklwd       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddwd         m0, m12, [r5 + 3 * mmsize]
    paddd           m6, m0
    pmaddwd         m0, m12, [r5 + 2 * mmsize]
    paddd           m8, m0
    pmaddwd         m0, m12, [r5 + 1 * mmsize]
    paddd           m10, m0
    pmaddwd         m12, [r5]
    movu            xm0, [r0 + r1 * 2]              ; m0 = row 14
    punpckhwd       xm2, xm13, xm0
    punpcklwd       xm13, xm0
    vinserti128     m13, m13, xm2, 1
    pmaddwd         m2, m13, [r5 + 3 * mmsize]
    paddd           m7, m2
    pmaddwd         m2, m13, [r5 + 2 * mmsize]
    paddd           m9, m2
    pmaddwd         m2, m13, [r5 + 1 * mmsize]
    paddd           m11, m2
    pmaddwd         m13, [r5]

%ifidn %1,sp
    paddd           m6, m14
    paddd           m7, m14
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m6, m7
    lea             r2, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m1, m4
    vextracti128    xm6, m4, 1
    movq            [r2], xm4
    movhps          [r2 + r3], xm4
    movq            [r2 + r3 * 2], xm6
    movhps          [r2 + r6], xm6
%else
    vpermq          m6, m6, 11011000b
    vpermq          m4, m4, 11011000b
    vextracti128    xm1, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r2], xm4
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm6
    movu            [r2 + r6], xm7
%endif

    movu            xm6, [r0 + r4]                  ; m6 = row 15
    punpckhwd       xm5, xm0, xm6
    punpcklwd       xm0, xm6
    vinserti128     m0, m0, xm5, 1
    pmaddwd         m5, m0, [r5 + 3 * mmsize]
    paddd           m8, m5
    pmaddwd         m5, m0, [r5 + 2 * mmsize]
    paddd           m10, m5
    pmaddwd         m5, m0, [r5 + 1 * mmsize]
    paddd           m12, m5
    pmaddwd         m0, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm2, [r0]                       ; m2 = row 16
    punpckhwd       xm3, xm6, xm2
    punpcklwd       xm6, xm2
    vinserti128     m6, m6, xm3, 1
    pmaddwd         m3, m6, [r5 + 3 * mmsize]
    paddd           m9, m3
    pmaddwd         m3, m6, [r5 + 2 * mmsize]
    paddd           m11, m3
    pmaddwd         m3, m6, [r5 + 1 * mmsize]
    paddd           m13, m3
    pmaddwd         m6, [r5]
    movu            xm3, [r0 + r1]                  ; m3 = row 17
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 3 * mmsize]
    paddd           m10, m4
    pmaddwd         m4, m2, [r5 + 2 * mmsize]
    paddd           m12, m4
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 18
    punpckhwd       xm2, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm2, 1
    pmaddwd         m2, m3, [r5 + 3 * mmsize]
    paddd           m11, m2
    pmaddwd         m2, m3, [r5 + 2 * mmsize]
    paddd           m13, m2
    pmaddwd         m3, [r5 + 1 * mmsize]
    paddd           m6, m3
    movu            xm2, [r0 + r4]                  ; m2 = row 19
    punpckhwd       xm7, xm4, xm2
    punpcklwd       xm4, xm2
    vinserti128     m4, m4, xm7, 1
    pmaddwd         m7, m4, [r5 + 3 * mmsize]
    paddd           m12, m7
    pmaddwd         m4, [r5 + 2 * mmsize]
    paddd           m0, m4
    lea             r0, [r0 + r1 * 4]
    movu            xm7, [r0]                       ; m7 = row 20
    punpckhwd       xm3, xm2, xm7
    punpcklwd       xm2, xm7
    vinserti128     m2, m2, xm3, 1
    pmaddwd         m3, m2, [r5 + 3 * mmsize]
    paddd           m13, m3
    pmaddwd         m2, [r5 + 2 * mmsize]
    paddd           m6, m2
    movu            xm3, [r0 + r1]                  ; m3 = row 21
    punpckhwd       xm2, xm7, xm3
    punpcklwd       xm7, xm3
    vinserti128     m7, m7, xm2, 1
    pmaddwd         m7, [r5 + 3 * mmsize]
    paddd           m0, m7
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 22
    punpckhwd       xm7, xm3, xm2
    punpcklwd       xm3, xm2
    vinserti128     m3, m3, xm7, 1
    pmaddwd         m3, [r5 + 3 * mmsize]
    paddd           m6, m3

%ifidn %1,sp
    paddd           m8, m14
    paddd           m9, m14
    paddd           m10, m14
    paddd           m11, m14
    paddd           m12, m14
    paddd           m13, m14
    paddd           m0, m14
    paddd           m6, m14
    psrad           m8, 12
    psrad           m9, 12
    psrad           m10, 12
    psrad           m11, 12
    psrad           m12, 12
    psrad           m13, 12
    psrad           m0, 12
    psrad           m6, 12
%else
    psrad           m8, 6
    psrad           m9, 6
    psrad           m10, 6
    psrad           m11, 6
    psrad           m12, 6
    psrad           m13, 6
    psrad           m0, 6
    psrad           m6, 6
%endif
    packssdw        m8, m9
    packssdw        m10, m11
    packssdw        m12, m13
    packssdw        m0, m6
    lea             r2, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m8, m10
    packuswb        m12, m0
    vpermd          m8, m1, m8
    vpermd          m12, m1, m12
    vextracti128    xm10, m8, 1
    vextracti128    xm0, m12, 1
    movq            [r2], xm8
    movhps          [r2 + r3], xm8
    movq            [r2 + r3 * 2], xm10
    movhps          [r2 + r6], xm10
    lea             r2, [r2 + r3 * 4]
    movq            [r2], xm12
    movhps          [r2 + r3], xm12
    movq            [r2 + r3 * 2], xm0
    movhps          [r2 + r6], xm0
%else
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm6, m0, 1
    movu            [r2], xm8
    movu            [r2 + r3], xm9
    movu            [r2 + r3 * 2], xm10
    movu            [r2 + r6], xm11
    lea             r2, [r2 + r3 * 4]
    movu            [r2], xm12
    movu            [r2 + r3], xm13
    movu            [r2 + r3 * 2], xm0
    movu            [r2 + r6], xm6
%endif

    lea             r2, [r2 + r3 * 4]
    sub             r0, r7
    dec             r8d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_S_AVX2_8xN sp, 16
FILTER_VER_LUMA_S_AVX2_8xN sp, 32
FILTER_VER_LUMA_S_AVX2_8xN ss, 16
FILTER_VER_LUMA_S_AVX2_8xN ss, 32

%macro PROCESS_LUMA_S_AVX2_W8_4R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r0, [r0 + r1 * 4]
    movu            xm4, [r0]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r0 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m4, [r5 + 1 * mmsize]
    paddd           m2, m4
    movu            xm6, [r0 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm4, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm4, 1
    pmaddwd         m4, m5, [r5 + 2 * mmsize]
    paddd           m1, m4
    pmaddwd         m5, [r5 + 1 * mmsize]
    paddd           m3, m5
    movu            xm4, [r0 + r4]                  ; m4 = row 7
    punpckhwd       xm5, xm6, xm4
    punpcklwd       xm6, xm4
    vinserti128     m6, m6, xm5, 1
    pmaddwd         m5, m6, [r5 + 3 * mmsize]
    paddd           m0, m5
    pmaddwd         m6, [r5 + 2 * mmsize]
    paddd           m2, m6
    lea             r0, [r0 + r1 * 4]
    movu            xm5, [r0]                       ; m5 = row 8
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 3 * mmsize]
    paddd           m1, m6
    pmaddwd         m4, [r5 + 2 * mmsize]
    paddd           m3, m4
    movu            xm6, [r0 + r1]                  ; m6 = row 9
    punpckhwd       xm4, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm4, 1
    pmaddwd         m5, [r5 + 3 * mmsize]
    paddd           m2, m5
    movu            xm4, [r0 + r1 * 2]              ; m4 = row 10
    punpckhwd       xm5, xm6, xm4
    punpcklwd       xm6, xm4
    vinserti128     m6, m6, xm5, 1
    pmaddwd         m6, [r5 + 3 * mmsize]
    paddd           m3, m6

%ifidn %1,sp
    paddd           m0, m7
    paddd           m1, m7
    paddd           m2, m7
    paddd           m3, m7
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
%ifidn %1,sp
    packuswb        m0, m2
    mova            m4, [interp8_hps_shuf]
    vpermd          m0, m4, m0
    vextracti128    xm2, m0, 1
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
%endif
%endmacro

%macro FILTER_VER_LUMA_S_AVX2_8x4 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_8x4, 4, 6, 8
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif

    PROCESS_LUMA_S_AVX2_W8_4R %1
    lea             r4, [r3 * 3]
%ifidn %1,sp
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r4], xm2
%else
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r4], xm3
%endif
    RET
%endmacro

FILTER_VER_LUMA_S_AVX2_8x4 sp
FILTER_VER_LUMA_S_AVX2_8x4 ss

%macro PROCESS_LUMA_AVX2_W8_16R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 2 * mmsize]
    paddd           m1, m7
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m3, m7
    pmaddwd         m5, [r5]
    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 3 * mmsize]
    paddd           m0, m8
    pmaddwd         m8, m6, [r5 + 2 * mmsize]
    paddd           m2, m8
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    paddd           m4, m8
    pmaddwd         m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhwd       xm9, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddwd         m9, m7, [r5 + 3 * mmsize]
    paddd           m1, m9
    pmaddwd         m9, m7, [r5 + 2 * mmsize]
    paddd           m3, m9
    pmaddwd         m9, m7, [r5 + 1 * mmsize]
    paddd           m5, m9
    pmaddwd         m7, [r5]
    movu            xm9, [r7 + r1]                  ; m9 = row 9
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m2, m10
    pmaddwd         m10, m8, [r5 + 2 * mmsize]
    paddd           m4, m10
    pmaddwd         m10, m8, [r5 + 1 * mmsize]
    paddd           m6, m10
    pmaddwd         m8, [r5]
    movu            xm10, [r7 + r1 * 2]             ; m10 = row 10
    punpckhwd       xm11, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddwd         m11, m9, [r5 + 3 * mmsize]
    paddd           m3, m11
    pmaddwd         m11, m9, [r5 + 2 * mmsize]
    paddd           m5, m11
    pmaddwd         m11, m9, [r5 + 1 * mmsize]
    paddd           m7, m11
    pmaddwd         m9, [r5]
    movu            xm11, [r7 + r4]                 ; m11 = row 11
    punpckhwd       xm12, xm10, xm11
    punpcklwd       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddwd         m12, m10, [r5 + 3 * mmsize]
    paddd           m4, m12
    pmaddwd         m12, m10, [r5 + 2 * mmsize]
    paddd           m6, m12
    pmaddwd         m12, m10, [r5 + 1 * mmsize]
    paddd           m8, m12
    pmaddwd         m10, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm12, [r7]                      ; m12 = row 12
    punpckhwd       xm13, xm11, xm12
    punpcklwd       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddwd         m13, m11, [r5 + 3 * mmsize]
    paddd           m5, m13
    pmaddwd         m13, m11, [r5 + 2 * mmsize]
    paddd           m7, m13
    pmaddwd         m13, m11, [r5 + 1 * mmsize]
    paddd           m9, m13
    pmaddwd         m11, [r5]

%ifidn %1,sp
    paddd           m0, m14
    paddd           m1, m14
    paddd           m2, m14
    paddd           m3, m14
    paddd           m4, m14
    paddd           m5, m14
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
    packssdw        m4, m5
%ifidn %1,sp
    packuswb        m0, m2
    mova            m5, [interp8_hps_shuf]
    vpermd          m0, m5, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif

    movu            xm13, [r7 + r1]                 ; m13 = row 13
    punpckhwd       xm0, xm12, xm13
    punpcklwd       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddwd         m0, m12, [r5 + 3 * mmsize]
    paddd           m6, m0
    pmaddwd         m0, m12, [r5 + 2 * mmsize]
    paddd           m8, m0
    pmaddwd         m0, m12, [r5 + 1 * mmsize]
    paddd           m10, m0
    pmaddwd         m12, [r5]
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 14
    punpckhwd       xm1, xm13, xm0
    punpcklwd       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddwd         m1, m13, [r5 + 3 * mmsize]
    paddd           m7, m1
    pmaddwd         m1, m13, [r5 + 2 * mmsize]
    paddd           m9, m1
    pmaddwd         m1, m13, [r5 + 1 * mmsize]
    paddd           m11, m1
    pmaddwd         m13, [r5]

%ifidn %1,sp
    paddd           m6, m14
    paddd           m7, m14
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m6, m7
    lea             r8, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m5, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm1, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%endif

    movu            xm1, [r7 + r4]                  ; m1 = row 15
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m2, m0, [r5 + 3 * mmsize]
    paddd           m8, m2
    pmaddwd         m2, m0, [r5 + 2 * mmsize]
    paddd           m10, m2
    pmaddwd         m2, m0, [r5 + 1 * mmsize]
    paddd           m12, m2
    pmaddwd         m0, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm2, [r7]                       ; m2 = row 16
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m3, m1, [r5 + 3 * mmsize]
    paddd           m9, m3
    pmaddwd         m3, m1, [r5 + 2 * mmsize]
    paddd           m11, m3
    pmaddwd         m3, m1, [r5 + 1 * mmsize]
    paddd           m13, m3
    pmaddwd         m1, [r5]
    movu            xm3, [r7 + r1]                  ; m3 = row 17
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 3 * mmsize]
    paddd           m10, m4
    pmaddwd         m4, m2, [r5 + 2 * mmsize]
    paddd           m12, m4
    pmaddwd         m2, [r5 + 1 * mmsize]
    paddd           m0, m2
    movu            xm4, [r7 + r1 * 2]              ; m4 = row 18
    punpckhwd       xm2, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm2, 1
    pmaddwd         m2, m3, [r5 + 3 * mmsize]
    paddd           m11, m2
    pmaddwd         m2, m3, [r5 + 2 * mmsize]
    paddd           m13, m2
    pmaddwd         m3, [r5 + 1 * mmsize]
    paddd           m1, m3
    movu            xm2, [r7 + r4]                  ; m2 = row 19
    punpckhwd       xm6, xm4, xm2
    punpcklwd       xm4, xm2
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 3 * mmsize]
    paddd           m12, m6
    pmaddwd         m4, [r5 + 2 * mmsize]
    paddd           m0, m4
    lea             r7, [r7 + r1 * 4]
    movu            xm6, [r7]                       ; m6 = row 20
    punpckhwd       xm7, xm2, xm6
    punpcklwd       xm2, xm6
    vinserti128     m2, m2, xm7, 1
    pmaddwd         m7, m2, [r5 + 3 * mmsize]
    paddd           m13, m7
    pmaddwd         m2, [r5 + 2 * mmsize]
    paddd           m1, m2
    movu            xm7, [r7 + r1]                  ; m7 = row 21
    punpckhwd       xm2, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm2, 1
    pmaddwd         m6, [r5 + 3 * mmsize]
    paddd           m0, m6
    movu            xm2, [r7 + r1 * 2]              ; m2 = row 22
    punpckhwd       xm3, xm7, xm2
    punpcklwd       xm7, xm2
    vinserti128     m7, m7, xm3, 1
    pmaddwd         m7, [r5 + 3 * mmsize]
    paddd           m1, m7

%ifidn %1,sp
    paddd           m8, m14
    paddd           m9, m14
    paddd           m10, m14
    paddd           m11, m14
    paddd           m12, m14
    paddd           m13, m14
    paddd           m0, m14
    paddd           m1, m14
    psrad           m8, 12
    psrad           m9, 12
    psrad           m10, 12
    psrad           m11, 12
    psrad           m12, 12
    psrad           m13, 12
    psrad           m0, 12
    psrad           m1, 12
%else
    psrad           m8, 6
    psrad           m9, 6
    psrad           m10, 6
    psrad           m11, 6
    psrad           m12, 6
    psrad           m13, 6
    psrad           m0, 6
    psrad           m1, 6
%endif
    packssdw        m8, m9
    packssdw        m10, m11
    packssdw        m12, m13
    packssdw        m0, m1
    lea             r8, [r8 + r3 * 4]

%ifidn %1,sp
    packuswb        m8, m10
    packuswb        m12, m0
    vpermd          m8, m5, m8
    vpermd          m12, m5, m12
    vextracti128    xm10, m8, 1
    vextracti128    xm0, m12, 1
    movq            [r8], xm8
    movhps          [r8 + r3], xm8
    movq            [r8 + r3 * 2], xm10
    movhps          [r8 + r6], xm10
    lea             r8, [r8 + r3 * 4]
    movq            [r8], xm12
    movhps          [r8 + r3], xm12
    movq            [r8 + r3 * 2], xm0
    movhps          [r8 + r6], xm0
%else
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vpermq          m12, m12, 11011000b
    vpermq          m0, m0, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    vextracti128    xm13, m12, 1
    vextracti128    xm1, m0, 1
    movu            [r8], xm8
    movu            [r8 + r3], xm9
    movu            [r8 + r3 * 2], xm10
    movu            [r8 + r6], xm11
    lea             r8, [r8 + r3 * 4]
    movu            [r8], xm12
    movu            [r8 + r3], xm13
    movu            [r8 + r3 * 2], xm0
    movu            [r8 + r6], xm1
%endif
%endmacro

%macro FILTER_VER_LUMA_AVX2_Nx16 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_%2x16, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, %2 / 8
.loopW:
    PROCESS_LUMA_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_Nx16 sp, 16
FILTER_VER_LUMA_AVX2_Nx16 sp, 32
FILTER_VER_LUMA_AVX2_Nx16 sp, 64
FILTER_VER_LUMA_AVX2_Nx16 ss, 16
FILTER_VER_LUMA_AVX2_Nx16 ss, 32
FILTER_VER_LUMA_AVX2_Nx16 ss, 64

%macro FILTER_VER_LUMA_AVX2_NxN 3
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%3_%1x%2, 4, 12, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4

%ifidn %3,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif

    lea             r6, [r3 * 3]
    lea             r11, [r1 * 4]
    mov             r9d, %2 / 16
.loopH:
    mov             r10d, %1 / 8
.loopW:
    PROCESS_LUMA_AVX2_W8_16R %3
%ifidn %3,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r10d
    jnz             .loopW
    sub             r7, r11
    lea             r0, [r7 - 2 * %1 + 16]
%ifidn %3,sp
    lea             r2, [r8 + r3 * 4 - %1 + 8]
%else
    lea             r2, [r8 + r3 * 4 - 2 * %1 + 16]
%endif
    dec             r9d
    jnz             .loopH
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_NxN 16, 32, sp
FILTER_VER_LUMA_AVX2_NxN 16, 64, sp
FILTER_VER_LUMA_AVX2_NxN 24, 32, sp
FILTER_VER_LUMA_AVX2_NxN 32, 32, sp
FILTER_VER_LUMA_AVX2_NxN 32, 64, sp
FILTER_VER_LUMA_AVX2_NxN 48, 64, sp
FILTER_VER_LUMA_AVX2_NxN 64, 32, sp
FILTER_VER_LUMA_AVX2_NxN 64, 48, sp
FILTER_VER_LUMA_AVX2_NxN 64, 64, sp
FILTER_VER_LUMA_AVX2_NxN 16, 32, ss
FILTER_VER_LUMA_AVX2_NxN 16, 64, ss
FILTER_VER_LUMA_AVX2_NxN 24, 32, ss
FILTER_VER_LUMA_AVX2_NxN 32, 32, ss
FILTER_VER_LUMA_AVX2_NxN 32, 64, ss
FILTER_VER_LUMA_AVX2_NxN 48, 64, ss
FILTER_VER_LUMA_AVX2_NxN 64, 32, ss
FILTER_VER_LUMA_AVX2_NxN 64, 48, ss
FILTER_VER_LUMA_AVX2_NxN 64, 64, ss

%macro FILTER_VER_LUMA_S_AVX2_12x16 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_12x16, 4, 9, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    PROCESS_LUMA_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    mova            m7, m14
    PROCESS_LUMA_AVX2_W4_16R %1
    RET
%endif
%endmacro

FILTER_VER_LUMA_S_AVX2_12x16 sp
FILTER_VER_LUMA_S_AVX2_12x16 ss

%macro FILTER_VER_LUMA_S_AVX2_16x12 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_16x12, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, 2
.loopW:
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 2 * mmsize]
    paddd           m1, m7
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m3, m7
    pmaddwd         m5, [r5]
    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 3 * mmsize]
    paddd           m0, m8
    pmaddwd         m8, m6, [r5 + 2 * mmsize]
    paddd           m2, m8
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    paddd           m4, m8
    pmaddwd         m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhwd       xm9, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddwd         m9, m7, [r5 + 3 * mmsize]
    paddd           m1, m9
    pmaddwd         m9, m7, [r5 + 2 * mmsize]
    paddd           m3, m9
    pmaddwd         m9, m7, [r5 + 1 * mmsize]
    paddd           m5, m9
    pmaddwd         m7, [r5]
    movu            xm9, [r7 + r1]                  ; m9 = row 9
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m2, m10
    pmaddwd         m10, m8, [r5 + 2 * mmsize]
    paddd           m4, m10
    pmaddwd         m10, m8, [r5 + 1 * mmsize]
    paddd           m6, m10
    pmaddwd         m8, [r5]
    movu            xm10, [r7 + r1 * 2]             ; m10 = row 10
    punpckhwd       xm11, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm11, 1
    pmaddwd         m11, m9, [r5 + 3 * mmsize]
    paddd           m3, m11
    pmaddwd         m11, m9, [r5 + 2 * mmsize]
    paddd           m5, m11
    pmaddwd         m11, m9, [r5 + 1 * mmsize]
    paddd           m7, m11
    pmaddwd         m9, [r5]
    movu            xm11, [r7 + r4]                 ; m11 = row 11
    punpckhwd       xm12, xm10, xm11
    punpcklwd       xm10, xm11
    vinserti128     m10, m10, xm12, 1
    pmaddwd         m12, m10, [r5 + 3 * mmsize]
    paddd           m4, m12
    pmaddwd         m12, m10, [r5 + 2 * mmsize]
    paddd           m6, m12
    pmaddwd         m12, m10, [r5 + 1 * mmsize]
    paddd           m8, m12
    pmaddwd         m10, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm12, [r7]                      ; m12 = row 12
    punpckhwd       xm13, xm11, xm12
    punpcklwd       xm11, xm12
    vinserti128     m11, m11, xm13, 1
    pmaddwd         m13, m11, [r5 + 3 * mmsize]
    paddd           m5, m13
    pmaddwd         m13, m11, [r5 + 2 * mmsize]
    paddd           m7, m13
    pmaddwd         m13, m11, [r5 + 1 * mmsize]
    paddd           m9, m13
    pmaddwd         m11, [r5]

%ifidn %1,sp
    paddd           m0, m14
    paddd           m1, m14
    paddd           m2, m14
    paddd           m3, m14
    paddd           m4, m14
    paddd           m5, m14
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
    packssdw        m4, m5

%ifidn %1,sp
    packuswb        m0, m2
    mova            m5, [interp8_hps_shuf]
    vpermd          m0, m5, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif

    movu            xm13, [r7 + r1]                 ; m13 = row 13
    punpckhwd       xm0, xm12, xm13
    punpcklwd       xm12, xm13
    vinserti128     m12, m12, xm0, 1
    pmaddwd         m0, m12, [r5 + 3 * mmsize]
    paddd           m6, m0
    pmaddwd         m0, m12, [r5 + 2 * mmsize]
    paddd           m8, m0
    pmaddwd         m12, [r5 + 1 * mmsize]
    paddd           m10, m12
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 14
    punpckhwd       xm1, xm13, xm0
    punpcklwd       xm13, xm0
    vinserti128     m13, m13, xm1, 1
    pmaddwd         m1, m13, [r5 + 3 * mmsize]
    paddd           m7, m1
    pmaddwd         m1, m13, [r5 + 2 * mmsize]
    paddd           m9, m1
    pmaddwd         m13, [r5 + 1 * mmsize]
    paddd           m11, m13

%ifidn %1,sp
    paddd           m6, m14
    paddd           m7, m14
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m6, m7
    lea             r8, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m5, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm1, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm1
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%endif

    movu            xm1, [r7 + r4]                  ; m1 = row 15
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m2, m0, [r5 + 3 * mmsize]
    paddd           m8, m2
    pmaddwd         m0, [r5 + 2 * mmsize]
    paddd           m10, m0
    lea             r7, [r7 + r1 * 4]
    movu            xm2, [r7]                       ; m2 = row 16
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m3, m1, [r5 + 3 * mmsize]
    paddd           m9, m3
    pmaddwd         m1, [r5 + 2 * mmsize]
    paddd           m11, m1
    movu            xm3, [r7 + r1]                  ; m3 = row 17
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m2, [r5 + 3 * mmsize]
    paddd           m10, m2
    movu            xm4, [r7 + r1 * 2]              ; m4 = row 18
    punpckhwd       xm2, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm2, 1
    pmaddwd         m3, [r5 + 3 * mmsize]
    paddd           m11, m3

%ifidn %1,sp
    paddd           m8, m14
    paddd           m9, m14
    paddd           m10, m14
    paddd           m11, m14
    psrad           m8, 12
    psrad           m9, 12
    psrad           m10, 12
    psrad           m11, 12
%else
    psrad           m8, 6
    psrad           m9, 6
    psrad           m10, 6
    psrad           m11, 6
%endif
    packssdw        m8, m9
    packssdw        m10, m11
    lea             r8, [r8 + r3 * 4]

%ifidn %1,sp
    packuswb        m8, m10
    vpermd          m8, m5, m8
    vextracti128    xm10, m8, 1
    movq            [r8], xm8
    movhps          [r8 + r3], xm8
    movq            [r8 + r3 * 2], xm10
    movhps          [r8 + r6], xm10
    add             r2, 8
%else
    vpermq          m8, m8, 11011000b
    vpermq          m10, m10, 11011000b
    vextracti128    xm9, m8, 1
    vextracti128    xm11, m10, 1
    movu            [r8], xm8
    movu            [r8 + r3], xm9
    movu            [r8 + r3 * 2], xm10
    movu            [r8 + r6], xm11
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_S_AVX2_16x12 sp
FILTER_VER_LUMA_S_AVX2_16x12 ss

%macro FILTER_VER_LUMA_S_AVX2_16x4 1
INIT_YMM avx2
cglobal interp_8tap_vert_%1_16x4, 4, 7, 8, 0 - gprsize
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m7, [pd_526336]
%else
    add             r3d, r3d
%endif
    mov             dword [rsp], 2
.loopW:
    PROCESS_LUMA_S_AVX2_W8_4R %1
    lea             r6, [r3 * 3]
%ifidn %1,sp
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
    add             r2, 8
%else
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
    add             r2, 16
%endif
    lea             r6, [8 * r1 - 16]
    sub             r0, r6
    dec             dword [rsp]
    jnz             .loopW
    RET
%endmacro

FILTER_VER_LUMA_S_AVX2_16x4 sp
FILTER_VER_LUMA_S_AVX2_16x4 ss

%macro PROCESS_LUMA_S_AVX2_W8_8R 1
    movu            xm0, [r0]                       ; m0 = row 0
    movu            xm1, [r0 + r1]                  ; m1 = row 1
    punpckhwd       xm2, xm0, xm1
    punpcklwd       xm0, xm1
    vinserti128     m0, m0, xm2, 1
    pmaddwd         m0, [r5]
    movu            xm2, [r0 + r1 * 2]              ; m2 = row 2
    punpckhwd       xm3, xm1, xm2
    punpcklwd       xm1, xm2
    vinserti128     m1, m1, xm3, 1
    pmaddwd         m1, [r5]
    movu            xm3, [r0 + r4]                  ; m3 = row 3
    punpckhwd       xm4, xm2, xm3
    punpcklwd       xm2, xm3
    vinserti128     m2, m2, xm4, 1
    pmaddwd         m4, m2, [r5 + 1 * mmsize]
    paddd           m0, m4
    pmaddwd         m2, [r5]
    lea             r7, [r0 + r1 * 4]
    movu            xm4, [r7]                       ; m4 = row 4
    punpckhwd       xm5, xm3, xm4
    punpcklwd       xm3, xm4
    vinserti128     m3, m3, xm5, 1
    pmaddwd         m5, m3, [r5 + 1 * mmsize]
    paddd           m1, m5
    pmaddwd         m3, [r5]
    movu            xm5, [r7 + r1]                  ; m5 = row 5
    punpckhwd       xm6, xm4, xm5
    punpcklwd       xm4, xm5
    vinserti128     m4, m4, xm6, 1
    pmaddwd         m6, m4, [r5 + 2 * mmsize]
    paddd           m0, m6
    pmaddwd         m6, m4, [r5 + 1 * mmsize]
    paddd           m2, m6
    pmaddwd         m4, [r5]
    movu            xm6, [r7 + r1 * 2]              ; m6 = row 6
    punpckhwd       xm7, xm5, xm6
    punpcklwd       xm5, xm6
    vinserti128     m5, m5, xm7, 1
    pmaddwd         m7, m5, [r5 + 2 * mmsize]
    paddd           m1, m7
    pmaddwd         m7, m5, [r5 + 1 * mmsize]
    paddd           m3, m7
    pmaddwd         m5, [r5]
    movu            xm7, [r7 + r4]                  ; m7 = row 7
    punpckhwd       xm8, xm6, xm7
    punpcklwd       xm6, xm7
    vinserti128     m6, m6, xm8, 1
    pmaddwd         m8, m6, [r5 + 3 * mmsize]
    paddd           m0, m8
    pmaddwd         m8, m6, [r5 + 2 * mmsize]
    paddd           m2, m8
    pmaddwd         m8, m6, [r5 + 1 * mmsize]
    paddd           m4, m8
    pmaddwd         m6, [r5]
    lea             r7, [r7 + r1 * 4]
    movu            xm8, [r7]                       ; m8 = row 8
    punpckhwd       xm9, xm7, xm8
    punpcklwd       xm7, xm8
    vinserti128     m7, m7, xm9, 1
    pmaddwd         m9, m7, [r5 + 3 * mmsize]
    paddd           m1, m9
    pmaddwd         m9, m7, [r5 + 2 * mmsize]
    paddd           m3, m9
    pmaddwd         m9, m7, [r5 + 1 * mmsize]
    paddd           m5, m9
    pmaddwd         m7, [r5]
    movu            xm9, [r7 + r1]                  ; m9 = row 9
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m2, m10
    pmaddwd         m10, m8, [r5 + 2 * mmsize]
    paddd           m4, m10
    pmaddwd         m8, [r5 + 1 * mmsize]
    paddd           m6, m8
    movu            xm10, [r7 + r1 * 2]             ; m10 = row 10
    punpckhwd       xm8, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm8, 1
    pmaddwd         m8, m9, [r5 + 3 * mmsize]
    paddd           m3, m8
    pmaddwd         m8, m9, [r5 + 2 * mmsize]
    paddd           m5, m8
    pmaddwd         m9, [r5 + 1 * mmsize]
    paddd           m7, m9
    movu            xm8, [r7 + r4]                  ; m8 = row 11
    punpckhwd       xm9, xm10, xm8
    punpcklwd       xm10, xm8
    vinserti128     m10, m10, xm9, 1
    pmaddwd         m9, m10, [r5 + 3 * mmsize]
    paddd           m4, m9
    pmaddwd         m10, [r5 + 2 * mmsize]
    paddd           m6, m10
    lea             r7, [r7 + r1 * 4]
    movu            xm9, [r7]                       ; m9 = row 12
    punpckhwd       xm10, xm8, xm9
    punpcklwd       xm8, xm9
    vinserti128     m8, m8, xm10, 1
    pmaddwd         m10, m8, [r5 + 3 * mmsize]
    paddd           m5, m10
    pmaddwd         m8, [r5 + 2 * mmsize]
    paddd           m7, m8

%ifidn %1,sp
    paddd           m0, m11
    paddd           m1, m11
    paddd           m2, m11
    paddd           m3, m11
    paddd           m4, m11
    paddd           m5, m11
    psrad           m0, 12
    psrad           m1, 12
    psrad           m2, 12
    psrad           m3, 12
    psrad           m4, 12
    psrad           m5, 12
%else
    psrad           m0, 6
    psrad           m1, 6
    psrad           m2, 6
    psrad           m3, 6
    psrad           m4, 6
    psrad           m5, 6
%endif
    packssdw        m0, m1
    packssdw        m2, m3
    packssdw        m4, m5

%ifidn %1,sp
    packuswb        m0, m2
    mova            m5, [interp8_hps_shuf]
    vpermd          m0, m5, m0
    vextracti128    xm2, m0, 1
    movq            [r2], xm0
    movhps          [r2 + r3], xm0
    movq            [r2 + r3 * 2], xm2
    movhps          [r2 + r6], xm2
%else
    vpermq          m0, m0, 11011000b
    vpermq          m2, m2, 11011000b
    vextracti128    xm1, m0, 1
    vextracti128    xm3, m2, 1
    movu            [r2], xm0
    movu            [r2 + r3], xm1
    movu            [r2 + r3 * 2], xm2
    movu            [r2 + r6], xm3
%endif

    movu            xm10, [r7 + r1]                 ; m10 = row 13
    punpckhwd       xm0, xm9, xm10
    punpcklwd       xm9, xm10
    vinserti128     m9, m9, xm0, 1
    pmaddwd         m9, [r5 + 3 * mmsize]
    paddd           m6, m9
    movu            xm0, [r7 + r1 * 2]              ; m0 = row 14
    punpckhwd       xm1, xm10, xm0
    punpcklwd       xm10, xm0
    vinserti128     m10, m10, xm1, 1
    pmaddwd         m10, [r5 + 3 * mmsize]
    paddd           m7, m10

%ifidn %1,sp
    paddd           m6, m11
    paddd           m7, m11
    psrad           m6, 12
    psrad           m7, 12
%else
    psrad           m6, 6
    psrad           m7, 6
%endif
    packssdw        m6, m7
    lea             r8, [r2 + r3 * 4]

%ifidn %1,sp
    packuswb        m4, m6
    vpermd          m4, m5, m4
    vextracti128    xm6, m4, 1
    movq            [r8], xm4
    movhps          [r8 + r3], xm4
    movq            [r8 + r3 * 2], xm6
    movhps          [r8 + r6], xm6
%else
    vpermq          m4, m4, 11011000b
    vpermq          m6, m6, 11011000b
    vextracti128    xm5, m4, 1
    vextracti128    xm7, m6, 1
    movu            [r8], xm4
    movu            [r8 + r3], xm5
    movu            [r8 + r3 * 2], xm6
    movu            [r8 + r6], xm7
%endif
%endmacro

%macro FILTER_VER_LUMA_AVX2_Nx8 2
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_%2x8, 4, 10, 12
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m11, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, %2 / 8
.loopW:
    PROCESS_LUMA_S_AVX2_W8_8R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    RET
%endif
%endmacro

FILTER_VER_LUMA_AVX2_Nx8 sp, 32
FILTER_VER_LUMA_AVX2_Nx8 sp, 16
FILTER_VER_LUMA_AVX2_Nx8 ss, 32
FILTER_VER_LUMA_AVX2_Nx8 ss, 16

%macro FILTER_VER_LUMA_S_AVX2_32x24 1
INIT_YMM avx2
%if ARCH_X86_64 == 1
cglobal interp_8tap_vert_%1_32x24, 4, 10, 15
    mov             r4d, r4m
    shl             r4d, 7
    add             r1d, r1d

%ifdef PIC
    lea             r5, [pw_LumaCoeffVer]
    add             r5, r4
%else
    lea             r5, [pw_LumaCoeffVer + r4]
%endif

    lea             r4, [r1 * 3]
    sub             r0, r4
%ifidn %1,sp
    mova            m14, [pd_526336]
%else
    add             r3d, r3d
%endif
    lea             r6, [r3 * 3]
    mov             r9d, 4
.loopW:
    PROCESS_LUMA_AVX2_W8_16R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loopW
    lea             r9, [r1 * 4]
    sub             r7, r9
    lea             r0, [r7 - 48]
%ifidn %1,sp
    lea             r2, [r8 + r3 * 4 - 24]
%else
    lea             r2, [r8 + r3 * 4 - 48]
%endif
    mova            m11, m14
    mov             r9d, 4
.loop:
    PROCESS_LUMA_S_AVX2_W8_8R %1
%ifidn %1,sp
    add             r2, 8
%else
    add             r2, 16
%endif
    add             r0, 16
    dec             r9d
    jnz             .loop
    RET
%endif
%endmacro

FILTER_VER_LUMA_S_AVX2_32x24 sp
FILTER_VER_LUMA_S_AVX2_32x24 ss

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_32x32(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------;
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_32x32, 4,7,6
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m5,           [pw_2000]
    mova               m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1
    mov                r6d,         32
    dec                r0
    test                r5d,      r5d
    je                 .loop
    sub                r0 ,         r1
    add                r6d ,        3

.loop
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 8]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           m5
    vpermq            m3,           m3,          11011000b
    movu             [r2],         m3

    vbroadcasti128    m3,           [r0 + 16]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 24]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           m5
    vpermq            m3,              m3,          11011000b
    movu             [r2 + 32],         m3

    add                r2,           r3
    add                r0,           r1
    dec               r6d
    jnz                .loop
   RET

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_16x16(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------;
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_16x16, 4,7,6
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m5,           [pw_2000]
    mova               m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1
    mov                r6d,         16
    dec                r0
    test                r5d,        r5d
    je                 .loop
    sub                r0 ,         r1
    add                r6d ,        3

.loop
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 8]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           m5
    vpermq            m3,           m3,          11011000b
    movu              [r2],         m3

    add                r2,          r3
    add                r0,          r1
    dec                r6d
    jnz                .loop
   RET

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_16xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PS_16xN_AVX2 2
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_%1x%2, 4,7,6
    mov                    r4d,        r4m
    mov                    r5d,        r5m
    add                    r3d,        r3d

%ifdef PIC
    lea                    r6,         [tab_ChromaCoeff]
    vpbroadcastd           m0,         [r6 + r4 * 4]
%else
    vpbroadcastd           m0,         [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128         m2,         [pw_1]
    vbroadcasti128         m5,         [pw_2000]
    mova                   m1,         [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1
    mov                    r6d,        %2
    dec                    r0
    test                   r5d,        r5d
    je                     .loop
    sub                    r0 ,        r1
    add                    r6d ,       3

.loop
    ; Row 0
    vbroadcasti128         m3,         [r0]
    pshufb                 m3,         m1
    pmaddubsw              m3,         m0
    pmaddwd                m3,         m2
    vbroadcasti128         m4,         [r0 + 8]
    pshufb                 m4,         m1
    pmaddubsw              m4,         m0
    pmaddwd                m4,         m2

    packssdw               m3,         m4
    psubw                  m3,         m5

    vpermq                 m3,         m3,          11011000b
    movu                   [r2],       m3

    add                    r2,         r3
    add                    r0,         r1
    dec                    r6d
    jnz                    .loop
    RET
%endmacro

    IPFILTER_CHROMA_PS_16xN_AVX2  16 , 32
    IPFILTER_CHROMA_PS_16xN_AVX2  16 , 12
    IPFILTER_CHROMA_PS_16xN_AVX2  16 , 8
    IPFILTER_CHROMA_PS_16xN_AVX2  16 , 4

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_32xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PS_32xN_AVX2 2
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_%1x%2, 4,7,6
    mov                r4d,          r4m
    mov                r5d,          r5m
    add                r3d,          r3d

%ifdef PIC
    lea                r6,           [tab_ChromaCoeff]
    vpbroadcastd       m0,           [r6 + r4 * 4]
%else
    vpbroadcastd       m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m5,           [pw_2000]
    mova               m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1
    mov                r6d,          %2
    dec                r0
    test               r5d,          r5d
    je                 .loop
    sub                r0 ,          r1
    add                r6d ,         3

.loop
    ; Row 0
    vbroadcasti128     m3,           [r0]
    pshufb             m3,           m1
    pmaddubsw          m3,           m0
    pmaddwd            m3,           m2
    vbroadcasti128     m4,           [r0 + 8]
    pshufb             m4,           m1
    pmaddubsw          m4,           m0
    pmaddwd            m4,           m2

    packssdw           m3,           m4
    psubw              m3,           m5

    vpermq             m3,           m3,          11011000b
    movu              [r2],          m3

    vbroadcasti128     m3,           [r0 + 16]
    pshufb             m3,           m1
    pmaddubsw          m3,           m0
    pmaddwd            m3,           m2
    vbroadcasti128     m4,           [r0 + 24]
    pshufb             m4,           m1
    pmaddubsw          m4,           m0
    pmaddwd            m4,           m2

    packssdw           m3,           m4
    psubw              m3,           m5

    vpermq             m3,           m3,          11011000b
    movu               [r2 + 32],    m3

    add                r2,           r3
    add                r0,           r1
    dec                r6d
    jnz                .loop
    RET
%endmacro

IPFILTER_CHROMA_PS_32xN_AVX2  32 , 16
IPFILTER_CHROMA_PS_32xN_AVX2  32 , 24
IPFILTER_CHROMA_PS_32xN_AVX2  32 , 8
;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_4x4(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_4x4, 4,7,5
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec                r0
    test                r5d,       r5d
    je                 .label
    sub                r0 , r1

.label
    ; Row 0-1
    movu              xm3,           [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 2-3
    lea               r0,           [r0 + r1 * 2]
    movu              xm4,           [r0]
    vinserti128       m4,           m4,      [r0 + r1],     1
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           [pw_2000]
    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4
    lea               r2,           [r2 + r3 * 2]
    movhps            [r2],         xm3
    movhps            [r2 + r3],    xm4

    test                r5d,        r5d
    jz                .end
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]

    ;Row 5-6
    movu              xm3,          [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 7
    lea               r0,           [r0 + r1 * 2]
    vbroadcasti128    m4,           [r0]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           [pw_2000]

    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4
    lea               r2,           [r2 + r3 * 2]
    movhps            [r2],         xm3
.end
   RET

cglobal interp_4tap_horiz_ps_4x2, 4,7,5
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec                r0
    test                r5d,       r5d
    je                 .label
    sub                r0 , r1

.label
    ; Row 0-1
    movu              xm3,           [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    packssdw          m3,           m3
    psubw             m3,           [pw_2000]
    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4

    test              r5d,          r5d
    jz                .end
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]

    ;Row 2-3
    movu              xm3,          [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 5
    lea               r0,           [r0 + r1 * 2]
    vbroadcasti128    m4,           [r0]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           [pw_2000]

    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4
    lea               r2,           [r2 + r3 * 2]
    movhps            [r2],         xm3
.end
   RET

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_4xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------;
%macro IPFILTER_CHROMA_PS_4xN_AVX2 2
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_%1x%2, 4,7,5
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    mov              r4,                %2
    dec              r0
    test             r5d,       r5d
    je               .loop
    sub              r0 ,               r1


.loop
    sub              r4d,           4
    ; Row 0-1
    movu              xm3,          [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 2-3
    lea               r0,           [r0 + r1 * 2]
    movu              xm4,          [r0]
    vinserti128       m4,           m4,      [r0 + r1],     1
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           [pw_2000]
    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4
    lea               r2,           [r2 + r3 * 2]
    movhps            [r2],         xm3
    movhps            [r2 + r3],    xm4

    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]

    test              r4d,          r4d
    jnz               .loop
    test                r5d,        r5d
    jz                .end

    ;Row 5-6
    movu              xm3,          [r0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 7
    lea               r0,           [r0 + r1 * 2]
    vbroadcasti128    m4,           [r0]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           [pw_2000]

    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movq              [r2+r3],      xm4
    lea               r2,           [r2 + r3 * 2]
    movhps            [r2],         xm3
.end
RET
%endmacro

    IPFILTER_CHROMA_PS_4xN_AVX2  4 , 8
    IPFILTER_CHROMA_PS_4xN_AVX2  4 , 16
;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_8x8(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------;
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_8x8, 4,7,6
    mov             r4d, r4m
    mov             r5d, r5m
    add             r3d, r3d

%ifdef PIC
    lea               r6,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r6 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,           [pw_1]
    vbroadcasti128     m5,           [pw_2000]
    mova               m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    mov                r6d,      4
    dec                r0
    test                r5d,     r5d
    je                 .loop
    sub                r0 ,      r1
    add                r6d ,     1

.loop
     dec               r6d
    ; Row 0
    vbroadcasti128    m3,           [r0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    psubw             m3,           m5

    vpermq            m3,           m3,          11011000b
    vextracti128      xm4,          m3,     1
    movu             [r2],         xm3
    movu             [r2 + r3],    xm4

    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]
    test               r6d,          r6d
    jnz               .loop
    test              r5d,         r5d
    je                .end

    ;Row 11
    vbroadcasti128    m3,           [r0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    packssdw          m3,           m3
    psubw             m3,           m5
    vpermq            m3,           m3,          11011000b
    movu             [r2],         xm3
.end
   RET

INIT_YMM avx2 
cglobal interp_4tap_horiz_pp_4x2, 4,6,4
    mov             r4d, r4m
%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128    m1,           [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table

    ; Row 0-1
    movu              xm2,          [r0 - 1]
    vinserti128       m2,           m2,      [r0 + r1 - 1],     1
    pshufb            m2,           m1
    pmaddubsw         m2,           m0
    pmaddwd           m2,           [pw_1]

    packssdw          m2,           m2
    pmulhrsw          m2,           [pw_512]
    vextracti128      xm3,          m2,     1
    packuswb          xm2,          xm3

    movd              [r2],         xm2
    pextrd            [r2+r3],      xm2,     2
    RET

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_pp_32xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PP_32xN_AVX2 2
INIT_YMM avx2
cglobal interp_4tap_horiz_pp_%1x%2, 4,6,7
    mov             r4d, r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]
    mova              m6,           [pw_512]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          %2

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    vbroadcasti128    m4,           [r0 + 16]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + 20]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    movu              [r2],         m3
    add               r2,           r3
    add               r0,           r1
    dec               r4d
    jnz               .loop
    RET
%endmacro

IPFILTER_CHROMA_PP_32xN_AVX2 32, 16
IPFILTER_CHROMA_PP_32xN_AVX2 32, 24
IPFILTER_CHROMA_PP_32xN_AVX2 32, 8

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_pp_8xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PP_8xN_AVX2 2
INIT_YMM avx2
cglobal interp_4tap_horiz_pp_%1x%2, 4,6,6
    mov               r4d,    r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    movu              m1,           [tab_Tm]
    vpbroadcastd      m2,           [pw_1]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    sub               r0,           1
    mov               r4d,          %2

.loop:
    sub               r4d,          4
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           [pw_512]
    lea               r0,           [r0 + r1 * 2]

    ; Row 2
    vbroadcasti128    m4,           [r0 ]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    ; Row 3
    vbroadcasti128    m5,           [r0 + r1]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           [pw_512]

    packuswb          m3,           m4
    mova              m5,           [interp_4tap_8x8_horiz_shuf]
    vpermd            m3,           m5,     m3
    vextracti128      xm4,          m3,     1
    movq              [r2],         xm3
    movhps            [r2 + r3],    xm3
    lea               r2,           [r2 + r3 * 2]
    movq              [r2],         xm4
    movhps            [r2 + r3],    xm4
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1*2]
    test              r4d,          r4d
    jnz               .loop
    RET
%endmacro

IPFILTER_CHROMA_PP_8xN_AVX2   8 , 16
IPFILTER_CHROMA_PP_8xN_AVX2   8 , 32
IPFILTER_CHROMA_PP_8xN_AVX2   8 , 4

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_pp_4xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PP_4xN_AVX2 2
INIT_YMM avx2 
cglobal interp_4tap_horiz_pp_%1x%2, 4,6,6
    mov             r4d, r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    vpbroadcastd      m2,           [pw_1]
    vbroadcasti128    m1,           [tab_Tm]
    mov               r4d,          %2
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec                r0

.loop
    sub               r4d,          4
    ; Row 0-1
    movu              xm3,          [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    vinserti128       m3,           m3,      [r0 + r1],     1
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2

    ; Row 2-3
    lea               r0,           [r0 + r1 * 2]
    movu              xm4,          [r0]                      ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    vinserti128       m4,           m4,      [r0 + r1],     1
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2

    packssdw          m3,           m4
    pmulhrsw          m3,           [pw_512]
    vextracti128      xm4,          m3,                     1
    packuswb          xm3,          xm4

    movd              [r2],         xm3
    pextrd            [r2+r3],      xm3,                    2
    lea               r2,           [r2 + r3 * 2]
    pextrd            [r2],         xm3,                    1
    pextrd            [r2+r3],      xm3,                    3

    lea               r0,           [r0 + r1 * 2]
    lea               r2,           [r2 + r3 * 2]
    test              r4d,          r4d
    jnz               .loop
    RET
%endmacro

IPFILTER_CHROMA_PP_4xN_AVX2  4 , 8
IPFILTER_CHROMA_PP_4xN_AVX2  4 , 16

%macro IPFILTER_LUMA_PS_32xN_AVX2 2
INIT_YMM avx2
cglobal interp_8tap_horiz_ps_%1x%2, 4, 7, 8
    mov                         r5d,               r5m
    mov                         r4d,               r4m
%ifdef PIC
    lea                         r6,                [tab_LumaCoeff]
    vpbroadcastq                m0,                [r6 + r4 * 8]
%else
    vpbroadcastq                m0,                [tab_LumaCoeff + r4 * 8]
%endif
    mova                        m6,                [tab_Lm + 32]
    mova                        m1,                [tab_Lm]
    mov                         r4d,                %2                           ;height
    add                         r3d,               r3d
    vbroadcasti128              m2,                [pw_1]
    mova                        m7,                [interp8_hps_shuf]

    ; register map
    ; m0      - interpolate coeff
    ; m1 , m6 - shuffle order table
    ; m2      - pw_1


    sub                         r0,                3
    test                        r5d,               r5d
    jz                          .label
    lea                         r6,                [r1 * 3]                     ; r8 = (N / 2 - 1) * srcStride
    sub                         r0,                r6
    add                         r4d,                7

.label
    lea                         r6,                 [pw_2000]
.loop
    ; Row 0
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m3,             m6           ; row 0 (col 4 to 7)
    pshufb                      m3,                m1                           ; shuffled based on the col order tab_Lm row 0 (col 0 to 3)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4


    vbroadcasti128              m4,                [r0 + 8]
    pshufb                      m5,                m4,            m6            ;row 0 (col 12 to 15)
    pshufb                      m4,                m1                           ;row 0 (col 8 to 11)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5

    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]

    movu                        [r2],              m3                          ;row 0

    vbroadcasti128              m3,                [r0 + 16]
    pshufb                      m4,                m3,             m6           ; row 0 (col 20 to 23)
    pshufb                      m3,                m1                           ; row 0 (col 16 to 19)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 24]
    pshufb                      m5,                m4,            m6            ;row 0 (col 28 to 31)
    pshufb                      m4,                m1                           ;row 0 (col 24 to 27)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5

    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]

    movu                        [r2 + 32],         m3                          ;row 0

    add                         r0,                r1
    add                         r2,                r3
    dec                         r4d
    jnz                         .loop
    RET
%endmacro

IPFILTER_LUMA_PS_32xN_AVX2 32 , 32
IPFILTER_LUMA_PS_32xN_AVX2 32 , 16
IPFILTER_LUMA_PS_32xN_AVX2 32 , 24
IPFILTER_LUMA_PS_32xN_AVX2 32 , 8
IPFILTER_LUMA_PS_32xN_AVX2 32 , 64

INIT_YMM avx2
cglobal interp_8tap_horiz_ps_48x64, 4, 7, 8
    mov                         r5d,               r5m
    mov                         r4d,               r4m
%ifdef PIC
    lea                         r6,                [tab_LumaCoeff]
    vpbroadcastq                m0,                [r6 + r4 * 8]
%else
    vpbroadcastq                m0,                [tab_LumaCoeff + r4 * 8]
%endif
    mova                        m6,                [tab_Lm + 32]
    mova                        m1,                [tab_Lm]
    mov                         r4d,               64                           ;height
    add                         r3d,               r3d
    vbroadcasti128              m2,                [pw_2000]
    mova                        m7,                [pw_1]

    ; register map
    ; m0      - interpolate coeff
    ; m1 , m6 - shuffle order table
    ; m2      - pw_2000

    sub                         r0,                3
    test                        r5d,               r5d
    jz                          .label
    lea                         r6,                [r1 * 3]                     ; r6 = (N / 2 - 1) * srcStride
    sub                         r0,                r6                           ; r0(src)-r6
    add                         r4d,                7                            ; blkheight += N - 1  (7 - 1 = 6 ; since the last one row not in loop)

.label
    lea                         r6,                [interp8_hps_shuf]
.loop
    ; Row 0
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m3,             m6           ; row 0 (col 4 to 7)
    pshufb                      m3,                m1                           ; shuffled based on the col order tab_Lm row 0 (col 0 to 3)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 8]
    pshufb                      m5,                m4,             m6            ;row 0 (col 12 to 15)
    pshufb                      m4,                m1                           ;row 0 (col 8 to 11)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m7
    pmaddwd                     m5,                m7
    packssdw                    m4,                m5
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4
    mova                        m5,                [r6]
    vpermd                      m3,                m5,             m3
    psubw                       m3,                m2
    movu                        [r2],              m3                          ;row 0

    vbroadcasti128              m3,                [r0 + 16]
    pshufb                      m4,                m3,             m6           ; row 0 (col 20 to 23)
    pshufb                      m3,                m1                           ; row 0 (col 16 to 19)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 24]
    pshufb                      m5,                m4,             m6            ;row 0 (col 28 to 31)
    pshufb                      m4,                m1                           ;row 0 (col 24 to 27)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m7
    pmaddwd                     m5,                m7
    packssdw                    m4,                m5
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4
    mova                        m5,                [r6]
    vpermd                      m3,                m5,               m3
    psubw                       m3,                m2
    movu                        [r2 + 32],         m3                          ;row 0

    vbroadcasti128              m3,                [r0 + 32]
    pshufb                      m4,                m3,             m6           ; row 0 (col 36 to 39)
    pshufb                      m3,                m1                           ; row 0 (col 32 to 35)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 40]
    pshufb                      m5,                m4,            m6            ;row 0 (col 44 to 47)
    pshufb                      m4,                m1                           ;row 0 (col 40 to 43)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m7
    pmaddwd                     m5,                m7
    packssdw                    m4,                m5
    pmaddwd                     m3,                m7
    pmaddwd                     m4,                m7
    packssdw                    m3,                m4
    mova                        m5,                [r6]
    vpermd                      m3,                m5,               m3
    psubw                       m3,                m2
    movu                        [r2 + 64],         m3                          ;row 0

    add                         r0,                r1
    add                         r2,                r3
    dec                         r4d
    jnz                         .loop
    RET

INIT_YMM avx2
cglobal interp_8tap_horiz_pp_24x32, 4,6,8
    sub               r0,         3
    mov               r4d,        r4m
%ifdef PIC
    lea               r5,         [tab_LumaCoeff]
    vpbroadcastd      m0,         [r5 + r4 * 8]
    vpbroadcastd      m1,         [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,         [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,         [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,         [tab_Tm + 16]
    vpbroadcastd      m7,         [pw_1]
    lea               r5,         [tab_Tm]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,        32
.loop:
    ; Row 0
    vbroadcasti128    m4,         [r0]                        ; [x E D C B A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,         m4,     m3
    pshufb            m4,         [r5]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7

    vbroadcasti128    m5,         [r0 + 8]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [r5]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7
    packssdw          m4,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,         [pw_512]

    vbroadcasti128    m2,         [r0 + 16]
    pshufb            m5,         m2,     m3
    pshufb            m2,         [r5]
    pmaddubsw         m2,         m0
    pmaddubsw         m5,         m1
    paddw             m2,         m5
    pmaddwd           m2,         m7

    packssdw          m2,         m2
    pmulhrsw          m2,         [pw_512]
    packuswb          m4,         m2
    vpermq            m4,         m4,     11011000b
    vextracti128      xm5,        m4,     1
    pshufd            xm4,        xm4,    11011000b
    pshufd            xm5,        xm5,    11011000b

    movu              [r2],       xm4
    movq              [r2 + 16],  xm5
    add               r0,         r1
    add               r2,         r3
    dec               r4d
    jnz               .loop
    RET

INIT_YMM avx2
cglobal interp_8tap_horiz_pp_12x16, 4,6,8
    sub               r0,        3
    mov               r4d,       r4m
%ifdef PIC
    lea               r5,        [tab_LumaCoeff]
    vpbroadcastd      m0,        [r5 + r4 * 8]
    vpbroadcastd      m1,        [r5 + r4 * 8 + 4]
%else
    vpbroadcastd      m0,         [tab_LumaCoeff + r4 * 8]
    vpbroadcastd      m1,         [tab_LumaCoeff + r4 * 8 + 4]
%endif
    movu              m3,         [tab_Tm + 16]
    vpbroadcastd      m7,         [pw_1]
    lea               r5,         [tab_Tm]

    ; register map
    ; m0 , m1 interpolate coeff
    ; m2 , m2  shuffle order table
    ; m7 - pw_1

    mov               r4d,        8
.loop:
    ; Row 0
    vbroadcasti128    m4,         [r0]                        ;first 8 element
    pshufb            m5,         m4,     m3
    pshufb            m4,         [r5]
    pmaddubsw         m4,         m0
    pmaddubsw         m5,         m1
    paddw             m4,         m5
    pmaddwd           m4,         m7

    vbroadcasti128    m5,         [r0 + 8]                    ; element 8 to 11
    pshufb            m6,         m5,     m3
    pshufb            m5,         [r5]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7

    packssdw          m4,         m5                          ; [17 16 15 14 07 06 05 04 13 12 11 10 03 02 01 00]
    pmulhrsw          m4,         [pw_512]

    ;Row 1
    vbroadcasti128    m2,         [r0 + r1]
    pshufb            m5,         m2,     m3
    pshufb            m2,         [r5]
    pmaddubsw         m2,         m0
    pmaddubsw         m5,         m1
    paddw             m2,         m5
    pmaddwd           m2,         m7

    vbroadcasti128    m5,         [r0 + r1 + 8]
    pshufb            m6,         m5,     m3
    pshufb            m5,         [r5]
    pmaddubsw         m5,         m0
    pmaddubsw         m6,         m1
    paddw             m5,         m6
    pmaddwd           m5,         m7

    packssdw          m2,         m5
    pmulhrsw          m2,         [pw_512]
    packuswb          m4,         m2
    vpermq            m4,         m4,     11011000b
    vextracti128      xm5,        m4,     1
    pshufd            xm4,        xm4,    11011000b
    pshufd            xm5,        xm5,    11011000b

    movq              [r2],       xm4
    pextrd            [r2+8],     xm4,    2
    movq              [r2 + r3],  xm5
    pextrd            [r2+r3+8],  xm5,    2
    lea               r0,         [r0 + r1 * 2]
    lea               r2,         [r2 + r3 * 2]
    dec               r4d
    jnz              .loop
    RET

;-------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_pp_16xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PP_16xN_AVX2 2
INIT_YMM avx2
cglobal interp_4tap_horiz_pp_%1x%2, 4, 6, 7
    mov               r4d,          r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m6,           [pw_512]
    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          %2/2

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + r1 + 4]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    vextracti128      xm4,          m3,       1
    movu              [r2],         xm3
    movu              [r2 + r3],    xm4
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]
    dec               r4d
    jnz               .loop
    RET
%endmacro

IPFILTER_CHROMA_PP_16xN_AVX2 16 , 8
IPFILTER_CHROMA_PP_16xN_AVX2 16 , 32
IPFILTER_CHROMA_PP_16xN_AVX2 16 , 12
IPFILTER_CHROMA_PP_16xN_AVX2 16 , 4

%macro IPFILTER_LUMA_PS_64xN_AVX2 1
INIT_YMM avx2
cglobal interp_8tap_horiz_ps_64x%1, 4, 7, 8
    mov                         r5d,               r5m
    mov                         r4d,               r4m
%ifdef PIC
    lea                         r6,                [tab_LumaCoeff]
    vpbroadcastq                m0,                [r6 + r4 * 8]
%else
    vpbroadcastq                m0,                [tab_LumaCoeff + r4 * 8]
%endif
    mova                        m6,                [tab_Lm + 32]
    mova                        m1,                [tab_Lm]
    mov                         r4d,               %1                           ;height
    add                         r3d,               r3d
    vbroadcasti128              m2,                [pw_1]
    mova                        m7,                [interp8_hps_shuf]

    ; register map
    ; m0      - interpolate coeff
    ; m1 , m6 - shuffle order table
    ; m2      - pw_2000

    sub                         r0,                3
    test                        r5d,               r5d
    jz                          .label
    lea                         r6,                [r1 * 3]
    sub                         r0,                r6                           ; r0(src)-r6
    add                         r4d,               7                            ; blkheight += N - 1

.label
    lea                         r6,                [pw_2000]
.loop
    ; Row 0
    vbroadcasti128              m3,                [r0]                         ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb                      m4,                m3,             m6           ; row 0 (col 4 to 7)
    pshufb                      m3,                m1                           ; shuffled based on the col order tab_Lm row 0 (col 0 to 3)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 8]
    pshufb                      m5,                m4,            m6            ;row 0 (col 12 to 15)
    pshufb                      m4,                m1                           ;row 0 (col 8 to 11)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]
    movu                        [r2],              m3                          ;row 0

    vbroadcasti128              m3,                [r0 + 16]
    pshufb                      m4,                m3,             m6           ; row 0 (col 20 to 23)
    pshufb                      m3,                m1                           ; row 0 (col 16 to 19)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 24]
    pshufb                      m5,                m4,            m6            ;row 0 (col 28 to 31)
    pshufb                      m4,                m1                           ;row 0 (col 24 to 27)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]
    movu                        [r2 + 32],         m3                          ;row 0

    vbroadcasti128              m3,                [r0 + 32]
    pshufb                      m4,                m3,             m6           ; row 0 (col 36 to 39)
    pshufb                      m3,                m1                           ; row 0 (col 32 to 35)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 40]
    pshufb                      m5,                m4,            m6            ;row 0 (col 44 to 47)
    pshufb                      m4,                m1                           ;row 0 (col 40 to 43)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]
    movu                        [r2 + 64],         m3                          ;row 0
    vbroadcasti128              m3,                [r0 + 48]
    pshufb                      m4,                m3,             m6           ; row 0 (col 52 to 55)
    pshufb                      m3,                m1                           ; row 0 (col 48 to 51)
    pmaddubsw                   m3,                m0
    pmaddubsw                   m4,                m0
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4

    vbroadcasti128              m4,                [r0 + 56]
    pshufb                      m5,                m4,            m6            ;row 0 (col 60 to 63)
    pshufb                      m4,                m1                           ;row 0 (col 56 to 59)
    pmaddubsw                   m4,                m0
    pmaddubsw                   m5,                m0
    pmaddwd                     m4,                m2
    pmaddwd                     m5,                m2
    packssdw                    m4,                m5
    pmaddwd                     m3,                m2
    pmaddwd                     m4,                m2
    packssdw                    m3,                m4
    vpermd                      m3,                m7,               m3
    psubw                       m3,                [r6]
    movu                        [r2 + 96],         m3                          ;row 0

    add                          r0,                r1
    add                          r2,                r3
    dec                          r4d
    jnz                         .loop
    RET
%endmacro

IPFILTER_LUMA_PS_64xN_AVX2 64
IPFILTER_LUMA_PS_64xN_AVX2 48
IPFILTER_LUMA_PS_64xN_AVX2 32
IPFILTER_LUMA_PS_64xN_AVX2 16

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_8xN(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------
%macro IPFILTER_CHROMA_PS_8xN_AVX2 1
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_8x%1, 4,7,6
    mov                r4d,             r4m
    mov                r5d,             r5m
    add                r3d,             r3d

%ifdef PIC
    lea                r6,              [tab_ChromaCoeff]
    vpbroadcastd       m0,              [r6 + r4 * 4]
%else
    vpbroadcastd       m0,              [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,              [pw_1]
    vbroadcasti128     m5,              [pw_2000]
    mova               m1,              [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    mov                r6d,             %1/2
    dec                r0
    test               r5d,             r5d
    jz                 .loop
    sub                r0 ,             r1
    inc                r6d

.loop
    ; Row 0
    vbroadcasti128     m3,              [r0]
    pshufb             m3,              m1
    pmaddubsw          m3,              m0
    pmaddwd            m3,              m2

    ; Row 1
    vbroadcasti128     m4,              [r0 + r1]
    pshufb             m4,              m1
    pmaddubsw          m4,              m0
    pmaddwd            m4,              m2
    packssdw           m3,              m4
    psubw              m3,              m5
    vpermq             m3,              m3,          11011000b
    vextracti128       xm4,             m3,          1
    movu               [r2],            xm3
    movu               [r2 + r3],       xm4

    lea                r2,              [r2 + r3 * 2]
    lea                r0,              [r0 + r1 * 2]
    dec                r6d
    jnz                .loop
    test               r5d,             r5d
    jz                 .end

    ;Row 11
    vbroadcasti128     m3,              [r0]
    pshufb             m3,              m1
    pmaddubsw          m3,              m0
    pmaddwd            m3,              m2
    packssdw           m3,              m3
    psubw              m3,              m5
    vpermq             m3,              m3,          11011000b
    movu               [r2],            xm3
.end
   RET
%endmacro

    IPFILTER_CHROMA_PS_8xN_AVX2  2
    IPFILTER_CHROMA_PS_8xN_AVX2  32
    IPFILTER_CHROMA_PS_8xN_AVX2  16
    IPFILTER_CHROMA_PS_8xN_AVX2  6
    IPFILTER_CHROMA_PS_8xN_AVX2  4

INIT_YMM avx2
cglobal interp_4tap_horiz_ps_2x4, 4, 7, 3
    mov                r4d,            r4m
    mov                r5d,            r5m
    add                r3d,            r3d
%ifdef PIC
    lea                r6,             [tab_ChromaCoeff]
    vpbroadcastd       m0,             [r6 + r4 * 4]
%else
    vpbroadcastd       m0,             [tab_ChromaCoeff + r4 * 4]
%endif

    mova               xm3,            [pw_2000]
    dec                r0
    test               r5d,            r5d
    jz                 .label
    sub                r0,             r1

.label
    lea                r6,             [r1 * 3]
    movq               xm1,            [r0]
    movhps             xm1,            [r0 + r1]
    movq               xm2,            [r0 + r1 * 2]
    movhps             xm2,            [r0 + r6]

    vinserti128        m1,             m1,          xm2,          1
    pshufb             m1,             [interp4_hps_shuf]
    pmaddubsw          m1,             m0
    pmaddwd            m1,             [pw_1]
    vextracti128       xm2,            m1,          1
    packssdw           xm1,            xm2
    psubw              xm1,            xm3

    lea                r4,             [r3 * 3]
    movd               [r2],           xm1
    pextrd             [r2 + r3],      xm1,         1
    pextrd             [r2 + r3 * 2],  xm1,         2
    pextrd             [r2 + r4],      xm1,         3

    test               r5d,            r5d
    jz                .end
    lea                r2,             [r2 + r3 * 4]
    lea                r0,             [r0 + r1 * 4]

    movq               xm1,            [r0]
    movhps             xm1,            [r0 + r1]
    movq               xm2,            [r0 + r1 * 2]
    vinserti128        m1,             m1,          xm2,          1
    pshufb             m1,             [interp4_hps_shuf]
    pmaddubsw          m1,             m0
    pmaddwd            m1,             [pw_1]
    vextracti128       xm2,            m1,          1
    packssdw           xm1,            xm2
    psubw              xm1,            xm3

    movd               [r2],           xm1
    pextrd             [r2 + r3],      xm1,         1
    pextrd             [r2 + r3 * 2],  xm1,         2
.end
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_ps_2x8, 4, 7, 7
    mov               r4d,           r4m
    mov               r5d,           r5m
    add               r3d,           r3d

%ifdef PIC
    lea               r6,            [tab_ChromaCoeff]
    vpbroadcastd      m0,            [r6 + r4 * 4]
%else
    vpbroadcastd      m0,            [tab_ChromaCoeff + r4 * 4]
%endif
    vbroadcasti128    m6,            [pw_2000]
    test              r5d,            r5d
    jz                .label
    sub               r0,             r1

.label
    mova              m4,            [interp4_hps_shuf]
    mova              m5,            [pw_1]
    dec               r0
    lea               r4,            [r1 * 3]
    movq              xm1,           [r0]                                   ;row 0
    movhps            xm1,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    movhps            xm2,           [r0 + r4]
    vinserti128       m1,            m1,          xm2,          1
    lea               r0,            [r0 + r1 * 4]
    movq              xm3,           [r0]
    movhps            xm3,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    movhps            xm2,           [r0 + r4]
    vinserti128       m3,            m3,          xm2,          1

    pshufb            m1,            m4
    pshufb            m3,            m4
    pmaddubsw         m1,            m0
    pmaddubsw         m3,            m0
    pmaddwd           m1,            m5
    pmaddwd           m3,            m5
    packssdw          m1,            m3
    psubw             m1,            m6

    lea               r4,            [r3 * 3]
    vextracti128      xm2,           m1,          1

    movd              [r2],          xm1
    pextrd            [r2 + r3],     xm1,         1
    movd              [r2 + r3 * 2], xm2
    pextrd            [r2 + r4],     xm2,         1
    lea               r2,            [r2 + r3 * 4]
    pextrd            [r2],          xm1,         2
    pextrd            [r2 + r3],     xm1,         3
    pextrd            [r2 + r3 * 2], xm2,         2
    pextrd            [r2 + r4],     xm2,         3
    test              r5d,            r5d
    jz                .end

    lea               r0,            [r0 + r1 * 4]
    lea               r2,            [r2 + r3 * 4]
    movq              xm1,           [r0]                                   ;row 0
    movhps            xm1,           [r0 + r1]
    movq              xm2,           [r0 + r1 * 2]
    vinserti128       m1,            m1,          xm2,          1
    pshufb            m1,            m4
    pmaddubsw         m1,            m0
    pmaddwd           m1,            m5
    packssdw          m1,            m1
    psubw             m1,            m6
    vextracti128      xm2,           m1,          1

    movd              [r2],          xm1
    pextrd            [r2 + r3],     xm1,         1
    movd              [r2 + r3 * 2], xm2
.end
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_12x16, 4, 6, 7
    mov               r4d,          r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m6,           [pw_512]
    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          8

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]                    ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    ; Row 1
    vbroadcasti128    m4,           [r0 + r1]                   ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + r1 + 4]               ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    vextracti128      xm4,          m3,       1
    movq              [r2],         xm3
    pextrd            [r2+8],       xm3,      2
    movq              [r2 + r3],    xm4
    pextrd            [r2 + r3 + 8],xm4,      2
    lea               r2,           [r2 + r3 * 2]
    lea               r0,           [r0 + r1 * 2]
    dec               r4d
    jnz               .loop
    RET

INIT_YMM avx2
cglobal interp_4tap_horiz_pp_24x32, 4,6,7
    mov              r4d,           r4m

%ifdef PIC
    lea               r5,           [tab_ChromaCoeff]
    vpbroadcastd      m0,           [r5 + r4 * 4]
%else
    vpbroadcastd      m0,           [tab_ChromaCoeff + r4 * 4]
%endif

    mova              m1,           [interp4_horiz_shuf1]
    vpbroadcastd      m2,           [pw_1]
    mova              m6,           [pw_512]
    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    dec               r0
    mov               r4d,          32

.loop:
    ; Row 0
    vbroadcasti128    m3,           [r0]                        ; [x x x x x A 9 8 7 6 5 4 3 2 1 0]
    pshufb            m3,           m1
    pmaddubsw         m3,           m0
    pmaddwd           m3,           m2
    vbroadcasti128    m4,           [r0 + 4]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    packssdw          m3,           m4
    pmulhrsw          m3,           m6

    vbroadcasti128    m4,           [r0 + 16]
    pshufb            m4,           m1
    pmaddubsw         m4,           m0
    pmaddwd           m4,           m2
    vbroadcasti128    m5,           [r0 + 20]
    pshufb            m5,           m1
    pmaddubsw         m5,           m0
    pmaddwd           m5,           m2
    packssdw          m4,           m5
    pmulhrsw          m4,           m6

    packuswb          m3,           m4
    vpermq            m3,           m3,      11011000b

    vextracti128      xm4,          m3,       1
    movu              [r2],         xm3
    movq              [r2 + 16],    xm4
    add               r2,           r3
    add               r0,           r1
    dec               r4d
    jnz               .loop
    RET

;-----------------------------------------------------------------------------------------------------------------------------
; void interp_4tap_horiz_ps_6x8(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;-----------------------------------------------------------------------------------------------------------------------------;
INIT_YMM avx2 
cglobal interp_4tap_horiz_ps_6x8, 4,7,6
    mov                r4d,            r4m
    mov                r5d,            r5m
    add                r3d,            r3d

%ifdef PIC
    lea                r6,             [tab_ChromaCoeff]
    vpbroadcastd       m0,             [r6 + r4 * 4]
%else
    vpbroadcastd       m0,             [tab_ChromaCoeff + r4 * 4]
%endif

    vbroadcasti128     m2,             [pw_1]
    vbroadcasti128     m5,             [pw_2000]
    mova               m1,             [tab_Tm]

    ; register map
    ; m0 - interpolate coeff
    ; m1 - shuffle order table
    ; m2 - constant word 1

    mov               r6d,             8/2
    dec               r0
    test              r5d,             r5d
    jz                .loop
    sub               r0 ,             r1
    inc               r6d

.loop
    ; Row 0
    vbroadcasti128    m3,              [r0]
    pshufb            m3,              m1
    pmaddubsw         m3,              m0
    pmaddwd           m3,              m2

    ; Row 1
    vbroadcasti128    m4,              [r0 + r1]
    pshufb            m4,              m1
    pmaddubsw         m4,              m0
    pmaddwd           m4,              m2
    packssdw          m3,              m4
    psubw             m3,              m5
    vpermq            m3,              m3,          11011000b
    vextracti128      xm4,             m3,          1
    movq              [r2],            xm3
    pextrd            [r2 + 8],        xm3,         2
    movq              [r2 + r3],       xm4
    pextrd            [r2 + r3 + 8],   xm4,         2
    lea               r2,              [r2 + r3 * 2]
    lea               r0,              [r0 + r1 * 2]
    dec               r6d
    jnz              .loop
    test              r5d,             r5d
    jz               .end

    ;Row 11
    vbroadcasti128    m3,              [r0]
    pshufb            m3,              m1
    pmaddubsw         m3,              m0
    pmaddwd           m3,              m2
    packssdw          m3,              m3
    psubw             m3,              m5
    vextracti128      xm4,             m3,          1
    movq              [r2],            xm3
    movd              [r2+8],          xm4
.end
    RET