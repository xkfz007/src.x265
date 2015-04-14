;*****************************************************************************
;* Copyright (C) 2013 x265 project
;*
;* Authors: Nabajit Deka <nabajit@multicorewareinc.com>
;*          Murugan Vairavel <murugan@multicorewareinc.com>
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

tab_c_32:         times 4 dd 32
tab_c_n32768:     times 4 dd -32768
tab_c_524800:     times 4 dd 524800
tab_c_n8192:      times 8 dw -8192

tab_Tm16:         db 0, 1, 2, 3, 4,  5,  6, 7, 2, 3, 4,  5, 6, 7, 8, 9

tab_ChromaCoeff:  dw  0, 64,  0,  0
                  dw -2, 58, 10, -2
                  dw -4, 54, 16, -2
                  dw -6, 46, 28, -4
                  dw -4, 36, 36, -4
                  dw -4, 28, 46, -6
                  dw -2, 16, 54, -4
                  dw -2, 10, 58, -2

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

tab_LumaCoeff:    dw   0, 0,  0,  64,  0,   0,  0,  0
                  dw  -1, 4, -10, 58,  17, -5,  1,  0
                  dw  -1, 4, -11, 40,  40, -11, 4, -1
                  dw   0, 1, -5,  17,  58, -10, 4, -1

tab_LumaCoeffV:   times 4 dw 0, 0
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

SECTION .text

cextern pd_32
cextern pw_pixel_max
cextern pd_n32768

;------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_4x4(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;------------------------------------------------------------------------------------------------------------
%macro FILTER_HOR_LUMA_W4 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4, 7, 8

    mov         r4d, r4m
    sub         r0, 6
    shl         r4d, 4
    add         r1, r1
    add         r3, r3

%ifdef PIC
    lea         r6, [tab_LumaCoeff]
    mova        m0, [r6 + r4]
%else
    mova        m0, [tab_LumaCoeff + r4]
%endif

%ifidn %3, pp 
    mova        m1, [pd_32]
    pxor        m6, m6
    mova        m7, [pw_pixel_max]
%else
    mova        m1, [pd_n32768]
%endif

    mov         r4d, %2
%ifidn %3, ps
    cmp         r5m, byte 0
    je          .loopH
    lea         r6, [r1 + 2 * r1]
    sub         r0, r6
    add         r4d, 7
%endif

.loopH:
    movu        m2, [r0]                     ; m2 = src[0-7]
    movu        m3, [r0 + 16]                ; m3 = src[8-15]

    pmaddwd     m4, m2, m0
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m3, m2, 6                    ; m3 = src[3-10]
    pmaddwd     m3, m0
    phaddd      m5, m3

    phaddd      m4, m5
    paddd       m4, m1
%ifidn %3, pp
    psrad       m4, 6
    packusdw    m4, m4
    CLIPW       m4, m6, m7
%else
    psrad       m4, 2
    packssdw    m4, m4
%endif

    movh        [r2], m4

    add         r0, r1
    add         r2, r3

    dec         r4d
    jnz         .loopH
    RET
%endmacro

;------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_4x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W4 4, 4, pp
FILTER_HOR_LUMA_W4 4, 8, pp
FILTER_HOR_LUMA_W4 4, 16, pp

;---------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_4x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;---------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W4 4, 4, ps
FILTER_HOR_LUMA_W4 4, 8, ps
FILTER_HOR_LUMA_W4 4, 16, ps

;------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;------------------------------------------------------------------------------------------------------------
%macro FILTER_HOR_LUMA_W8 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4, 7, 8

    add         r1, r1
    add         r3, r3
    mov         r4d, r4m
    sub         r0, 6
    shl         r4d, 4

%ifdef PIC
    lea         r6, [tab_LumaCoeff]
    mova        m0, [r6 + r4]
%else
    mova        m0, [tab_LumaCoeff + r4]
%endif

%ifidn %3, pp 
    mova        m1, [pd_32]
    pxor        m7, m7
%else
    mova        m1, [pd_n32768]
%endif

    mov         r4d, %2
%ifidn %3, ps
    cmp         r5m, byte 0
    je          .loopH
    lea         r6, [r1 + 2 * r1]
    sub         r0, r6
    add         r4d, 7
%endif

.loopH:
    movu        m2, [r0]                     ; m2 = src[0-7]
    movu        m3, [r0 + 16]                ; m3 = src[8-15]

    pmaddwd     m4, m2, m0
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m3, m2, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m3, m2, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m3, m2, 14                   ; m3 = src[7-14]
    pmaddwd     m3, m0
    phaddd      m6, m3
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    CLIPW       m4, m7, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif

    movu        [r2], m4

    add         r0, r1
    add         r2, r3

    dec         r4d
    jnz         .loopH
    RET
%endmacro

;------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_8x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W8 8, 4, pp
FILTER_HOR_LUMA_W8 8, 8, pp
FILTER_HOR_LUMA_W8 8, 16, pp
FILTER_HOR_LUMA_W8 8, 32, pp

;---------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_8x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;---------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W8 8, 4, ps
FILTER_HOR_LUMA_W8 8, 8, ps
FILTER_HOR_LUMA_W8 8, 16, ps
FILTER_HOR_LUMA_W8 8, 32, ps

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_HOR_LUMA_W12 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4, 7, 8

    add         r1, r1
    add         r3, r3
    mov         r4d, r4m
    sub         r0, 6
    shl         r4d, 4

%ifdef PIC
    lea         r6, [tab_LumaCoeff]
    mova        m0, [r6 + r4]
%else
    mova        m0, [tab_LumaCoeff + r4]
%endif
%ifidn %3, pp 
    mova        m1, [pd_32]
%else
    mova        m1, [pd_n32768]
%endif

    mov         r4d, %2
%ifidn %3, ps
    cmp         r5m, byte 0
    je          .loopH
    lea         r6, [r1 + 2 * r1]
    sub         r0, r6
    add         r4d, 7
%endif

.loopH:
    movu        m2, [r0]                     ; m2 = src[0-7]
    movu        m3, [r0 + 16]                ; m3 = src[8-15]

    pmaddwd     m4, m2, m0
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m3, m2, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m3, m2, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m7, m3, m2, 14               ; m2 = src[7-14]
    pmaddwd     m7, m0
    phaddd      m6, m7
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif

    movu        [r2], m4

    movu        m2, [r0 + 32]                ; m2 = src[16-23]

    pmaddwd     m4, m3, m0                   ; m3 = src[8-15]
    palignr     m5, m2, m3, 2                ; m5 = src[9-16]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m2, m3, 4                ; m5 = src[10-17]
    pmaddwd     m5, m0
    palignr     m2, m3, 6                    ; m2 = src[11-18]
    pmaddwd     m2, m0
    phaddd      m5, m2
    phaddd      m4, m5
    paddd       m4, m1
%ifidn %3, pp 
    psrad       m4, 6
    packusdw    m4, m4
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    packssdw    m4, m4
%endif

    movh        [r2 + 16], m4

    add         r0, r1
    add         r2, r3

    dec         r4d
    jnz         .loopH
    RET
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_12x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W12 12, 16, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_12x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W12 12, 16, ps

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_HOR_LUMA_W16 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4, 7, 8

    add         r1, r1
    add         r3, r3
    mov         r4d, r4m
    sub         r0, 6
    shl         r4d, 4

%ifdef PIC
    lea         r6, [tab_LumaCoeff]
    mova        m0, [r6 + r4]
%else
    mova        m0, [tab_LumaCoeff + r4]
%endif

%ifidn %3, pp 
    mova        m1, [pd_32]
%else
    mova        m1, [pd_n32768]
%endif

    mov         r4d, %2
%ifidn %3, ps
    cmp         r5m, byte 0
    je          .loopH
    lea         r6, [r1 + 2 * r1]
    sub         r0, r6
    add         r4d, 7
%endif

.loopH:
%assign x 0
%rep %1 / 16
    movu        m2, [r0 + x]                 ; m2 = src[0-7]
    movu        m3, [r0 + 16 + x]            ; m3 = src[8-15]

    pmaddwd     m4, m2, m0
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m3, m2, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m3, m2, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m7, m3, m2, 14               ; m2 = src[7-14]
    pmaddwd     m7, m0
    phaddd      m6, m7
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif
    movu        [r2 + x], m4

    movu        m2, [r0 + 32 + x]            ; m2 = src[16-23]

    pmaddwd     m4, m3, m0                   ; m3 = src[8-15]
    palignr     m5, m2, m3, 2                ; m5 = src[9-16]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m2, m3, 4                ; m5 = src[10-17]
    pmaddwd     m5, m0
    palignr     m6, m2, m3, 6                ; m6 = src[11-18]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m2, m3, 8                ; m5 = src[12-19]
    pmaddwd     m5, m0
    palignr     m6, m2, m3, 10               ; m6 = src[13-20]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m2, m3, 12               ; m6 = src[14-21]
    pmaddwd     m6, m0
    palignr     m2, m3, 14                   ; m3 = src[15-22]
    pmaddwd     m2, m0
    phaddd      m6, m2
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif
    movu        [r2 + 16 + x], m4

%assign x x+32
%endrep

    add         r0, r1
    add         r2, r3

    dec         r4d
    jnz         .loopH
    RET
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_16x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 16, 4, pp
FILTER_HOR_LUMA_W16 16, 8, pp
FILTER_HOR_LUMA_W16 16, 12, pp
FILTER_HOR_LUMA_W16 16, 16, pp
FILTER_HOR_LUMA_W16 16, 32, pp
FILTER_HOR_LUMA_W16 16, 64, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_16x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 16, 4, ps
FILTER_HOR_LUMA_W16 16, 8, ps
FILTER_HOR_LUMA_W16 16, 12, ps
FILTER_HOR_LUMA_W16 16, 16, ps
FILTER_HOR_LUMA_W16 16, 32, ps
FILTER_HOR_LUMA_W16 16, 64, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_32x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 32, 8, pp
FILTER_HOR_LUMA_W16 32, 16, pp
FILTER_HOR_LUMA_W16 32, 24, pp
FILTER_HOR_LUMA_W16 32, 32, pp
FILTER_HOR_LUMA_W16 32, 64, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_32x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 32, 8, ps
FILTER_HOR_LUMA_W16 32, 16, ps
FILTER_HOR_LUMA_W16 32, 24, ps
FILTER_HOR_LUMA_W16 32, 32, ps
FILTER_HOR_LUMA_W16 32, 64, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_48x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 48, 64, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_48x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 48, 64, ps

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_64x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 64, 16, pp
FILTER_HOR_LUMA_W16 64, 32, pp
FILTER_HOR_LUMA_W16 64, 48, pp
FILTER_HOR_LUMA_W16 64, 64, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_64x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W16 64, 16, ps
FILTER_HOR_LUMA_W16 64, 32, ps
FILTER_HOR_LUMA_W16 64, 48, ps
FILTER_HOR_LUMA_W16 64, 64, ps

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_HOR_LUMA_W24 3
INIT_XMM sse4
cglobal interp_8tap_horiz_%3_%1x%2, 4, 7, 8

    add         r1, r1
    add         r3, r3
    mov         r4d, r4m
    sub         r0, 6
    shl         r4d, 4

%ifdef PIC
    lea         r6, [tab_LumaCoeff]
    mova        m0, [r6 + r4]
%else
    mova        m0, [tab_LumaCoeff + r4]
%endif
%ifidn %3, pp 
    mova        m1, [pd_32]
%else
    mova        m1, [pd_n32768]
%endif

    mov         r4d, %2
%ifidn %3, ps
    cmp         r5m, byte 0
    je          .loopH
    lea         r6, [r1 + 2 * r1]
    sub         r0, r6
    add         r4d, 7
%endif

.loopH:
    movu        m2, [r0]                     ; m2 = src[0-7]
    movu        m3, [r0 + 16]                ; m3 = src[8-15]

    pmaddwd     m4, m2, m0
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m3, m2, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m3, m2, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m7, m3, m2, 14               ; m7 = src[7-14]
    pmaddwd     m7, m0
    phaddd      m6, m7
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif
    movu        [r2], m4

    movu        m2, [r0 + 32]                ; m2 = src[16-23]

    pmaddwd     m4, m3, m0                   ; m3 = src[8-15]
    palignr     m5, m2, m3, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m2, m3, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m2, m3, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m2, m3, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m2, m3, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m2, m3, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m7, m2, m3, 14               ; m7 = src[7-14]
    pmaddwd     m7, m0
    phaddd      m6, m7
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif
    movu        [r2 + 16], m4

    movu        m3, [r0 + 48]                ; m3 = src[24-31]

    pmaddwd     m4, m2, m0                   ; m2 = src[16-23]
    palignr     m5, m3, m2, 2                ; m5 = src[1-8]
    pmaddwd     m5, m0
    phaddd      m4, m5

    palignr     m5, m3, m2, 4                ; m5 = src[2-9]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 6                ; m6 = src[3-10]
    pmaddwd     m6, m0
    phaddd      m5, m6
    phaddd      m4, m5
    paddd       m4, m1

    palignr     m5, m3, m2, 8                ; m5 = src[4-11]
    pmaddwd     m5, m0
    palignr     m6, m3, m2, 10               ; m6 = src[5-12]
    pmaddwd     m6, m0
    phaddd      m5, m6

    palignr     m6, m3, m2, 12               ; m6 = src[6-13]
    pmaddwd     m6, m0
    palignr     m7, m3, m2, 14               ; m7 = src[7-14]
    pmaddwd     m7, m0
    phaddd      m6, m7
    phaddd      m5, m6
    paddd       m5, m1
%ifidn %3, pp 
    psrad       m4, 6
    psrad       m5, 6
    packusdw    m4, m5
    pxor        m5, m5
    CLIPW       m4, m5, [pw_pixel_max]
%else
    psrad       m4, 2
    psrad       m5, 2
    packssdw    m4, m5
%endif
    movu        [r2 + 32], m4

    add         r0, r1
    add         r2, r3

    dec         r4d
    jnz         .loopH
    RET
%endmacro

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_pp_24x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx
;-------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W24 24, 32, pp

;----------------------------------------------------------------------------------------------------------------------------
; void interp_8tap_horiz_ps_24x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx, int isRowExt)
;----------------------------------------------------------------------------------------------------------------------------
FILTER_HOR_LUMA_W24 24, 32, ps

%macro FILTER_W2_2 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + r1]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1
%ifidn %1, pp
    psrad       m3,         6
    packusdw    m3,         m3
    CLIPW       m3,         m7,    m6
%else
    psrad       m3,         2
    packssdw    m3,         m3
%endif
    movd        [r2],       m3
    pextrd      [r2 + r3],  m3, 1
%endmacro

%macro FILTER_W4_2 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + r1]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + r1 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m7,    m6
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + r3],  m3
%endmacro

;-----------------------------------------------------------------------------
; void interp_4tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------
%macro FILTER_CHROMA_H 6
INIT_XMM sse4
cglobal interp_4tap_horiz_%3_%1x%2, 4, %4, %5

    add         r3,       r3
    add         r1,       r1
    sub         r0,       2
    mov         r4d,      r4m
    add         r4d,      r4d

%ifdef PIC
    lea         r%6,      [tab_ChromaCoeff]
    movh        m0,       [r%6 + r4 * 4]
%else
    movh        m0,       [tab_ChromaCoeff + r4 * 4]
%endif

    punpcklqdq  m0,       m0
    mova        m2,       [tab_Tm16]

%ifidn %3, ps
    mova        m1,       [tab_c_n32768]
    cmp         r5m, byte 0
    je          .skip
    sub         r0, r1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0

    %if %1 == 4
        movu        m4,         [r0 + 4]
        pshufb      m4,         m4, m2
        pmaddwd     m4,         m0
        phaddd      m3,         m4
    %else
        phaddd      m3,         m3
    %endif

    paddd       m3,         m1
    psrad       m3,         2
    packssdw    m3,         m3

    %if %1 == 2
        movd        [r2],       m3
    %else
        movh        [r2],       m3
    %endif

    add         r0, r1
    add         r2, r3
    FILTER_W%1_2 %3
    lea         r0,       [r0 + 2 * r1]
    lea         r2,       [r2 + 2 * r3]

.skip:

%else     ;%ifidn %3, ps
    pxor        m7,       m7
    mova        m6,       [pw_pixel_max]
    mova        m1,       [tab_c_32]
%endif    ;%ifidn %3, ps

    FILTER_W%1_2 %3

%rep (%2/2) - 1
    lea         r0,       [r0 + 2 * r1]
    lea         r2,       [r2 + 2 * r3]
    FILTER_W%1_2 %3
%endrep

RET
%endmacro

FILTER_CHROMA_H 2, 4, pp, 6, 8, 5
FILTER_CHROMA_H 2, 8, pp, 6, 8, 5
FILTER_CHROMA_H 4, 2, pp, 6, 8, 5
FILTER_CHROMA_H 4, 4, pp, 6, 8, 5
FILTER_CHROMA_H 4, 8, pp, 6, 8, 5
FILTER_CHROMA_H 4, 16, pp, 6, 8, 5

FILTER_CHROMA_H 2, 4, ps, 7, 5, 6
FILTER_CHROMA_H 2, 8, ps, 7, 5, 6
FILTER_CHROMA_H 4, 2, ps, 7, 6, 6
FILTER_CHROMA_H 4, 4, ps, 7, 6, 6
FILTER_CHROMA_H 4, 8, ps, 7, 6, 6
FILTER_CHROMA_H 4, 16, ps, 7, 6, 6

FILTER_CHROMA_H 2, 16, pp, 6, 8, 5
FILTER_CHROMA_H 4, 32, pp, 6, 8, 5
FILTER_CHROMA_H 2, 16, ps, 7, 5, 6
FILTER_CHROMA_H 4, 32, ps, 7, 6, 6


%macro FILTER_W6_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m4,         [r0 + 8]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m4,         m4
    paddd       m4,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m4,         6
    packusdw    m3,         m4
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m4,         2
    packssdw    m3,         m4
%endif
    movh        [r2],       m3
    pextrd      [r2 + 8],   m3, 2
%endmacro

cglobal chroma_filter_pp_6x1_internal
    FILTER_W6_1 pp
    ret

cglobal chroma_filter_ps_6x1_internal
    FILTER_W6_1 ps
    ret

%macro FILTER_W8_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + 8],   m3
%endmacro

cglobal chroma_filter_pp_8x1_internal
    FILTER_W8_1 pp
    ret

cglobal chroma_filter_ps_8x1_internal
    FILTER_W8_1 ps
    ret

%macro FILTER_W12_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + 8],   m3

    movu        m3,         [r0 + 16]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 20]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

%ifidn %1, pp
    psrad       m3,         6
    packusdw    m3,         m3
    CLIPW       m3,         m6, m7
%else
    psrad       m3,         2
    packssdw    m3,         m3
%endif
    movh        [r2 + 16],  m3
%endmacro

cglobal chroma_filter_pp_12x1_internal
    FILTER_W12_1 pp
    ret

cglobal chroma_filter_ps_12x1_internal
    FILTER_W12_1 ps
    ret

%macro FILTER_W16_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + 8],   m3

    movu        m3,         [r0 + 16]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 20]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 24]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 28]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 16],  m3
    movhps      [r2 + 24],  m3
%endmacro

cglobal chroma_filter_pp_16x1_internal
    FILTER_W16_1 pp
    ret

cglobal chroma_filter_ps_16x1_internal
    FILTER_W16_1 ps
    ret

%macro FILTER_W24_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + 8],   m3

    movu        m3,         [r0 + 16]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 20]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 24]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 28]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 16],  m3
    movhps      [r2 + 24],  m3

    movu        m3,         [r0 + 32]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 36]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 40]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 44]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 32],  m3
    movhps      [r2 + 40],  m3
%endmacro

cglobal chroma_filter_pp_24x1_internal
    FILTER_W24_1 pp
    ret

cglobal chroma_filter_ps_24x1_internal
    FILTER_W24_1 ps
    ret

%macro FILTER_W32_1 1
    movu        m3,         [r0]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2],       m3
    movhps      [r2 + 8],   m3

    movu        m3,         [r0 + 16]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 20]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 24]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 28]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 16],  m3
    movhps      [r2 + 24],  m3

    movu        m3,         [r0 + 32]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 36]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 40]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 44]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 32],  m3
    movhps      [r2 + 40],  m3

    movu        m3,         [r0 + 48]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + 52]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + 56]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + 60]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + 48],  m3
    movhps      [r2 + 56],  m3
%endmacro

cglobal chroma_filter_pp_32x1_internal
    FILTER_W32_1 pp
    ret

cglobal chroma_filter_ps_32x1_internal
    FILTER_W32_1 ps
    ret

%macro FILTER_W8o_1 2
    movu        m3,         [r0 + %2]
    pshufb      m3,         m3, m2
    pmaddwd     m3,         m0
    movu        m4,         [r0 + %2 + 4]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m3,         m4
    paddd       m3,         m1

    movu        m5,         [r0 + %2 + 8]
    pshufb      m5,         m5, m2
    pmaddwd     m5,         m0
    movu        m4,         [r0 + %2 + 12]
    pshufb      m4,         m4, m2
    pmaddwd     m4,         m0
    phaddd      m5,         m4
    paddd       m5,         m1
%ifidn %1, pp
    psrad       m3,         6
    psrad       m5,         6
    packusdw    m3,         m5
    CLIPW       m3,         m6,    m7
%else
    psrad       m3,         2
    psrad       m5,         2
    packssdw    m3,         m5
%endif
    movh        [r2 + %2],       m3
    movhps      [r2 + %2 + 8],   m3
%endmacro

%macro FILTER_W48_1 1
    FILTER_W8o_1 %1, 0
    FILTER_W8o_1 %1, 16
    FILTER_W8o_1 %1, 32
    FILTER_W8o_1 %1, 48
    FILTER_W8o_1 %1, 64
    FILTER_W8o_1 %1, 80
%endmacro

cglobal chroma_filter_pp_48x1_internal
    FILTER_W48_1 pp
    ret

cglobal chroma_filter_ps_48x1_internal
    FILTER_W48_1 ps
    ret

%macro FILTER_W64_1 1
    FILTER_W8o_1 %1, 0
    FILTER_W8o_1 %1, 16
    FILTER_W8o_1 %1, 32
    FILTER_W8o_1 %1, 48
    FILTER_W8o_1 %1, 64
    FILTER_W8o_1 %1, 80
    FILTER_W8o_1 %1, 96
    FILTER_W8o_1 %1, 112
%endmacro

cglobal chroma_filter_pp_64x1_internal
    FILTER_W64_1 pp
    ret

cglobal chroma_filter_ps_64x1_internal
    FILTER_W64_1 ps
    ret


;-----------------------------------------------------------------------------
; void interp_4tap_horiz_%3_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------

INIT_XMM sse4
%macro IPFILTER_CHROMA 6
cglobal interp_4tap_horiz_%3_%1x%2, 4, %5, %6

    add         r3,        r3
    add         r1,        r1
    sub         r0,         2
    mov         r4d,        r4m
    add         r4d,        r4d

%ifdef PIC
    lea         r%4,       [tab_ChromaCoeff]
    movh        m0,       [r%4 + r4 * 4]
%else
    movh        m0,       [tab_ChromaCoeff + r4 * 4]
%endif

    punpcklqdq  m0,       m0
    mova        m2,       [tab_Tm16]

%ifidn %3, ps
    mova        m1,       [tab_c_n32768]
    cmp         r5m, byte 0
    je          .skip
    sub         r0, r1
    call chroma_filter_%3_%1x1_internal
    add         r0, r1
    add         r2, r3
    call chroma_filter_%3_%1x1_internal
    add         r0, r1
    add         r2, r3
    call chroma_filter_%3_%1x1_internal
    add         r0, r1
    add         r2, r3
.skip:
%else
    mova        m1,         [tab_c_32]
    pxor        m6,         m6
    mova        m7,         [pw_pixel_max]
%endif

    call chroma_filter_%3_%1x1_internal
%rep %2 - 1
    add         r0,       r1
    add         r2,       r3
    call chroma_filter_%3_%1x1_internal
%endrep
RET
%endmacro
IPFILTER_CHROMA 6, 8, pp, 5, 6, 8
IPFILTER_CHROMA 8, 2, pp, 5, 6, 8
IPFILTER_CHROMA 8, 4, pp, 5, 6, 8
IPFILTER_CHROMA 8, 6, pp, 5, 6, 8
IPFILTER_CHROMA 8, 8, pp, 5, 6, 8
IPFILTER_CHROMA 8, 16, pp, 5, 6, 8
IPFILTER_CHROMA 8, 32, pp, 5, 6, 8
IPFILTER_CHROMA 12, 16, pp, 5, 6, 8
IPFILTER_CHROMA 16, 4, pp, 5, 6, 8
IPFILTER_CHROMA 16, 8, pp, 5, 6, 8
IPFILTER_CHROMA 16, 12, pp, 5, 6, 8
IPFILTER_CHROMA 16, 16, pp, 5, 6, 8
IPFILTER_CHROMA 16, 32, pp, 5, 6, 8
IPFILTER_CHROMA 24, 32, pp, 5, 6, 8
IPFILTER_CHROMA 32, 8, pp, 5, 6, 8
IPFILTER_CHROMA 32, 16, pp, 5, 6, 8
IPFILTER_CHROMA 32, 24, pp, 5, 6, 8
IPFILTER_CHROMA 32, 32, pp, 5, 6, 8

IPFILTER_CHROMA 6, 8, ps, 6, 7, 6
IPFILTER_CHROMA 8, 2, ps, 6, 7, 6
IPFILTER_CHROMA 8, 4, ps, 6, 7, 6
IPFILTER_CHROMA 8, 6, ps, 6, 7, 6
IPFILTER_CHROMA 8, 8, ps, 6, 7, 6
IPFILTER_CHROMA 8, 16, ps, 6, 7, 6
IPFILTER_CHROMA 8, 32, ps, 6, 7, 6
IPFILTER_CHROMA 12, 16, ps, 6, 7, 6
IPFILTER_CHROMA 16, 4, ps, 6, 7, 6
IPFILTER_CHROMA 16, 8, ps, 6, 7, 6
IPFILTER_CHROMA 16, 12, ps, 6, 7, 6
IPFILTER_CHROMA 16, 16, ps, 6, 7, 6
IPFILTER_CHROMA 16, 32, ps, 6, 7, 6
IPFILTER_CHROMA 24, 32, ps, 6, 7, 6
IPFILTER_CHROMA 32, 8, ps, 6, 7, 6
IPFILTER_CHROMA 32, 16, ps, 6, 7, 6
IPFILTER_CHROMA 32, 24, ps, 6, 7, 6
IPFILTER_CHROMA 32, 32, ps, 6, 7, 6

IPFILTER_CHROMA 6, 16, pp, 5, 6, 8
IPFILTER_CHROMA 8, 12, pp, 5, 6, 8
IPFILTER_CHROMA 8, 64, pp, 5, 6, 8
IPFILTER_CHROMA 12, 32, pp, 5, 6, 8
IPFILTER_CHROMA 16, 24, pp, 5, 6, 8
IPFILTER_CHROMA 16, 64, pp, 5, 6, 8
IPFILTER_CHROMA 24, 64, pp, 5, 6, 8
IPFILTER_CHROMA 32, 48, pp, 5, 6, 8
IPFILTER_CHROMA 32, 64, pp, 5, 6, 8
IPFILTER_CHROMA 6, 16, ps, 6, 7, 6
IPFILTER_CHROMA 8, 12, ps, 6, 7, 6
IPFILTER_CHROMA 8, 64, ps, 6, 7, 6
IPFILTER_CHROMA 12, 32, ps, 6, 7, 6
IPFILTER_CHROMA 16, 24, ps, 6, 7, 6
IPFILTER_CHROMA 16, 64, ps, 6, 7, 6
IPFILTER_CHROMA 24, 64, ps, 6, 7, 6
IPFILTER_CHROMA 32, 48, ps, 6, 7, 6
IPFILTER_CHROMA 32, 64, ps, 6, 7, 6

IPFILTER_CHROMA 48, 64, pp, 5, 6, 8
IPFILTER_CHROMA 64, 48, pp, 5, 6, 8
IPFILTER_CHROMA 64, 64, pp, 5, 6, 8
IPFILTER_CHROMA 64, 32, pp, 5, 6, 8
IPFILTER_CHROMA 64, 16, pp, 5, 6, 8
IPFILTER_CHROMA 48, 64, ps, 6, 7, 6
IPFILTER_CHROMA 64, 48, ps, 6, 7, 6
IPFILTER_CHROMA 64, 64, ps, 6, 7, 6
IPFILTER_CHROMA 64, 32, ps, 6, 7, 6
IPFILTER_CHROMA 64, 16, ps, 6, 7, 6


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

;-----------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_%3_%1x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-----------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_SS 4
INIT_XMM sse2
cglobal interp_4tap_vert_%3_%1x%2, 5, 7, %4 ,0-gprsize

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

%ifnidn %3, ss
    %ifnidn %3, ps
        mova      m7, [pw_pixel_max]
        %ifidn %3, pp
            mova      m6, [tab_c_32]
        %else
            mova      m6, [tab_c_524800]
        %endif
    %else
        mova      m6, [tab_c_n32768]
    %endif
%endif

.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_CHROMA_SP_W4_4R

%ifidn %3, ss
    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3
%elifidn %3, ps
    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6
    psrad     m0, 2
    psrad     m1, 2
    psrad     m2, 2
    psrad     m3, 2

    packssdw  m0, m1
    packssdw  m2, m3
%else
    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6
    %ifidn %3, pp
        psrad     m0, 6
        psrad     m1, 6
        psrad     m2, 6
        psrad     m3, 6
    %else
        psrad     m0, 10
        psrad     m1, 10
        psrad     m2, 10
        psrad     m3, 10
    %endif
    packssdw  m0, m1
    packssdw  m2, m3
    pxor      m5, m5
    CLIPW2    m0, m2, m5, m7
%endif

    movh      [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movh      [r5], m2
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

    FILTER_VER_CHROMA_SS 4, 4, ss, 6
    FILTER_VER_CHROMA_SS 4, 8, ss, 6
    FILTER_VER_CHROMA_SS 16, 16, ss, 6
    FILTER_VER_CHROMA_SS 16, 8, ss, 6
    FILTER_VER_CHROMA_SS 16, 12, ss, 6
    FILTER_VER_CHROMA_SS 12, 16, ss, 6
    FILTER_VER_CHROMA_SS 16, 4, ss, 6
    FILTER_VER_CHROMA_SS 4, 16, ss, 6
    FILTER_VER_CHROMA_SS 32, 32, ss, 6
    FILTER_VER_CHROMA_SS 32, 16, ss, 6
    FILTER_VER_CHROMA_SS 16, 32, ss, 6
    FILTER_VER_CHROMA_SS 32, 24, ss, 6
    FILTER_VER_CHROMA_SS 24, 32, ss, 6
    FILTER_VER_CHROMA_SS 32, 8, ss, 6

    FILTER_VER_CHROMA_SS 4, 4, ps, 7
    FILTER_VER_CHROMA_SS 4, 8, ps, 7
    FILTER_VER_CHROMA_SS 16, 16, ps, 7
    FILTER_VER_CHROMA_SS 16, 8, ps, 7
    FILTER_VER_CHROMA_SS 16, 12, ps, 7
    FILTER_VER_CHROMA_SS 12, 16, ps, 7
    FILTER_VER_CHROMA_SS 16, 4, ps, 7
    FILTER_VER_CHROMA_SS 4, 16, ps, 7
    FILTER_VER_CHROMA_SS 32, 32, ps, 7
    FILTER_VER_CHROMA_SS 32, 16, ps, 7
    FILTER_VER_CHROMA_SS 16, 32, ps, 7
    FILTER_VER_CHROMA_SS 32, 24, ps, 7
    FILTER_VER_CHROMA_SS 24, 32, ps, 7
    FILTER_VER_CHROMA_SS 32, 8, ps, 7

    FILTER_VER_CHROMA_SS 4, 4, sp, 8
    FILTER_VER_CHROMA_SS 4, 8, sp, 8
    FILTER_VER_CHROMA_SS 16, 16, sp, 8
    FILTER_VER_CHROMA_SS 16, 8, sp, 8
    FILTER_VER_CHROMA_SS 16, 12, sp, 8
    FILTER_VER_CHROMA_SS 12, 16, sp, 8
    FILTER_VER_CHROMA_SS 16, 4, sp, 8
    FILTER_VER_CHROMA_SS 4, 16, sp, 8
    FILTER_VER_CHROMA_SS 32, 32, sp, 8
    FILTER_VER_CHROMA_SS 32, 16, sp, 8
    FILTER_VER_CHROMA_SS 16, 32, sp, 8
    FILTER_VER_CHROMA_SS 32, 24, sp, 8
    FILTER_VER_CHROMA_SS 24, 32, sp, 8
    FILTER_VER_CHROMA_SS 32, 8, sp, 8

    FILTER_VER_CHROMA_SS 4, 4, pp, 8
    FILTER_VER_CHROMA_SS 4, 8, pp, 8
    FILTER_VER_CHROMA_SS 16, 16, pp, 8
    FILTER_VER_CHROMA_SS 16, 8, pp, 8
    FILTER_VER_CHROMA_SS 16, 12, pp, 8
    FILTER_VER_CHROMA_SS 12, 16, pp, 8
    FILTER_VER_CHROMA_SS 16, 4, pp, 8
    FILTER_VER_CHROMA_SS 4, 16, pp, 8
    FILTER_VER_CHROMA_SS 32, 32, pp, 8
    FILTER_VER_CHROMA_SS 32, 16, pp, 8
    FILTER_VER_CHROMA_SS 16, 32, pp, 8
    FILTER_VER_CHROMA_SS 32, 24, pp, 8
    FILTER_VER_CHROMA_SS 24, 32, pp, 8
    FILTER_VER_CHROMA_SS 32, 8, pp, 8


    FILTER_VER_CHROMA_SS 16, 24, ss, 6
    FILTER_VER_CHROMA_SS 12, 32, ss, 6
    FILTER_VER_CHROMA_SS 4, 32, ss, 6
    FILTER_VER_CHROMA_SS 32, 64, ss, 6
    FILTER_VER_CHROMA_SS 16, 64, ss, 6
    FILTER_VER_CHROMA_SS 32, 48, ss, 6
    FILTER_VER_CHROMA_SS 24, 64, ss, 6

    FILTER_VER_CHROMA_SS 16, 24, ps, 7
    FILTER_VER_CHROMA_SS 12, 32, ps, 7
    FILTER_VER_CHROMA_SS 4, 32, ps, 7
    FILTER_VER_CHROMA_SS 32, 64, ps, 7
    FILTER_VER_CHROMA_SS 16, 64, ps, 7
    FILTER_VER_CHROMA_SS 32, 48, ps, 7
    FILTER_VER_CHROMA_SS 24, 64, ps, 7

    FILTER_VER_CHROMA_SS 16, 24, sp, 8
    FILTER_VER_CHROMA_SS 12, 32, sp, 8
    FILTER_VER_CHROMA_SS 4, 32, sp, 8
    FILTER_VER_CHROMA_SS 32, 64, sp, 8
    FILTER_VER_CHROMA_SS 16, 64, sp, 8
    FILTER_VER_CHROMA_SS 32, 48, sp, 8
    FILTER_VER_CHROMA_SS 24, 64, sp, 8

    FILTER_VER_CHROMA_SS 16, 24, pp, 8
    FILTER_VER_CHROMA_SS 12, 32, pp, 8
    FILTER_VER_CHROMA_SS 4, 32, pp, 8
    FILTER_VER_CHROMA_SS 32, 64, pp, 8
    FILTER_VER_CHROMA_SS 16, 64, pp, 8
    FILTER_VER_CHROMA_SS 32, 48, pp, 8
    FILTER_VER_CHROMA_SS 24, 64, pp, 8


    FILTER_VER_CHROMA_SS 48, 64, ss, 6
    FILTER_VER_CHROMA_SS 64, 48, ss, 6
    FILTER_VER_CHROMA_SS 64, 64, ss, 6
    FILTER_VER_CHROMA_SS 64, 32, ss, 6
    FILTER_VER_CHROMA_SS 64, 16, ss, 6

    FILTER_VER_CHROMA_SS 48, 64, ps, 7
    FILTER_VER_CHROMA_SS 64, 48, ps, 7
    FILTER_VER_CHROMA_SS 64, 64, ps, 7
    FILTER_VER_CHROMA_SS 64, 32, ps, 7
    FILTER_VER_CHROMA_SS 64, 16, ps, 7

    FILTER_VER_CHROMA_SS 48, 64, sp, 8
    FILTER_VER_CHROMA_SS 64, 48, sp, 8
    FILTER_VER_CHROMA_SS 64, 64, sp, 8
    FILTER_VER_CHROMA_SS 64, 32, sp, 8
    FILTER_VER_CHROMA_SS 64, 16, sp, 8

    FILTER_VER_CHROMA_SS 48, 64, pp, 8
    FILTER_VER_CHROMA_SS 64, 48, pp, 8
    FILTER_VER_CHROMA_SS 64, 64, pp, 8
    FILTER_VER_CHROMA_SS 64, 32, pp, 8
    FILTER_VER_CHROMA_SS 64, 16, pp, 8


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

;---------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_%2_2x%1(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_W2 3
INIT_XMM sse4
cglobal interp_4tap_vert_%2_2x%1, 5, 6, %3

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

    mov       r4d, (%1/4)
%ifnidn %2, ss
    %ifnidn %2, ps
        pxor      m7, m7
        mova      m6, [pw_pixel_max]
        %ifidn %2, pp
            mova      m5, [tab_c_32]
        %else
            mova      m5, [tab_c_524800]
        %endif
    %else
        mova      m5, [tab_c_n32768]
    %endif
%endif

.loopH:
    PROCESS_CHROMA_SP_W2_4R r5
%ifidn %2, ss
    psrad     m0, 6
    psrad     m2, 6
    packssdw  m0, m2
%elifidn %2, ps
    paddd     m0, m5
    paddd     m2, m5
    psrad     m0, 2
    psrad     m2, 2
    packssdw  m0, m2
%else
    paddd     m0, m5
    paddd     m2, m5
    %ifidn %2, pp
        psrad     m0, 6
        psrad     m2, 6
    %else
        psrad     m0, 10
        psrad     m2, 10
    %endif
    packusdw  m0, m2
    CLIPW     m0, m7,    m6
%endif

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

FILTER_VER_CHROMA_W2 4, ss, 5
FILTER_VER_CHROMA_W2 8, ss, 5

FILTER_VER_CHROMA_W2 4, pp, 8
FILTER_VER_CHROMA_W2 8, pp, 8

FILTER_VER_CHROMA_W2 4, ps, 6
FILTER_VER_CHROMA_W2 8, ps, 6

FILTER_VER_CHROMA_W2 4, sp, 8
FILTER_VER_CHROMA_W2 8, sp, 8

FILTER_VER_CHROMA_W2 16, ss, 5
FILTER_VER_CHROMA_W2 16, pp, 8
FILTER_VER_CHROMA_W2 16, ps, 6
FILTER_VER_CHROMA_W2 16, sp, 8


;---------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_%1_4x2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_W4 3
INIT_XMM sse4
cglobal interp_4tap_vert_%2_4x%1, 5, 6, %3

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

%ifnidn %2, 2
    mov        r4d, %1/2
%endif

%ifnidn %2, ss
    %ifnidn %2, ps
        pxor      m6, m6
        mova      m5, [pw_pixel_max]
        %ifidn %2, pp
            mova      m4, [tab_c_32]
        %else
            mova      m4, [tab_c_524800]
        %endif
    %else
        mova      m4, [tab_c_n32768]
    %endif
%endif

%ifnidn %2, 2
.loop:
%endif

    movh       m0, [r0]
    movh       m1, [r0 + r1]
    punpcklwd  m0, m1                          ;m0=[0 1]
    pmaddwd    m0, [r5 + 0 *16]                ;m0=[0+1]  Row1

    lea        r0, [r0 + 2 * r1]
    movh       m2, [r0]
    punpcklwd  m1, m2                          ;m1=[1 2]
    pmaddwd    m1, [r5 + 0 *16]                ;m1=[1+2]  Row2

    movh       m3, [r0 + r1]
    punpcklwd  m2, m3                          ;m4=[2 3]
    pmaddwd    m2, [r5 + 1 * 16]
    paddd      m0, m2                          ;m0=[0+1+2+3]  Row1 done

    movh       m2, [r0 + 2 * r1]
    punpcklwd  m3, m2                          ;m5=[3 4]
    pmaddwd    m3, [r5 + 1 * 16]
    paddd      m1, m3                          ;m1=[1+2+3+4]  Row2 done

%ifidn %2, ss
    psrad     m0, 6
    psrad     m1, 6
    packssdw  m0, m1
%elifidn %2, ps
    paddd     m0, m4
    paddd     m1, m4
    psrad     m0, 2
    psrad     m1, 2
    packssdw  m0, m1
%else
    paddd     m0, m4
    paddd     m1, m4
    %ifidn %2, pp
        psrad     m0, 6
        psrad     m1, 6
    %else
        psrad     m0, 10
        psrad     m1, 10
    %endif
    packusdw  m0, m1
    CLIPW     m0, m6,    m5
%endif

    movh       [r2], m0
    movhps     [r2 + r3], m0

%ifnidn %2, 2
    lea        r2, [r2 + r3 * 2]
    dec        r4d
    jnz        .loop
%endif

    RET
%endmacro

FILTER_VER_CHROMA_W4 2, ss, 4
FILTER_VER_CHROMA_W4 2, pp, 7
FILTER_VER_CHROMA_W4 2, ps, 5
FILTER_VER_CHROMA_W4 2, sp, 7

FILTER_VER_CHROMA_W4 4, ss, 4
FILTER_VER_CHROMA_W4 4, pp, 7
FILTER_VER_CHROMA_W4 4, ps, 5
FILTER_VER_CHROMA_W4 4, sp, 7

;-------------------------------------------------------------------------------------------------------------------
; void interp_4tap_vertical_%1_6x8(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_W6 3
INIT_XMM sse4
cglobal interp_4tap_vert_%2_6x%1, 5, 7, %3

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

    mov       r4d, %1/4

%ifnidn %2, ss
    %ifnidn %2, ps
        mova      m7, [pw_pixel_max]
        %ifidn %2, pp
            mova      m6, [tab_c_32]
        %else
            mova      m6, [tab_c_524800]
        %endif
    %else
        mova      m6, [tab_c_n32768]
    %endif
%endif

.loopH:
    PROCESS_CHROMA_SP_W4_4R

%ifidn %2, ss
    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3
%elifidn %2, ps
    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6
    psrad     m0, 2
    psrad     m1, 2
    psrad     m2, 2
    psrad     m3, 2

    packssdw  m0, m1
    packssdw  m2, m3
%else
    paddd     m0, m6
    paddd     m1, m6
    paddd     m2, m6
    paddd     m3, m6
    %ifidn %2, pp
        psrad     m0, 6
        psrad     m1, 6
        psrad     m2, 6
        psrad     m3, 6
    %else
        psrad     m0, 10
        psrad     m1, 10
        psrad     m2, 10
        psrad     m3, 10
    %endif
    packssdw  m0, m1
    packssdw  m2, m3
    pxor      m5, m5
    CLIPW2    m0, m2, m5, m7
%endif

    movh      [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movh      [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [4 * r1 - 2 * 4]
    sub       r0, r5
    add       r2, 2 * 4

    PROCESS_CHROMA_SP_W2_4R r6

%ifidn %2, ss
    psrad     m0, 6
    psrad     m2, 6
    packssdw  m0, m2
%elifidn %2, ps
    paddd     m0, m6
    paddd     m2, m6
    psrad     m0, 2
    psrad     m2, 2
    packssdw  m0, m2
%else
    paddd     m0, m6
    paddd     m2, m6
    %ifidn %2, pp
        psrad     m0, 6
        psrad     m2, 6
    %else
        psrad     m0, 10
        psrad     m2, 10
    %endif
    packusdw  m0, m2
    CLIPW     m0, m5,    m7
%endif

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

FILTER_VER_CHROMA_W6 8, ss, 6
FILTER_VER_CHROMA_W6 8, ps, 7
FILTER_VER_CHROMA_W6 8, sp, 8
FILTER_VER_CHROMA_W6 8, pp, 8

FILTER_VER_CHROMA_W6 16, ss, 6
FILTER_VER_CHROMA_W6 16, ps, 7
FILTER_VER_CHROMA_W6 16, sp, 8
FILTER_VER_CHROMA_W6 16, pp, 8

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

;----------------------------------------------------------------------------------------------------------------
; void interp_4tap_vert_%3_%1x%2(int16_t *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;----------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_CHROMA_W8 4
INIT_XMM sse2
cglobal interp_4tap_vert_%3_%1x%2, 5, 6, %4

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

%ifidn %3, pp
    mova      m7, [tab_c_32]
%elifidn %3, sp
    mova      m7, [tab_c_524800]
%elifidn %3, ps
    mova      m7, [tab_c_n32768]
%endif

.loopH:
    PROCESS_CHROMA_SP_W8_2R

%ifidn %3, ss
    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3
%elifidn %3, ps
    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7
    psrad     m0, 2
    psrad     m1, 2
    psrad     m2, 2
    psrad     m3, 2

    packssdw  m0, m1
    packssdw  m2, m3
%else
    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7
    %ifidn %3, pp
        psrad     m0, 6
        psrad     m1, 6
        psrad     m2, 6
        psrad     m3, 6
    %else
        psrad     m0, 10
        psrad     m1, 10
        psrad     m2, 10
        psrad     m3, 10
    %endif
    packssdw  m0, m1
    packssdw  m2, m3
    pxor      m5, m5
    mova      m6, [pw_pixel_max]
    CLIPW2    m0, m2, m5, m6
%endif

    movu      [r2], m0
    movu      [r2 + r3], m2

    lea       r2, [r2 + 2 * r3]

    dec       r4d
    jnz       .loopH

    RET
%endmacro

FILTER_VER_CHROMA_W8 8, 2, ss, 7
FILTER_VER_CHROMA_W8 8, 4, ss, 7
FILTER_VER_CHROMA_W8 8, 6, ss, 7
FILTER_VER_CHROMA_W8 8, 8, ss, 7
FILTER_VER_CHROMA_W8 8, 16, ss, 7
FILTER_VER_CHROMA_W8 8, 32, ss, 7

FILTER_VER_CHROMA_W8 8, 2, sp, 8
FILTER_VER_CHROMA_W8 8, 4, sp, 8
FILTER_VER_CHROMA_W8 8, 6, sp, 8
FILTER_VER_CHROMA_W8 8, 8, sp, 8
FILTER_VER_CHROMA_W8 8, 16, sp, 8
FILTER_VER_CHROMA_W8 8, 32, sp, 8

FILTER_VER_CHROMA_W8 8, 2, ps, 8
FILTER_VER_CHROMA_W8 8, 4, ps, 8
FILTER_VER_CHROMA_W8 8, 6, ps, 8
FILTER_VER_CHROMA_W8 8, 8, ps, 8
FILTER_VER_CHROMA_W8 8, 16, ps, 8
FILTER_VER_CHROMA_W8 8, 32, ps, 8

FILTER_VER_CHROMA_W8 8, 2, pp, 8
FILTER_VER_CHROMA_W8 8, 4, pp, 8
FILTER_VER_CHROMA_W8 8, 6, pp, 8
FILTER_VER_CHROMA_W8 8, 8, pp, 8
FILTER_VER_CHROMA_W8 8, 16, pp, 8
FILTER_VER_CHROMA_W8 8, 32, pp, 8

FILTER_VER_CHROMA_W8 8, 12, ss, 7
FILTER_VER_CHROMA_W8 8, 64, ss, 7
FILTER_VER_CHROMA_W8 8, 12, sp, 8
FILTER_VER_CHROMA_W8 8, 64, sp, 8
FILTER_VER_CHROMA_W8 8, 12, ps, 8
FILTER_VER_CHROMA_W8 8, 64, ps, 8
FILTER_VER_CHROMA_W8 8, 12, pp, 8
FILTER_VER_CHROMA_W8 8, 64, pp, 8


INIT_XMM sse2
cglobal chroma_p2s, 3, 7, 3

    ; load width and height
    mov         r3d, r3m
    mov         r4d, r4m
    add         r1, r1

    ; load constant
    mova        m2, [tab_c_n8192]

.loopH:

    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5 * 2]

    movu        m0, [r6]
    psllw       m0, 4
    paddw       m0, m2

    movu        m1, [r6 + r1]
    psllw       m1, 4
    paddw       m1, m2

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

%macro PROCESS_LUMA_VER_W4_4R 0
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
; void interp_8tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_PP 2
INIT_XMM sse4
cglobal interp_8tap_vert_pp_%1x%2, 5, 7, 8 ,0-gprsize

    add       r1d, r1d
    add       r3d, r3d
    lea       r5, [r1 + 2 * r1]
    sub       r0, r5
    shl       r4d, 6

%ifdef PIC
    lea       r5, [tab_LumaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffV + r4]
%endif

    mova      m7, [pd_32]

    mov       dword [rsp], %2/4
.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_LUMA_VER_W4_4R

    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7

    psrad     m0, 6
    psrad     m1, 6
    psrad     m2, 6
    psrad     m3, 6

    packssdw  m0, m1
    packssdw  m2, m3

    pxor      m1, m1
    CLIPW2    m0, m2, m1, [pw_pixel_max]

    movh      [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movh      [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [8 * r1 - 2 * 4]
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

;-------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_pp_%1x%2(pixel *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;-------------------------------------------------------------------------------------------------------------
    FILTER_VER_LUMA_PP 4, 4
    FILTER_VER_LUMA_PP 8, 8
    FILTER_VER_LUMA_PP 8, 4
    FILTER_VER_LUMA_PP 4, 8
    FILTER_VER_LUMA_PP 16, 16
    FILTER_VER_LUMA_PP 16, 8
    FILTER_VER_LUMA_PP 8, 16
    FILTER_VER_LUMA_PP 16, 12
    FILTER_VER_LUMA_PP 12, 16
    FILTER_VER_LUMA_PP 16, 4
    FILTER_VER_LUMA_PP 4, 16
    FILTER_VER_LUMA_PP 32, 32
    FILTER_VER_LUMA_PP 32, 16
    FILTER_VER_LUMA_PP 16, 32
    FILTER_VER_LUMA_PP 32, 24
    FILTER_VER_LUMA_PP 24, 32
    FILTER_VER_LUMA_PP 32, 8
    FILTER_VER_LUMA_PP 8, 32
    FILTER_VER_LUMA_PP 64, 64
    FILTER_VER_LUMA_PP 64, 32
    FILTER_VER_LUMA_PP 32, 64
    FILTER_VER_LUMA_PP 64, 48
    FILTER_VER_LUMA_PP 48, 64
    FILTER_VER_LUMA_PP 64, 16
    FILTER_VER_LUMA_PP 16, 64

;---------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_%1x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_PS 2
INIT_XMM sse4
cglobal interp_8tap_vert_ps_%1x%2, 5, 7, 8 ,0-gprsize

    add       r1d, r1d
    add       r3d, r3d
    lea       r5, [r1 + 2 * r1]
    sub       r0, r5
    shl       r4d, 6

%ifdef PIC
    lea       r5, [tab_LumaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffV + r4]
%endif

    mova      m7, [pd_n32768]

    mov       dword [rsp], %2/4
.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_LUMA_VER_W4_4R

    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7

    psrad     m0, 2
    psrad     m1, 2
    psrad     m2, 2
    psrad     m3, 2

    packssdw  m0, m1
    packssdw  m2, m3

    movh      [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movh      [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [8 * r1 - 2 * 4]
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

;---------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_ps_%1x%2(pixel *src, intptr_t srcStride, int16_t *dst, intptr_t dstStride, int coeffIdx)
;---------------------------------------------------------------------------------------------------------------
    FILTER_VER_LUMA_PS 4, 4
    FILTER_VER_LUMA_PS 8, 8
    FILTER_VER_LUMA_PS 8, 4
    FILTER_VER_LUMA_PS 4, 8
    FILTER_VER_LUMA_PS 16, 16
    FILTER_VER_LUMA_PS 16, 8
    FILTER_VER_LUMA_PS 8, 16
    FILTER_VER_LUMA_PS 16, 12
    FILTER_VER_LUMA_PS 12, 16
    FILTER_VER_LUMA_PS 16, 4
    FILTER_VER_LUMA_PS 4, 16
    FILTER_VER_LUMA_PS 32, 32
    FILTER_VER_LUMA_PS 32, 16
    FILTER_VER_LUMA_PS 16, 32
    FILTER_VER_LUMA_PS 32, 24
    FILTER_VER_LUMA_PS 24, 32
    FILTER_VER_LUMA_PS 32, 8
    FILTER_VER_LUMA_PS 8, 32
    FILTER_VER_LUMA_PS 64, 64
    FILTER_VER_LUMA_PS 64, 32
    FILTER_VER_LUMA_PS 32, 64
    FILTER_VER_LUMA_PS 64, 48
    FILTER_VER_LUMA_PS 48, 64
    FILTER_VER_LUMA_PS 64, 16
    FILTER_VER_LUMA_PS 16, 64

;--------------------------------------------------------------------------------------------------------------
; void interp_8tap_vert_sp_%1x%2(int16_t *src, intptr_t srcStride, pixel *dst, intptr_t dstStride, int coeffIdx)
;--------------------------------------------------------------------------------------------------------------
%macro FILTER_VER_LUMA_SP 2
INIT_XMM sse4
cglobal interp_8tap_vert_sp_%1x%2, 5, 7, 8 ,0-gprsize

    add       r1d, r1d
    add       r3d, r3d
    lea       r5, [r1 + 2 * r1]
    sub       r0, r5
    shl       r4d, 6

%ifdef PIC
    lea       r5, [tab_LumaCoeffV]
    lea       r6, [r5 + r4]
%else
    lea       r6, [tab_LumaCoeffV + r4]
%endif

    mova      m7, [tab_c_524800]

    mov       dword [rsp], %2/4
.loopH:
    mov       r4d, (%1/4)
.loopW:
    PROCESS_LUMA_VER_W4_4R

    paddd     m0, m7
    paddd     m1, m7
    paddd     m2, m7
    paddd     m3, m7

    psrad     m0, 10
    psrad     m1, 10
    psrad     m2, 10
    psrad     m3, 10

    packssdw  m0, m1
    packssdw  m2, m3

    pxor      m1, m1
    CLIPW2    m0, m2, m1, [pw_pixel_max]

    movh      [r2], m0
    movhps    [r2 + r3], m0
    lea       r5, [r2 + 2 * r3]
    movh      [r5], m2
    movhps    [r5 + r3], m2

    lea       r5, [8 * r1 - 2 * 4]
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
    PROCESS_LUMA_VER_W4_4R

    psrad      m0, 6
    psrad      m1, 6
    packssdw   m0, m1
    movlps     [r2], m0
    movhps     [r2 + r3], m0

    psrad      m2, 6
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

;--------------------------------------------------------------------------------------------------
; void filterConvertPelToShort(pixel *src, intptr_t srcStride, int16_t *dst, int width, int height)
;--------------------------------------------------------------------------------------------------
INIT_XMM sse2
cglobal luma_p2s, 3, 7, 5

    add         r1, r1

    ; load width and height
    mov         r3d, r3m
    mov         r4d, r4m

    ; load constant
    mova        m4, [tab_c_n8192]

.loopH:

    xor         r5d, r5d
.loopW:
    lea         r6, [r0 + r5 * 2]

    movu        m0, [r6]
    psllw       m0, 4
    paddw       m0, m4

    movu        m1, [r6 + r1]
    psllw       m1, 4
    paddw       m1, m4

    movu        m2, [r6 + r1 * 2]
    psllw       m2, 4
    paddw       m2, m4

    lea         r6, [r6 + r1 * 2]
    movu        m3, [r6 + r1]
    psllw       m3, 4
    paddw       m3, m4

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
