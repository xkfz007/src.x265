;*****************************************************************************
;* Copyright (C) 2013 x265 project
;*
;* Authors: Praveen Kumar Tiwari <praveen@multicorewareinc.com>
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

SECTION .text

cextern pw_pixel_max

;-----------------------------------------------------------------------------
; void pixel_add_ps_4x4(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_4x4, 6, 6, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m1,     [pw_pixel_max]
    pxor    m0,     m0
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
    movh    m2,     [r2]
    movhps  m2,     [r2 + r4]
    movh    m3,     [r3]
    movhps  m3,     [r3 + r5]
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]
    movh    m4,     [r2]
    movhps  m4,     [r2 + r4]
    movh    m5,     [r3]
    movhps  m5,     [r3 + r5]

    paddw   m2,     m3
    paddw   m4,     m5
    CLIPW2  m2, m4, m0, m1

    movh    [r0],       m2
    movhps  [r0 + r1],  m2
    lea     r0,     [r0 + r1 * 2]
    movh    [r0],       m4
    movhps  [r0 + r1],  m4

    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_4x4, 6, 6, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    add         r5,         r5
    pmovzxbw    m0,         [r2]
    pmovzxbw    m2,         [r2 + r4]
    movh        m1,         [r3]
    movh        m3,         [r3 + r5]
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]
    pmovzxbw    m4,         [r2]
    pmovzxbw    m6,         [r2 + r4]
    movh        m5,         [r3]
    movh        m7,         [r3 + r5]

    paddw       m0,         m1
    paddw       m2,         m3
    paddw       m4,         m5
    paddw       m6,         m7
    packuswb    m0,         m0
    packuswb    m2,         m2
    packuswb    m4,         m4
    packuswb    m6,         m6

    movd        [r0],       m0
    movd        [r0 + r1],  m2
    lea         r0,         [r0 + r1 * 2]
    movd        [r0],       m4
    movd        [r0 + r1],  m6

    RET
%endif


;-----------------------------------------------------------------------------
; void pixel_add_ps_4x%2(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%macro PIXEL_ADD_PS_W4_H4 2
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_4x%2, 6, 7, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m1,     [pw_pixel_max]
    pxor    m0,     m0
    mov     r6d,    %2/4
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
.loop:
    movh    m2,     [r2]
    movhps  m2,     [r2 + r4]
    movh    m3,     [r3]
    movhps  m3,     [r3 + r5]
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]
    movh    m4,     [r2]
    movhps  m4,     [r2 + r4]
    movh    m5,     [r3]
    movhps  m5,     [r3 + r5]
    dec     r6d
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m2,     m3
    paddw   m4,     m5
    CLIPW2  m2, m4, m0, m1

    movh    [r0],       m2
    movhps  [r0 + r1],  m2
    lea     r0,     [r0 + r1 * 2]
    movh    [r0],       m4
    movhps  [r0 + r1],  m4
    lea     r0,     [r0 + r1 * 2]

    jnz     .loop
    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_4x%2, 6, 7, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    mov         r6d,        %2/4
    add         r5,         r5
.loop:
    pmovzxbw    m0,         [r2]
    pmovzxbw    m2,         [r2 + r4]
    movh        m1,         [r3]
    movh        m3,         [r3 + r5]
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]
    pmovzxbw    m4,         [r2]
    pmovzxbw    m6,         [r2 + r4]
    movh        m5,         [r3]
    movh        m7,         [r3 + r5]
    dec         r6d
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m1
    paddw       m2,         m3
    paddw       m4,         m5
    paddw       m6,         m7
    packuswb    m0,         m0
    packuswb    m2,         m2
    packuswb    m4,         m4
    packuswb    m6,         m6

    movd        [r0],       m0
    movd        [r0 + r1],  m2
    lea         r0,         [r0 + r1 * 2]
    movd        [r0],       m4
    movd        [r0 + r1],  m6
    lea         r0,         [r0 + r1 * 2]

    jnz         .loop
    RET
%endif
%endmacro

PIXEL_ADD_PS_W4_H4   4,  8


;-----------------------------------------------------------------------------
; void pixel_add_ps_8x%2(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%macro PIXEL_ADD_PS_W8_H4 2
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_8x%2, 6, 7, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m5,     [pw_pixel_max]
    pxor    m4,     m4
    mov     r6d,    %2/4
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
.loop:
    movu    m0,     [r2]
    movu    m2,     [r2 + r4]
    movu    m1,     [r3]
    movu    m3,     [r3 + r5]
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + r1],  m2

    movu    m0,     [r2]
    movu    m2,     [r2 + r4]
    movu    m1,     [r3]
    movu    m3,     [r3 + r5]
    dec     r6d
    lea     r0,     [r0 + r1 * 2]
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + r1],  m2
    lea     r0,     [r0 + r1 * 2]

    jnz     .loop
    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_8x%2, 6, 7, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    mov         r6d,        %2/4
    add         r5,         r5
.loop:
    pmovzxbw    m0,         [r2]
    pmovzxbw    m2,         [r2 + r4]
    movu        m1,         [r3]
    movu        m3,         [r3 + r5]
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]
    pmovzxbw    m4,         [r2]
    pmovzxbw    m6,         [r2 + r4]
    movu        m5,         [r3]
    movu        m7,         [r3 + r5]
    dec         r6d
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m1
    paddw       m2,         m3
    paddw       m4,         m5
    paddw       m6,         m7
    packuswb    m0,         m0
    packuswb    m2,         m2
    packuswb    m4,         m4
    packuswb    m6,         m6

    movh        [r0],       m0
    movh        [r0 + r1],  m2
    lea         r0,         [r0 + r1 * 2]
    movh        [r0],       m4
    movh        [r0 + r1],  m6
    lea         r0,         [r0 + r1 * 2]

    jnz         .loop
    RET
%endif
%endmacro

PIXEL_ADD_PS_W8_H4 8,  8
PIXEL_ADD_PS_W8_H4 8, 16


;-----------------------------------------------------------------------------
; void pixel_add_ps_16x%2(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%macro PIXEL_ADD_PS_W16_H4 2
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_16x%2, 6, 7, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m5,     [pw_pixel_max]
    pxor    m4,     m4
    mov     r6d,    %2/4
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
.loop:
    movu    m0,     [r2]
    movu    m2,     [r2 + 16]
    movu    m1,     [r3]
    movu    m3,     [r3 + 16]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + 16],  m2

    movu    m0,     [r2 + r4]
    movu    m2,     [r2 + r4 + 16]
    movu    m1,     [r3 + r5]
    movu    m3,     [r3 + r5 + 16]
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1],      m0
    movu    [r0 + r1 + 16], m2

    movu    m0,     [r2]
    movu    m2,     [r2 + 16]
    movu    m1,     [r3]
    movu    m3,     [r3 + 16]
    lea     r0,     [r0 + r1 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + 16],  m2

    movu    m0,     [r2 + r4]
    movu    m2,     [r2 + r4 + 16]
    movu    m1,     [r3 + r5]
    movu    m3,     [r3 + r5 + 16]
    dec     r6d
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1],      m0
    movu    [r0 + r1 + 16], m2
    lea     r0,     [r0 + r1 * 2]

    jnz     .loop
    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_16x%2, 6, 7, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    mov         r6d,        %2/4
    add         r5,         r5
.loop:
    pmovzxbw    m0,         [r2]
    pmovzxbw    m1,         [r2 + 8]
    pmovzxbw    m4,         [r2 + r4]
    pmovzxbw    m5,         [r2 + r4 + 8]
    movu        m2,         [r3]
    movu        m3,         [r3 + 16]
    movu        m6,         [r3 + r5]
    movu        m7,         [r3 + r5 + 16]
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m2
    paddw       m1,         m3
    paddw       m4,         m6
    paddw       m5,         m7
    packuswb    m0,         m1
    packuswb    m4,         m5

    movu        [r0],       m0
    movu        [r0 + r1],  m4

    pmovzxbw    m0,         [r2]
    pmovzxbw    m1,         [r2 + 8]
    pmovzxbw    m4,         [r2 + r4]
    pmovzxbw    m5,         [r2 + r4 + 8]
    movu        m2,         [r3]
    movu        m3,         [r3 + 16]
    movu        m6,         [r3 + r5]
    movu        m7,         [r3 + r5 + 16]
    dec         r6d
    lea         r0,         [r0 + r1 * 2]
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m2
    paddw       m1,         m3
    paddw       m4,         m6
    paddw       m5,         m7
    packuswb    m0,         m1
    packuswb    m4,         m5

    movu        [r0],       m0
    movu        [r0 + r1],  m4
    lea         r0,         [r0 + r1 * 2]

    jnz         .loop
    RET
%endif
%endmacro

PIXEL_ADD_PS_W16_H4 16, 16
PIXEL_ADD_PS_W16_H4 16, 32


;-----------------------------------------------------------------------------
; void pixel_add_ps_32x%2(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%macro PIXEL_ADD_PS_W32_H2 2
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_32x%2, 6, 7, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m5,     [pw_pixel_max]
    pxor    m4,     m4
    mov     r6d,    %2/2
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
.loop:
    movu    m0,     [r2]
    movu    m2,     [r2 + 16]
    movu    m1,     [r3]
    movu    m3,     [r3 + 16]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + 16],  m2

    movu    m0,     [r2 + 32]
    movu    m2,     [r2 + 48]
    movu    m1,     [r3 + 32]
    movu    m3,     [r3 + 48]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + 32],  m0
    movu    [r0 + 48],  m2

    movu    m0,     [r2 + r4]
    movu    m2,     [r2 + r4 + 16]
    movu    m1,     [r3 + r5]
    movu    m3,     [r3 + r5 + 16]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1],      m0
    movu    [r0 + r1 + 16], m2

    movu    m0,     [r2 + r4 + 32]
    movu    m2,     [r2 + r4 + 48]
    movu    m1,     [r3 + r5 + 32]
    movu    m3,     [r3 + r5 + 48]
    dec     r6d
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1 + 32], m0
    movu    [r0 + r1 + 48], m2
    lea     r0,     [r0 + r1 * 2]

    jnz     .loop
    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_32x%2, 6, 7, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    mov         r6d,        %2/2
    add         r5,         r5
.loop:
    pmovzxbw    m0,         [r2]
    pmovzxbw    m1,         [r2 + 8]
    pmovzxbw    m2,         [r2 + 16]
    pmovzxbw    m3,         [r2 + 24]
    movu        m4,         [r3]
    movu        m5,         [r3 + 16]
    movu        m6,         [r3 + 32]
    movu        m7,         [r3 + 48]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0],       m0
    movu        [r0 + 16],  m2

    pmovzxbw    m0,         [r2 + r4]
    pmovzxbw    m1,         [r2 + r4 + 8]
    pmovzxbw    m2,         [r2 + r4 + 16]
    pmovzxbw    m3,         [r2 + r4 + 24]
    movu        m4,         [r3 + r5]
    movu        m5,         [r3 + r5 + 16]
    movu        m6,         [r3 + r5 + 32]
    movu        m7,         [r3 + r5 + 48]
    dec         r6d
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0 + r1],      m0
    movu        [r0 + r1 + 16], m2
    lea         r0,         [r0 + r1 * 2]

    jnz         .loop
    RET
%endif
%endmacro

PIXEL_ADD_PS_W32_H2 32, 32
PIXEL_ADD_PS_W32_H2 32, 64


;-----------------------------------------------------------------------------
; void pixel_add_ps_64x%2(pixel *dest, intptr_t destride, pixel *src0, int16_t *scr1, intptr_t srcStride0, intptr_t srcStride1)
;-----------------------------------------------------------------------------
%macro PIXEL_ADD_PS_W64_H2 2
%if HIGH_BIT_DEPTH
INIT_XMM sse2
cglobal pixel_add_ps_64x%2, 6, 7, 6, dest, destride, src0, scr1, srcStride0, srcStride1
    mova    m5,     [pw_pixel_max]
    pxor    m4,     m4
    mov     r6d,    %2/2
    add     r4,     r4
    add     r5,     r5
    add     r1,     r1
.loop:
    movu    m0,     [r2]
    movu    m2,     [r2 + 16]
    movu    m1,     [r3]
    movu    m3,     [r3 + 16]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0],       m0
    movu    [r0 + 16],  m2

    movu    m0,     [r2 + 32]
    movu    m2,     [r2 + 48]
    movu    m1,     [r3 + 32]
    movu    m3,     [r3 + 48]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + 32],  m0
    movu    [r0 + 48],  m2

    movu    m0,     [r2 + 64]
    movu    m2,     [r2 + 80]
    movu    m1,     [r3 + 64]
    movu    m3,     [r3 + 80]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + 64],  m0
    movu    [r0 + 80],  m2

    movu    m0,     [r2 + 96]
    movu    m2,     [r2 + 112]
    movu    m1,     [r3 + 96]
    movu    m3,     [r3 + 112]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + 96],  m0
    movu    [r0 + 112], m2

    movu    m0,     [r2 + r4]
    movu    m2,     [r2 + r4 + 16]
    movu    m1,     [r3 + r5]
    movu    m3,     [r3 + r5 + 16]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1],      m0
    movu    [r0 + r1 + 16], m2

    movu    m0,     [r2 + r4 + 32]
    movu    m2,     [r2 + r4 + 48]
    movu    m1,     [r3 + r5 + 32]
    movu    m3,     [r3 + r5 + 48]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1 + 32], m0
    movu    [r0 + r1 + 48], m2

    movu    m0,     [r2 + r4 + 64]
    movu    m2,     [r2 + r4 + 80]
    movu    m1,     [r3 + r5 + 64]
    movu    m3,     [r3 + r5 + 80]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1 + 64], m0
    movu    [r0 + r1 + 80], m2

    movu    m0,     [r2 + r4 + 96]
    movu    m2,     [r2 + r4 + 112]
    movu    m1,     [r3 + r5 + 96]
    movu    m3,     [r3 + r5 + 112]
    dec     r6d
    lea     r2,     [r2 + r4 * 2]
    lea     r3,     [r3 + r5 * 2]

    paddw   m0,     m1
    paddw   m2,     m3
    CLIPW2  m0, m2, m4, m5

    movu    [r0 + r1 + 96],     m0
    movu    [r0 + r1 + 112],    m2
    lea     r0,     [r0 + r1 * 2]

    jnz     .loop
    RET
%else
INIT_XMM sse4
cglobal pixel_add_ps_64x%2, 6, 7, 8, dest, destride, src0, scr1, srcStride0, srcStride1
    mov         r6d,        %2/2
    add         r5,         r5
.loop:
    pmovzxbw    m0,         [r2]
    pmovzxbw    m1,         [r2 + 8]
    pmovzxbw    m2,         [r2 + 16]
    pmovzxbw    m3,         [r2 + 24]
    movu        m4,         [r3]
    movu        m5,         [r3 + 16]
    movu        m6,         [r3 + 32]
    movu        m7,         [r3 + 48]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0],       m0
    movu        [r0 + 16],  m2

    pmovzxbw    m0,         [r2 + 32]
    pmovzxbw    m1,         [r2 + 40]
    pmovzxbw    m2,         [r2 + 48]
    pmovzxbw    m3,         [r2 + 56]
    movu        m4,         [r3 + 64]
    movu        m5,         [r3 + 80]
    movu        m6,         [r3 + 96]
    movu        m7,         [r3 + 112]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0 + 32],  m0
    movu        [r0 + 48],  m2

    pmovzxbw    m0,         [r2 + r4]
    pmovzxbw    m1,         [r2 + r4 + 8]
    pmovzxbw    m2,         [r2 + r4 + 16]
    pmovzxbw    m3,         [r2 + r4 + 24]
    movu        m4,         [r3 + r5]
    movu        m5,         [r3 + r5 + 16]
    movu        m6,         [r3 + r5 + 32]
    movu        m7,         [r3 + r5 + 48]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0 + r1],      m0
    movu        [r0 + r1 + 16], m2

    pmovzxbw    m0,         [r2 + r4 + 32]
    pmovzxbw    m1,         [r2 + r4 + 40]
    pmovzxbw    m2,         [r2 + r4 + 48]
    pmovzxbw    m3,         [r2 + r4 + 56]
    movu        m4,         [r3 + r5 + 64]
    movu        m5,         [r3 + r5 + 80]
    movu        m6,         [r3 + r5 + 96]
    movu        m7,         [r3 + r5 + 112]
    dec         r6d
    lea         r2,         [r2 + r4 * 2]
    lea         r3,         [r3 + r5 * 2]

    paddw       m0,         m4
    paddw       m1,         m5
    paddw       m2,         m6
    paddw       m3,         m7
    packuswb    m0,         m1
    packuswb    m2,         m3

    movu        [r0 + r1 + 32], m0
    movu        [r0 + r1 + 48], m2
    lea         r0,         [r0 + r1 * 2]

    jnz         .loop
    RET
%endif
%endmacro

PIXEL_ADD_PS_W64_H2 64, 64
