;*****************************************************************************
;* Copyright (C) 2013 x265 project
;*
;* Authors: Min Chen <chenm001@163.com>
;*          Praveen Kumar Tiwari <praveen@multicorewareinc.com>
;*          Nabajit Deka <nabajit@multicorewareinc.com>
;*          Dnyaneshwar Gorade <dnyaneshwar@multicorewareinc.com>
;*          Murugan Vairavel <murugan@multicorewareinc.com>
;*          Yuvaraj Venkatesh <yuvaraj@multicorewareinc.com>
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

SECTION_RODATA 32
pb_31:      times 32 db 31
pb_15:      times 32 db 15
pb_movemask_32:  times 32 db 0x00
                 times 32 db 0xFF

SECTION .text
cextern pb_1
cextern pb_128
cextern pb_2
cextern pw_2
cextern pb_movemask


;============================================================================================================
; void saoCuOrgE0(pixel * rec, int8_t * offsetEo, int lcuWidth, int8_t* signLeft, intptr_t stride)
;============================================================================================================
INIT_XMM sse4
cglobal saoCuOrgE0, 5, 5, 8, rec, offsetEo, lcuWidth, signLeft, stride

    mov         r4d, r4m
    mova        m4,  [pb_128]                ; m4 = [80]
    pxor        m5,  m5                      ; m5 = 0
    movu        m6,  [r1]                    ; m6 = offsetEo

    movzx       r1d, byte [r3]
    inc         r3
    neg         r1b
    movd        m0, r1d
    lea         r1, [r0 + r4]
    mov         r4d, r2d

.loop:
    movu        m7, [r0]                    ; m7 = rec[x]
    movu        m2, [r0 + 1]                ; m2 = rec[x+1]

    pxor        m1, m7, m4
    pxor        m3, m2, m4
    pcmpgtb     m2, m1, m3
    pcmpgtb     m3, m1
    pand        m2, [pb_1]
    por         m2, m3

    pslldq      m3, m2, 1
    por         m3, m0

    psignb      m3, m4                      ; m3 = signLeft
    pxor        m0, m0
    palignr     m0, m2, 15
    paddb       m2, m3
    paddb       m2, [pb_2]                  ; m2 = uiEdgeType
    pshufb      m3, m6, m2
    pmovzxbw    m2, m7                      ; rec
    punpckhbw   m7, m5
    pmovsxbw    m1, m3                      ; offsetEo
    punpckhbw   m3, m3
    psraw       m3, 8
    paddw       m2, m1
    paddw       m7, m3
    packuswb    m2, m7
    movu        [r0], m2

    add         r0q, 16
    sub         r2d, 16
    jnz        .loop

    movzx       r3d, byte [r3]
    neg         r3b
    movd        m0, r3d
.loopH:
    movu        m7, [r1]                    ; m7 = rec[x]
    movu        m2, [r1 + 1]                ; m2 = rec[x+1]

    pxor        m1, m7, m4
    pxor        m3, m2, m4
    pcmpgtb     m2, m1, m3
    pcmpgtb     m3, m1
    pand        m2, [pb_1]
    por         m2, m3

    pslldq      m3, m2, 1
    por         m3, m0

    psignb      m3, m4                      ; m3 = signLeft
    pxor        m0, m0
    palignr     m0, m2, 15
    paddb       m2, m3
    paddb       m2, [pb_2]                  ; m2 = uiEdgeType
    pshufb      m3, m6, m2
    pmovzxbw    m2, m7                      ; rec
    punpckhbw   m7, m5
    pmovsxbw    m1, m3                      ; offsetEo
    punpckhbw   m3, m3
    psraw       m3, 8
    paddw       m2, m1
    paddw       m7, m3
    packuswb    m2, m7
    movu        [r1], m2

    add         r1q, 16
    sub         r4d, 16
    jnz        .loopH
    RET

INIT_YMM avx2
cglobal saoCuOrgE0, 5, 5, 7, rec, offsetEo, lcuWidth, signLeft, stride

    mov                 r4d,        r4m
    vbroadcasti128      m4,         [pb_128]                   ; m4 = [80]
    vbroadcasti128      m6,         [r1]                       ; m6 = offsetEo
    movzx               r1d,        byte [r3]
    neg                 r1b
    movd                xm0,        r1d
    movzx               r1d,        byte [r3 + 1]
    neg                 r1b
    movd                xm1,        r1d
    vinserti128         m0,         m0,        xm1,           1

.loop:
    movu                xm5,        [r0]                       ; xm5 = rec[x]
    movu                xm2,        [r0 + 1]                   ; xm2 = rec[x + 1]
    vinserti128         m5,         m5,        [r0 + r4],     1
    vinserti128         m2,         m2,        [r0 + r4 + 1], 1

    pxor                m1,         m5,        m4
    pxor                m3,         m2,        m4
    pcmpgtb             m2,         m1,        m3
    pcmpgtb             m3,         m1
    pand                m2,         [pb_1]
    por                 m2,         m3

    pslldq              m3,         m2,        1
    por                 m3,         m0

    psignb              m3,         m4                         ; m3 = signLeft
    pxor                m0,         m0
    palignr             m0,         m2,        15
    paddb               m2,         m3
    paddb               m2,         [pb_2]                     ; m2 = uiEdgeType
    pshufb              m3,         m6,        m2
    pmovzxbw            m2,         xm5                        ; rec
    vextracti128        xm5,        m5,        1
    pmovzxbw            m5,         xm5
    pmovsxbw            m1,         xm3                        ; offsetEo
    vextracti128        xm3,        m3,        1
    pmovsxbw            m3,         xm3
    paddw               m2,         m1
    paddw               m5,         m3
    packuswb            m2,         m5
    vpermq              m2,         m2,        11011000b
    movu                [r0],       xm2
    vextracti128        [r0 + r4],  m2,        1

    add                 r0q,        16
    sub                 r2d,        16
    jnz                 .loop
    RET

;==================================================================================================
; void saoCuOrgE1(pixel *pRec, int8_t *m_iUpBuff1, int8_t *m_iOffsetEo, Int iStride, Int iLcuWidth)
;==================================================================================================
INIT_XMM sse4
cglobal saoCuOrgE1, 3, 5, 8, pRec, m_iUpBuff1, m_iOffsetEo, iStride, iLcuWidth
    mov         r3d, r3m
    mov         r4d, r4m
    pxor        m0,    m0                      ; m0 = 0
    mova        m6,    [pb_2]                  ; m6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    mova        m7,    [pb_128]
    shr         r4d,   4
.loop
    movu        m1,    [r0]                    ; m1 = pRec[x]
    movu        m2,    [r0 + r3]               ; m2 = pRec[x + iStride]

    pxor        m3,    m1,    m7
    pxor        m4,    m2,    m7
    pcmpgtb     m2,    m3,    m4
    pcmpgtb     m4,    m3
    pand        m2,    [pb_1]
    por         m2,    m4

    movu        m3,    [r1]                    ; m3 = m_iUpBuff1

    paddb       m3,    m2
    paddb       m3,    m6

    movu        m4,    [r2]                    ; m4 = m_iOffsetEo
    pshufb      m5,    m4,    m3

    psubb       m3,    m0,    m2
    movu        [r1],  m3

    pmovzxbw    m2,    m1
    punpckhbw   m1,    m0
    pmovsxbw    m3,    m5
    punpckhbw   m5,    m5
    psraw       m5,    8

    paddw       m2,    m3
    paddw       m1,    m5
    packuswb    m2,    m1
    movu        [r0],  m2

    add         r0,    16
    add         r1,    16
    dec         r4d
    jnz         .loop
    RET

INIT_YMM avx2
cglobal saoCuOrgE1, 3, 5, 8, pRec, m_iUpBuff1, m_iOffsetEo, iStride, iLcuWidth
    mov           r3d,    r3m
    mov           r4d,    r4m
    movu          xm0,    [r2]                    ; xm0 = m_iOffsetEo
    mova          xm6,    [pb_2]                  ; xm6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    mova          xm7,    [pb_128]
    shr           r4d,    4
.loop
    movu          xm1,    [r0]                    ; xm1 = pRec[x]
    movu          xm2,    [r0 + r3]               ; xm2 = pRec[x + iStride]

    pxor          xm3,    xm1,    xm7
    pxor          xm4,    xm2,    xm7
    pcmpgtb       xm2,    xm3,    xm4
    pcmpgtb       xm4,    xm3
    pand          xm2,    [pb_1]
    por           xm2,    xm4

    movu          xm3,    [r1]                    ; xm3 = m_iUpBuff1

    paddb         xm3,    xm2
    paddb         xm3,    xm6

    pshufb        xm5,    xm0,    xm3
    pxor          xm4,    xm4
    psubb         xm3,    xm4,    xm2
    movu          [r1],   xm3

    pmovzxbw      m2,     xm1
    pmovsxbw      m3,     xm5

    paddw         m2,     m3
    vextracti128  xm3,    m2,     1
    packuswb      xm2,    xm3
    movu          [r0],   xm2

    add           r0,     16
    add           r1,     16
    dec           r4d
    jnz           .loop
    RET

;========================================================================================================
; void saoCuOrgE1_2Rows(pixel *pRec, int8_t *m_iUpBuff1, int8_t *m_iOffsetEo, Int iStride, Int iLcuWidth)
;========================================================================================================
INIT_XMM sse4
cglobal saoCuOrgE1_2Rows, 3, 5, 8, pRec, m_iUpBuff1, m_iOffsetEo, iStride, iLcuWidth
    mov         r3d,        r3m
    mov         r4d,        r4m
    pxor        m0,         m0                      ; m0 = 0
    mova        m7,         [pb_128]
    shr         r4d,        4
.loop
    movu        m1,         [r0]                    ; m1 = pRec[x]
    movu        m2,         [r0 + r3]               ; m2 = pRec[x + iStride]

    pxor        m3,         m1,         m7
    pxor        m4,         m2,         m7
    pcmpgtb     m6,         m3,         m4
    pcmpgtb     m5,         m4,         m3
    pand        m6,         [pb_1]
    por         m6,         m5

    movu        m5,         [r0 + r3 * 2]
    pxor        m3,         m5,         m7
    pcmpgtb     m5,         m4,         m3
    pcmpgtb     m3,         m4
    pand        m5,         [pb_1]
    por         m5,         m3

    movu        m3,         [r1]                    ; m3 = m_iUpBuff1
    paddb       m3,         m6
    paddb       m3,         [pb_2]

    movu        m4,         [r2]                    ; m4 = m_iOffsetEo
    pshufb      m4,         m3

    psubb       m3,         m0,         m6
    movu        [r1],       m3

    pmovzxbw    m6,         m1
    punpckhbw   m1,         m0
    pmovsxbw    m3,         m4
    punpckhbw   m4,         m4
    psraw       m4,         8

    paddw       m6,         m3
    paddw       m1,         m4
    packuswb    m6,         m1
    movu        [r0],       m6

    movu        m3,         [r1]                    ; m3 = m_iUpBuff1
    paddb       m3,         m5
    paddb       m3,         [pb_2]

    movu        m4,         [r2]                    ; m4 = m_iOffsetEo
    pshufb      m4,         m3
    psubb       m3,         m0,         m5
    movu        [r1],       m3

    pmovzxbw    m5,         m2
    punpckhbw   m2,         m0
    pmovsxbw    m3,         m4
    punpckhbw   m4,         m4
    psraw       m4,         8

    paddw       m5,         m3
    paddw       m2,         m4
    packuswb    m5,         m2
    movu        [r0 + r3],  m5

    add         r0,         16
    add         r1,         16
    dec         r4d
    jnz         .loop
    RET

INIT_YMM avx2
cglobal saoCuOrgE1_2Rows, 3, 5, 7, pRec, m_iUpBuff1, m_iOffsetEo, iStride, iLcuWidth
    mov             r3d,        r3m
    mov             r4d,        r4m
    pxor            m0,         m0                           ; m0 = 0
    vbroadcasti128  m5,         [pb_128]
    vbroadcasti128  m6,         [r2]                         ; m6 = m_iOffsetEo
    shr             r4d,        4
.loop
    movu            xm1,        [r0]                         ; m1 = pRec[x]
    movu            xm2,        [r0 + r3]                    ; m2 = pRec[x + iStride]
    vinserti128     m1,         m1,       xm2,            1
    vinserti128     m2,         m2,       [r0 + r3 * 2],  1

    pxor            m3,         m1,       m5
    pxor            m4,         m2,       m5
    pcmpgtb         m2,         m3,       m4
    pcmpgtb         m4,         m3
    pand            m2,         [pb_1]
    por             m2,         m4

    movu            xm3,        [r1]                         ; xm3 = m_iUpBuff
    psubb           m4,         m0,       m2
    vinserti128     m3,         m3,       xm4,            1
    paddb           m3,         m2
    paddb           m3,         [pb_2]
    pshufb          m2,         m6,       m3
    vextracti128    [r1],       m4,       1

    pmovzxbw        m4,         xm1
    vextracti128    xm3,        m1,       1
    pmovzxbw        m3,         xm3
    pmovsxbw        m1,         xm2
    vextracti128    xm2,        m2,       1
    pmovsxbw        m2,         xm2

    paddw           m4,         m1
    paddw           m3,         m2
    packuswb        m4,         m3
    vpermq          m4,         m4,       11011000b
    movu            [r0],       xm4
    vextracti128    [r0 + r3],  m4,       1

    add             r0,         16
    add             r1,         16
    dec             r4d
    jnz             .loop
    RET

;======================================================================================================================================================
; void saoCuOrgE2(pixel * rec, int8_t * bufft, int8_t * buff1, int8_t * offsetEo, int lcuWidth, intptr_t stride)
;======================================================================================================================================================
INIT_XMM sse4
cglobal saoCuOrgE2, 5, 6, 8, rec, bufft, buff1, offsetEo, lcuWidth
    mov         r4d,   r4m
    mov         r5d,   r5m
    pxor        m0,    m0                      ; m0 = 0
    mova        m6,    [pb_2]                  ; m6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    mova        m7,    [pb_128]
    inc         r1
    movh        m5,    [r0 + r4]
    movhps      m5,    [r1 + r4]

.loop
    movu        m1,    [r0]                    ; m1 = rec[x]
    movu        m2,    [r0 + r5 + 1]           ; m2 = rec[x + stride + 1]
    pxor        m3,    m1,    m7
    pxor        m4,    m2,    m7
    pcmpgtb     m2,    m3,    m4
    pcmpgtb     m4,    m3
    pand        m2,    [pb_1]
    por         m2,    m4
    movu        m3,    [r2]                    ; m3 = buff1

    paddb       m3,    m2
    paddb       m3,    m6                      ; m3 = edgeType

    movu        m4,    [r3]                    ; m4 = offsetEo
    pshufb      m4,    m3

    psubb       m3,    m0,    m2
    movu        [r1],  m3

    pmovzxbw    m2,    m1
    punpckhbw   m1,    m0
    pmovsxbw    m3,    m4
    punpckhbw   m4,    m4
    psraw       m4,    8

    paddw       m2,    m3
    paddw       m1,    m4
    packuswb    m2,    m1
    movu        [r0],  m2

    add         r0,    16
    add         r1,    16
    add         r2,    16
    sub         r4,    16
    jg          .loop

    movh        [r0 + r4], m5
    movhps      [r1 + r4], m5
    RET

INIT_YMM avx2
cglobal saoCuOrgE2, 5, 6, 7, rec, bufft, buff1, offsetEo, lcuWidth
    mov            r4d,   r4m
    mov            r5d,   r5m
    pxor           xm0,   xm0                     ; xm0 = 0
    mova           xm5,   [pb_128]
    inc            r1
    movq           xm6,   [r0 + r4]
    movhps         xm6,   [r1 + r4]

    movu           xm1,   [r0]                    ; xm1 = rec[x]
    movu           xm2,   [r0 + r5 + 1]           ; xm2 = rec[x + stride + 1]
    pxor           xm3,   xm1,   xm5
    pxor           xm4,   xm2,   xm5
    pcmpgtb        xm2,   xm3,   xm4
    pcmpgtb        xm4,   xm3
    pand           xm2,   [pb_1]
    por            xm2,   xm4
    movu           xm3,   [r2]                    ; xm3 = buff1

    paddb          xm3,   xm2
    paddb          xm3,   [pb_2]                  ; xm3 = edgeType

    movu           xm4,   [r3]                    ; xm4 = offsetEo
    pshufb         xm4,   xm3

    psubb          xm3,   xm0,   xm2
    movu           [r1],  xm3

    pmovzxbw       m2,    xm1
    pmovsxbw       m3,    xm4

    paddw          m2,    m3
    vextracti128   xm3,   m2,    1
    packuswb       xm2,   xm3
    movu           [r0],  xm2

    movq           [r0 + r4], xm6
    movhps         [r1 + r4], xm6
    RET

INIT_YMM avx2
cglobal saoCuOrgE2_32, 5, 6, 8, rec, bufft, buff1, offsetEo, lcuWidth
    mov             r4d,   r4m
    mov             r5d,   r5m
    pxor            m0,    m0                      ; m0 = 0
    vbroadcasti128  m7,    [pb_128]
    vbroadcasti128  m5,    [r3]                    ; m5 = offsetEo
    inc             r1
    movq            xm6,   [r0 + r4]
    movhps          xm6,   [r1 + r4]

.loop:
    movu            m1,    [r0]                    ; m1 = rec[x]
    movu            m2,    [r0 + r5 + 1]           ; m2 = rec[x + stride + 1]
    pxor            m3,    m1,    m7
    pxor            m4,    m2,    m7
    pcmpgtb         m2,    m3,    m4
    pcmpgtb         m4,    m3
    pand            m2,    [pb_1]
    por             m2,    m4
    movu            m3,    [r2]                    ; m3 = buff1

    paddb           m3,    m2
    paddb           m3,    [pb_2]                  ; m3 = edgeType

    pshufb          m4,    m5,    m3

    psubb           m3,    m0,    m2
    movu            [r1],  m3

    pmovzxbw        m2,    xm1
    vextracti128    xm1,   m1,    1
    pmovzxbw        m1,    xm1
    pmovsxbw        m3,    xm4
    vextracti128    xm4,   m4,    1
    pmovsxbw        m4,    xm4

    paddw           m2,    m3
    paddw           m1,    m4
    packuswb        m2,    m1
    vpermq          m2,    m2,    11011000b
    movu            [r0],  m2

    add             r0,    32
    add             r1,    32
    add             r2,    32
    sub             r4,    32
    jg              .loop

    movq            [r0 + r4], xm6
    movhps          [r1 + r4], xm6
    RET

;=======================================================================================================
;void saoCuOrgE3(pixel *rec, int8_t *upBuff1, int8_t *m_offsetEo, intptr_t stride, int startX, int endX)
;=======================================================================================================
INIT_XMM sse4
cglobal saoCuOrgE3, 3,6,8
    mov             r3d, r3m
    mov             r4d, r4m
    mov             r5d, r5m

    ; save latest 2 pixels for case startX=1 or left_endX=15
    movh            m7, [r0 + r5]
    movhps          m7, [r1 + r5 - 1]

    ; move to startX+1
    inc             r4d
    add             r0, r4
    add             r1, r4
    sub             r5d, r4d
    pxor            m0, m0                      ; m0 = 0
    movu            m6, [pb_2]                  ; m6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

.loop:
    movu            m1, [r0]                    ; m1 = pRec[x]
    movu            m2, [r0 + r3]               ; m2 = pRec[x + iStride]

    psubusb         m3, m2, m1
    psubusb         m4, m1, m2
    pcmpeqb         m3, m0
    pcmpeqb         m4, m0
    pcmpeqb         m2, m1

    pabsb           m3, m3
    por             m4, m3
    pandn           m2, m4                      ; m2 = iSignDown

    movu            m3, [r1]                    ; m3 = m_iUpBuff1

    paddb           m3, m2
    paddb           m3, m6                      ; m3 = uiEdgeType

    movu            m4, [r2]                    ; m4 = m_iOffsetEo
    pshufb          m5, m4, m3

    psubb           m3, m0, m2
    movu            [r1 - 1], m3

    pmovzxbw        m2, m1
    punpckhbw       m1, m0
    pmovsxbw        m3, m5
    punpckhbw       m5, m5
    psraw           m5, 8

    paddw           m2, m3
    paddw           m1, m5
    packuswb        m2, m1
    movu            [r0], m2

    add             r0, 16
    add             r1, 16

    sub             r5, 16
    jg             .loop

    ; restore last pixels (up to 2)
    movh            [r0 + r5], m7
    movhps          [r1 + r5 - 1], m7
    RET

INIT_YMM avx2
cglobal saoCuOrgE3, 3, 6, 8
    mov             r3d,  r3m
    mov             r4d,  r4m
    mov             r5d,  r5m

    ; save latest 2 pixels for case startX=1 or left_endX=15
    movq            xm7,  [r0 + r5]
    movhps          xm7,  [r1 + r5 - 1]

    ; move to startX+1
    inc             r4d
    add             r0,   r4
    add             r1,   r4
    sub             r5d,  r4d
    pxor            xm0,  xm0                     ; xm0 = 0
    mova            xm6,  [pb_2]                  ; xm6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    movu            xm5,  [r2]                    ; xm5 = m_iOffsetEo

.loop:
    movu            xm1,  [r0]                    ; xm1 = pRec[x]
    movu            xm2,  [r0 + r3]               ; xm2 = pRec[x + iStride]

    psubusb         xm3,  xm2,  xm1
    psubusb         xm4,  xm1,  xm2
    pcmpeqb         xm3,  xm0
    pcmpeqb         xm4,  xm0
    pcmpeqb         xm2,  xm1

    pabsb           xm3,  xm3
    por             xm4,  xm3
    pandn           xm2,  xm4                     ; xm2 = iSignDown

    movu            xm3,  [r1]                    ; xm3 = m_iUpBuff1

    paddb           xm3,  xm2
    paddb           xm3,  xm6                     ; xm3 = uiEdgeType

    pshufb          xm4,  xm5,  xm3

    psubb           xm3,  xm0,  xm2
    movu            [r1 - 1],   xm3

    pmovzxbw        m2,   xm1
    pmovsxbw        m3,   xm4

    paddw           m2,   m3
    vextracti128    xm3,  m2,   1
    packuswb        xm2,  xm3
    movu            [r0], xm2

    add             r0,   16
    add             r1,   16

    sub             r5,   16
    jg             .loop

    ; restore last pixels (up to 2)
    movq            [r0 + r5],     xm7
    movhps          [r1 + r5 - 1], xm7
    RET

INIT_YMM avx2
cglobal saoCuOrgE3_32, 3, 6, 8
    mov             r3d,  r3m
    mov             r4d,  r4m
    mov             r5d,  r5m

    ; save latest 2 pixels for case startX=1 or left_endX=15
    movq            xm7,  [r0 + r5]
    movhps          xm7,  [r1 + r5 - 1]

    ; move to startX+1
    inc             r4d
    add             r0,   r4
    add             r1,   r4
    sub             r5d,  r4d
    pxor            m0,   m0                      ; m0 = 0
    mova            m6,   [pb_2]                  ; m6 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    vbroadcasti128  m5,   [r2]                    ; m5 = m_iOffsetEo

.loop:
    movu            m1,   [r0]                    ; m1 = pRec[x]
    movu            m2,   [r0 + r3]               ; m2 = pRec[x + iStride]

    psubusb         m3,   m2,   m1
    psubusb         m4,   m1,   m2
    pcmpeqb         m3,   m0
    pcmpeqb         m4,   m0
    pcmpeqb         m2,   m1

    pabsb           m3,   m3
    por             m4,   m3
    pandn           m2,   m4                      ; m2 = iSignDown

    movu            m3,   [r1]                    ; m3 = m_iUpBuff1

    paddb           m3,   m2
    paddb           m3,   m6                      ; m3 = uiEdgeType

    pshufb          m4,   m5,   m3

    psubb           m3,   m0,   m2
    movu            [r1 - 1],   m3

    pmovzxbw        m2,   xm1
    vextracti128    xm1,  m1,   1
    pmovzxbw        m1,   xm1
    pmovsxbw        m3,   xm4
    vextracti128    xm4,  m4,   1
    pmovsxbw        m4,   xm4

    paddw           m2,   m3
    paddw           m1,   m4
    packuswb        m2,   m1
    vpermq          m2,   m2,   11011000b
    movu            [r0], m2

    add             r0,   32
    add             r1,   32
    sub             r5,   32
    jg             .loop

    ; restore last pixels (up to 2)
    movq            [r0 + r5],     xm7
    movhps          [r1 + r5 - 1], xm7
    RET

;=====================================================================================
; void saoCuOrgB0(pixel* rec, const pixel* offset, int lcuWidth, int lcuHeight, int stride)
;=====================================================================================
INIT_XMM sse4
cglobal saoCuOrgB0, 4, 7, 8

    mov         r3d, r3m
    mov         r4d, r4m

    shr         r2d, 4
    movu        m3, [r1 + 0]      ; offset[0-15]
    movu        m4, [r1 + 16]     ; offset[16-31]
    pxor        m7, m7            ; m7 =[0]
.loopH
    mov         r5d, r2d
    xor         r6,  r6

.loopW
    movu        m2, [r0 + r6]     ; m0 = [rec]
    psrlw       m1, m2, 3
    pand        m1, [pb_31]       ; m1 = [index]
    pcmpgtb     m0, m1, [pb_15]   ; m2 = [mask]

    pshufb      m6, m3, m1
    pshufb      m5, m4, m1

    pblendvb    m6, m6, m5, m0    ; NOTE: don't use 3 parameters style, x264 macro have some bug!

    pmovzxbw    m1, m2            ; rec
    punpckhbw   m2, m7

    pmovsxbw    m0, m6            ; offset
    punpckhbw   m6, m6
    psraw       m6, 8

    paddw       m1, m0
    paddw       m2, m6
    packuswb    m1, m2

    movu        [r0 + r6], m1
    add         r6d, 16
    dec         r5d
    jnz         .loopW

    lea         r0, [r0 + r4]

    dec         r3d
    jnz         .loopH
    RET

INIT_YMM avx2
cglobal saoCuOrgB0, 4, 7, 8

    mov             r3d,        r3m
    mov             r4d,        r4m
    mova            m7,         [pb_31]
    vbroadcasti128  m3,         [r1 + 0]            ; offset[0-15]
    vbroadcasti128  m4,         [r1 + 16]           ; offset[16-31]
    lea             r6,         [r4 * 2]
    sub             r6d,        r2d
    shr             r2d,        4
    mov             r1d,        r3d
    shr             r3d,        1
.loopH
    mov             r5d,        r2d
.loopW
    movu            xm2,        [r0]                ; m2 = [rec]
    vinserti128     m2,         m2,  [r0 + r4],  1
    psrlw           m1,         m2,  3
    pand            m1,         m7                  ; m1 = [index]
    pcmpgtb         m0,         m1,  [pb_15]        ; m0 = [mask]

    pshufb          m6,         m3,  m1
    pshufb          m5,         m4,  m1

    pblendvb        m6,         m6,  m5,  m0        ; NOTE: don't use 3 parameters style, x264 macro have some bug!

    pmovzxbw        m1,         xm2                 ; rec
    vextracti128    xm2,        m2,  1
    pmovzxbw        m2,         xm2
    pmovsxbw        m0,         xm6                 ; offset
    vextracti128    xm6,        m6,  1
    pmovsxbw        m6,         xm6

    paddw           m1,         m0
    paddw           m2,         m6
    packuswb        m1,         m2
    vpermq          m1,         m1,  11011000b

    movu            [r0],       xm1
    vextracti128    [r0 + r4],  m1,  1
    add             r0,         16
    dec             r5d
    jnz             .loopW

    add             r0,         r6
    dec             r3d
    jnz             .loopH
    test            r1b,        1
    jz              .end
    mov             r5d,        r2d
.loopW1
    movu            xm2,        [r0]                ; m2 = [rec]
    psrlw           xm1,        xm2, 3
    pand            xm1,        xm7                 ; m1 = [index]
    pcmpgtb         xm0,        xm1, [pb_15]        ; m0 = [mask]

    pshufb          xm6,        xm3, xm1
    pshufb          xm5,        xm4, xm1

    pblendvb        xm6,        xm6, xm5, xm0       ; NOTE: don't use 3 parameters style, x264 macro have some bug!

    pmovzxbw        m1,         xm2                 ; rec
    pmovsxbw        m0,         xm6                 ; offset

    paddw           m1,         m0
    vextracti128    xm0,        m1,  1
    packuswb        xm1,        xm0

    movu            [r0],       xm1
    add             r0,         16
    dec             r5d
    jnz             .loopW1
.end
    RET

;============================================================================================================
; void calSign(int8_t *dst, const Pixel *src1, const Pixel *src2, const int width)
;============================================================================================================
INIT_XMM sse4
cglobal calSign, 4,5,6
    mova        m0,     [pb_128]
    mova        m1,     [pb_1]

    sub         r1,     r0
    sub         r2,     r0

    mov         r4d,    r3d
    shr         r3d,    4
    jz         .next
.loop:
    movu        m2,     [r0 + r1]            ; m2 = pRec[x]
    movu        m3,     [r0 + r2]            ; m3 = pTmpU[x]
    pxor        m4,     m2,     m0
    pxor        m3,     m0
    pcmpgtb     m5,     m4,     m3
    pcmpgtb     m3,     m4
    pand        m5,     m1
    por         m5,     m3
    movu        [r0],   m5

    add         r0,     16
    dec         r3d
    jnz        .loop

    ; process partial
.next:
    and         r4d, 15
    jz         .end

    movu        m2,     [r0 + r1]            ; m2 = pRec[x]
    movu        m3,     [r0 + r2]            ; m3 = pTmpU[x]
    pxor        m4,     m2,     m0
    pxor        m3,     m0
    pcmpgtb     m5,     m4,     m3
    pcmpgtb     m3,     m4
    pand        m5,     m1
    por         m5,     m3

    lea         r3,     [pb_movemask + 16]
    sub         r3,     r4
    movu        xmm0,   [r3]
    movu        m3,     [r0]
    pblendvb    m5,     m5,     m3,     xmm0
    movu        [r0],   m5

.end:
    RET

INIT_YMM avx2
cglobal calSign, 4, 5, 6
    vbroadcasti128  m0,     [pb_128]
    mova            m1,     [pb_1]

    sub             r1,     r0
    sub             r2,     r0

    mov             r4d,    r3d
    shr             r3d,    5
    jz              .next
.loop:
    movu            m2,     [r0 + r1]            ; m2 = pRec[x]
    movu            m3,     [r0 + r2]            ; m3 = pTmpU[x]
    pxor            m4,     m2,     m0
    pxor            m3,     m0
    pcmpgtb         m5,     m4,     m3
    pcmpgtb         m3,     m4
    pand            m5,     m1
    por             m5,     m3
    movu            [r0],   m5

    add             r0,     mmsize
    dec             r3d
    jnz             .loop

    ; process partial
.next:
    and             r4d,    31
    jz              .end

    movu            m2,     [r0 + r1]            ; m2 = pRec[x]
    movu            m3,     [r0 + r2]            ; m3 = pTmpU[x]
    pxor            m4,     m2,     m0
    pxor            m3,     m0
    pcmpgtb         m5,     m4,     m3
    pcmpgtb         m3,     m4
    pand            m5,     m1
    por             m5,     m3

    lea             r3,     [pb_movemask_32 + 32]
    sub             r3,     r4
    movu            m0,     [r3]
    movu            m3,     [r0]
    pblendvb        m5,     m5,     m3,     m0
    movu            [r0],   m5

.end:
    RET
