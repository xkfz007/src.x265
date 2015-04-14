/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com.
 *****************************************************************************/

#include "primitives.h"
#include "x265.h"

/* The #if logic here must match the file lists in CMakeLists.txt */
#if X265_ARCH_X86
#if defined(__INTEL_COMPILER)
#define HAVE_SSE3
#define HAVE_SSSE3
#define HAVE_SSE4
#define HAVE_AVX2
#elif defined(__GNUC__)
#if __clang__ || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 3)
#define HAVE_SSE3
#define HAVE_SSSE3
#define HAVE_SSE4
#endif
#if __clang__ || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 7)
#define HAVE_AVX2
#endif
#elif defined(_MSC_VER)
#define HAVE_SSE3
#define HAVE_SSSE3
#define HAVE_SSE4
#if _MSC_VER >= 1700 // VC11
#define HAVE_AVX2
#endif
#endif // compiler checks
#endif // if X265_ARCH_X86

namespace x265 {
// private x265 namespace

void Setup_Vec_DCTPrimitives_sse3(EncoderPrimitives&);
void Setup_Vec_DCTPrimitives_ssse3(EncoderPrimitives&);
void Setup_Vec_DCTPrimitives_sse41(EncoderPrimitives&);

/* Use primitives for the best available vector architecture */
void Setup_Instrinsic_Primitives(EncoderPrimitives &p, int cpuMask)
{
#ifdef HAVE_SSE3
    if (cpuMask & X265_CPU_SSE3)
    {
        Setup_Vec_DCTPrimitives_sse3(p);
    }
#endif
#ifdef HAVE_SSSE3
    if (cpuMask & X265_CPU_SSSE3)
    {
        Setup_Vec_DCTPrimitives_ssse3(p);
    }
#endif
#ifdef HAVE_SSE4
    if (cpuMask & X265_CPU_SSE4)
    {
        Setup_Vec_DCTPrimitives_sse41(p);
    }
#endif
    (void)p;
    (void)cpuMask;
}
}
