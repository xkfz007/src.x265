// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include <stdint.h>

#ifndef __PERF_H__
#define __PERF_H__

#ifdef __cplusplus
extern "C"
{
#endif

// Benchmark switch: 0 (disabled), 1 (program only), 2 (program and segment).
#define PERF_SWITCH 2

// True if the samples must be printed.
#define PERF_SEG_PRINT 0

// Segment count shift factor to avoid overflows in the sample array. Set this
// to 32 when the number of segment executions is unknown.
#define PERF_SHIFT_FACTOR 32

// Filtering mask bits used.
#define PERF_FILTER_MASK 0x0

// Number of cycles spent in the program without segment profiling.
#define PERF_BASE_CYCLE 0

// Number of microseconds/cycles at the beginning of the program/segment.
extern uint64_t perf_prog_usec_start;
extern uint64_t perf_prog_cycle_start;
extern uint64_t perf_seg_cycle_start;

// Filtering mask. A sample is only collected when the mask is zero.
extern uint32_t perf_filter_mask;

// Number of times the measured segment was executed.
extern uint64_t perf_seg_count;

// Total number of cycles spent in the segment.
extern uint64_t perf_seg_cycles;

// Array of segment samples (timestamp delta for each segment execution).
extern uint64_t perf_seg_samples[32768];

// Must be called on program entry and exit.
void PERF_PROG_ENTER();
void PERF_PROG_LEAVE();

// Set a bit of the filter mask on the measured control path.
static inline void PERF_FILTER_ENTER(uint32_t mask)
{
    #if PERF_SWITCH == 2
    perf_filter_mask &= ~mask;
    #endif
}

static inline void PERF_FILTER_LEAVE(uint32_t mask)
{
    #if PERF_SWITCH == 2
    perf_filter_mask |= mask;
    #endif
}

// Benchmark the measured segment.
static inline void PERF_BENCH_ENTER()
{
    #if PERF_SWITCH == 2
    __asm__ __volatile__(
        "rdtscp\n\t"
        "shl        $32, %%rdx\n\t"
        "add        %%rdx, %%rax\n\t"
        : "=a"(perf_seg_cycle_start)
        :
        : "rcx", "rdx"
        );
    #endif
}

static inline void PERF_BENCH_LEAVE()
{
    #if PERF_SWITCH == 2
    #define STR1(X) #X
    #define STR(X) STR1(X)
    __asm__ __volatile__(
        "rdtscp\n\t"                                        // Get the processor time stamp.
        "shl        $32, %%rdx\n\t"
        "add        %%rax, %%rdx\n\t"
        "sub        %4, %%rdx\n\t"                          // Subtract the start time stamp.
        "mov        %0, %%rcx\n\t"                          // Get the segment count.
        "shr        $"STR(PERF_SHIFT_FACTOR)", %%rcx\n\t"   // Shift the segment count.
        "lea        %1, %%rax\n\t"                          // Get the sample location.
        "movnti     %%rdx, (%%rax,%%rcx,8)\n\t"             // Store the sample (non-temporal).
        "mov        %2, %%ecx\n\t"                          // Load the filter mask.
        "xor        %%eax, %%eax\n\t"                       // Set a flag if the mask is null.
        "test       %%ecx, %%ecx\n\t"
        "sete       %%al\n\t"
        "add        %%rax, %0\n\t"                          // Increment the segment count.
        "sub        $1, %%rax\n\t"                          // Mask the delta with the flag.
        "not        %%rax\n\t"
        "and        %%rax, %%rdx\n\t"
        "add        %%rdx, %3\n\t"                          // Add the delta.
        :
        : "m"(perf_seg_count), "m"(*perf_seg_samples), "m"(perf_filter_mask), "m"(perf_seg_cycles), "m" (perf_seg_cycle_start)
        : "rax", "rcx", "rdx"
        );
    #undef STR1
    #undef STR
    #endif
}

#ifdef __cplusplus
}
#endif

#endif

