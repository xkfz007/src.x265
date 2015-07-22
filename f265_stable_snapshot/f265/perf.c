// Copyright (c) 2014, VANTRIX CORPORATION. All rights reserved. See LICENSE.txt
// for the full license text.

#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <sched.h>
#include <errno.h>
#include "perf.h"

uint64_t perf_prog_usec_start = 0;
uint64_t perf_prog_cycle_start = 0;
uint64_t perf_seg_cycle_start = 0;
uint32_t perf_filter_mask = PERF_FILTER_MASK;
uint64_t perf_seg_count = 0;
uint64_t perf_seg_cycles = 0;
// This padding is required to avoid flushing the cache line containing the
// variables above when the samples are written with movnti.
uint8_t perf_padding[64];
uint64_t perf_seg_samples[32768];
uint8_t perf_padding2[64];

// Median computation. Return 0 if there are no elements in the array. This
// Quickselect routine is based on the algorithm described in "Numerical recipes
// in C", Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
// 0-521-43108-5. This code by Nicolas Devillard - 1998. Public domain.
uint64_t get_median(uint64_t arr[], int n)
{
    #define ELEM_SWAP(a,b) { uint64_t t=(a);(a)=(b);(b)=t; }
    int low, high;
    int median;
    int middle, ll, hh;

    if (!n) return 0; /* No element. */
    low = 0; high = n-1; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median];

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]);
            return arr[median];
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]);
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]);
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]);

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]);

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]);
        do hh--; while (arr[hh]  > arr[low]);

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]);
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]);

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
    #undef ELEM_SWAP
}

uint64_t perf_get_date()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (uint64_t)tv.tv_sec*1000000 + (uint64_t)tv.tv_usec;
}

void PERF_PROG_ENTER()
{
    #if PERF_SWITCH
    // Set the processor affinity.
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(1, &mask);
    if (sched_setaffinity(0, sizeof(mask), &mask))
    {
        perror("Cannot set processor affinity");
        exit(1);
    }

    // Set the real-time priority, if possible.
    struct sched_param sp;
    memset(&sp, 0, sizeof(sp));
    sp.sched_priority = 99;
    if (sched_setscheduler(getpid(), SCHED_FIFO, &sp))
    {
        // Commented out, since this happens frequently.
        #if 0
        perror("Cannot set real-time priority");
        #endif
    }

    // Zero the samples to avoid page faults.
    memset(perf_seg_samples, 0, sizeof(perf_seg_samples));

    // Get the start time.
    perf_prog_usec_start = perf_get_date();
    __asm__ __volatile__(
        "rdtscp\n\t"
        "shl        $32, %%rdx\n\t"
        "add        %%rdx, %%rax\n\t"
        : "=a"(perf_prog_cycle_start)
        :
        : "rcx", "rdx"
        );
   #endif
}

void PERF_PROG_LEAVE()
{
    #if PERF_SWITCH
    // Get the end time.
    uint64_t perf_prog_usec_end = perf_get_date();
    uint64_t perf_prog_cycle_end;
    __asm__ __volatile__(
        "sfence\n\t" // Flush non-temporal writes.
        "rdtscp\n\t"
        "shl        $32, %%rdx\n\t"
        "add        %%rdx, %%rax\n\t"
        : "=a"(perf_prog_cycle_end)
        :
        : "rcx", "rdx"
        );

    uint64_t perf_prog_usec_delta = perf_prog_usec_end  - perf_prog_usec_start;
    double perf_prog_sec_delta = ((float)perf_prog_usec_delta )/1000000.0;
    uint64_t perl_prog_cycle_delta = perf_prog_cycle_end - perf_prog_cycle_start;
    double perf_prog_overhead = (float)perl_prog_cycle_delta / (float)PERF_BASE_CYCLE;
    double perf_seg_ratio = (float)perf_seg_cycles/(float)perl_prog_cycle_delta;
    double perf_seg_sec_delta = perf_prog_sec_delta * perf_seg_ratio;
    uint32_t shift_factor = (perf_seg_count < 32768) ? 0 : (1+(uint32_t)log2f(perf_seg_count>>15));
    uint32_t nb_samples = (perf_seg_count>>PERF_SHIFT_FACTOR);
    uint64_t seg_min = (uint64_t)(-1llu), seg_max = 0, seg_sum = 0, seg_mean = 0;

    // Compute the sample metrics. The median computation is done in place so
    // the samples must be copied.
    for (uint32_t i = 0; i < nb_samples; i++)
    {
        uint64_t s = perf_seg_samples[i];
        if (s < seg_min) seg_min = s;
        if (s > seg_max) seg_max = s;
        seg_sum += s;
    }
    if (!nb_samples) seg_min = 0;
    else seg_mean = seg_sum / nb_samples;
    uint64_t *mem = (uint64_t*)malloc(sizeof(perf_seg_samples));
    if (mem == NULL) { printf("OOM.\n"); exit(1); }
    memcpy(mem, perf_seg_samples, sizeof(perf_seg_samples));
    uint64_t seg_med = get_median(mem, nb_samples);
    free(mem);
    double seg_sum_ratio = perf_seg_cycles ? (float)(seg_sum<<PERF_SHIFT_FACTOR) / (float)perf_seg_cycles : 0;

    #define ulli unsigned long long int
    #if PERF_SEG_PRINT
    printf("Samples:\n");
    for (uint32_t i = 0; i < nb_samples; i++)
        printf("%llu ", (ulli)perf_seg_samples[i]);
    printf("\n");
    #endif
    printf("Program time:     %.2f seconds, %llu cycles, overhead %.2f.\n",
           perf_prog_sec_delta, (ulli)perl_prog_cycle_delta, perf_prog_overhead);
    printf("Segment time:     %.2f seconds, %llu cycles, ratio %.2f%%.\n",
           perf_seg_sec_delta, (ulli)perf_seg_cycles, perf_seg_ratio*100.0);
    printf("Segment count:    exec %llu, samples %d, shift factor %d (use %d).\n",
           (ulli)perf_seg_count, nb_samples, PERF_SHIFT_FACTOR, shift_factor);
    printf("Segment cycles:   min %llu, median %llu, mean %llu, max %llu, correlation %.2f.\n",
           (ulli)seg_min, (ulli)seg_med, (ulli)seg_mean, (ulli)seg_max, seg_sum_ratio);
    #undef ulli
    #endif
}

