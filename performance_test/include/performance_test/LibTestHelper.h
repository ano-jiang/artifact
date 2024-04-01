#ifdef __INTEL_COMPILER
#include "mathimf.h"
#else
#define _GNU_SOURCE
#include "math.h"
#endif

#if defined(__x86_64__) || defined(__ppc64__) || defined(__i386__)
#include <x86intrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

int additionallyIgnoreThisInput(float x);

void RunTest(char *FileName) {
    unsigned long time_total = 0, count = 0;
    unsigned long time_t1, time_t2;
    unsigned long accum;
    double res;
    unsigned int dummy;

    for (count = 0x0; count < 0x100000000; count += 0x1) {
        struct timeval start, end;

        float_x xbase;
        xbase.x = count;
        float x = xbase.f;
#ifdef IGNORE_SPECIAL_INPUT
        if (additionallyIgnoreThisInput(x)) continue;
#endif

#ifdef __GNUC__
    #if defined(__x86_64__) || defined(__ppc64__) || defined(__i386__)
        // x64平台的代码
        do {
            time_t1 = __rdtscp(&dummy);
            res = __ELEM__(x);
            time_t2 = __rdtscp(&dummy);
        } while (time_t1 >= time_t2);
        if (res == 0.0) accum += 1;
        time_total += (time_t2 - time_t1);

    #else
        gettimeofday(&start, NULL);
        res = __ELEM__(x);
        gettimeofday(&end, NULL);
        time_total+=((end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec);
        if (res == 0.0) accum += 1;
    #endif
#endif

    }
    printf("%lu\n", time_total);
}
