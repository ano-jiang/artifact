#define _GNU_SOURCE

#include <stdio.h>

#include "prlibm.h"
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>


// MAXVAL = 0111 1111 0111 1111 1111 1111 1111 1111 11
#define MAXVAL 3.40282361850336062550457001444955389952e+38
// MAXm1VAL = 0111 1111 0111 1111 1111 1111 1111 1111 10
#define MAXm1VAL 3.40282356779733661637539395458142568448e+38

double RoundToFloat34RO(double val) {
    // Take care of special cases
    if (val == 0.0) return 0.0;
    if (val == 1.0 / 0.0) return val;
    if (val == -1.0 / 0.0) return val;
    if (val > 0) {
        if (val < ldexp(1.0, -150)) return ldexp(1.0, -151);
        if (val > MAXm1VAL) return MAXVAL;

    } else {
        if (val > ldexp(-1.0, -150)) return ldexp(-1.0, -151);
        if (val < -MAXm1VAL) return -MAXVAL;
    }

    // At this point we have AT LEAST 2 precision bits. Guaranteed 1+ mantissa bit.
    double_x dx;
    dx.d = val;

    // Get exp value of val
    int exp;
    double frac = frexp(val, &exp);
    exp--;

    // with exp value, figure out how many precision/mantissa bits i get to have
    // # precision = min(26, max(2, 152 + exp))
    // # mantissa  = min(25, max(1, 151 + exp))
    long numMantissa = 151L + exp;
    if (numMantissa < 1) numMantissa = 1;
    if (numMantissa > 25) numMantissa = 25;

    unsigned long sticky = dx.x & ((1LU << (52LU - (unsigned long) numMantissa)) - 1LU);
    dx.x -= sticky;
    if (sticky != 0LU) dx.x |= 1LU << (52LU - (unsigned long) numMantissa);

    return dx.d;
}

void RunTest(char *FileName) {

    struct stat s;
    unsigned long prlibm_fast_wrongResult = 0;

    float x;
    double prlibm_fast_res;


    int fd = open(FileName, O_RDONLY);

    // Get Size of oracle file
    int status = fstat(fd, &s);
    size_t file_size = s.st_size;

    // MMap oracle file
    double *oracle = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (oracle == MAP_FAILED) {
        printf("mmap %s failed: %s\n", FileName, strerror(errno));
        exit(0);
    }
#pragma omp parallel for
    for (unsigned long count = 0x0; count < 0x100000000; count++) {
        float_x xbase;
        xbase.x = count;
        x = xbase.f;


        prlibm_fast_res = RoundToFloat34RO(__PRLIBM_FAST_ELEM__(x));
        double oracleResult = oracle[count];

        // Now check if the two values are exactly the same
        if ((oracle[count] != prlibm_fast_res) &&
            (oracleResult == oracleResult || prlibm_fast_res == prlibm_fast_res)) {
            printf("prlibm result is %a\n", prlibm_fast_res);
            printf("oracle_result is %a\n", oracleResult);
            printf("the input is %x, float value is %a\n", xbase.x, xbase.f);
#pragma omp atomic
            prlibm_fast_wrongResult++;
        }

        if (count % 10000000 == 0) {
            printf("count = %lu (%lu)\r", count, prlibm_fast_wrongResult);
            fflush(stdout);
        }
    }

    // Un-mmap oracle file
    munmap(oracle, file_size);
    close(fd);

    printf("Wrong results: \n");
    printf("prlibm-latest   wrong result: %lu\n", prlibm_fast_wrongResult);
}
