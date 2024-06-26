#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <x86intrin.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include<omp.h>

typedef union {
    double d;
    unsigned long x;
} doubleX;

typedef union {
    float f;
    unsigned int x;
} floatX;

void RunTest(char* FileName) {

    struct stat s;
    unsigned long count = 0;
    unsigned long wrongResult = 0;
    unsigned long totalTime = 0;
    unsigned int t1, t2, dummy;

    float x, res;
    floatX xbase;

    int fd = open(FileName, O_RDONLY);

    // Get Size of oracle file
    int status = fstat(fd, &s);
    size_t file_size = s.st_size;

    // MMap oracle file
    float* oracle = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (oracle == MAP_FAILED) {
        printf("mmap %s failed: %s\n", FileName, strerror(errno));
        exit(0);
    }
#pragma omp parallel for
    for (count = 0x0; count < 0x100000000; count++) {
        xbase.x = count;
        x = xbase.f;
        res = __PRLIBM_FAST_ELEM__(x);
        float oracleResult = oracle[count];

        // Check whether they're both NaN. Since NaN != NaN, we need this
        // check
        if (oracleResult != oracleResult && res != res) continue;

        // Now check if the two values are exactly the same
        if (oracle[count] != res) {
#pragma omp atomic
            wrongResult++;
        }
    }

    // Un-mmap oracle file
    munmap(oracle, file_size);
    close(fd);
    printf("glibc double wrong result: %lu\n", wrongResult);
}
