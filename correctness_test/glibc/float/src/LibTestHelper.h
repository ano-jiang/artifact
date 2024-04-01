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


void RunTest(char *FileName) {

    struct stat s;
    unsigned long prlibm_fast_wrongResult = 0;

    float x,prlibm_fast_res;
//    double ;


    int fd = open(FileName, O_RDONLY);

    // Get Size of oracle file
    int status = fstat(fd, &s);
    size_t file_size = s.st_size;

    // MMap oracle file
    float *oracle = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (oracle == MAP_FAILED) {
        printf("mmap %s failed: %s\n", FileName, strerror(errno));
        exit(0);
    }
#pragma omp parallel for
    for (unsigned long count = 0x0; count < 0x100000000; count++) {
        float_x xbase;
        xbase.x = count;
        x = xbase.f;


        prlibm_fast_res = __PRLIBM_FAST_ELEM__(x);
        float oracleResult = oracle[count];

        // Now check if the two values are exactly the same
        if ((oracle[count] != prlibm_fast_res) &&
            (oracleResult == oracleResult || prlibm_fast_res == prlibm_fast_res)) {
//            printf("glibc result is %a\n", prlibm_fast_res);
//            printf("oracle_result is %a\n", oracleResult);
//            printf("the input is %x, float value is %a\n", xbase.x, xbase.f);
#pragma omp atomic
            prlibm_fast_wrongResult++;
        }

//        if (count % 10000000 == 0) {
//            printf("count = %lu (%lu)\r", count, prlibm_fast_wrongResult);
//            fflush(stdout);
//        }
    }

    // Un-mmap oracle file
    munmap(oracle, file_size);
    close(fd);

    printf("Wrong results: \n");
    printf("glibc float wrong result: %lu\n", prlibm_fast_wrongResult);
}
