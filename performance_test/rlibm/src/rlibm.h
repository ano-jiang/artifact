#ifndef PRLIBM_RLIBM_H
#define PRLIBM_RLIBM_H

#include <stdint.h>

typedef union {
    double d;
    uint64_t x;
} double_x;

typedef union {
    float f;
    uint32_t x;
} float_x;
double rlibm_log2_estrin_fma(float);


double rlibm_log10_estrin_fma(float);


double rlibm_log_estrin_fma(float);


double rlibm_exp2_estrin_fma(float);


double rlibm_exp10_estrin_fma(float);


double rlibm_exp_estrin_fma(float);

#endif //PRLIBM_RLIBM_H
