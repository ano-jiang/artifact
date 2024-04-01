#ifndef __PRLIBM_LIBM_H__
#define __PRLIBM_LIBM_H__

#include <stdint.h>

typedef union {
  double d;
  uint64_t x;
} double_x;

typedef union {
  float f;
  uint32_t x;
} float_x;

double prlibm_log2(float);


double prlibm_log10(float);


double prlibm_log(float);


double prlibm_exp2(float);


double prlibm_exp10(float);


double prlibm_exp(float);


#endif
