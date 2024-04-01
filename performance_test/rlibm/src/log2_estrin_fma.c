#define __ELEM__ rlibm_log2_estrin_fma
#include "rlibm.h"
#include "performance_test/LibTestHelper.h"

int additionallyIgnoreThisInput(float x) {
  float_x fix, fit, spec;
  fix.f = x;
  int m = 0;
  
  if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
    if ((fix.x & 0x7FFFFFFF) == 0) return 1;
    if (fix.x > 0x7FFFFFFF) return 1;
    if (fix.x >= 0x7F800000) return 1;
    int exp;
    spec.f = frexpf(fix.f, &exp);
    if (spec.x == 0x3f000000) return 1;
  }
  
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  
  if (fix.x == 0) return 1;
  
  return 0;
}


int main(int argc, char** argv) {
  RunTest("log2_estrin_fma.txt");
  return 0;
}
