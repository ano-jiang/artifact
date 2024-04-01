#define __ELEM__ rlibm_exp2_estrin_fma
#include "rlibm.h"
#include "performance_test/LibTestHelper.h"

int additionallyIgnoreThisInput(float x) {
  float_x fx;
  fx.f = x;
  
  // Take care of special cases
  if ((fx.x & 0x7FFFFFFF) == 0) return 1;
  if (fx.x <= 876128826)  return 1;
  if (1124073472 <= fx.x && fx.x <= 3015223867)  return 1;
  if (fx.x >= 3272998913)  return 1;
  
  // Perform range reduction
  double xp = x * 64;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  
  int M = N1 / 64;
  int J = N2;
  double R = x - N * 0.015625;
  
  if (R == 0.0 && N2 == 0)  return 1;
  
  return 0;
}


int main(int argc, char** argv) {
  RunTest("exp2_estrin_fma.txt");
  return 0;
}
