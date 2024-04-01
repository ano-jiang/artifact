#include "prlibm.h"
#include "log2.h"
#include <math.h>

double prlibm_log2(float x) {
    float_x fix, fit, spec;
    fix.f = x;
    int m = 0;

    if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
        if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
            fix.x = 0xFF800000;
            return fix.f;
        }

        if (fix.x > 0x7FFFFFFF) { // Log(-val) = NaN
            return (x - x) / 0;

        }

        if (fix.x >= 0x7F800000) {
            return x + x;
        }

        // Special case when we have denormal input and exact result
        int exp;
        spec.f = frexpf(fix.f, &exp);
        if (spec.x == 0x3f000000) return (double) (exp - 1);

        fix.f *= 8.388608e+06;
        m -= 23;
    }

    m += fix.x >> 23;
    m -= 127;
    fix.x &= 0x007FFFFF;

    if (fix.x == 0) {
        return (double) m;
    }

    fix.x |= 0x3F800000;

    fit.x = fix.x & 0x007F0000;
    int FIndex = fit.x >> 16;
    fit.x |= 0x3F800000;

    double f = fix.f - fit.f;
    f *= __log_oneByF[FIndex];

    double y;

    // FMA coefficients
    double coeffs[] = {
//            0x1.71547652b83f7p+0,
//            -0x1.7154765134475p-1,
//            0x1.ec709d3094c1bp-2,
//            -0x1.71591d35e2e7ap-2,
//            0x1.284a8df1d9342p-2
//            0x1.71547652b8284p+0,
//            -0x1.715476514b4ebp-1,
//            0x1.ec709da508be9p-2,
//            -0x1.71590d8d35cc5p-2,
//            0x1.282e1b8bbc1d4p-2
            0x1.71547652bcde3p+0,
            -0x1.7154769679dd8p-1,
            0x1.ec7198ec61291p-2,
            -0x1.72033bee9c2d6p-2,
            0x1.4f082e01903edp-2
//            1.4426950408890186761112772728665731847286224365234375000000000000000000000000e+00,
//            -7.2134752026808135472180083525017835199832916259765625000000000000000000000000e-01,
//            4.8089833840384849095173080968379508703947067260742187500000000000000000000000e-01,
//            -3.6069150582691344997243731995695270597934722900390625000000000000000000000000e-01,
//            2.8934690273881724653648461753618903458118438720703125000000000000000000000000e-01

    };

#if 0
Polynomial: y=0x1.71547652bcde3p+0 x^(1) + -0x1.7154769679dd8p-1 x^(2) + 0x1.ec7198ec61291p-2 x^(3) +
            -0x1.72033bee9c2d6p-2 x^(4) + 0x1.4f082e01903edp-2 x^(5)
Polynomial: y=0x1.71547652b8284p+0 x^(1) + -0x1.715476514b4ebp-1 x^(2) + 0x1.ec709da508be9p-2 x^(3) +
            -0x1.71590d8d35cc5p-2 x^(4) + 0x1.282e1b8bbc1d4p-2 x^(5)
  Polynomial: y=0x1.71547652b83f7p+0 x^(1) + -0x1.7154765134475p-1 x^(2) + 0x1.ec709d3094c1bp-2 x^(3) +
            -0x1.71591d35e2e7ap-2 x^(4) + 0x1.284a8df1d9342p-2 x^(5)
      p0: number of violated intervals: 0, total iterations=41
  p0: VIOLATING INPUTS BELOW THRESHOLD:
  p0: starting input is -0x1.fc05f417d05f4p-8
  Polynomial: y=0x1.71547652b83f7p+0 x^(1) + -0x1.7154765134475p-1 x^(2) + 0x1.ec709d3094c1bp-2 x^(3) +
          -0x1.71591d35e2e7ap-2 x^(4) + 0x1.284a8df1d9342p-2 x^(5)
    Polynomial: y=0x1.71547652b83f7p+0 x^(1) + -0x1.7154765134475p-1 x^(2) + 0x1.ec709d3094c1bp-2 x^(3) +
          -0x1.71591d35e2e7ap-2 x^(4) + 0x1.284a8df1d9342p-2 x^(5)
#endif

    double xsquare = f * f;
    double temp1 = fma(f, coeffs[4], coeffs[3]);
    double temp2 = fma(f, coeffs[2], coeffs[1]);
    double temp3 = fma(xsquare, temp1, temp2);
    double temp4 = xsquare * temp3;
    y = fma(f, coeffs[0], temp4);
    y += __log2_lut[FIndex];
    y += m;

    return y;
}
