#include "prlibm.h"
#include "log.h"
#include <math.h>

#define LN2HIGH 0.69314718055994528622676398299518041312694549560546875

double prlibm_log(float x) {
    float_x fix, fit;
    fix.f = x;
    int m = 0;

    if (x == 1.0) return 0.0;

    if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
        if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infinity
            fix.x = 0xFF800000;
            return fix.f;
        }

        if (fix.x > 0x7FFFFFFF) return (x - x) / 0; // Log(-val) = NaN
        if (fix.x >= 0x7F800000) return x + x;
        fix.f *= 8.388608e+06;
        m -= 23;
    }

    m += fix.x >> 23;
    m -= 127;
    fix.x &= 0x007FFFFF;
    fix.x |= 0x3F800000;

    fit.x = fix.x & 0x007F0000;
    int FIndex = fit.x >> 16;
    fit.x |= 0x3F800000;

    double f = fix.f - fit.f;
    f *= log_oneByF[FIndex];

    // Find the index of polynomial coefficients
    double_x dX;
    dX.d = f;
    double y;
#if 0
    p0: number of violated intervals: 4, total iterations=55
p0: VIOLATING INPUTS BELOW THRESHOLD:
p0: starting input is 0x0p+0
p0: violated_input is 0x1.67f6db6db6db7p-10, lb is 0x1.67b7a57460101p-10, ub is 0x1.67b7cfcd9f5p-10
p0: violated_input is 0x1.d6f7e432f7e44p-10, lb is 0x1.d68bb6f7c2101p-10, ub is 0x1.d68be43409dp-10
p0: violated_input is 0x1.fbd2361d2361ep-10, lb is 0x1.fb54746536101p-10, ub is 0x1.fb54885bc7dp-10
p0: violated_input is 0x1.57497b425ed09p-9, lb is 0x1.56d6991a29e81p-9, ub is 0x1.56d69cf495b7fp-9
Polynomial: y=0x1.fffffffffffd1p-1 x^(1) + -0x1.ffffffaa6e1ffp-2 x^(2) + 0x1.5553b37049aecp-2 x^(3) +
        -0x1.fb889a1fc999ep-3 x^(4)
#endif
    if (f < 0x1.5a8f8d28ac42fp-9) {
        if (f == 0x1.67f6db6db6db7p-10) {
            y = 0x1.67b7a57460101p-10;
        } else if (f == 0x1.d6f7e432f7e44p-10) {
            y = 0x1.d68bb6f7c2101p-10;
        } else if (f == 0x1.fbd2361d2361ep-10) {
            y = 0x1.fb54746536101p-10;
        } else if (f == 0x1.57497b425ed09p-9) {
            y = 0x1.56d6991a29e81p-9;
        } else {

            double coeffs[4] = {
                    0x1.fffffffffffd1p-1,
                    -0x1.ffffffaa6e1ffp-2,
                    0x1.5553b37049aecp-2,
                    -0x1.fb889a1fc999ep-3
            };

            double temp1 = fma(f, coeffs[2], coeffs[1]);
            double xsquare = f * f;
            double temp2 = fma(xsquare, coeffs[3], temp1);
            double temp3 = xsquare * temp2;
            y = fma(f, coeffs[0], temp3);
        }
    } else {
#if 0
        p1: number of violated intervals: 6, total iterations=22
p1: VIOLATING INPUTS BELOW THRESHOLD:
p1: starting input is 0x1.5a8f8d28ac42fp-9
p1: violated_input is 0x1.78d3dcb08d3ddp-9, lb is 0x1.78496eb2b2bc1p-9, ub is 0x1.784974d38489fp-9
p1: violated_input is 0x1.e8a1fd1b7af01p-9, lb is 0x1.e7b9668c47041p-9, ub is 0x1.e7b974fcfd0bfp-9
p1: violated_input is 0x1.3155555555555p-8, lb is 0x1.309fcf6432ca1p-8, ub is 0x1.309fcfa21506fp-8
p1: violated_input is 0x1.740a7ac29eb0ap-8, lb is 0x1.72fd098c2a761p-8, ub is 0x1.72fd2857b4e9fp-8
p1: violated_input is 0x1.a486d6f63aa04p-8, lb is 0x1.a32ee8debeb9bp-8, ub is 0x1.a32eea0a93b95p-8
p1: violated_input is 0x1.f9fcp-8, lb is 0x1.f80a850000001p-8, ub is 0x1.f80a85fffffffp-8
Polynomial: y=0x1.ffffffff5adf6p-1 x^(1) + -0x1.fffffb33f4724p-2 x^(2) + 0x1.554ed3225f9f3p-2 x^(3) + -0x1.f869f77cf2e1cp-3 x^(4)
#endif

        if (f == 0x1.78d3dcb08d3ddp-9) {
            y = 0x1.78496eb2b2bc1p-9;
        } else if (f == 0x1.e8a1fd1b7af01p-9) {
            y = 0x1.e7b9668c47041p-9;
        } else if (f == 0x1.3155555555555p-8) {
            y = 0x1.309fcf6432ca1p-8;
        } else if (f == 0x1.740a7ac29eb0ap-8) {
            y = 0x1.72fd098c2a761p-8;
        } else if (f == 0x1.a486d6f63aa04p-8) {
            y = 0x1.a32ee8debeb9bp-8;
        } else if (f == 0x1.f9fcp-8) {
            y = 0x1.f80a850000001p-8;
        } else {
            double coeffs[4] = {
                    0x1.ffffffff5adf6p-1,
                    -0x1.fffffb33f4724p-2,
                    0x1.554ed3225f9f3p-2,
                    -0x1.f869f77cf2e1cp-3
            };


            double temp1 = fma(f, coeffs[2], coeffs[1]);
            double xsquare = f * f;
            double temp2 = fma(xsquare, coeffs[3], temp1);
            double temp3 = xsquare * temp2;
            y = fma(f, coeffs[0], temp3);
        }
    }

    y += ln_lutHIGH[FIndex];
    y += m * LN2HIGH;

    return y;
}
