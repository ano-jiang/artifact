#include <math.h>
#include "prlibm.h"
#include "exp2.h"

// Estrin+FMA version of the piecewise polynomial for exp

double prlibm_exp(float x) {
    float_x fx;
    fx.f = x;

    // Take care of special cases
    if ((fx.x & 0x7FFFFFFF) == 0) return 1.0;

    if (fx.x <= 872415231) {
        if (fx.x <= 864026623) return 1.0000000298023223876953125;
        return 1.0000000894069671630859375;
    }

    if (1118925336 <= fx.x && fx.x <= 3011510272) {
        if (fx.x < 0x80000000) {
            if (fx.x < 0x7F800000) return 3.40282361850336062550457001444955389952e+38;
            if (fx.x == 0x7F800000) return 1.0 / 0.0;
            return 0.0 / 0.0;
        }

        // negative counterpart
        if (fx.x <= 3003121664) return 0.99999998509883880615234375;

        return 0.99999995529651641845703125;
    }

    if (fx.x >= 3268407733) {
        if (fx.x == 0xFF800000) return 0.0;
        if (fx.x < 0xFF800000) return ldexp(1.0, -151);
        return 0.0 / 0.0;
    }

    // Perform range reduction
    double xp = x * 92.332482616893656768297660164535045623779296875;
    int N = (int) xp;
    int N2 = N % 64;
    if (N2 < 0) N2 += 64;
    int N1 = N - N2;

    int M = N1 / 64;
    int J = N2;
    double R = x - N *
                   0.01083042469624914509729318723429969395510852336883544921875;

    double_x dX;
    dX.d = R;
    double y;
#if 0
    Trying to generate a polynomial  with 5 terms
p0: number of violated intervals: 3, total iterations=277
p0: VIOLATING INPUTS BELOW THRESHOLD:
p0: starting input is -0x1.62e42f490ap-7
p0: violated_input is -0x1.f925ff514p-9, lb is 0x1.fe07d2e08aadp-1, ub is 0x1.fe07e367cfc55p-1
p0: violated_input is -0x1.e08d3ep-9, lb is 0x1.fe20540000001p-1, ub is 0x1.fe2054fffffffp-1
p0: violated_input is -0x1.c57f87ab4d8p-9, lb is 0x1.fe3b49143b1b9p-1, ub is 0x1.fe3b4a7683067p-1
Polynomial: y=0x1.0000000000004p+0 x^(0) + 0x1.fffffffffa5adp-1 x^(1) + 0x1.ffffffaabedfp-2 x^(2) +
        0x1.5554722d18f4bp-3 x^(3) + 0x1.53a5b8cdcb223p-5 x^(4)
    p0: number of violated intervals: 1, total iterations=995
p0: VIOLATING INPUTS BELOW THRESHOLD:
p0: starting input is -0x1.62e42f490ap-7
p0: violated_input is -0x1.f925ff514p-9, lb is 0x1.fe07d2e08aadp-1, ub is 0x1.fe07e367cfc55p-1
Polynomial: y=0x1.0000000000009p+0 x^(0) + 0x1.00000000078fp+0 x^(1) + 0x1.0000002934215p-1 x^(2) + 0x1.5556024e1aecbp-3 x^(3)
        + 0x1.568896be3321bp-5 x^(4) + 0x1.d21fa89ec68p-7 x^(5)
#endif
    if (R < -0x1.9e76dcp-24) {
        if (R == -0x1.f925ff514p-9) {
            y = 0x1.fe07d2e08aadp-1;
        } else if(R == -0x1.e08d3ep-9){
            y = 0x1.fe20540000001p-1;
        }else if(R ==-0x1.c57f87ab4d8p-9){
            y=0x1.fe3b49143b1b9p-1;
        } else {
            double coeffs[] = {
//                    0x1.0000000000009p+0,
//                    0x1.00000000078fp+0,
//                    0x1.0000002934215p-1,
//                    0x1.5556024e1aecbp-3,
//                    0x1.568896be3321bp-5,
//                    0x1.d21fa89ec68p-7
                    0x1.0000000000004p+0,
                    0x1.fffffffffa5adp-1,
                    0x1.ffffffaabedfp-2,
                    0x1.5554722d18f4bp-3,
                    0x1.53a5b8cdcb223p-5
            };

//            double xsquare = R * R;
//            double temp1 = fma(R, coeffs[1], coeffs[0]);
//            double temp2 = fma(R, coeffs[5], coeffs[4]);
//            double temp3 = fma(R, coeffs[3], coeffs[2]);
//
//            double temp4 = fma(xsquare, temp2, temp3);
//            y = fma(xsquare, temp4, temp1);
            double xsquare = R * R;
            double temp1 = fma(R, coeffs[1], coeffs[0]);
            double temp2 = fma(R, coeffs[3], coeffs[2]);
            double temp3 = fma(xsquare, coeffs[4], temp2);
            y = fma(xsquare, temp3, temp1);
        }
    } else {
#if 0
        p1: number of violated intervals: 0, total iterations=57
p1: VIOLATING INPUTS BELOW THRESHOLD:
p1: starting input is -0x1.9e76dcp-24
Polynomial: y=0x1.ffffffffffffcp-1 x^(0) + 0x1.00000000021d8p+0 x^(1) + 0x1.ffffffed036cp-2 x^(2) + 0x1.5555730aea46ep-3
        x^(3) + 0x1.55289a7207dffp-5 x^(4) + 0x1.2d4e385563f7cp-7 x^(5)
#endif
        double coeffs[] = {
                0x1.ffffffffffffcp-1,
                0x1.00000000021d8p+0,
                0x1.ffffffed036cp-2,
                0x1.5555730aea46ep-3,
                0x1.55289a7207dffp-5,
                0x1.2d4e385563f7cp-7
        };
        double xsquare = R * R;
        double temp1 = fma(R, coeffs[1], coeffs[0]);
        double temp2 = fma(R, coeffs[5], coeffs[4]);
        double temp3 = fma(R, coeffs[3], coeffs[2]);

        double temp4 = fma(xsquare, temp2, temp3);
        y = fma(xsquare, temp4, temp1);

    }

    // Perform output compensation
    y *= ldexp(exp2JBy64[J], M);
    return y;
}
