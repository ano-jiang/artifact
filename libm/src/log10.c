#include "prlibm.h"
#include "log10.h"
#include <math.h>

#define LOG102HIGH 0.30102999566398114250631579125183634459972381591796875
#define LOG102LOW  5.27074231034726570126349709198449199648263806413338306011695522101945243775844573974609375e-17

double prlibm_log10(float x) {
//    float_x fix, fit;
//    fix.f = x;
//    int m = 0;
//
//    if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
//        if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
//            fix.x = 0xFF800000;
//            return fix.f;
//        }
//
//        if (fix.x > 0x7FFFFFFF) { // Log(-val) = NaN
//            return (x - x) / 0;
//
//        }
//
//        if (fix.x >= 0x7F800000) {
//            return x + x;
//        }
//
//        fix.f *= 8.388608e+06;
//        m -= 23;
//    }
//
//    switch (fix.x) {
//        case 0x3f800000 :
//            return 0.0;
//        case 0x41200000 :
//            return 1.0;
//        case 0x42c80000 :
//            return 2.0;
//        case 0x447a0000 :
//            return 3.0;
//        case 0x461c4000 :
//            return 4.0;
//        case 0x47c35000 :
//            return 5.0;
//        case 0x49742400 :
//            return 6.0;
//        case 0x4b189680 :
//            return 7.0;
//        case 0x4cbebc20 :
//            return 8.0;
//        case 0x4e6e6b28 :
//            return 9.0;
//        case 0x501502f9 :
//            return 10.0;
//    }
//
//    m += fix.x >> 23;
//    m -= 127;
//    fix.x &= 0x007FFFFF;
//    fix.x |= 0x3F800000;
//
//    fit.x = fix.x & 0x007F0000;
//    int FIndex = fit.x >> 16;
//    fit.x |= 0x3F800000;
//
//    double f = fix.f - fit.f;
//    f *= log_oneByF[FIndex];
//
//    double y;
//
//    if (f < 0x1.5a8f8d28ac42fp-9) {
//        if (f < 0x1.5a91111111111p-10) {
//            // 1st sub-domain
//#if 0
//Polynomial: y=0x1.bcb7b152566b9p-2 x^(1) + -0x1.bcb7ae3f54b59p-3 x^(2) + 0x1.286a6ae2f8166p-3 x^(3) +
//        -0x1.88240065ec912p-4 x^(4)
//
//Polynomial: y=0x1.bcb7b152566b9p-2 x^(1) + -0x1.bcb7ae3f54b59p-3 x^(2) + 0x1.286a6ae2f8166p-3 x^(3) +
//        -0x1.88240065ec912p-4 x^(4)
//            p0: number of violated intervals: 0, total iterations=19
//p0: VIOLATING INPUTS BELOW THRESHOLD:
//p0: starting input is 0x0p+0
//Polynomial: y=0x1.bcb7b152566b9p-2 x^(1) + -0x1.bcb7ae3f54b59p-3 x^(2) + 0x1.286a6ae2f8166p-3 x^(3) +
//        -0x1.88240065ec912p-4 x^(4)
//#endif
//            double coeffs[] = {
//                    0x1.bcb7b152566b9p-2,
//                    -0x1.bcb7ae3f54b59p-3,
//                    0x1.286a6ae2f8166p-3,
//                    -0x1.88240065ec912p-4
//            };
//
//            double temp1 = fma(f, coeffs[2], coeffs[1]);
//            double xsquare = f * f;
//            double temp2 = fma(xsquare, coeffs[3], temp1);
//            double temp3 = xsquare * temp2;
//            y = fma(f, coeffs[0], temp3);
//
//        } else {
//            // 2nd sub-doomain
//#if 0
//            number of violated intervals: 1, total iterations=27
//VIOLATING INPUTS BELOW THRESHOLD:
//starting input is 0x1.5a91111111111p-10
//violated_input is 0x1.1fddb0d3224f3p-9, lb is 0x1.f386956531508p-11, ub is 0x1.f386ad8d4abbdp-11
//Polynomial: y=0x1.bcb7b1516bf32p-2 x^(1) + -0x1.bcb7a4ab5429ep-3 x^(2) + 0x1.28608cd5626aap-3 x^(3)
//        + -0x1.994d5e748761p-4 x^(4)
//
//            p1: number of violated intervals: 1, total iterations=23
//p1: VIOLATING INPUTS BELOW THRESHOLD:
//p1: starting input is 0x1.5a91111111111p-10
//p1: violated_input is 0x1.1fddb0d3224f3p-9, lb is 0x1.f386956531508p-11, ub is 0x1.f386ad8d4abbdp-11
//Polynomial: y=0x1.bcb7b1516bf32p-2 x^(1) + -0x1.bcb7a4ab5429ep-3 x^(2) + 0x1.28608cd5626aap-3 x^(3) +
//        -0x1.994d5e748761p-4 x^(4)
//            p1: number of violated intervals: 1, total iterations=55
//p1: VIOLATING INPUTS BELOW THRESHOLD:
//p1: starting input is 0x1.5a91111111111p-10
//p1: violated_input is 0x1.1fddb0d3224f3p-9, lb is 0x1.f386956531508p-11, ub is 0x1.f386ad8d4abbdp-11
//Polynomial: y=0x1.bcb7b1516bf32p-2 x^(1) + -0x1.bcb7a4ab5429ep-3 x^(2) + 0x1.28608cd5626aap-3 x^(3) +
//        -0x1.994d5e748761p-4 x^(4)
//#endif
//            if (f == 0x1.1fddb0d3224f3p-9) {
//                y = 0x1.f386956531508p-11;
//            } else {
//
//                double coeffs[] = {
//                        0x1.bcb7b1516bf32p-2,
//                        -0x1.bcb7a4ab5429ep-3,
//                        0x1.28608cd5626aap-3,
//                        -0x1.994d5e748761p-4
//                };
//
//                double temp1 = fma(f, coeffs[2], coeffs[1]);
//                double xsquare = f * f;
//                double temp2 = fma(xsquare, coeffs[3], temp1);
//                double temp3 = xsquare * temp2;
//                y = fma(f, coeffs[0], temp3);
//
//            }
//        }
//    } else {
//        if (f < 0x1.03652e52e52e6p-8) {
//            // 3rd sub-domain
//#if 0
//            total iterations=816
//VIOLATING INPUTS BELOW THRESHOLD:
//starting input is 0x1.5a8f8d28ac42fp-9
//violated_input is 0x1.8ff099fc267fp-9, lb is 0x1.5adab7bb93889p-10, ub is 0x1.5adabbecb58d9p-10
//violated_input is 0x1.bde34a2b10bf6p-9, lb is 0x1.82a2d0cd6ee1p-10, ub is 0x1.82a2d501e19b2p-10
//Polynomial: y=0x1.bcb7b1514e819p-2 x^(1) + -0x1.bcb7a939b55c8p-3 x^(2) + 0x1.2870443390d58p-3 x^(3) +
//        -0x1.b2a2fa1d5534dp-4 x^(4)
//
//            p2: number of violated intervals: 3, total iterations=152
//p2: VIOLATING INPUTS BELOW THRESHOLD:
//p2: starting input is 0x1.5a8f8d28ac42fp-9
//p2: violated_input is 0x1.fbf5f5f5f5f5fp-9, lb is 0x1.b85b3cbbd9d74p-10, ub is 0x1.b85b3cbfd9c72p-10
//p2: violated_input is 0x1.bde34a2b10bf6p-9, lb is 0x1.82a2d0cd6ee1p-10, ub is 0x1.82a2d501e19b2p-10
//p2: violated_input is 0x1.8ff099fc267fp-9, lb is 0x1.5adab7bb93889p-10, ub is 0x1.5adabbecb58d9p-10
//Polynomial: y=0x1.bcb7b14f2e0aap-2 x^(1) + -0x1.bcb79a2b5a2aep-3 x^(2) + 0x1.285eb11241f0ep-3 x^(3) +
//        -0x1.a518081cdab1fp-4 x^(4)
//            p2: number of violated intervals: 3, total iterations=332
//p2: VIOLATING INPUTS BELOW THRESHOLD:
//p2: starting input is 0x1.5a8f8d28ac42fp-9
//p2: violated_input is 0x1.8ff099fc267fp-9, lb is 0x1.5adab7bb93889p-10, ub is 0x1.5adabbecb58d9p-10
//p2: violated_input is 0x1.bde34a2b10bf6p-9, lb is 0x1.82a2d0cd6ee1p-10, ub is 0x1.82a2d501e19b2p-10
//p2: violated_input is 0x1.fbf5f5f5f5f5fp-9, lb is 0x1.b85b3cbbd9d74p-10, ub is 0x1.b85b3cbfd9c72p-10
//Polynomial: y=0x1.bcb7b14f2e0aap-2 x^(1) + -0x1.bcb79a2b5a2aep-3 x^(2) + 0x1.285eb11241f0ep-3 x^(3) +
//        -0x1.a518081cdab1fp-4 x^(4)
//#endif
//            if (f == 0x1.8ff099fc267fp-9) {
//                y = 0x1.5adab7bb93889p-10;
//            } else if (f == 0x1.bde34a2b10bf6p-9) {
//                y = 0x1.82a2d0cd6ee1p-10;
//            } else {
//                double coeffs[] = {
//                        0x1.bcb7b1514e819p-2,
//                        -0x1.bcb7a939b55c8p-3,
//                        0x1.2870443390d58p-3,
//                        -0x1.b2a2fa1d5534dp-4
//                };
//
//                double temp1 = fma(f, coeffs[2], coeffs[1]);
//                double xsquare = f * f;
//                double temp2 = fma(xsquare, coeffs[3], temp1);
//                double temp3 = xsquare * temp2;
//                y = fma(f, coeffs[0], temp3);
//            }
//        } else {
//            // 4th sub-domain
//#if 0
//            number of violated intervals: 2, total iterations=92
//VIOLATING INPUTS BELOW THRESHOLD:
//starting input is 0x1.03652e52e52e6p-8
//violated_input is 0x1.4212f684bda13p-8, lb is 0x1.1710980c1c9efp-9, ub is 0x1.171098627003bp-9
//violated_input is 0x1.95e4d9364d937p-8, lb is 0x1.5f77b475327f9p-9, ub is 0x1.5f77bc53b21bbp-9
//Polynomial: y=0x1.bcb7b151dd50ap-2 x^(1) + -0x1.bcb7ad0365d8cp-3 x^(2) + 0x1.28749ea54c589p-3 x^(3) +
//        -0x1.b60309bf9d835p-4 x^(4)
//
//            p3: number of violated intervals: 2, total iterations=101
//p3: VIOLATING INPUTS BELOW THRESHOLD:
//p3: starting input is 0x1.03652e52e52e6p-8
//p3: violated_input is 0x1.4212f684bda13p-8, lb is 0x1.1710980c1c9efp-9, ub is 0x1.171098627003bp-9
//p3: violated_input is 0x1.95e4d9364d937p-8, lb is 0x1.5f77b475327f9p-9, ub is 0x1.5f77bc53b21bbp-9
//Polynomial: y=0x1.bcb7b151dd50ap-2 x^(1) + -0x1.bcb7ad0365d8cp-3 x^(2) + 0x1.28749ea54c589p-3 x^(3) +
//        -0x1.b60309bf9d835p-4 x^(4)
//
//            p3: number of violated intervals: 2, total iterations=16
//p3: VIOLATING INPUTS BELOW THRESHOLD:
//p3: starting input is 0x1.03652e52e52e6p-8
//p3: violated_input is 0x1.2af84a062b2e5p-8, lb is 0x1.03176656347c4p-9, ub is 0x1.03176698a5a14p-9
//p3: violated_input is 0x1.4212f684bda13p-8, lb is 0x1.1710980c1c9efp-9, ub is 0x1.171098627003bp-9
//Polynomial: y=0x1.bcb7b150e96dep-2 x^(1) + -0x1.bcb7a8b2a108p-3 x^(2) + 0x1.2871731216007p-3 x^(3) +
//        -0x1.b47f8ef897553p-4 x^(4)
//#endif
//            if (f == 0x1.4212f684bda13p-8) {
//                y = 0x1.1710980c1c9efp-9;
//            } else if (f == 0x1.95e4d9364d937p-8) {
//                y = 0x1.5f77b475327f9p-9;
//            } else {
//                double coeffs[] = {
//                        0x1.bcb7b151dd50ap-2,
//                        -0x1.bcb7ad0365d8cp-3,
//                        0x1.28749ea54c589p-3,
//                        -0x1.b60309bf9d835p-4
//                };
//
//                double temp1 = fma(f, coeffs[2], coeffs[1]);
//                double xsquare = f * f;
//                double temp2 = fma(xsquare, coeffs[3], temp1);
//                double temp3 = xsquare * temp2;
//                y = fma(f, coeffs[0], temp3);
//            }
//        }
//    }
//
//    y += m * LOG102LOW;
//    y += log10_lut[FIndex];
//    y += m * LOG102HIGH;
//
//    return y;
    float_x fix, fit;
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

        fix.f *= 8.388608e+06;
        m -= 23;
    }

    switch (fix.x) {
        case 0x3f800000 :
            return 0.0;
        case 0x41200000 :
            return 1.0;
        case 0x42c80000 :
            return 2.0;
        case 0x447a0000 :
            return 3.0;
        case 0x461c4000 :
            return 4.0;
        case 0x47c35000 :
            return 5.0;
        case 0x49742400 :
            return 6.0;
        case 0x4b189680 :
            return 7.0;
        case 0x4cbebc20 :
            return 8.0;
        case 0x4e6e6b28 :
            return 9.0;
        case 0x501502f9 :
            return 10.0;
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

    double y;

    if (f < 0x1.5a8f8d28ac42fp-9) {
        if (f < 0x1.5a91111111111p-10) {
            // 1st sub-domain
#if 0
            p0: number of violated intervals: 0, total iterations=19
p0: VIOLATING INPUTS BELOW THRESHOLD:
p0: starting input is 0x0p+0
Polynomial: y=0x1.bcb7b152566b9p-2 x^(1) + -0x1.bcb7ae3f54b59p-3 x^(2) + 0x1.286a6ae2f8166p-3 x^(3) +
        -0x1.88240065ec912p-4 x^(4)
#endif
            double coeffs[] = {
                    0x1.bcb7b152566b9p-2,
                    -0x1.bcb7ae3f54b59p-3,
                    0x1.286a6ae2f8166p-3,
                    -0x1.88240065ec912p-4
            };

            double temp1 = fma(f, coeffs[2], coeffs[1]);
            double xsquare = f * f;
            double temp2 = fma(xsquare, coeffs[3], temp1);
            double temp3 = xsquare * temp2;
            y = fma(f, coeffs[0], temp3);

        } else {
            // 2nd sub-doomain
#if 0
            p1: number of violated intervals: 1, total iterations=55
p1: VIOLATING INPUTS BELOW THRESHOLD:
p1: starting input is 0x1.5a91111111111p-10
p1: violated_input is 0x1.1fddb0d3224f3p-9, lb is 0x1.f386956531508p-11, ub is 0x1.f386ad8d4abbdp-11
Polynomial: y=0x1.bcb7b1516bf32p-2 x^(1) + -0x1.bcb7a4ab5429ep-3 x^(2) + 0x1.28608cd5626aap-3 x^(3) +
        -0x1.994d5e748761p-4 x^(4)
#endif
            if (f == 0x1.1fddb0d3224f3p-9) {
                y = 0x1.f386956531508p-11;
            } else {

                double coeffs[] = {
                        0x1.bcb7b1516bf32p-2,
                        -0x1.bcb7a4ab5429ep-3,
                        0x1.28608cd5626aap-3,
                        -0x1.994d5e748761p-4
                };

                double temp1 = fma(f, coeffs[2], coeffs[1]);
                double xsquare = f * f;
                double temp2 = fma(xsquare, coeffs[3], temp1);
                double temp3 = xsquare * temp2;
                y = fma(f, coeffs[0], temp3);

            }
        }
    } else {
        if (f < 0x1.03652e52e52e6p-8) {
            // 3rd sub-domain
#if 0
            p2: number of violated intervals: 3, total iterations=332
p2: VIOLATING INPUTS BELOW THRESHOLD:
p2: starting input is 0x1.5a8f8d28ac42fp-9
p2: violated_input is 0x1.8ff099fc267fp-9, lb is 0x1.5adab7bb93889p-10, ub is 0x1.5adabbecb58d9p-10
p2: violated_input is 0x1.bde34a2b10bf6p-9, lb is 0x1.82a2d0cd6ee1p-10, ub is 0x1.82a2d501e19b2p-10
p2: violated_input is 0x1.fbf5f5f5f5f5fp-9, lb is 0x1.b85b3cbbd9d74p-10, ub is 0x1.b85b3cbfd9c72p-10
Polynomial: y=0x1.bcb7b14f2e0aap-2 x^(1) + -0x1.bcb79a2b5a2aep-3 x^(2) + 0x1.285eb11241f0ep-3 x^(3) +
        -0x1.a518081cdab1fp-4 x^(4)
#endif
            if (f == 0x1.8ff099fc267fp-9) {
                y = 0x1.5adab7bb93889p-10;
            } else if (f == 0x1.bde34a2b10bf6p-9) {
                y = 0x1.82a2d0cd6ee1p-10;
            } else if(f==0x1.fbf5f5f5f5f5fp-9){
                y=0x1.b85b3cbbd9d74p-10;
            } else {
                double coeffs[] = {
                        0x1.bcb7b14f2e0aap-2,
                        -0x1.bcb79a2b5a2aep-3,
                        0x1.285eb11241f0ep-3,
                        -0x1.a518081cdab1fp-4
                };

                double temp1 = fma(f, coeffs[2], coeffs[1]);
                double xsquare = f * f;
                double temp2 = fma(xsquare, coeffs[3], temp1);
                double temp3 = xsquare * temp2;
                y = fma(f, coeffs[0], temp3);
            }
        } else {
            // 4th sub-domain
#if 0
            p3: number of violated intervals: 2, total iterations=16
p3: VIOLATING INPUTS BELOW THRESHOLD:
p3: starting input is 0x1.03652e52e52e6p-8
p3: violated_input is 0x1.2af84a062b2e5p-8, lb is 0x1.03176656347c4p-9, ub is 0x1.03176698a5a14p-9
p3: violated_input is 0x1.4212f684bda13p-8, lb is 0x1.1710980c1c9efp-9, ub is 0x1.171098627003bp-9
Polynomial: y=0x1.bcb7b150e96dep-2 x^(1) + -0x1.bcb7a8b2a108p-3 x^(2) + 0x1.2871731216007p-3 x^(3) +
        -0x1.b47f8ef897553p-4 x^(4)
#endif
            if (f == 0x1.2af84a062b2e5p-8) {
                y = 0x1.03176656347c4p-9;
            } else if (f == 0x1.4212f684bda13p-8) {
                y = 0x1.1710980c1c9efp-9;
            } else {
                double coeffs[] = {
                        0x1.bcb7b150e96dep-2,
                        -0x1.bcb7a8b2a108p-3,
                        0x1.2871731216007p-3,
                        -0x1.b47f8ef897553p-4
                };

                double temp1 = fma(f, coeffs[2], coeffs[1]);
                double xsquare = f * f;
                double temp2 = fma(xsquare, coeffs[3], temp1);
                double temp3 = xsquare * temp2;
                y = fma(f, coeffs[0], temp3);
            }
        }
    }

    y += m * LOG102LOW;
    y += log10_lut[FIndex];
    y += m * LOG102HIGH;

    return y;
}
