#define MPFR_PREC 200
#include "mpfr.h"
#include "prlibm.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

mpfr_t mval;

float MpfrCalculate(float x);
float rlibmTest(float x);

float ToMPFR(mpfr_t mval) {
    double retVal = mpfr_get_d(mval, MPFR_RNDN);

    if (retVal == 0) return 0.0f;
    if (retVal != retVal) {
        return retVal;
    }

    if (mpfr_cmp_d(mval, pow(2, -150)) <= 0 && mpfr_cmp_d(mval, -pow(2, -150)) >= 0) {
        return 0.0f;
    }

    long exp;
    double fr = mpfr_get_d_2exp(&exp, mval, MPFR_RNDN);
    fr *= 2;
    exp--;

    if (mpfr_cmp_d(mval, 0.0) > 0) {
        if (mpfr_cmp_d(mval, 1.5 * pow(2, -149)) < 0) return pow(2, -149);
        if (mpfr_cmp_d(mval, pow(2, -148)) < 0) return pow(2, -148);

    } else {
        if (mpfr_cmp_d(mval, -1.5 * pow(2, -149)) > 0) return -pow(2, -149);
        if (mpfr_cmp_d(mval, -pow(2, -148)) > 0) return -pow(2, -148);
    }

    if (exp >= -148 && exp <= -127) {
        int prec = 150 + exp;
        mpfr_t r;
        mpfr_init2(r, prec);
        mpfr_set(r, mval, MPFR_RNDN);
        retVal = mpfr_get_d(r, MPFR_RNDN);
        mpfr_clear(r);
        return retVal;
    } else {
        mpfr_t r;
        mpfr_init2(r, 24);
        mpfr_set(r, mval, MPFR_RNDN);
        retVal = mpfr_get_d(r, MPFR_RNDN);
        mpfr_clear(r);
        return retVal;
    }
}

unsigned int m_ulpf(float x, float y) {
    if (x == 0)
        x = 0; // -0 == 0
    if (y == 0)
        y = 0; // -0 == 0

    if (x != x && y != y) return 0;

    if (x != x)
        return UINT_MAX - 1; // Maximum error
    if (y != y)
        return UINT_MAX - 1; // Maximum error

    int xx = *((int *)&x);
    xx = xx < 0 ? INT_MAX - xx : xx;

    int yy = *((int *)&y);
    yy = yy < 0 ? INT_MAX - yy : yy;
    return xx >= yy ? xx - yy : yy - xx;
}

void RunCorrectnessTest(char const* FunctionName, char* resFileName) {
    int printcount = 0;
    mpfr_init2(mval, MPFR_PREC);

    unsigned long wrongRlibmCount = 0;
    unsigned long count = 0;

    unsigned int maxUlpRlibm = 0;
    float maxXRlibm = 0;

//    FILE* f = fopen(resFileName, "w");

    float x;
    float_x xbase;
//    for (count = 0x0; count < 0x100000000; count++) {
        xbase.x = 0x3f7e3274;
        x = xbase.f;

        float bmy = MpfrCalculate(x);
        float bres = rlibmTest(x);

        // if bres is nan and bmy is nan, continue
//        if (bres != bres && bmy != bmy) continue;

        // Otherwise check if the output is correct
        if (bres != bmy) {
            float_x y;
            y.f=bres;
            if (printcount < 10) {
                printf("count = %lu\n", count);
                printf("bres = %x\n", y.x);
                y.f=bmy;
                printf("bmy  = %x\n", y.x);
                printcount++;
            }
            wrongRlibmCount++;
            unsigned int error = m_ulpf(bres, bmy);
            if (error > maxUlpRlibm) {
                maxUlpRlibm = error;
                maxXRlibm = x;
            }
        }
//    }

    printf("%s TEST RESULT:\n", FunctionName);
    if (wrongRlibmCount == 0) {
        printf("RLIBM returns correct result for all inputs\n");
    } else {
        printf("RLIBM: Found %lu/%lu inputs with wrong result\n", wrongRlibmCount, count);
        printf("RLIBM's max ulp error: %u\n", maxUlpRlibm);
        printf("RLIBM max ulp error at: %.100e\n", maxXRlibm);
    }

//    fclose(f);
    mpfr_clear(mval);
}

float MpfrCalculate(float x) {
    mpfr_set_d(mval, (float)x, MPFR_RNDN);
    mpfr_log2(mval, mval, MPFR_RNDN);
    return ToMPFR(mval);
}

float rlibmTest(float x) {
    float result = prlibm_log2(x);
    return result;
}

int main(int argc, char** argv) {

//    if (argc != 2) {
//        printf("Usage: %s <output file>\n", argv[0]);
//        exit(0);
//    }

    RunCorrectnessTest("Float - Log2", argv[1]);
    return 0;
}