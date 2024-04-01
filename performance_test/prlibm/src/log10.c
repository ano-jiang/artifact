#define __ELEM__ prlibm_log10

#include "prlibm.h"
#include "performance_test/LibTestHelper.h"

int additionallyIgnoreThisInput(float x) {
    float_x fix, fit;
    fix.f = x;
    int m = 0;

    if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
        if ((fix.x & 0x7FFFFFFF) == 0) return 1;

        if (fix.x > 0x7FFFFFFF) return 1;

        if (fix.x >= 0x7F800000) return 1;
    }

    switch (fix.x) {
        case 0x3f800000 :
            return 1;
        case 0x41200000 :
            return 1;
        case 0x42c80000 :
            return 1;
        case 0x447a0000 :
            return 1;
        case 0x461c4000 :
            return 1;
        case 0x47c35000 :
            return 1;
        case 0x49742400 :
            return 1;
        case 0x4b189680 :
            return 1;
        case 0x4cbebc20 :
            return 1;
        case 0x4e6e6b28 :
            return 1;
        case 0x501502f9 :
            return 1;
    }

    return 0;
}

int main(int argc, char **argv) {
    RunTest("./prlibm/logs/log10.txt");
    return 0;
}
