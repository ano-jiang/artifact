#define __ELEM__ exp10

#include "prlibm.h"
#include "performance_test/LibTestHelper.h"

int additionallyIgnoreThisInput(float x) {
    float_x fx;
    fx.f = x;

    // Take care of special cases
    if ((fx.x & 0x7FFFFFFF) == 0) return 1;

    if (fx.x <= 861821911) return 1;

    if (1109008539 <= fx.x && fx.x <= 3000916953) return 1;

    if (fx.x >= 3258228278) return 1;

    // If x == 0.0, 1.0, 2.0, ..., 10.0, then it's also special case
    switch (fx.x) {
        case 0x00000000:
        case 0x80000000:
        case 0x3f800000:
        case 0x40000000:
        case 0x40400000:
        case 0x40800000:
        case 0x40a00000:
        case 0x40c00000:
        case 0x40e00000:
        case 0x41000000:
        case 0x41100000:
        case 0x41200000:
            return 1;
    }

    return 0;
}

int main(int argc, char **argv) {
    RunTest("./prlibm/logs/exp10.txt");
    return 0;
}
