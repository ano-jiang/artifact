#define __ELEM__ exp

#include "prlibm.h"
#include "performance_test/LibTestHelper.h"

int additionallyIgnoreThisInput(float x) {
    float_x fx;
    fx.f = x;

    // Take care of special cases
    if ((fx.x & 0x7FFFFFFF) == 0) return 1;
    if (fx.x <= 872415231) return 1;
    if (1118925336 <= fx.x && fx.x <= 3011510272) return 1;
    if (fx.x >= 3268407733) return 1;

    return 0;
}


int main(int argc, char **argv) {
    RunTest("./prlibm/logs/exp.txt");
    return 0;
}
