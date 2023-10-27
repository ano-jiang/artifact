set -e

rm -rf cmake-build-release/CMakeCache.txt
source env.sh
cmake -B cmake-build-release/
cmake --build cmake-build-release -- VERBOSE=1 
