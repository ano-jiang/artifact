find_package(OpenMP)
set(MPI_EXECUTABLE_SUFFIX ".mpich")
find_package(MPI)

add_executable(polygen-estrin-fma-naive piecewise-small-estrin-fma.cpp rlibm-fast.h)

target_link_libraries(polygen-estrin-fma-naive $ENV{SOPLEX_LIB} gmp z m OpenMP::OpenMP_CXX MPI::MPI_CXX)
target_include_directories(polygen-estrin-fma-naive PRIVATE $ENV{SOPLEX_INCLUDE})
target_compile_options(polygen-estrin-fma-naive PRIVATE -std=c++17 -O2 -g -DTBB_USE_DEBUG=ON)
