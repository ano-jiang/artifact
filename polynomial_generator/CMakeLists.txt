find_package(OpenMP)
set(MPI_EXECUTABLE_SUFFIX ".mpich")
find_package(MPI)

find_package(Threads REQUIRED)

add_executable(polygen-estrin-fma piecewise-small-estrin-fma.cpp rlibm-fast.h)
target_link_libraries(polygen-estrin-fma $ENV{SOPLEX_LIB} gmp z m OpenMP::OpenMP_CXX MPI::MPI_CXX)
set_target_properties(polygen-estrin-fma PROPERTIES
                    CXX_STANDARD 17
                    CXX_STANDARD_REQUIRED ON
                    CXX_EXTENSIONS OFF)
target_include_directories(polygen-estrin-fma PRIVATE $ENV{SOPLEX_INCLUDE})
target_compile_options(polygen-estrin-fma PRIVATE -O2 -g -DTBB_USE_DEBUG=ON)
