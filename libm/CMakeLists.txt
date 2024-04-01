file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.c include/*.h)

#set(CMAKE_C_STANDARD 11)
add_library(prlibm STATIC ${srcs})

target_compile_options(prlibm PRIVATE -std=c11 -march=native -O3)
target_include_directories(prlibm PUBLIC include/libm)

target_link_libraries(prlibm PUBLIC m)
