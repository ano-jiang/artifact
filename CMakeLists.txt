cmake_minimum_required(VERSION 3.21)
project(prlibm C CXX)
add_subdirectory(polynomial_generator)
add_subdirectory(libm)
add_subdirectory(correctness_test)
add_subdirectory(performance_test)
add_subdirectory(prlibm-naive)
