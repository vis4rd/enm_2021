cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(IMN_LAB9)
set(CMAKE_COLOR_MAKEFILE ON)

add_executable(imn_lab9)

target_sources(imn_lab9 PRIVATE "main.cpp")
target_compile_options(imn_lab9 PRIVATE "-Wall" "-Wextra" "-pedantic" "-std=c++17")
target_link_libraries(imn_lab9 PRIVATE gsl gslcblas)
