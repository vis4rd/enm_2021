cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(IMN_LAB5)

set(CMAKE_COLOR_MAKEFILE True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(imn_lab5.out)

target_compile_options(imn_lab5.out PRIVATE "-Wall" "-Wextra" "-pedantic" "-std=c++20")
target_sources(imn_lab5.out PRIVATE ${SOURCES} main.cpp)
target_include_directories(imn_lab5.out PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(imn_lab5.out PRIVATE m)
