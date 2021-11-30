cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(IMN_LAB6_LIB C)

set(CMAKE_COLOR_MAKEFILE True)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.c")

add_library(imn_lab6_lib STATIC)

target_compile_options(imn_lab6_lib
    PRIVATE "-Wall"
            "-Wextra"
            "-pedantic"
            "-std=c11"
    PUBLIC
            "-Wno-unused-parameter"
            "-Wno-unused-variable"
            "-Wno-format-extra-args"
            "-Wno-unused-but-set-variable"
            "-Wno-maybe-uninitialized"
)
target_sources(imn_lab6_lib PRIVATE ${SOURCES})
target_include_directories(imn_lab6_lib INTERFACE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(imn_lab6_lib PRIVATE m)
