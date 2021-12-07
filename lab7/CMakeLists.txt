cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(IMN_LAB7 C)

set(CMAKE_COLOR_MAKEFILE True)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.c")

add_executable(imn_lab7)

target_compile_options(imn_lab7
    PRIVATE "-Wall"
            "-Wextra"
            "-pedantic"
            "-std=c11"
)
target_sources(imn_lab7 PRIVATE ${SOURCES} main.c)
target_include_directories(imn_lab7 PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(imn_lab7 PRIVATE m)
