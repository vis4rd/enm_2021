cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(IMN_LAB6)

set(CMAKE_COLOR_MAKEFILE True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_subdirectory("${PROJECT_SOURCE_DIR}/lib")

add_executable(imn_lab6.out)

target_compile_options(imn_lab6.out PRIVATE "-Wall" "-Wextra" "-pedantic")
target_compile_features(imn_lab6.out PRIVATE "cxx_std_20")
target_sources(imn_lab6.out PRIVATE ${SOURCES} main.cpp)
target_include_directories(imn_lab6.out PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(imn_lab6.out PRIVATE m imn_lab6_lib)
