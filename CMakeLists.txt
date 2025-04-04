cmake_minimum_required(VERSION 3.10)
project(pamc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing
enable_testing()
include(GoogleTest)

# Find required libraries
find_package(GSL REQUIRED)
find_package(GTest REQUIRED)

# Define paths
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/tests)

# Compiler optimizations for g++ and Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
endif()

# Collect source files automatically
file(GLOB_RECURSE SRC_FILES "${SRC_DIR}/*.cpp")

# Collect all include directories recursively
file(GLOB_RECURSE HEADER_DIRS LIST_DIRECTORIES true ${INCLUDE_DIR}/*)
set(INCLUDE_DIRS ${INCLUDE_DIR})
foreach(dir ${HEADER_DIRS})
    if(IS_DIRECTORY ${dir})
        list(APPEND INCLUDE_DIRS ${dir})
    endif()
endforeach()

# Define the main library
add_library(pamc ${SRC_FILES})
target_include_directories(pamc PUBLIC ${INCLUDE_DIRS})
target_compile_options(pamc PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)

# Define the test executable
file(GLOB_RECURSE TEST_SOURCES "${TEST_DIR}/*.cpp")
add_executable(pamcTest ${TEST_SOURCES})
target_link_libraries(pamcTest PRIVATE pamc GTest::GTest GTest::Main GSL::gsl GSL::gslcblas)

# Discover and register tests
gtest_discover_tests(pamcTest)