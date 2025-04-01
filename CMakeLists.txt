cmake_minimum_required(VERSION 3.10)
project(PAMC-Project VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required libraries
find_package(GSL REQUIRED CONFIG)
find_package(GTest REQUIRED)

# Define paths
set(SRC_DIR src)
set(INCLUDE_DIR include)
set(TEST_DIR tests)

# Compiler optimizations for g++ and Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
endif()

# Collect source files automatically
file(GLOB_RECURSE PAMC_SOURCES "${SRC_DIR}/*.cpp")

# Define the main library
add_library(PAMC ${PAMC_SOURCES})
target_include_directories(PAMC PUBLIC ${INCLUDE_DIR})
target_compile_options(PAMC PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)

# Enable testing
enable_testing()
include(GoogleTest)

# Define the test executable
file(GLOB_RECURSE TEST_SOURCES "${TEST_DIR}/*.cpp")
add_executable(PAMCTest ${TEST_SOURCES})
target_link_libraries(PAMCTest PRIVATE PAMC GTest::GTest GTest::Main GSL::gsl GSL::gslcblas)

# Discover and register tests
gtest_discover_tests(PAMCTest)
