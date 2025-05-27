cmake_minimum_required(VERSION 3.10)
project(pamc LANGUAGES CXX)

# Enable debug info and disable optimizations for debugging
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address -fno-omit-frame-pointer")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# Enable and configure testing
enable_testing()
find_package(GTest REQUIRED)

file(GLOB TEST_FILES ${CMAKE_SOURCE_DIR}/tests/*.cpp)

foreach(test_src ${TEST_FILES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE GTest::gtest_main)
    gtest_discover_tests(${test_name})
endforeach()
