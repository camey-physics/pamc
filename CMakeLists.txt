cmake_minimum_required(VERSION 3.10)
project(pamc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require GSL
find_package(GSL REQUIRED)

# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")

# Enable common warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# Enable testing
enable_testing()
find_package(GTest REQUIRED)

# Glob all test files (recursively through subdirs like tests/models)
file(GLOB_RECURSE TEST_FILES ${CMAKE_SOURCE_DIR}/tests/*.cpp)

# Source files for model implementations
set(MODEL_SOURCES
  ${CMAKE_SOURCE_DIR}/src/models/IsingModel.cpp
  ${CMAKE_SOURCE_DIR}/src/models/Ising3DHelpers.cpp
  ${CMAKE_SOURCE_DIR}/src/models/EAModel3DHelpers.cpp 
)

foreach(test_src ${TEST_FILES})
  get_filename_component(test_name ${test_src} NAME_WE)
  add_executable(${test_name} ${test_src} ${MODEL_SOURCES})
  target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${test_name} PRIVATE GTest::gtest_main GSL::gsl)
  gtest_discover_tests(${test_name})
endforeach()

# Build examples
add_executable(run_ising examples/run_ising.cpp ${MODEL_SOURCES})
target_include_directories(run_ising PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(run_ising PRIVATE GSL::gsl)

add_executable(run_3D_EA examples/run_3D_EA.cpp ${MODEL_SOURCES})
target_include_directories(run_3D_EA PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(run_3D_EA PRIVATE GSL::gsl)