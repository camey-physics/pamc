cmake_minimum_required(VERSION 3.10)
project(pamc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require GSL
find_package(GSL REQUIRED)

# Require OpenMP
find_package(OpenMP REQUIRED)

# Enable testing
enable_testing()
find_package(GTest REQUIRED)

# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")

# Enable common warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# Common includes and libraries
include_directories(${CMAKE_SOURCE_DIR}/include)
set(COMMON_LIBS GSL::gsl OpenMP::OpenMP_CXX)

# Source files for model implementations
set(MODEL_SOURCES
  ${CMAKE_SOURCE_DIR}/src/models/IsingModel.cpp
  ${CMAKE_SOURCE_DIR}/src/models/Ising3DHelpers.cpp
  ${CMAKE_SOURCE_DIR}/src/models/EAModel3DHelpers.cpp 
)

# Build examples
add_executable(run_ising examples/run_ising.cpp ${MODEL_SOURCES})
target_link_libraries(run_ising PRIVATE ${COMMON_LIBS})

add_executable(run_3D_EA examples/run_3D_EA.cpp ${MODEL_SOURCES})
target_link_libraries(run_3D_EA PRIVATE ${COMMON_LIBS})

# Glob all test files (recursively)
file(GLOB_RECURSE TEST_FILES ${CMAKE_SOURCE_DIR}/tests/*.cpp)

foreach(test_src ${TEST_FILES})
  get_filename_component(test_name ${test_src} NAME_WE)
  add_executable(${test_name} ${test_src} ${MODEL_SOURCES})
  target_link_libraries(${test_name} PRIVATE GTest::gtest_main ${COMMON_LIBS})
  gtest_discover_tests(${test_name})
endforeach()
