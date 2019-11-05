# - Basic CMake setup
# Check version meets our requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(DigitalFilter)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Need ROOT library separately
find_package(ROOT REQUIRED)

# find the headers
include_directories(.)

# Build a dynamic library from our sources
add_library(DFilter SHARED DFilter.h DFilter.cpp)

# Build a testing executable separately for file output
include_directories(. test ${ROOT_INCLUDE_DIRS})
add_executable(testdf test/testdfilter.cpp)
target_link_libraries(testdf PUBLIC DFilter ${ROOT_LIBRARIES})

# Build the testing code, tell CTest about it
enable_testing()
set(CMAKE_CXX_STANDARD 11)

# main catch library compile only once
add_library(Catch test/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/test)

# add test executable
include_directories(. test)
add_executable(filtertest test/filtertest.cpp)
target_link_libraries(filtertest PUBLIC Catch DFilter)

# CTest comand
add_test(NAME filtertest
  COMMAND filtertest -s
)

