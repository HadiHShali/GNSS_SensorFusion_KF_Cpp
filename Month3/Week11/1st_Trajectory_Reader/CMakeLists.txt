cmake_minimum_required(VERSION 3.15)
project(TrajReader LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
add_executable(test_trajectory_reader
    apps/test_trajectory_reader.cpp
    src/trajectory_reader.cpp
    src/utils.cpp)

