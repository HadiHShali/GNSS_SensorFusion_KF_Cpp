cmake_minimum_required(VERSION 3.15)
project(KFClassTest LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
# Library
# These 3 lines set up a two-stage build — compile once as a library, then link into an executable.
add_library(kf_gps STATIC src/kalman_filter_gps.cpp)
add_executable(test_kf_static apps/test_kf_static.cpp)
target_link_libraries(test_kf_static PRIVATE kf_gps)
