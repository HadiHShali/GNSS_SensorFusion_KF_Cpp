cmake_minimum_required(VERSION 3.15)
project(KFClass LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
# Library
add_library(kf_gps STATIC src/kalman_filter_gps.cpp)
# We'll add a test executable in apps/ tomorrow (Day 4)