# CMakeLists.txt for Week 6 Day 3
cmake_minimum_required(VERSION 3.15)
project(LeastSquaresSolver LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
# Tell CMake where Eigen lives
include_directories("C:/eigen-3.4.0")
 
add_executable(ls_solver apps/ls_solver.cpp)