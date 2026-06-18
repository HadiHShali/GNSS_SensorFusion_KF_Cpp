cmake_minimum_required(VERSION 3.15)
project(RobustSolver LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
 
add_executable(robust_solver apps/robust_solver_main.cpp src/gps_ls_solver.cpp src/utils.cpp)
