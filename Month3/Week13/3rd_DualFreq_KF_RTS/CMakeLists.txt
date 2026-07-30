cmake_minimum_required(VERSION 3.15)
project(DualFreqSpp_KF_RTS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)

# Library: KF class (with BiasEstimator, rtsSmooth) + trajectory CSV reader + utils
add_library(kf_gps      STATIC src/kalman_filter_gps.cpp)
add_library(traj_reader STATIC src/trajectory_reader.cpp)
add_library(utils_lib   STATIC src/utils.cpp)

# Executable
add_executable(DualFreqSpp_KF_RTS apps/DualFreqSpp_KF_RTS.cpp)
target_link_libraries(DualFreqSpp_KF_RTS PRIVATE
    kf_gps traj_reader utils_lib)