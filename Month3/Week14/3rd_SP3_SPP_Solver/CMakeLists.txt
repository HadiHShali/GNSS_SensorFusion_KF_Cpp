cmake_minimum_required(VERSION 3.15)
project(SppPreciseProducts LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
add_library(satellites_lib STATIC src/satellites.cpp)
add_library(solver_lib     STATIC src/gps_ls_solver.cpp)
add_library(utils_lib      STATIC src/utils.cpp)
add_library(rinex_obs_lib  STATIC src/rinex_obs.cpp)
add_library(sp3_lib        STATIC src/sp3_reader.cpp)
add_library(clk_lib        STATIC src/clk_reader.cpp)
 
add_executable(spp_processor_precise apps/spp_processor_precise.cpp)
target_link_libraries(spp_processor_precise PRIVATE
    satellites_lib solver_lib utils_lib rinex_obs_lib sp3_lib clk_lib)
