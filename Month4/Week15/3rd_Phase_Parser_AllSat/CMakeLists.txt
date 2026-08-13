cmake_minimum_required(VERSION 3.15)
project(PhaseParserAllSat LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
add_library(rinex_obs_lib STATIC src/rinex_obs.cpp)
 
add_executable(test_phase_parser_allSat apps/test_phase_parser_allSat.cpp)
target_link_libraries(test_phase_parser_allSat PRIVATE rinex_obs_lib)

