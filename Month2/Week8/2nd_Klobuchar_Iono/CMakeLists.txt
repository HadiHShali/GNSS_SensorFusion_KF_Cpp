cmake_minimum_required(VERSION 3.15)
project(KlobucharIono LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen library (for any future code using matrices)
include_directories("C:/eigen-3.4.0")

# Your project headers
include_directories(include)

# Visual Studio needs this so M_PI works in <cmath>
add_definitions(-D_USE_MATH_DEFINES)

# ── Today's test program: standalone Klobuchar test ──
add_executable(test_klobuchar
    apps/test_klobuchar.cpp
    src/satellites.cpp
)