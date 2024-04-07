cmake_minimum_required(VERSION 3.0.0)
project(N_Body_Simulation VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_FLAGS "-Wall -O3 -fopenmp")

add_executable(N_Body_Simulation src/main.cpp src/oct_tree.cpp)