cmake_minimum_required(VERSION 3.17)
project(SharedGenerator)

set(CMAKE_CXX_STANDARD 14)
include_directories("C:/Xilinx/Vivado/2019.2/include")

add_executable(SharedGenerator main.cpp lognum.cpp lognum.h ConversionEngine.cpp)

SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")