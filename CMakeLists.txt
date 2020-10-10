cmake_minimum_required(VERSION 3.17)
project(SharedGenerator)

set(CMAKE_CXX_STANDARD 14)

add_executable(SharedGenerator main.cpp)

SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")