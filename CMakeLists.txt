cmake_minimum_required(VERSION 3.0)
add_subdirectory(Lpopc)
add_subdirectory(Lpopc/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)