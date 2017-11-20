cmake_minimum_required(VERSION 3.3)

project(trnsEigen)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11   -march=x86-64 -Wunused-result")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopt-info-vec")

set(SOURCE_FILES TestEigen.cpp)
add_executable(trnsEigen ${SOURCE_FILES})#
target_link_libraries(trnsEigen
        openblas
        )