cmake_minimum_required(VERSION 2.7)
project(SymTriang)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-DMYBLAS -std=gnu++11" )

include_directories(${MKL_INC} "../../symGen"  )

set(SOURCE_FILES main.cpp)
add_executable(SymTriang ${SOURCE_FILES})

target_link_libraries(SymTriang   )