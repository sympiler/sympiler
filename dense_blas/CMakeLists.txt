cmake_minimum_required( VERSION 3.0)
project(sym_dense_blas)

include_directories(includes )

add_library (sym_dense_blas STATIC BLAS.cpp)