cmake_minimum_required( VERSION 2.5)
project(dense_blas_library)

include_directories(includes )

add_library (dense_blas_library STATIC BLAS.cpp)