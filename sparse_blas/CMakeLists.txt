cmake_minimum_required( VERSION 2.5)
project(sparse_blas_library)

include_directories(includes
        ${BLAS_INC}
        ${SPARSE_UTIL_INC}
        ${MKL_INC}
        )

add_library (sparse_blas_library STATIC
        spmv.cpp
        sptrsv.cpp
        add.cpp
        spic0.cpp
        spilu0.cpp
        sptrsv_lt.cpp
        spldl.cpp
        other_operations.cpp
        diagonal_solver.cpp
        )
target_link_libraries(sparse_blas_library
        ${ALL_LIBS}
        )