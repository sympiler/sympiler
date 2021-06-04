cmake_minimum_required( VERSION 2.5)
project(sym_sparse_blas)

include_directories(includes
        includes/parsy
        ${BLAS_INC}
        ${SPARSE_EXTN_INC}
        ${SPARSE_UTIL_INC}
        ${MKL_INC}
        )
link_directories(
        ${METIS_LIBRARY}
)
file(GLOB_RECURSE PARSY_SRC parsy/*.cpp)



if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
add_library (sym_sparse_blas STATIC
        spmv.cpp
        sptrsv.cpp
        add.cpp
        spic0.cpp
        spilu0.cpp
        sptrsv_lt.cpp
        spldl.cpp
        other_operations.cpp
        diagonal_solver.cpp
        ${PARSY_SRC}
        )

#add_dependencies(sym_sparse_blas metis)
target_link_libraries(sym_sparse_blas
        ${ALL_LIBS}
        )
