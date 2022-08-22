cmake_minimum_required( VERSION 3.0)
project(sym_sparse_blas)

include_directories(
        #include
        #include/parsy
        ${BLAS_INC}
       # ${SPARSE_EXTN_INC}
       # ${SPARSE_UTIL_INC}
        )

if(MKL_FOUND)
file(GLOB_RECURSE PARSY_SRC parsy/*.cpp)
endif()

file(GLOB_RECURSE LINSOLVE_SRC linear_solver/*.cpp)

set(BLAS_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas/BLAS.cpp" )

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
        spcholesky.cpp
        other_operations.cpp
        diagonal_solver.cpp
        ${BLAS_SRC_FILES}
        ${PARSY_SRC}
        ${LINSOLVE_SRC}
        )

#add_dependencies(sym_sparse_blas metis)
target_link_libraries(sym_sparse_blas PRIVATE
        ${ALL_LIBS}
        )
if(OpenMP_FOUND)
    target_link_libraries(sym_sparse_blas PRIVATE OpenMP::OpenMP_CXX)
endif()

#[[if(SYMPILER_USE_BLAS_MKL)
    target_link_libraries(sym_sparse_blas PRIVATE mkl::mkl)
    target_compile_definitions(sym_sparse_blas PRIVATE "MKL_BLAS")
else]]
if(SYMPILER_USE_BLAS_OpenBLAS)
    target_link_libraries(sym_sparse_blas PRIVATE OpenBLAS::OpenBLAS)
    target_compile_definitions(sym_sparse_blas PRIVATE "OPENBLAS")
endif()
