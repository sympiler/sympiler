include_directories(
        ${BLAS_INC}
        ${SPARSE_LIB_INC}
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${CODE_FUSION_INC}
        ${MKL_INC})

link_directories(
        ${SPARSE_LIB_LIB}
        ${SPARSE_UTIL_LIB}
        ${MKL_LIB}
        ${MKL_LIB2}
        ${PAPI_LIB}
        ${PAPI_LIB2})


add_executable(Catch_tests ${SOURCE_FILES}
        spmv_tests.cpp
        sptrsv_tests.cpp
        transpose_tests.cpp
        )

target_link_libraries(Catch_tests
        ${ALL_LIBS}
        sym_sparse_blas
        )