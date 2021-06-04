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

set(SOURCE_FILES
        ${SPARSE_LIB_SRC}/spmv.cpp
        ${SPARSE_LIB_SRC}/sptrsv.cpp
        ${SPARSE_LIB_SRC}/add.cpp
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/BCSCMatrix.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${BLAS_SRC}/BLAS.cpp
        )

add_executable(Catch_tests ${SOURCE_FILES}
        spmv_tests.cpp
        sptrsv_tests.cpp
        transpose_tests.cpp
        )

target_link_libraries(Catch_tests
        ${ALL_LIBS})