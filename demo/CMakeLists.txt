cmake_minimum_required(VERSION 3.0)



include_directories(includes
        #${SPARSE_UTIL_INC}
        #${SPARSE_EXTN_INC}
        ${SPARSE_LIB_INC}
        ${LBC_INC}
        )

link_directories(
        ${CMAKE_CURRENT_LIST_DIR}/../lib/
)

project(sptrsv_demo)
set(SPTRSVD_SOURCE_FILES
        ${METIS_SOURCE_FILES}
       # ${ALL_SOURCE_FILES}
        FusionDemo.cpp
        )

add_executable(sptrsv_demo
        ${SPTRSVD_SOURCE_FILES}
        sptrsv_demo.cpp )


target_link_libraries(sptrsv_demo
        ${ALL_LIBS}
        sym_sparse_blas
        )