cmake_minimum_required(VERSION 3.15)



include_directories(includes
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


if(OpenMP_FOUND)
    target_link_libraries(sptrsv_demo PRIVATE
            ${ALL_LIBS}
            sym_sparse_blas
            OpenMP::OpenMP_CXX
            )
else()
    target_link_libraries(sptrsv_demo PRIVATE
            ${ALL_LIBS}
            sym_sparse_blas
            )

endif()
