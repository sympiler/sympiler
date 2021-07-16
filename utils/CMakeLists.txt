cmake_minimum_required( VERSION 3.0)
project(sympiler_utils)

include_directories(includes
        external/includes/
        )


file(GLOB_RECURSE EXT_SRC external/*.cpp external/*.c)

#message(" ===> ${METIS_SOURCE_FILES}")
add_library (sympiler_utils STATIC
        BCSCMatrix.cpp
        lfactor_creation.cpp
        sparse_inspector.cpp
        sparse_io.cpp
        test_utils.cpp
        ${METIS_SOURCE_FILES}
        ${EXT_SRC}

#        external/mmio.cpp #TODO: remove this and use SMP

        )
#target_link_libraries(sympiler_utils)