cmake_minimum_required( VERSION 2.5)
project(sympiler_utils)

include_directories(includes
        external/includes/
        )


file(GLOB_RECURSE EXT_SRC external/*.cpp external/*.c)

add_library (sympiler_utils STATIC
        BCSCMatrix.cpp
        etree.cpp
        lfactor_creation.cpp
        sparse_inspector.cpp
        sparse_io.cpp
        sparse_utilities.cpp
        test_utils.cpp
        utils.cpp
        ${EXT_SRC}

#        external/mmio.cpp #TODO: remove this and use SMP
#        external/amd_order.c
#        external/amd_1.c
#        external/amd_2.c

        )
#target_link_libraries(sympiler_utils)