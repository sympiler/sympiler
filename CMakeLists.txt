cmake_minimum_required( VERSION 3.15)


list(APPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake
        "${CMAKE_CURRENT_LIST_DIR}/cmake/third_party"
        )
include(DownloadProject)

project(sympiler_project)

set(SYMPILER_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(SYMPILER_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)

option(SYMPILER_USE_TBB "Use TBB for threading" OFF)

# Defining paths
set(SPARSE_LIB_INC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/includes")
set(SPARSE_LIB_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(SPARSE_LIB_SRC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas")
AUX_SOURCE_DIRECTORY(${SPARSE_LIB_SRC} SPARSE_LIB_SRC_FILES)

set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/utils")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")
#[[
set(SPARSE_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/utils/include")
set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/utils")
#AUX_SOURCE_DIRECTORY(${SPARSE_UTIL_SRC} SPARSE_UTIL_SRC_FILES)
set(SPARSE_UTIL_SRC_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/BCSCMatrix.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/etree.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_inspector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_io.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_utilities.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/lfactor_creation.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/utils.cpp")

set(SPARSE_EXTN_INC "${SPARSE_UTIL_SRC}/external/include")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")
AUX_SOURCE_DIRECTORY(${SPARSE_EXTN_SRC} SPARSE_EXTN_SRC_FILES)
]]



set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/include")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/src/")
set(LBC_LIB "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/lib/")
set(LBC_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/")
#AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)


#AUX_SOURCE_DIRECTORY(${BLAS_SRC} BLAS_SRC_FILES)
#set(BLAS_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/dense_blas/BLAS.cpp" )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)



# Find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    if(APPLE) #TODO: set path for brew install libomp
        set(OpenMP_CXX_INCLUDE_DIRS  "/usr/local/include/")
        include_directories(
                ${OpenMP_CXX_INCLUDE_DIRS}
        )
        link_libraries(
                OpenMP::OpenMP_CXX
        )
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmacosx-version-min=10.6")
    endif()
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP)


endif()

if(SYMPILER_USE_TBB)
    add_definitions(-DUSE_TBB)
    message(STATUS "TBB threading is enabled.")
endif()
#add_definitions( -DPROFILE)
#add_definitions( -DLOAD)

include(metis)
include(tbb)

find_package(MKL OPTIONAL_COMPONENTS)
if(MKL_FOUND)
    add_definitions(-DMKL) # add openblas
else()
    message(STATUS "Intel MKL not found, skipping supernodal Cholesky lib")
endif()

### Find PAPI for Profiling
### Can be cloned directly from https://bitbucket.org/icl/papi.git
#if(EXISTS "$ENV{PAPIROOT}")
#    if(EXISTS "$ENV{PAPIROOT}/src")
#        set(PAPI_LIB "$ENV{PAPIROOT}/src")
#        set(PAPI_LIB2 "${PAPI_LIB}/libpfm4/lib/")
#    elseif(EXISTS "$ENV{PAPIROOT}/lib")
#        set(PAPI_LIB "$ENV{PAPIROOT}/lib")
#    endif()
#    if(EXISTS "$ENV{PAPIROOT}/include")
#        set(PAPI_INC "$ENV{PAPIROOT}/include")
#    else()
#        set(PAPI_INC "$ENV{PAPIROOT}/src")
#    endif()
#
#    add_definitions(-DPAPI)
#    set(PAPI_LIB_NAME "papi")
#    set(PAPI_LIB_NAME2 "pfm")
#    message(STATUS "PAPI found")
#else()
#    message(STATUS "PAPI not found, building without PAPI.")
#    set(PAPI_LIB_NAME "")
#endif()

include_directories(
        ${MKL_INCLUDE_DIR}
        ${BLAS_INC}
        ${LBC_INC}
        include
        include/parsy/
        include/linear_solver/
)


set(ALL_LIBS
        ${MKL_LIBRARIES}
        ${PAPI_LIB_NAME}
        ${PAPI_LIB_NAME2}
        metis::metis
        tbb_static
        lbc
        )

set(ALL_SOURCE_FILES
        ${BLAS_SRC_FILES}
        ${METIS_SOURCE_FILES})

set(CMAKE_CXX_STANDARD 11)
if(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11 -w -m64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")
endif()

Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")



if(METIS_FOUND)
    add_definitions(-DMETIS)
    add_subdirectory(demo)
endif()

add_subdirectory(aggregation)
add_subdirectory(sparse_blas)
add_subdirectory(Catch_tests)
add_subdirectory(code_gen)
#add_subdirectory(utils)  ## TODO: rmeove the folder
#add_subdirectory(dense_blas)
