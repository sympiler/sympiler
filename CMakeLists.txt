cmake_minimum_required( VERSION 3.16)


list(APPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake
        "${CMAKE_CURRENT_LIST_DIR}/cmake/third_party"
        )
include(DownloadProject)

project(sympiler_project)

#----------------------------[ global configuration ]--------------------------#

set(CMAKE_CXX_STANDARD 11)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(SYMPILER_IS_TOPLEVEL OFF)
else()
    set(SYMPILER_IS_TOPLEVEL ON)
endif()

if(UNIX )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)


#----------------------------------[ options ]---------------------------------#

option(SYMPILER_BUILD_EXAMPLES "Build SYMPILER Examples" ${SYMPILER_IS_TOPLEVEL})
option(SYMPILER_BUILD_DOCS "Build SYMPILER Documentation" ${SYMPILER_IS_TOPLEVEL})


set(SYMPILER_BLAS_BACKEND "OpenBLAS" CACHE STRING "BLAS implementation for SYMPILER to use")
set_property(CACHE SYMPILER_BLAS_BACKEND PROPERTY STRINGS MKL OpenBLAS Apple)

if(${SYMPILER_BLAS_BACKEND} STREQUAL "MKL")
    set(SYMPILER_USE_BLAS_MKL      ON)
    set(SYMPILER_USE_BLAS_OpenBLAS OFF)
    set(SYMPILER_USE_BLAS_Apple OFF)
elseif(${SYMPILER_BLAS_BACKEND} STREQUAL "OpenBLAS")
    set(SYMPILER_USE_BLAS_MKL      OFF)
    set(SYMPILER_USE_BLAS_OpenBLAS ON)
    set(SYMPILER_USE_BLAS_Apple OFF)
elseif(${SYMPILER_BLAS_BACKEND} STREQUAL "Apple")
    set(SYMPILER_USE_BLAS_MKL      OFF)
    set(SYMPILER_USE_BLAS_OpenBLAS OFF)
    set(SYMPILER_USE_BLAS_Apple ON)
else()
    message(FATAL_ERROR "unrecognized value for `SYMPILER_BLAS_BACKEND` option: '${SYMPILER_BLAS_BACKEND}'")
endif()



set(SYMPILER_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(SYMPILER_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)

option(SYMPILER_USE_TBB "Use TBB for threading" OFF)
option(Sympiler_USE_SCOTCH  "Use SCOTCH" OFF)


# Defining paths
set(SPARSE_LIB_INC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/includes")
set(SPARSE_LIB_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(SPARSE_LIB_SRC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas")
AUX_SOURCE_DIRECTORY(${SPARSE_LIB_SRC} SPARSE_LIB_SRC_FILES)

set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/utils")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")



set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/include")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/src/")
set(LBC_LIB "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/lib/")
set(LBC_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/aggregation/")



#--------------------------------[ third party ]-------------------------------#
if(SYMPILER_USE_BLAS_MKL AND NOT TARGET MKL::MKL)
    #message(FATAL_ERROR "should not go there")
    set(MKL_INTERFACE "lp64")
    find_package(MKL CONFIG)

    if(MKL_FOUND)
        add_definitions(-DMKL)
        add_definitions(-DCBLAS)
    else()
        message(FATAL_ERROR "Intel MKL not found, please use OpenBLAS")
    endif()
elseif(SYMPILER_USE_BLAS_Apple)
    set(BLA_VENDOR Apple)
    find_package(BLAS)
    add_definitions(-DAPPLEBLAS)
    add_definitions(-DCBLAS)
elseif(SYMPILER_USE_BLAS_OpenBLAS)
    #set(openblas_WITHOUT_LAPACK ON) # for win support
    include(openblas)
    add_definitions(-DOPENBLAS)
    add_definitions(-DCBLAS)
endif()


# Find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    #message(FATAL_ERROR "don't use openmp")
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
    add_definitions(-DENABLE_OPENMP)
endif()

if(SYMPILER_USE_TBB)
    add_definitions(-DUSE_TBB)
    include(tbb)
    set(ALL_LIBS
            tbb_static)
    message(STATUS "TBB threading is enabled.")
endif()


include(metis)


if(Sympiler_USE_SCOTCH)
    #include(scotch)  # TODO: fixed this cmake part

    add_definitions(-DSCOTCH)
    link_directories(
            ${SCOTCH_ROOT}/build/lib/
    )
    set(ALL_LIBS
            ${ALL_LIBS}
            ${SCOTCH_ROOT}/build/lib/libscotch.a
            ${SCOTCH_ROOT}/build/lib/libscotcherr.a
            )
    include_directories(
            ${SCOTCH_ROOT}/build/src/include/
    )
    message(STATUS "SCOTCH is enabled from ${SCOTCH_ROOT}")
endif()



### Find PAPI for Profiling
#add_definitions( -DPROFILE)
#add_definitions( -DLOAD)
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
        #${MKL_INCLUDE_DIR}
        #${BLAS_INC}
        ${LBC_INC}
#        $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
        include
        include/sympiler/parsy/
        include/sympiler/linear_solver/
)


set(ALL_LIBS
        ${ALL_LIBS}
       # ${MKL_LIBRARIES}
        ${PAPI_LIB_NAME}
        ${PAPI_LIB_NAME2}
        metis::metis

        aggregation
        )



set(ALL_SOURCE_FILES
        ${BLAS_SRC_FILES}
        ${METIS_SOURCE_FILES})

set(CMAKE_CXX_STANDARD 11)
if(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11 -w -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")
endif()
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11 -w -m64")
endif()

Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")



if(TARGET metis::metis )
    add_definitions(-DMETIS)
    add_subdirectory(demo)
endif()

add_subdirectory(aggregation)
add_subdirectory(sparse_blas)
add_subdirectory(Catch_tests)
add_subdirectory(code_gen)

