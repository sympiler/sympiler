cmake_minimum_required( VERSION 2.5)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DownloadProject)

set(SYMPILER_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(SYMPILER_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)

# Defining paths
set(SPARSE_LIB_INC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/includes")
set(SPARSE_LIB_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(SPARSE_LIB_SRC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas")
AUX_SOURCE_DIRECTORY(${SPARSE_LIB_SRC} SPARSE_LIB_SRC_FILES)

set(SPARSE_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/utils/includes")
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

set(SPARSE_EXTN_INC "${SPARSE_UTIL_SRC}/external/includes")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")
AUX_SOURCE_DIRECTORY(${SPARSE_EXTN_SRC} SPARSE_EXTN_SRC_FILES)



set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/lbc/src/includes")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/lbc/src/")
AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)


set(BLAS_INC "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas/includes")
set(BLAS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas")
#AUX_SOURCE_DIRECTORY(${BLAS_SRC} BLAS_SRC_FILES)
set(BLAS_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas/BLAS.cpp" )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)



# Find OpenMP
if(APPLE)
    #    SET(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
    #    SET(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
    #    SET(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
    #    SET(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

#find_package(OPENMP)
find_package(OpenMP)

if (OPENMP_FOUND)
    #include_directories("${OPENMP_INCLUDES}")
    #include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
    link_directories("${OPENMP_LIBRARIES}")
    #link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DOPENMP)
endif()

#add_definitions( -DPROFILE)
#add_definitions( -DLOAD)


### Find Metis if available

#if(EXISTS "$ENV{METISINC}" AND EXISTS "$ENV{METISLIB}")
#    set(METIS_INC "$ENV{METISINC}")
#    set(METIS_LIB "$ENV{METISLIB}")
#    add_definitions(-DMETIS)
#    set(METIS_LIB_NAME "metis")
#    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
#    set(METIS 1)
#    message(STATUS "METIS found")
#elseif(EXISTS "$ENV{METISROOT}")
#    set(METIS_INC "$ENV{METISROOT}/include")
#    set(METIS_LIB "$ENV{METISROOT}/libmetis")
#    if(NOT EXISTS ${METIS_LIB})
#        set(METIS_LIB "$ENV{METISROOT}/lib")
#    endif()
#    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
#    add_definitions(-DMETIS)
#    set(METIS 1)
#    set(METIS_LIB_NAME "metis")
#    message(STATUS "METIS found")
#else()
#    message(STATUS "METIS not found, building without metis.")
#    set(METIS_LIB_NAME "")
#    set(METIS 0)
#endif()
#find_package(METIS OPTIONAL_COMPONENTS)


find_package(METIS OPTIONAL_COMPONENTS)
# set(METIS_LIBRARY "")
 set(METIS_INCLUDE_DIR "")
#ELSE()
#    set(ignoreMe "${METIS_ROOT_PATH}")
#ENDIF()

if(NOT METIS_VERSION_OK)
    #    find_program(MAKE_EXE NAMES gmake nmake make)
    #    set(metis_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/metis")
    #    set(metis_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${metis_INSTALL_DIR})
    #set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_LIST_DIR}/external)
    # set(GKLIB_PATH ${SYMPILER_EXTERNAL}/metis/metis-src/GKLIB/ CACHE PATH "d")
    download_project(PROJ metis
            GIT_REPOSITORY https://github.com/cheshmi/METIS.git  # https://github.com/KarypisLab/METIS.git
            # BUILD_IN_SOURCE 1
            #PREFIX  ${metis_BINARY_DIR}
            #     URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
            PREFIX  ${SYMPILER_EXTERNAL}/metis/
            #SOURCE_DIR ${SYMPILER_EXTERNAL}/metis/metis-src/
            BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/metis/
            UPDATE_DISCONNECTED 1
            #CMAKE_ARGS  -DGKLIB_PATH=${metis_SOURCE_DIR}/GKlib/
            # BUILD_COMMAND ${MAKE_EXE}
            )
    add_subdirectory(${metis_SOURCE_DIR})
    set(METIS_INCLUDE_DIR ${metis_SOURCE_DIR}/include)
    set(METIS_INCLUDES ${metis_SOURCE_DIR}/include/)
    set(METIS_LIBRARY ${SYMPILER_ROOT_DIR}/lib/)
    set(METIS_LIBRARIES ${SYMPILER_ROOT_DIR}/lib/)
endif()


find_package(MKL REQUIRED)
add_definitions(-DMKL) # add openblas

### Find Intel MKL if available
#if(EXISTS "$ENV{MKLROOT}" OR EXISTS "$ENV{MKL_ROOT}")
#    if(EXISTS "$ENV{MKLROOT}")
#        set(MKLPATH "$ENV{MKLROOT}")
#    else()
#        set(MKLPATH "$ENV{MKL_ROOT}")
#    endif()
#
#    if(APPLE)
#        #find_package(OpenMP REQUIRED)
#        set(MKL_INC "${MKLPATH}/include")
#        set(MKL_LIB "${MKLPATH}/lib")
#        set(MKL_LIB2 "${MKLPATH}/../compiler/lib")
#    endif()
#    if(UNIX AND NOT APPLE)
#        set(MKL_INC "${MKLPATH}/include")
#        set(MKL_LIB "${MKLPATH}/lib/intel64")
#        # for iomp lib
#        set(MKL_LIB2 "${MKLPATH}/../compiler/lib/intel64")
#    endif()
#    add_definitions(-DMKL)
#    set(MKL 1)
#    set(MKL_LIB_NAME "mkl_intel_lp64"
#            "pthread"
#            "mkl_intel_thread"
#            "iomp5"
#            "mkl_sequential"
#            "mkl_core"
#            #"mkl_sequential"
#                     )
##                     "mkl_intel_ilp64")
#    message(STATUS "Intel MKL found")
#else()
#    message(STATUS "Intel MKL not found, building without MKL.")
#    set(MKL_LIB_NAME "")
#    set(MKL 0)
#endif()

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
        ${METIS_INCLUDE_DIR}
        ${METIS_INCLUDES}
)


set(ALL_LIBS
        ${MKL_LIBRARIES}
        ${PAPI_LIB_NAME}
        ${PAPI_LIB_NAME2}
        ${METIS_LIB_NAME}
        )


set(ALL_SOURCE_FILES
        ${SPARSE_LIB_SRC_FILES}
        ${SPARSE_UTIL_SRC_FILES}
        ${SPARSE_EXTN_SRC_FILES}
    #    ${APPS_SRC_FILES} # TODO:some of them are mixed with blas, should be cleaned
        ${LBC_SRC_FILES}
        ${BLAS_SRC_FILES}
        ${METIS_SOURCE_FILES})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11 -w -m64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")




if(EXISTS "${METIS_INCLUDES}")
    add_definitions(-DMETIS)
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
    add_subdirectory(demo)
endif()

add_subdirectory(sparse_blas)
add_subdirectory(lbc)
add_subdirectory(Catch_tests)
add_subdirectory(code_gen)
add_subdirectory(utils)
add_subdirectory(dense_blas)
