cmake_minimum_required( VERSION 2.5)

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



set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/lbc/includes")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/lbc")
AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)


set(BLAS_INC "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas/includes")
set(BLAS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas")
#AUX_SOURCE_DIRECTORY(${BLAS_SRC} BLAS_SRC_FILES)
set(BLAS_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/dense_blas/BLAS.cpp" )

### for intel omp this should be -openmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_definitions( -DPROFILE)
#add_definitions( -DLOAD)


### Find Metis if available

if(EXISTS "$ENV{METISINC}" AND EXISTS "$ENV{METISLIB}") 
    set(METIS_INC "$ENV{METISINC}")
    set(METIS_LIB "$ENV{METISLIB}")
    add_definitions(-DMETIS)
    set(METIS_LIB_NAME "metis")
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
    set(METIS 1)
    message(STATUS "METIS found")
elseif(EXISTS "$ENV{METISROOT}")
    set(METIS_INC "$ENV{METISROOT}/include")
    set(METIS_LIB "$ENV{METISROOT}/libmetis")
    if(NOT EXISTS ${METIS_LIB})
        set(METIS_LIB "$ENV{METISROOT}/lib")
    endif()
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
    add_definitions(-DMETIS)
    set(METIS 1)
    set(METIS_LIB_NAME "metis")
    message(STATUS "METIS found")
else()
    message(STATUS "METIS not found, building without metis.")
    set(METIS_LIB_NAME "")
    set(METIS 0)
endif()

### Find Intel MKL if available
if(EXISTS "$ENV{MKLROOT}" OR EXISTS "$ENV{MKL_ROOT}")
    if(EXISTS "$ENV{MKLROOT}")
        set(MKLPATH "$ENV{MKLROOT}")
    else()
        set(MKLPATH "$ENV{MKL_ROOT}")
    endif()

    if(APPLE)
        #find_package(OpenMP REQUIRED)
        set(MKL_INC "${MKLPATH}/include")
        set(MKL_LIB "${MKLPATH}/lib")
        set(MKL_LIB2 "${MKLPATH}/../compiler/lib")
    endif()
    if(UNIX AND NOT APPLE)
        set(MKL_INC "${MKLPATH}/include")
        set(MKL_LIB "${MKLPATH}/lib/intel64")
        # for iomp lib
        set(MKL_LIB2 "${MKLPATH}/../compiler/lib/intel64")
    endif()
    add_definitions(-DMKL)
    set(MKL 1)
    set(MKL_LIB_NAME "mkl_intel_lp64"
            "pthread"
            "mkl_intel_thread"
            "iomp5"
            "mkl_sequential"
            "mkl_core"
            #"mkl_sequential"
                     )
#                     "mkl_intel_ilp64")
    message(STATUS "Intel MKL found")
else()
    message(STATUS "Intel MKL not found, building without MKL.")
    set(MKL_LIB_NAME "")
    set(MKL 0)
endif()

### Find PAPI for Profiling
### Can be cloned directly from https://bitbucket.org/icl/papi.git
if(EXISTS "$ENV{PAPIROOT}")
    if(EXISTS "$ENV{PAPIROOT}/src")
        set(PAPI_LIB "$ENV{PAPIROOT}/src")
        set(PAPI_LIB2 "${PAPI_LIB}/libpfm4/lib/")
    elseif(EXISTS "$ENV{PAPIROOT}/lib")
        set(PAPI_LIB "$ENV{PAPIROOT}/lib")
    endif()
    if(EXISTS "$ENV{PAPIROOT}/include")
        set(PAPI_INC "$ENV{PAPIROOT}/include")
    else()
        set(PAPI_INC "$ENV{PAPIROOT}/src")
    endif()

    add_definitions(-DPAPI)
    set(PAPI_LIB_NAME "papi")
    set(PAPI_LIB_NAME2 "pfm")
    message(STATUS "PAPI found")
else()
    message(STATUS "PAPI not found, building without PAPI.")
    set(PAPI_LIB_NAME "")
endif()

set(ALL_LIBS
        ${MKL_LIB_NAME}
        ${PAPI_LIB_NAME}
        ${PAPI_LIB_NAME2}
        ${METIS_LIB_NAME})

set(ALL_SOURCE_FILES
        ${SPARSE_LIB_SRC_FILES}
        ${SPARSE_UTIL_SRC_FILES}
        ${SPARSE_EXTN_SRC_FILES}
        ${CODE_FUSION_SRC_FILES}
        ${DATA_FUSION_SRC_FILES}
    #    ${APPS_SRC_FILES} # TODO:some of them are mixed with blas, should be cleaned
        ${LBC_SRC_FILES}
        ${BLAS_SRC_FILES}
        ${METIS_SOURCE_FILES})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -m64")
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
