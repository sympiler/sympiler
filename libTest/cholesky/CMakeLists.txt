cmake_minimum_required(VERSION 3.1)
project(CholOpenMP)
#TODO use env instead
set(MKL_INC "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/include")
set(MKL_LIB "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64")
set(MMM "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64/")

set(SUITE_INC "/home/kazem/programs/SuiteSparse/include")
set(SUITE_LIB "/home/kazem/programs/SuiteSparse/lib")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL -Wall -std=gnu++11  -m64  ")

set(SOURCE_FILES main.cpp PB_Cholesky.h MyBLAS.h Ordering.h Inspection_Block.h Etree.h def.h PostOrder.h DFS.h Reach.h Inspection_Prune.h ColumnCount.h Transpose.h Util.h LSparsity.h SparseUtils.h colCounts.h Inspection_BlockC.h AMD.h)

include_directories(${MKL_INC} ${SUITE_INC})
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB})
add_executable(CholOpenMP ${SOURCE_FILES})
target_link_libraries(CholOpenMP
        #mkl_blacs_intelmpi_ilp64
        #mkl_blacs_intelmpi_lp64
        #mkl_blacs_openmpi_ilp64
    #    mkl_blacs_openmpi_lp64
#        mkl_blacs_sgimpt_ilp64
   #     mkl_blacs_sgimpt_lp64
  #      mkl_def
        #mkl_gf_lp64

 #       mkl_blas95_ilp64
 #       mkl_cdft_core

        #mkl_intel_ilp64
        #gomp
        iomp5
        mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        #mkl_blacs_openmpi_lp64
        #iomp
#[[        mkl_blacs_openmpi_lp64
        mkl_blacs_openmpi_ilp64
        mkl_pgi_thread]]
        #mkl_tbb_thread
        #Libs from CHOLMOD for ordering
        amd
        cholmod
        suitesparseconfig


        #openblas
        )