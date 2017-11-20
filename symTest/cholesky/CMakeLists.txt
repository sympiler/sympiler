cmake_minimum_required(VERSION 3.1)
project(symChol)
#TODO use env instead




set(SOURCE_FILES main.cpp PB_Cholesky.h
        Etree.h def.h PostOrder.h DFS.h Reach.h Inspection_Prune.h ColumnCount.h
        Transpose.h LSparsity.h SparseUtils.h
        Inspection_BlockC.h)

include_directories(../../symGen ${MKL_INC} ${SUITE_INC})
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB})
add_executable(symChol ${SOURCE_FILES})
target_link_libraries(symChol
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
        #iomp5
        mkl_sequential
        mkl_intel_lp64
        mkl_core
        mkl_gnu_thread
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