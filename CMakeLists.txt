cmake_minimum_required(VERSION 2.7)

set(MKL_INC "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/include")
set(MKL_LIB "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64")
set(MMM "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64/")

set(SUITE_INC "/home/kazem/programs/SuiteSparse/include")
set(SUITE_LIB "/home/kazem/programs/SuiteSparse/lib")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DCHOL -DMKL -DVERIFY -DCOLDCACHE  -Wall -std=gnu++11  -m64  ")

add_subdirectory(sympiler)
add_subdirectory(symTest)
add_subdirectory(libTest)
