#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET scotch::scotch)
    return()
endif()

    message(STATUS "Third-party (external): creating target 'scotch::scotch'")

    include(FetchContent)
    FetchContent_Declare(
            scotch
            GIT_REPOSITORY https://gitlab.inria.fr/scotch/scotch.git# git@git.corp.adobe.com:CTL-third-party/metis
            GIT_TAG       10931343abc47eb43565087697538ac3919289b9
    )

    FetchContent_GetProperties(scotch)
    if(NOT scotch_POPULATED)
        FetchContent_Populate(scotch)
    endif()

    list(APPEND CMAKE_MODULE_PATH ${scotch_SOURCE_DIR}/cmake/Modules)
    set(BUILD_PTSCOTCH OFF CACHE BOOL " " FORCE)
    set(BISON_EXECUTABLE /usr/local/Cellar/bison/3.8.2/bin/bison)
    add_subdirectory(${scotch_SOURCE_DIR} ${scotch_SOURCE_DIR}/build EXCLUDE_FROM_ALL)



# Create metis target
#    file(GLOB INC_FILES
#            "${scotch_SOURCE_DIR}/src/*.h"
#            )
#    file(GLOB SRC_FILES
#            "${scotch_SOURCE_DIR}/src/*.c"
#            )
##    list(REMOVE_ITEM SRC_FILES "${metis_SOURCE_DIR}/GKlib/gkregex.c")
#
#    add_library(scotch SHARED ${INC_FILES} ${SRC_FILES})
#    add_library(scotch::scotch ALIAS scotch)
#
##[[
#    if(MSVC)
#        target_compile_definitions(metis PUBLIC USE_GKREGEX)
#        target_compile_definitions(metis PUBLIC "__thread=__declspec(thread)")
#    endif()
#]]
#
#    target_include_directories(scotch PRIVATE "${scotch_SOURCE_DIR}/src/esmumps")
#    target_include_directories(scotch PRIVATE "${scotch_SOURCE_DIR}/src/libscotch")
#    target_include_directories(scotch PRIVATE "${scotch_SOURCE_DIR}/src/misc")
#    target_include_directories(scotch PRIVATE "${scotch_SOURCE_DIR}/src/scotch")
#
#    include(GNUInstallDirs)
#    target_include_directories(metis SYSTEM PUBLIC
#            "$<BUILD_INTERFACE:${metis_SOURCE_DIR}/include>"
#            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
#            )

    set_target_properties(scotch PROPERTIES FOLDER third_party)

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
            "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_compile_options(scotch PRIVATE
                "-Wno-unused-variable"
                "-Wno-sometimes-uninitialized"
                "-Wno-absolute-value"
                "-Wno-shadow"
                )
    elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        target_compile_options(scotch PRIVATE
                "-w" # Disallow all warnings from metis.
                )
    elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
        target_compile_options(scotch PRIVATE
                "/w" # Disable all warnings from metis!
                )
    endif()

    # Install rules
#[[    set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME scotch)
    set(METIS_FOUND TRUE)
    install(DIRECTORY ${scotch_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(TARGETS scotch EXPORT Scotch_Targets)
    install(EXPORT Scotch_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch NAMESPACE scotch::)]]
