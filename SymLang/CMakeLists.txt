cmake_minimum_required(VERSION 3.15)

project(SymLang)

set(CMAKE_CXX_STANDARD 20)

set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src")
set(INC_DIR "${CMAKE_CURRENT_LIST_DIR}/include")

file(GLOB_RECURSE PROJECT_SOURCES "${SRC_DIR}" "${INC_DIR}" *.cpp *.h)

message("Configuring SymLang Library with Source: ${PROJECT_SOURCES}")

add_library(SymLang SHARED ${PROJECT_SOURCES})
target_include_directories(SymLang PRIVATE "${INC_DIR}")

add_executable(SymLang-Driver main.cpp)
target_link_libraries(SymLang-Driver SymLang)
target_include_directories(SymLang-Driver PUBLIC "${INC_DIR}")

# =============================================================================
#   Add custom target dox_doxygen to generate doxygen documentation for the
#   nasoq lib.
#
#   Dependencies: doxygen, graphviz.
# =============================================================================
find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    add_custom_target( SymLang_doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)