
if(TARGET tbb::tbb)
    return()
endif()


    message(STATUS "Third-party (internal): creating target 'tbb::tbb'")

    include(FetchContent)
    FetchContent_Declare(
            tbb
            GIT_REPOSITORY https://github.com/wjakob/tbb.git
            GIT_TAG        344fa84f34089681732a54f5def93a30a3056ab9
    )

    FetchContent_GetProperties(tbb)
    if(NOT tbb_POPULATED)
        FetchContent_Populate(tbb)
    endif()

    set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
    set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
#message("=== ${tbb_SOURCE_DIR}")
    add_subdirectory(${tbb_SOURCE_DIR} build EXCLUDE_FROM_ALL)
#    add_library(tbb::tbb ALIAS tbb)

include_directories(
        ${tbb_SOURCE_DIR}/include
)



#    set_target_properties(tbb PROPERTIES FOLDER third_party)



    # Install rules
#    set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME tbb)
#    install(DIRECTORY ${tbb_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#    install(TARGETS tbb EXPORT Tbb_Targets)
#    install(EXPORT TBB_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tbb NAMESPACE tbb::)


