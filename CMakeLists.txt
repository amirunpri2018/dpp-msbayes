cmake_minimum_required (VERSION 2.6)

# get local cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/abacus/cmake/Modules")

project("dpp-msbayes" C)

enable_testing()

set(PROJECT_MAJOR_VERSION 0)
set(PROJECT_MINOR_VERSION 3)
set(PROJECT_PATCH_VERSION 0)
set(PROJECT_VERSION
    ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION})

set (STATIC_LINKING OFF CACHE BOOL "build with static links")

# configure a header file to pass CMake settings to source code
configure_file (
    "${PROJECT_SOURCE_DIR}/src/dpp_msbayes_config.h.in"
    "${PROJECT_SOURCE_DIR}/src/dpp_msbayes_config.h"
    )

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(src)
add_subdirectory(src/abacus)
add_subdirectory(src/old-model)
add_subdirectory(src/hashtab)

find_package (GSL)
if (NOT GSL_FOUND)
    message (FATAL_ERROR
        "GSL was not found. Please install GSL and try again."
        )
endif(NOT GSL_FOUND)

