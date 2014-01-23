cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(MBLIB Fortran)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

include(ConfigParentSettings)

if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules)
endif()

include_directories(
${CMAKE_Fortran_MODULE_DIRECTORY}
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

include(ConfigCompilerFlags)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/mb.F90
    ${PROJECT_SOURCE_DIR}/src/mb_precision.F90
    ${PROJECT_SOURCE_DIR}/src/mb_variables.F90
    ${PROJECT_SOURCE_DIR}/src/mb_io.F90
    ${PROJECT_SOURCE_DIR}/src/mb_utilities.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/integral_interface.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/dalton_interface.F
    )

add_library(mblib ${SOURCES})

install(TARGETS mblib ARCHIVE DESTINATION lib)
