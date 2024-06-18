#
# (C) Copyright 2024 MetOffice Crown Copyright
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
#
# MetOffice-JEDI Integration Test against UFO
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ufo-ci VERSION 1.0.0 LANGUAGES C CXX Fortran)

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

set(ENABLE_MPI ON CACHE BOOL "Compile with MPI")
set(ENABLE_OMP ON CACHE BOOL "Compile with OpenMP")
option(ENABLE_IODA_DATA "Obtain ioda test data from ioda-data repository (vs tarball)" ON)
option(ENABLE_UFO_DATA "Obtain ufo test data from ufo-data repository (vs tarball)" ON)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jedicmake")
if(NOT DEFINED jedicmake_DIR)
  set(jedicmake_DIR "${CMAKE_BINARY_DIR}/jedicmake")
endif()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/oops" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/vader" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/saber" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ioda" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ufo-data" )
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ufo" )
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/opsinputs" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/monio" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/nemo-feedback" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/orca-jedi")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ops-um-jedi")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jjtests")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lfric-jedi")

enable_testing()
