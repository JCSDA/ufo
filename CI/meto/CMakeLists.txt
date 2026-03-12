# (C) Crown Copyright 2025 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(
  metoffice-jedi-ci
  VERSION 1.0.0
  LANGUAGES C CXX Fortran
  DESCRIPTION "Met Office JEDI Integration Test"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_MPI ON CACHE BOOL "Compile with MPI")
set(ENABLE_OMP ON CACHE BOOL "Compile with OpenMP")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jedicmake")
if(NOT DEFINED jedicmake_DIR)
  set(jedicmake_DIR "${CMAKE_BINARY_DIR}/jedicmake")
endif()

set(
  include_names
  "$ENV{REPO}" "opsinputs" "orca-jedi" "ops-um-jedi" "jjtests" "lfric-jedi")
if(DEFINED ENV{REPO_TEST_DATA})
  list(APPEND include_names "$ENV{REPO_TEST_DATA}")
endif()
foreach(
  name
  "oops"
  "jedi-model-data"
  "vader"
  "saber"
  "ioda-data"
  "ioda"
  "ufo-data"
  "ufo"
  "opsinputs"
  "monio"
  "nemo-feedback"
  "orca-jedi"
  "ops-um-jedi"
  "jjtests"
  "lfric-jedi"
)
  list(FIND include_names ${name} IS_INCLUDE_NAME)
  if(${IS_INCLUDE_NAME} GREATER -1)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/${name}")
  else()
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/${name}" EXCLUDE_FROM_ALL)
  endif()
endforeach()

enable_testing()
