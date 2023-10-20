# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for utils

#
ufo_add_test( NAME    test_ufo_utils_metoffice_sort
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestMetOfficeSort.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/utils_metoffice_sort.yaml"
              MPI     1
              LIBS    ufo
              LABELS  utils metoffice
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    test_ufo_utils_variable_name_map
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variable_name_map.yaml"
              MPI     1
              LIBS    ufo
              LABELS  utils variablenamemap
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
