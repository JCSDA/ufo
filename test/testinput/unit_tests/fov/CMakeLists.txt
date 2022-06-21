# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for fov

#
ufo_add_test( NAME    test_ufo_fov_amsua
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestFieldOfView.cc ../../../ufo/FieldOfView.h
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n19_fov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests fov
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_fov_ssmis
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestFieldOfView.cc ../../../ufo/FieldOfView.h
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_f17_fov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests fov
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
