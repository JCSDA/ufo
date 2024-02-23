# (C) Copyright 2024 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for fov

# Test FOV footprint and antenna power for instrument amsua_n19
ufo_add_test( NAME    test_ufo_fov_amsua
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestFieldOfView.cc ../../../ufo/FieldOfView.h
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n19_fov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests fov
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
# Test FOV footprint and antenna power for instrument ssmis_f17
ufo_add_test( NAME    test_ufo_fov_ssmis
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestFieldOfView.cc ../../../ufo/FieldOfView.h
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_f17_fov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests fov
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
# Test helper methods to sample and reduce over fields of view
ufo_add_test( NAME    test_ufo_sample_and_reduce_over_fov
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestSampleAndReduceOverFieldOfView.cc
                      ../../../ufo/SampleAndReduceOverFieldOfView.h
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sample_and_reduce_over_fov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests fov
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
