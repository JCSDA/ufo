# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Radio Acoustic Sounding Systems (RASS)


# H(x) test
ecbuild_add_test( TARGET  test_ufo_instrument_rass_gfs_HofX
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/rass_gfs_HofX.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     1
                  LABELS  rass HofX instrument
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

