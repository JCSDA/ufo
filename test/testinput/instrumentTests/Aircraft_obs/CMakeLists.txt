# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Aircraft observations (winds, temperature, some with moisture)


# H(x) test
ecbuild_add_test( TARGET  test_ufo_instrument_aircraft_gfs_HofX
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft_gfs_HofX.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LABELS  aircraft instrument HofX
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data )

