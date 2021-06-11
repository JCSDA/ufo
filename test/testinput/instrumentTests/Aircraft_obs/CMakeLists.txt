# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Aircraft observations (winds, temperature, some with moisture)


# H(x) test
ufo_add_test(     NAME    instrument_aircraft_gfs_HofX
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft_gfs_HofX.yaml"
                  LIBS    ufo
                  LABELS  aircraft instrument HofX
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                  TEST_DEPENDS ufo_get_ufo_test_data )

# H(x) with QC
#ufo_add_test(     TARGET  test_ufo_aircraft_gfs_qc
#                  TIER 1
#                  ECBUILD
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
#                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft_gfs_qc.yaml"
#                  LABELS  aircraft instrument QC
#                  ENVIRONMENT OOPS_TRAPFPE=1
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
#                  LIBS    ufo
#                  TEST_DEPENDS test_ufo_aircraft_gfs_HofX )

