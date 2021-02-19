# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Scatterometer surface winds


# H(x) test
ecbuild_add_test( TARGET  test_ufo_scatwinds_gfs_HofX
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwinds_gfs_HofX.yaml"
                  LABELS  scatwinds HofX instrument
                  ENVIRONMENT OOPS_TRAPFPE=1
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data )

# H(x) with QC
#ecbuild_add_test( TARGET  test_ufo_scatwinds_gfs_qc
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
#                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwinds_gfs_qc.yaml"
#                  LABELS  scatwinds QC instrument
#                  ENVIRONMENT OOPS_TRAPFPE=1
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
#                  LIBS    ufo
#                  TEST_DEPENDS test_ufo_scatwinds_HofX )

