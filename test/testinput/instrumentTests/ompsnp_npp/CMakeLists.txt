# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Ozone measurements/profiles from NPP satellite.


# H(x) test
ufo_add_test( NAME    instrument_ompsnp_npp_gfs_HofX
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_gfs_HofX.yaml"
              LIBS    ufo
              LABELS  ozone instrument HofX
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# H(x) with QC
ufo_add_test( NAME    instrument_ompsnp_npp_gfs_qc
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_gfs_qc.yaml"
              LIBS    ufo
              LABELS  ozone instrument QC
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
