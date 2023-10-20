# (C) Copyright 2020-2023 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for rawinsonde and dropsonde observations

# H(x) test
ufo_add_test( NAME    instrument_sonde_gfs_HofX
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sonde_gfs_HofX.yaml"
              MPI     1
              LIBS    ufo
              LABELS  sonde HofX instrument
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# H(x) with QC
ufo_add_test( NAME    instrument_sonde_gfs_qc
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sonde_gfs_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  sonde QC instrument
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# H(x) with QC (GMAO setup)
ufo_add_test( NAME    instrument_sonde_geos_qc
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sonde_geos_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  sonde GEOS UV QC instrument
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
