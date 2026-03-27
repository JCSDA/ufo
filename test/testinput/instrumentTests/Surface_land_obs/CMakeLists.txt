# (C) Copyright 2020-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for Surface-based land observations (METAR, SYNOP, mesonet)


# H(x) with QC
ufo_add_test( NAME    instrument_sfcLand_geos_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sfcLand_geos_qc.yaml"
              MPI     1
              LABELS  sfcLand QC instrument
              LIBS    ufo
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

