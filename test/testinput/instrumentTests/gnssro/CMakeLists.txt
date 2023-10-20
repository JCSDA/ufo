# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for GNSSRO


# H(x) with QC

ufo_add_test( NAME    instrument_gnssro_gfs_qc
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/GnssroBnd_gfs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gnssro instrument HofX
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
