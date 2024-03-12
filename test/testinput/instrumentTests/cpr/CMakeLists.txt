# (C) Copyright 2020-2023 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for CPR


if( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
    # H(x) test
    ufo_add_test( NAME  instrument_cpr_cloudsat_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cpr_cloudsat_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm cpr cloudsat HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
