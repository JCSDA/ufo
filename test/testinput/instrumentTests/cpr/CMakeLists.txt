# (C) Copyright 2020-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for CPR


if( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
    # H(x) test
    ufo_add_test( NAME  instrument_cpr_cloudsat_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cpr_cloudsat_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm cpr cloudsat HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
    # TL/AD test
    ufo_add_test( NAME  instrument_cpr_cloudsat_gfs_TLAD
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cpr_cloudsat_gfs_TLAD.yaml"
        MPI     1
        LABELS  crtm cpr cloudsat linopr instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
    # ObsDiagnostics test: reflectivity Jacobians (active-sensor K-matrix path)
    ufo_add_test( NAME  instrument_cpr_cloudsat_gfs_obsdiag_jacobian
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cpr_cloudsat_gfs_obsdiag_jacobian.yaml"
        MPI     1
        LABELS  crtm cpr cloudsat obsdiag instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
