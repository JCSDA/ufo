# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for SEVIRI


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm seviri m11 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS get_crtm_coeffs )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm seviri m11 bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm seviri m11 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
