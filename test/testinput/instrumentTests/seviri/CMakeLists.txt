# (C) Copyright 2020-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for SEVIRI


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm seviri m11 HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS get_crtm_coeffs )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm seviri m11 bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_seviri_m11_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm seviri m11 qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
