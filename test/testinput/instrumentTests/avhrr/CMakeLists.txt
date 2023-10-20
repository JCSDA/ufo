# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for AVHRR


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_avhrr_metop-a_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_metop-a_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm avhrr metop-a HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_avhrr_n18_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_n18_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm avhrr n18 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_avhrr_metop-a_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_metop-a_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm avhrr metop-a HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_avhrr_n18_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_n18_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm avhrr n18 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_avhrr_metop-a_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_metop-a_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm avhrr metop-a qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_avhrr_n18_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/avhrr_n18_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm avhrr n18 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
