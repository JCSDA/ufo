# (C) Copyright 2020-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for IASI

if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_iasi_metop-a_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-a_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm iasi metop-a HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_iasi_metop-b_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-b_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm iasi metop-b HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_iasi_metop-a_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-a_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm iasi metop-a bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_iasi_metop-b_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-b_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm iasi metop-b bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_iasi_metop-a_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-a_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm iasi metop-a qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_iasi_metop-b_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_metop-b_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm iasi metop-b qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
