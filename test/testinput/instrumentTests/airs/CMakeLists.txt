# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for AIRS


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm airs aqua HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm airs aqua bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/airs/airs_aqua_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm airs aqua qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC (obs error input via namelist)
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_qc_obin
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/airs/airs_aqua_gfs_HofX_qc_obin.yaml"
        MPI     1
        LABELS  crtm airs aqua qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
