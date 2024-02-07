# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for ATMS


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_atms_n20_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_n20_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm atms n20 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_atms_npp_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_npp_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm atms npp HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_atms_n20_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_n20_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm atms n20 bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_atms_npp_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_npp_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm atms npp bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC in Skylab
    ufo_add_test( NAME    test_ufo_atms_skylab_filters
        TIER    1
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_skylab_filters.yaml"
        MPI     1
        LIBS    ufo
        LABELS  crtm filters
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )
 
    # H(x) with QC
    ufo_add_test( NAME  instrument_atms_n20_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/atms/atms_n20_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm atms n20 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_atms_npp_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/atms/atms_npp_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm atms npp qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
