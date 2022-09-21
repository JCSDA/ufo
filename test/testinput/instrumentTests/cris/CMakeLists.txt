# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for CrIS

if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX.yaml"
        MPI     1
        LABELS  cris n20 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX.yaml"
        MPI     1
        LABELS  cris npp HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  cris-fsr n20 bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ufo_add_test( NAME    test_ufo_cris-fsr_n20_geos_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_geos_qc_filters.yaml"
        MPI     1
        LABELS  cris-fsr n20 geos bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  cris-fsr npp bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    # H(x) with QC
    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  cris-fsr n20 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  cris-fsr npp qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo
        TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )
endif( crtm_FOUND )
