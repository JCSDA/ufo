# (C) Copyright 2020-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for CrIS

if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm cris n20 HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm cris npp HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm cris-fsr n20 bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm cris-fsr npp bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_cris-fsr_n20_geos_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_geos_HofX_qc.yaml"
        MPI     1
        LABELS  crtm cris-fsr n20 qc instrument geos tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_cris-fsr_n20_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_n20_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm cris-fsr n20 qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_cris-fsr_npp_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris-fsr_npp_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm cris-fsr npp qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
