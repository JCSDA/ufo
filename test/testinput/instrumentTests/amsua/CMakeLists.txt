# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for AMSU-A


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_amsua_aqua_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_aqua_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua aqua HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-a_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-a_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua metop-a HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-b_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-b_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua metop-b HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-c_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-c_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua metop-c HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n15_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n15_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua n15 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n18_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n18_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua n18 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n19_gfs_HofX
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n19_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm amsua n19 HofX instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_amsua_aqua_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_aqua_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua aqua bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-a_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-a_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua metop-a bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-b_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-b_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua metop-b bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-c_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-c_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua metop-c bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n15_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n15_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua n15 bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

#   ufo_add_test( NAME  instrument_amsua_n18_gfs_HofX_bc
#       TIER    2
#       ECBUILD
#       COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
#       ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n18_gfs_HofX_bc.yaml"
#       MPI     1
#       LABELS  crtm amsua n18 bc instrument
#       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
#       LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n19_gfs_HofX_bc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n19_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm amsua n19 bc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_amsua_aqua_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_aqua_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua aqua qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-a_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/amsua/amsua_metop-a_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua metop-a qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-b_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-b_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua metop-b qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_metop-c_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-c_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua metop-c qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n15_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/amsua/amsua_n15_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua n15 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n18_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/amsua/amsua_n18_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua n18 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    ufo_add_test( NAME  instrument_amsua_n19_gfs_HofX_qc
        TIER    2
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CRTM_CMAKE_TESTINPUT}/instrumentTests/amsua/amsua_n19_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm amsua n19 qc instrument
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
