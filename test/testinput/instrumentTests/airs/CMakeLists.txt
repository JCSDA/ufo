# (C) Copyright 2020-2026 UCAR.
#

# Instrument tests for AIRS


if( crtm_FOUND )
    # H(x) test
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX.yaml"
        MPI     1
        LABELS  crtm airs aqua HofX instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with Bias Correction (bc)
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_bc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX_bc.yaml"
        MPI     1
        LABELS  crtm airs aqua bc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_qc
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX_qc.yaml"
        MPI     1
        LABELS  crtm airs aqua qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )

    # H(x) with QC (obs error input via namelist)
    ufo_add_test( NAME  instrument_airs_aqua_gfs_HofX_qc_obin
        ECBUILD
        COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
        ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_aqua_gfs_HofX_qc_obin.yaml"
        MPI     1
        LABELS  crtm airs aqua qc instrument tier2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
        LIBS    ufo )
endif( crtm_FOUND )
