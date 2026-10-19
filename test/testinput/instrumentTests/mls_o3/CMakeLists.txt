# (C) Copyright 2020-2026 UCAR.
#

# Instrument tests for MLS O3


# H(x) test for MLS O3 v5 NRT
ufo_add_test( NAME  instrument_mls_o3_nrt_v5_hofx
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_nrt_v5.yaml"
    MPI     1
    LABELS  mls o3 nrt v5 HofX instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )

# Obs filter test for MLS O3 v5 NRT
ufo_add_test( NAME  instrument_mls_o3_nrt_v5_filters
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_nrt_v5_qc.yaml"
    MPI     1
    LABELS  mls o3 nrt v5 filters instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )

# TL/AD test for MLS O3 v5 NRT
ufo_add_test( NAME  instrument_mls_o3_nrt_v5_linopr
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_nrt_v5.yaml"
    MPI     1
    LABELS  mls o3 nrt v5 linopr instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )

# H(x) test for MLS O3 v6 RES
ufo_add_test( NAME  instrument_mls_o3_res_v6_hofx
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_res_v6.yaml"
    MPI     1
    LABELS  mls o3 res v6 HofX instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )

# Obs filter test for MLS O3 v6 RES
ufo_add_test( NAME  instrument_mls_o3_res_v6_filters
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_res_v6_qc.yaml"
    MPI     1
    LABELS  mls o3 res v6 filters instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )

# TL/AD test for MLS O3 v6 RES
ufo_add_test( NAME  instrument_mls_o3_res_v6_linopr
    ECBUILD
    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
    ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mls_o3_res_v6.yaml"
    MPI     1
    LABELS  mls o3 res v6 linopr instrument tier1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
    LIBS    ufo )
