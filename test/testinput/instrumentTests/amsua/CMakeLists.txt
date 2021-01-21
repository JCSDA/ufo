# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for AMSU-A


# H(x) test
if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_aqua_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_aqua_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua aqua HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_metop-a_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-a_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua metop-a HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_metop-b_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-b_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua metop-b HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_metop-c_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_metop-c_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua metop-c HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_n15_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n15_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua n15 HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_n18_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n18_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua n18 HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_amsua_n19_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_n19_gfs_HofX.yaml"
                      MPI     1
                      LABELS  amsua n19 HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )
endif( crtm_FOUND )
