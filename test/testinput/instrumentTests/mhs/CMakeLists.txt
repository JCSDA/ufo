# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for MHS


# H(x) test
if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_instrument_mhs_metop-b_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_metop-b_gfs_HofX.yaml"
                      MPI     1
                      LABELS  mhs metop-b HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_mhs_metop-c_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_metop-c_gfs_HofX.yaml"
                      MPI     1
                      LABELS  mhs metop-c HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_instrument_mhs_n19_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_n19_gfs_HofX.yaml"
                      MPI     1
                      LABELS  mhs n19 HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data )
endif( crtm_FOUND )