# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for SEVIRI


if( crtm_FOUND )
    # H(x) test
    ecbuild_add_test( TARGET  test_ufo_instrument_seviri_m11_gfs_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX.yaml"
                      MPI     1
                      LABELS  seviri m11 HofX instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data get_crtm_coeffs )

    # H(x) with Bias Correction (bc)
    ecbuild_add_test( TARGET  test_ufo_instrument_seviri_m11_gfs_bc_HofX
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_HofX_bc.yaml"
                      MPI     1
                      LABELS  seviri m11 HofX_bc instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data get_crtm_coeffs )

    # H(x) with QC
    ecbuild_add_test( TARGET  test_ufo_instrument_seviri_m11_gfs_qc_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_m11_gfs_qc_filters.yaml"
                      MPI     1
                      LABELS  seviri m11 qc_filters instrument
                      ENVIRONMENT OOPS_TRAPFPE=1
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                      LIBS    ufo
                      TEST_DEPENDS ufo_get_ioda_test_data ufo_get_ufo_test_data ufo_get_crtm_test_data get_crtm_coeffs )
endif( crtm_FOUND )
