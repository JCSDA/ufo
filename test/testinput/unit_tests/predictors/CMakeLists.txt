# (C) Copyright 2022 UCAR.
# (C) Crown Copyright 2024, the Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for predictors

#
if( crtm_FOUND )
    ufo_add_test( NAME    test_ufo_amsua_crtm_bc_channelnobc_geos
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc_channelnobc_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_amsua_crtm_bc_geos
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_amsua_crtm_bc_obsoperator
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_amsua_crtm_bc_tlad
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc.yaml"
                  MPI     2
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    test_ufo_bias_coeff
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_coeff_increment
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsBiasIncrement.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_covariance_amsua
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsBiasCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_covariance_details
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsBiasCovarianceDetails.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_coeff_readwrite
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestUfoObsBias.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_coeff_increment_readwrite
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestUfoObsBiasIncrement.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_linear_op
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestLinearObsBiasOperator.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_linear_op.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_bias_single_channel
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_single_channel.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_constant_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/constant_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_legendre_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/legendre_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_interpolate_data_from_file_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/interpolate_data_from_file_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_orbital_angle_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/orbital_angle_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_satellite_selector_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satellite_selector_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( crtm_FOUND )
    ufo_add_test( NAME    test_ufo_ssmis_f17_gfs_backgroundcheck_bc
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CRTM_CMAKE_TESTINPUT}/unit_tests/predictors/ssmis_f17_gfs_backgroundcheck_bc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_ssmis_f17_gfs_backgroundcheck_nbc
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_f17_gfs_backgroundcheck_nbc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    test_ufo_thickness_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/thickness_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    test_ufo_read_bias_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/read_bias_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    test_ufo_obs_metadata_predictor
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obs_metadata_predictor.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
