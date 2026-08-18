# (C) Copyright 2026 UCAR.
# (C) Crown Copyright 2024, the Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for predictors

#
if( crtm_FOUND )
    ufo_add_test( NAME    amsua_crtm_bc_channelnobc_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc_channelnobc_geos.yaml"
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_crtm_bc_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc_geos.yaml"
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_crtm_bc_obsoperator
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc.yaml"
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_crtm_bc_tlad
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_bc.yaml"
                  MPI     2
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    bias_coeff
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_coeff_increment
              ECBUILD
              SOURCES ../../../mains/TestObsBiasIncrement.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_covariance_amsua
              ECBUILD
              SOURCES ../../../mains/TestObsBiasCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_covariance_aircraft
              ECBUILD
              COMMAND ufo_bias_covariance_amsua
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov_byrecord.yaml"
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_covariance_details
              ECBUILD
              SOURCES ../../../mains/TestObsBiasCovarianceDetails.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_covariance_details_byrecord
              ECBUILD
              COMMAND ufo_bias_covariance_details
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov_byrecord.yaml"
              LABELS  unit_tests predictors
              TEST_DEPENDS ufo_bias_covariance_aircraft
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_covariance_details_byrecord_2
              ECBUILD
              COMMAND ufo_bias_covariance_details
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov_byrecord.yaml"
              MPI     2
              LABELS  unit_tests predictors
              TEST_DEPENDS ufo_bias_covariance_details_byrecord
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_coeff_readwrite
              ECBUILD
              SOURCES ../../../mains/TestUfoObsBias.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_coeff_readwrite_byrecord
              ECBUILD
              COMMAND ufo_bias_coeff_readwrite
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff_cov_byrecord.yaml"
              MPI     2
              LABELS  unit_tests predictors
              TEST_DEPENDS ufo_bias_covariance_details_byrecord_2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_coeff_increment_readwrite
              ECBUILD
              SOURCES ../../../mains/TestUfoObsBiasIncrement.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_coeff.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_linear_op
              ECBUILD
              SOURCES ../../../mains/TestLinearObsBiasOperator.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_linear_op.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_linear_op_tlad
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_linear_op_tlad.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
# Test for application of bias prediction by record
ufo_add_test( NAME    aircraft_bc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft_bc.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    aircraft_bc_mpi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft_bc_mpi.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    bias_single_channel
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/bias_single_channel.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    constant_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/constant_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    cosine_of_latitude_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cosine_of_latitude_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    legendre_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/legendre_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    interpolate_data_from_file_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/interpolate_data_from_file_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    orbital_angle_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/orbital_angle_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    satellite_selector_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satellite_selector_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    sine_of_latitude_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sine_of_latitude_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( crtm_FOUND )
    ufo_add_test( NAME    ssmis_f17_gfs_backgroundcheck_bc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_f17_gfs_backgroundcheck_bc.yaml"
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    ssmis_f17_gfs_backgroundcheck_nbc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_f17_gfs_backgroundcheck_nbc.yaml"
                  LIBS    ufo
                  LABELS  crtm unit_tests predictors
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    SymmCldImpact_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/SymmCldImpact_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    thickness_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/thickness_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    read_bias_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/read_bias_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    obs_metadata_predictor
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obs_metadata_predictor.yaml"
              LIBS    ufo
              LABELS  unit_tests predictors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
