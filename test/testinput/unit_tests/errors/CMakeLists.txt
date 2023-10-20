# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for errors

#
ufo_add_test( NAME    test_ufo_metoffice_error_matrices
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestMetOfficeRadianceErrorMatrices.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/metoffice_radiance_error_matrices.yaml"
              MPI     1
              LIBS    ufo
              LABELS  errors metoffice
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_obserrorcrossvarcorr
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsErrorCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorcrossvarcorr.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_obserrorcrossvarcov
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsErrorCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorcrossvarcov.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    test_ufo_obserrorwithingroupcorr
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsErrorCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorwithingroupcorr.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Test the diagnostics application
ufo_add_test( NAME    test_ufo_obserrorwithingroupcorrdiags
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/ufo_obserrorcov_diags.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorwithingroupcorrdiags.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    test_ufo_obserrordiagonal
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestObsErrorCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrordiagonal.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    test_ufo_reconditioning
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestReconditioning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reconditioning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    test_ufo_reconditioning_mpi
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestReconditioning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reconditioning.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
