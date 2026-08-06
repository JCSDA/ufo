# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for errors


ufo_add_test( NAME    metoffice_error_matrices
              ECBUILD
              SOURCES ../../../mains/TestMetOfficeRadianceErrorMatrices.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/metoffice_radiance_error_matrices.yaml"
              MPI     1
              LIBS    ufo
              LABELS  errors metoffice
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrorcrossvarcorr
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsError.x
              DEPENDS test_ObsError.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorcrossvarcorr.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrorcrossvarcov
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsError.x
              DEPENDS test_ObsError.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorcrossvarcov.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrorwithingroupcorr
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsError.x
              DEPENDS test_ObsError.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorwithingroupcorr.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrordiffusion
              ECBUILD
              SOURCES ../../../mains/TestObsErrorCovariance.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrordiffusion.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserror_dirac_diffusion
              ECBUILD
              SOURCES ../../../mains/TestObsErrorDiffusion.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/diracDiffusion.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Test the diagnostics application
ufo_add_test( NAME    obserrorwithingroupcorrdiags
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/ufo_obserrorcov_diags.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrorwithingroupcorrdiags.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrordiagonal
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsError.x
              DEPENDS test_ObsError.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrordiagonal.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    obserrordiagonal_inv_gamma
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsError.x
              DEPENDS test_ObsError.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserrordiagonal_inv_gamma.yaml"
              MPI     4
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    reconditioning
              ECBUILD
              SOURCES ../../../mains/TestReconditioning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reconditioning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    reconditioning_mpi
              ECBUILD
              SOURCES ../../../mains/TestReconditioning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reconditioning.yaml"
              MPI     2
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

ufo_add_test( NAME    reconditioning_within_group
              TIER    1
              ECBUILD
              SOURCES ../../../mains/TestReconditioning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reconditioning_withinGroup.yaml"
              MPI     1
              LIBS    ufo
              LABELS  errors
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
