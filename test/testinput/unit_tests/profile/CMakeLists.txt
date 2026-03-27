# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for profile

#
ufo_add_test( NAME    conventionalprofileprocessing_average_pressure_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_pressure_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_pressure_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_pressure_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_relativehumidity_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_relativehumidity_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_relativehumidity_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_relativehumidity_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_single_level_profile_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_single_level_profile_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_single_level_profile_MPI_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_single_level_profile_MPI_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_temperature_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_temperature_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_temperature_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_temperature_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_windspeed_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_windspeed_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_average_windspeed_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_average_windspeed_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_modobs_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_modobs_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_modobs_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_modobs_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_modobs_separatefilters_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_repobs_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_repobs_differenterrors_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_differenterrors_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_repobs_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_repobs_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_repobs_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_bkgqc_repobs_separatefilters_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_monolithicfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_monolithicfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_monolithicfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_monolithicfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_obsfunction_fillaveragedprofiledata
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_obsfunction_fillaveragedprofiledata.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_obsfunction_profileaverageobs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_obsfunction_profileaverageobs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_opr_average
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_opr_average.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_linopr_average
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_linopr_average.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_pressure_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_pressure_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_pressure_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_pressure_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_RH_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_RH_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_RH_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_RH_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_separatefilters
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_separatefilters.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_separatefilters
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_separatefilters.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_sondeflags_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_sondeflags_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_sondeflags_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_sondeflags_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_UInterp_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_UInterp_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_UInterp_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_UInterp_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_UInterpAlternative_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_UInterpAlternative_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_UInterpAlternative_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_unittest_oneprofileMPI
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittest_oneprofileMPI.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_unittests
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_unittests
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_unittests_profiledataholder
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_profiledataholder.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_unittests_vertaverage_ascending
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_vertaverage_ascending.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_unittests_vertaverage
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_vertaverage.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_winproflags_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_winproflags_obsfilter
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_winproflags_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_MPI_winproflags_OPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_wrongOPScomparison
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_wrongOPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    conventionalprofileprocessing_fewerobsthanranks
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_fewerobsthanranks.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
