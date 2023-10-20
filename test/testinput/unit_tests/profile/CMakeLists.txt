# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for profile

#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_pressure_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_pressure_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_pressure_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_pressure_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_relativehumidity_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_relativehumidity_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_relativehumidity_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_relativehumidity_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_single_level_profile_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_single_level_profile_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_temperature_obsfilte
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_temperature_obsfilte
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_temperature_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_temperature_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_windspeed_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_windspeed_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_average_windspeed_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_average_windspeed_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_modobs_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_modobs_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_separatefilters_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_repobs_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_repobs_differenterrors_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_differenterrors_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_repobs_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_separatefilters_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_monolithicfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_monolithicfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_monolithicfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_monolithicfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_obsfunction_fillaveragedprofiledata
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_obsfunction_fillaveragedprofiledata.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_obsfunction_profileaverageobs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_obsfunction_profileaverageobs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_opr_average
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_opr_average.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_linopr_average
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_linopr_average.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_pressure_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_pressure_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_pressure_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_pressure_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_pressure_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_RH_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_RH_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_RH_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_RH_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_RH_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_separatefilters
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_separatefilters.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_separatefilters
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_separatefilters.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_sondeflags_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_sondeflags_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_sondeflags_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_sondeflags_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_UInterp_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_UInterp_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_UInterp_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_UInterp_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_UInterpAlternative_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_UInterpAlternative_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_UInterpAlternative_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_UInterpAlternative_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_unittest_oneprofileMPI
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittest_oneprofileMPI.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_unittests
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_unittests
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_unittests_profiledataholder
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_profiledataholder.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_unittests_vertaverage_ascending
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_vertaverage_ascending.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_unittests_vertaverage
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_unittests_vertaverage.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_winproflags_obsfilter
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_obsfilter.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_winproflags_obsfilter
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_obsfilter.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_winproflags_OPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_MPI_winproflags_OPScomparison
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_conventionalprofileprocessing_wrongOPScomparison
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/conventionalprofileprocessing_wrongOPScomparison.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests profile
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
