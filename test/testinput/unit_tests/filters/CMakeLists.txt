# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for filters

add_subdirectory(actions)
add_subdirectory(obsfunctions)


ufo_add_test( NAME    amsua_seaice_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_seaice_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    assign_error_using_where
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/assign_error_using_where.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    ensemblestatistics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_EnsembleObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ensemblestatistics.yaml"
              MPI     8
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    ensemblestatistics_channels
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_EnsembleObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ensemblestatistics_channels.yaml"
              MPI     8
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    filter_stage
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/filter_stage.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    filter_stage_mpi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/filter_stage.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobendmetoffice_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobendmetoffice_qc_profile
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_qc_profile.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobendmetoffice_obserror
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_obserror.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrorefmetoffice_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobndnbam_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobndnbam_qc_geos
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_qc_geos.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssrobndnbam_qc_action_obserr_inflation
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_qc_action_obserr_inflation.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    gnssro_obs_error
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssro_obs_error.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    gnssro_super_refraction_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssro_super_refraction_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    metar_qc_filters
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/metar_qc_filters.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    obserror_assign_unittests
              ECBUILD
              SOURCES ../../../mains/TestObsErrorAssign.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obserror_assign_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    obsfilterdata
              ECBUILD
              SOURCES ../../../mains/TestObsFilterData.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsfilterdata.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    ompsnp_npp_qc_L127
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_qc_L127.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_MPI2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_MPI3
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data.yaml"
              MPI     3
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_inefficient_distribution
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data_inefficient_distribution.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_MPI2_rank0
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data_MPI2_rank0.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_time_window_selection
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data_time_window_selection.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    print_filter_data_time_window_selection_MPI2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/print_filter_data_time_window_selection.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    processwhere
              ECBUILD
              SOURCES ../../../mains/TestProcessWhere.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/processwhere.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_acceptlist
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_acceptlist.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_backgroundcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_backgroundcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_backgroundcheck_ensemblestatistics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_EnsembleObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_backgroundcheck_ensemblestatistics.yaml"
              MPI     8
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_bayesianbackgroundqcflags
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_bayesianbackgroundqcflags.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_bayesian_background_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_bayesian_background_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_boundscheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_boundscheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_creatediagnosticflags
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_creatediagnosticflags.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    function_copy_flags_from_extended_to_original_space
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_copy_flags_from_ext_to_orig_space.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_defer_to_post
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_defer_to_post.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_derivative_dpdt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_derivative_dpdt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_derivative_dxdt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_derivative_dxdt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_differencecheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_differencecheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_finalcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_finalcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_mean
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_mean.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_mean_error
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_mean_error.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_median
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_median.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_ops.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_records_are_single_obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_records_are_single_obs.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_satrad
              ECBUILD
              SOURCES ../../../mains/TestGaussianThinning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_satrad.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_gauss_thinning_unittests
              ECBUILD
              SOURCES ../../../mains/TestGaussianThinning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_gauss_thinning_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_grosserrorwholereport
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_grosserrorwholereport.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_historycheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_historycheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_historycheck_mpi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_historycheck.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_historycheck_unittests
              ECBUILD
              SOURCES ../../../mains/TestHistoryCheck.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_historycheck_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_buddy_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_buddy_check.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_buddy_check_unittests
              ECBUILD
              SOURCES ../../../mains/TestMetOfficeBuddyCheck.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_buddy_check_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_buddy_pair_finder
              ECBUILD
              SOURCES ../../../mains/TestMetOfficeBuddyPairFinder.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_buddy_pair_finder.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_pressure_consistency_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_pressure_consistency_check.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_modelbestfitpressure
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_modelbestfitpressure.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_modelobthreshold
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_modelobthreshold.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_multichannel_thinning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_multichannel_thinning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_oceancolor_preqc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_oceancolor_preqc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( ${gsw_FOUND} )
    ufo_add_test( NAME    qc_ocean_vertical_stability_check
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_ocean_vertical_stability_check.yaml"
                  MPI     4
                  LIBS    ufo
                  LABELS  filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${gsw_FOUND} )
#
ufo_add_test( NAME    qc_poisson_disk_thinning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_poisson_disk_thinning.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_poisson_disk_thinning_parallel
              ECBUILD
              SOURCES ../../../mains/TestParallelPoissonDiskThinning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_poisson_disk_thinning_parallel.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_poisson_disk_thinning_unittests
              ECBUILD
              SOURCES ../../../mains/TestPoissonDiskThinning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_poisson_disk_thinning_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_duplicate_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_duplicate_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_met_office_duplicate_check_parallel
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_met_office_duplicate_check.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_preqc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_preqc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_preqc_halo
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_preqc_halo.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_processamvqi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_processamvqi.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_average_obs_to_mod_levels
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_profile_average_obs_to_mod_levels.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_profile_background_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_profile_background_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_profile_fewobs_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_profile_fewobs_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_profile_unflagobs_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_profile_unflagobs_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_spike_and_step_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_spike_and_step_check.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_stuckcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_stuckcheck.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_stuckcheck_unittests
              ECBUILD
              SOURCES ../../../mains/TestStuckCheck.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_stuckcheck_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_thinning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_thinning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_temporal_thinning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_temporal_thinning.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_temporal_thinning_records_are_single_obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_temporal_thinning_records_are_single_obs.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_temporal_thinning_unittests
              ECBUILD
              SOURCES ../../../mains/TestTemporalThinning.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_temporal_thinning_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_trackcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_trackcheck.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_trackcheckship
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_trackcheckship.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_trackcheck_unittests
              ECBUILD
              SOURCES ../../../mains/TestTrackCheck.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_trackcheck_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_trackcheckship_unittests
              ECBUILD
              SOURCES ../../../mains/TestTrackCheckShip.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_trackcheckship_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_trackcheckship_mainloop_unittests
              ECBUILD
              SOURCES ../../../mains/TestTrackCheckShipMainLoop.cc
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_trackcheckship_mainloop_unittests.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    qc_variableassignment
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_variableassignment.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    satname
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satname.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    satwind_inversion_correction
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satwind_inversion_correction.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    iasi_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    superob
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/superob.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    superob_4PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/superob.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    superob_radar
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/superob_radar.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    superob_radar_4PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/superob_radar.yaml"
              MPI     4
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_parameter_substitution
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/parameter_substitution.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    polygon_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/polygon_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Filters that require CRTM in ObsOperator

if( crtm_FOUND )
#
    ufo_add_test( NAME    airs_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsr2_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsr2_qc_halo
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_qc_halo.yaml"
                  MPI     4
                  LIBS    ufo
                  LABELS  crtm filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_allsky_gfs_gsi_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_allsky_gfs_gsi_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_allsky_gsi_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_allsky_gsi_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_qc_clwretmw
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_qc_clwretmw.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_qc_filters_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_qc_filters_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsua_qc_miss_val
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_qc_miss_val.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    atms_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    atms_n20_qc_filters_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_n20_qc_filters_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    cris_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    cris_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    cris_qc_land
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_qc_land.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsr2_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_skylab_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
                  TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data ) 
#
    ufo_add_test( NAME    tropics_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropics_skylab_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    gmi_skylab_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_skylab_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    gmi_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    gmi_qc_filters_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_qc_filters_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    iasi_qc_filters
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_qc_filters.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    mhs_qc_filters_geos
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_qc_filters_geos.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    qc_flags_true
              	  ECBUILD
              	  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              	  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/abi_g16_qc_flags.yaml"
              	  MPI     1
              	  LIBS    ufo
              	  LABELS  filters
		  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND )

# Filters that require RTTOV in ObsOperator

if( ${rttov_FOUND} )
    ufo_add_test( NAME    amsr2_rttov_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_rttov_qc${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    amsr2_rttov_ops_qc_rttovonedvarcheck
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_rttov_ops_qc_rttovonedvarcheck${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    atms_rttov_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_qc${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    atms_rttov_ops_qc_rttovonedvarcheck
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_ops_qc_rttovonedvarcheck${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
# DAR: I cannot make this test work
    ufo_add_test( NAME    iasi_rttov_ops_qc_rttovonedvarcheck
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops_qc_rttovonedvarcheck${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     2
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters ci_oneapi_disable
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    mtgirs_reconrad_obsop_check
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mtgirs_reconrad_obsop_check${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters ci_oneapi_disable
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    ssmis_rttov_ops_qc_rttovonedvarcheck
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_rttov_ops_qc_rttovonedvarcheck${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atovs_rttovonedvar_multiplatform
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttovonedvar_multiplatform${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atovs_rttovonedvar_transmittance
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttovonedvar_transmittance${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atovs_rttovonedvar_transmittance_multiplatform
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttovonedvar_transmittance_multiplatform${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${rttov_FOUND} )


# Filters that require ROPP in ObsOperator

if( ${ropp-ufo_FOUND} )

    ufo_add_test( NAME    gnssrobndropp1d_qc
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp1d_qc.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  ropp-ufo filters
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${ropp-ufo_FOUND} )

ufo_add_test( NAME    preprocess
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/ufo_preprocess.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/preprocess.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
