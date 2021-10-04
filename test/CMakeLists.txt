# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include macros for creating links and symlinks
include( ufo_functions )
include( git_functions )

# Create Data directory for test input config and symlink all files
list( APPEND ufo_test_input
  testinput/abi_ahi_crtm.yaml
  testinput/adt.yaml
  testinput/aircraft.yaml
  testinput/airs_crtm.yaml
  testinput/airs_qc_filters.yaml
  testinput/amsr2_rttov_qc.yaml
  testinput/amsr2_crtm.yaml
  testinput/amsr2_clw_ret.yaml
  testinput/amsr2_qc.yaml
  testinput/amsr2_qc_halo.yaml
  testinput/amsr2_rttov_ops_qc_rttovonedvarcheck.yaml
  testinput/amsua_crtm.yaml
  testinput/amsua_crtm_bc.yaml
  testinput/amsua_qc.yaml
  testinput/amsua_qc_clwretmw.yaml
  testinput/amsua_qc_filters.yaml
  testinput/amsua_qc_miss_val.yaml
  testinput/amsua_allsky_gsi_qc.yaml
  testinput/amsua_seaice_qc.yaml
  testinput/aod_crtm.yaml
  testinput/aod_luts.yaml
  testinput/aod_metoffice.yaml
  testinput/atms_crtm.yaml
  testinput/atms_qc_filters.yaml
  testinput/atms_rttov_ops_qc_rttovonedvarcheck.yaml
  testinput/atms_rttov_ops.yaml
  testinput/atms_rttov_qc.yaml
  testinput/amsua_rttovcpp.yaml
  testinput/background_error_vert_interp.yaml
  testinput/background_error_identity.yaml
  testinput/bias_coeff.yaml
  testinput/bias_coeff_cov.yaml
  testinput/bias_linear_op.yaml
  testinput/bias_single_channel.yaml
  testinput/categorical.yaml
  testinput/coolskin.yaml
  testinput/composite.yaml
  testinput/constant_predictor.yaml
  testinput/cris_crtm.yaml
  testinput/cris_qc.yaml
  testinput/cris_qc_filters.yaml
  testinput/cris_qc_land.yaml
  testinput/empty.yaml
  testinput/function_bennartzscatindex.yaml
  testinput/function_cloudcost.yaml
  testinput/function_clouddetect_iasi.yaml
  testinput/function_clouddetect_airs.yaml
  testinput/function_clouddetect_cris.yaml
  testinput/function_clwmatchidx.yaml
  testinput/function_clwmatchidx_atms.yaml
  testinput/function_clwret.yaml
  testinput/function_clwret_atms.yaml
  testinput/function_clwretmean.yaml
  testinput/function_clwretmean_atms.yaml
  testinput/function_clwret_hofx.yaml
  testinput/function_clwret_hofx_atms.yaml
  testinput/function_clwret_obsval.yaml
  testinput/function_clwret_obsval_atms.yaml
  testinput/function_errfgrosschk.yaml
  testinput/function_errfjsfc.yaml
  testinput/function_errflat.yaml
  testinput/function_errftopo.yaml
  testinput/function_errftransmittop.yaml
  testinput/function_errfwavenum.yaml
  testinput/function_hydrometeorchk.yaml
  testinput/function_hydrometeorchk_atms.yaml
  testinput/function_bgddepartureanomaly.yaml
  testinput/function_nsstret.yaml
  testinput/function_obserrorinflatesfcp.yaml
  testinput/function_obserrorfactorconv.yaml
  testinput/function_obserrmean.yaml
  testinput/function_obserrmean_abi.yaml
  testinput/function_obserrmean_mhs.yaml
  testinput/function_obserrmean_atms.yaml
  testinput/function_reperr.yaml
  testinput/function_setsurfacetype.yaml
  testinput/function_setsurfacetype_ssmis.yaml
  testinput/function_satwind_indiv_errors.yaml
  testinput/function_scatret.yaml
  testinput/function_scatret_atms.yaml
  testinput/function_errfsdoei.yaml
  testinput/function_errfsdoei_atms.yaml
  testinput/function_velocity.yaml
  testinput/function_obserrorstepwiselinear.yaml
  testinput/gmi_crtm.yaml
  testinput/gmi_clw_ret.yaml
  testinput/gmi_qc_filters.yaml
  testinput/gmi_scan_angles.yaml
  testinput/genericprof.yaml
  testinput/geovals.yaml
  testinput/geos_aod.yaml
  testinput/aod_extinction.yaml
  testinput/geovals_spec.yaml
  testinput/gnssrobendmetoffice.yaml
  testinput/gnssrobendmetoffice_qc.yaml
  testinput/gnssrobendmetoffice_obserror.yaml
  testinput/gnssrobendmetoffice_nopseudo.yaml
  testinput/gnssrobendmetoffice_qc.yaml
  testinput/gnssrorefmetoffice.yaml
  testinput/gnssrobndropp1d.yaml
  testinput/gnssrobndropp1d_qc.yaml
  testinput/gnssrobndropp2d.yaml
  testinput/gnssrobndnbam.yaml
  testinput/gnssro_obs_error.yaml
  testinput/gnssro_obs_error_halo.yaml
  testinput/gnssro_obs_error_ropp.yaml
  testinput/gnssrobndnbam_backgroundcheck_qc.yaml
  testinput/gnssrobndnbam_super_refraction.yaml
  testinput/gnssrobndnbam_qc_action_obserr_inflation.yaml
  testinput/gnssroref.yaml
  testinput/gnssrobndnbam_0obs.yaml
  testinput/gnssrobndnbam_1obs.yaml
  testinput/gnssro_domain_check.yaml
  testinput/gome_metop-a.yaml
  testinput/gome_metop-a_flipz.yaml
  testinput/groundgnssmetoffice.yaml
  testinput/groundgnssropp.yaml
  testinput/gsisfcmodel.yaml
  testinput/hirs4_crtm.yaml
  testinput/iasi_crtm.yaml
  testinput/iasi_qc.yaml
  testinput/iasi_qc_filters.yaml
  testinput/identity.yaml
  testinput/interpolate_data_from_file_predictor.yaml
  testinput/legendre_predictor.yaml
  testinput/locations.yaml
  testinput/metar_qc_filters.yaml
  testinput/mhs_crtm.yaml
  testinput/obsdiag_background_error_identity.yaml
  testinput/obsdiag_background_error_upper_air_and_surface_variables.yaml
  testinput/obsdiag_background_error_vert_interp_air_pressure.yaml
  testinput/obsdiag_background_error_vert_interp_height.yaml
  testinput/obsdiag_crtm_airs_jacobian.yaml
  testinput/obsdiag_crtm_airs_optics.yaml
  testinput/obsdiag_crtm_amsua_jacobian.yaml
  testinput/obsdiag_crtm_amsua_optics.yaml
  testinput/obsdiag_crtm_atms_jacobian.yaml
  testinput/obsdiag_crtm_atms_optics.yaml
  testinput/obsdiag_crtm_cris_jacobian.yaml
  testinput/obsdiag_crtm_cris_optics.yaml
  testinput/obsdiag_crtm_iasi_jacobian.yaml
  testinput/obsdiag_crtm_iasi_optics.yaml
  testinput/obserror_assign_unittests.yaml
  testinput/obserrordiagonal.yaml
  testinput/obserrorcrossvarcorr.yaml
  testinput/obserrorcrossvarcov.yaml
  testinput/obsfilterdata.yaml
  testinput/omi_aura.yaml
  testinput/omi_aura_flipz.yaml
  testinput/ompsnp_npp.yaml
  testinput/ompsnp_npp_flipz.yaml
  testinput/ompsnp_npp_L127.yaml
  testinput/ompsnp_npp_qc_L127.yaml
  testinput/ompstc_npp.yaml
  testinput/ompstc_npp_flipz.yaml
  testinput/ompslp_npp.yaml
  testinput/operatorutils.yaml
  testinput/orbital_angle_predictor.yaml
  testinput/processwhere.yaml
  testinput/parameters.yaml
  testinput/parameters_older_eckit.yaml
  testinput/primitive_variables.yaml
  testinput/metoffice_radiance_error_matrices.yaml
  testinput/qc_actions.yaml
  testinput/qc_backgroundcheck.yaml
  testinput/qc_bayesianbackgroundqcflags.yaml
  testinput/qc_bayesian_background_check.yaml
  testinput/qc_boundscheck.yaml
  testinput/qc_velocitycheck.yaml
  testinput/qc_defer_to_post.yaml
  testinput/qc_derivative_dpdt.yaml
  testinput/qc_derivative_dxdt.yaml
  testinput/qc_differencecheck.yaml
  testinput/qc_finalcheck.yaml
  testinput/qc_gauss_thinning.yaml
  testinput/qc_gauss_thinning_ops.yaml
  testinput/qc_gauss_thinning_unittests.yaml
  testinput/qc_grosserrorwholereport.yaml
  testinput/qc_historycheck_unittests.yaml
  testinput/qc_met_office_buddy_check.yaml
  testinput/qc_met_office_buddy_check_unittests.yaml
  testinput/qc_met_office_buddy_pair_finder.yaml
  testinput/qc_modelbestfitpressure.yaml
  testinput/qc_modelobthreshold.yaml
  testinput/qc_multichannel_thinning.yaml
  testinput/qc_oceancolor_preqc.yaml
  testinput/qc_performaction.yaml
  testinput/qc_poisson_disk_thinning.yaml
  testinput/qc_poisson_disk_thinning_unittests.yaml
  testinput/qc_poisson_disk_thinning_parallel.yaml
  testinput/qc_preqc.yaml
  testinput/qc_preqc_halo.yaml
  testinput/qc_processamvqi.yaml
  testinput/qc_profile_background_check.yaml
  testinput/qc_profile_fewobs_check.yaml
  testinput/qc_stuckcheck.yaml
  testinput/qc_stuckcheck_unittests.yaml
  testinput/qc_thinning.yaml
  testinput/qc_temporal_thinning.yaml
  testinput/qc_temporal_thinning_unittests.yaml
  testinput/qc_trackcheck.yaml
  testinput/qc_trackcheck_unittests.yaml
  testinput/qc_trackcheckship.yaml
  testinput/qc_trackcheckship_mainloop_unittests.yaml
  testinput/qc_trackcheckship_unittests.yaml
  testinput/qc_variableassignment.yaml
  testinput/qc_acceptlist.yaml
  testinput/radiosonde.yaml
  testinput/radialvelocity.yaml
  testinput/reflectivity.yaml
  testinput/satellite_selector_predictor.yaml
  testinput/satname.yaml
  testinput/sattcwv.yaml
  testinput/satwind.yaml
  testinput/satwind_inversion_correction.yaml
  testinput/satwind_metoffice.yaml
  testinput/sbuv2_n19.yaml
  testinput/sbuv2_n19_flipz.yaml
  testinput/sbuv2_n19_L127.yaml
  testinput/scatwind_neutral_metoffice.yaml
  testinput/seaicefrac.yaml
  testinput/seaicethick.yaml
  testinput/chleuzintegr.yaml
  testinput/sea_surface_temp.yaml
  testinput/seviri_crtm.yaml
  testinput/smap_crtm.yaml
  testinput/sndrd1-4_crtm.yaml
  testinput/sfcpcorrected.yaml
  testinput/thickness_predictor.yaml
  testinput/conventionalprofileprocessing_monolithicfilter.yaml
  testinput/conventionalprofileprocessing_OPScomparison.yaml
  testinput/conventionalprofileprocessing_wrongOPScomparison.yaml
  testinput/conventionalprofileprocessing_separatefilters.yaml
  testinput/conventionalprofileprocessing_unittests.yaml
  testinput/conventionalprofileprocessing_unittest_oneprofileMPI.yaml
  testinput/conventionalprofileprocessing_unittests_vertaverage.yaml
  testinput/conventionalprofileprocessing_unittests_vertaverage_ascending.yaml
  testinput/conventionalprofileprocessing_unittests_vertinterp.yaml
  testinput/conventionalprofileprocessing_unittests_profiledataholder.yaml
  testinput/conventionalprofileprocessing_RH_obsfilter.yaml
  testinput/conventionalprofileprocessing_RH_OPScomparison.yaml
  testinput/conventionalprofileprocessing_UInterp_obsfilter.yaml
  testinput/conventionalprofileprocessing_UInterp_OPScomparison.yaml
  testinput/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml
  testinput/conventionalprofileprocessing_UInterpAlternative_OPScomparison.yaml
  testinput/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml
  testinput/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml
  testinput/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml
  testinput/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml
  testinput/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml
  testinput/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml
  testinput/conventionalprofileprocessing_sondeflags_obsfilter.yaml
  testinput/conventionalprofileprocessing_sondeflags_OPScomparison.yaml
  testinput/conventionalprofileprocessing_winproflags_obsfilter.yaml
  testinput/conventionalprofileprocessing_winproflags_OPScomparison.yaml
  testinput/conventionalprofileprocessing_pressure_obsfilter.yaml
  testinput/conventionalprofileprocessing_pressure_OPScomparison.yaml
  testinput/conventionalprofileprocessing_average_pressure_obsfilter.yaml
  testinput/conventionalprofileprocessing_average_pressure_OPScomparison.yaml
  testinput/conventionalprofileprocessing_average_temperature_obsfilter.yaml
  testinput/conventionalprofileprocessing_average_temperature_OPScomparison.yaml
  testinput/conventionalprofileprocessing_average_windspeed_obsfilter.yaml
  testinput/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml
  testinput/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml
  testinput/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml
  testinput/conventionalprofileprocessing_opr_average.yaml
  testinput/ssmis_f17_gfs_backgroundcheck_bc.yaml
  testinput/ssmis_f17_gfs_backgroundcheck_nbc.yaml
  testinput/timeoper.yaml
  testinput/tprof.yaml
  testinput/tropomi_no2.yaml
  testinput/variables.yaml
  testinput/windprof.yaml
  testinput/unit_tests/composite_GsiSfcModel_SfcPCorrected.yaml
  testinput/unit_tests/composite_vertInterp_SfcPCorrected.yaml
  testinput/unit_tests/function_clwret_ssmis.yaml
  testinput/unit_tests/function_conditional.yaml
  testinput/unit_tests/function_drawvaluefromfile.yaml
  testinput/unit_tests/function_tropo_estimate.yaml
  testinput/unit_tests/function_lam_domaincheck.yaml
  testinput/unit_tests/function_lam_domaincheck_circle.yaml
  testinput/unit_tests/function_ModelSurfaceHeightAdjustMarineWinds.yaml
  testinput/unit_tests/function_ModelSurfaceHeightAdjustRelativeHumidity.yaml
  testinput/unit_tests/function_ModelSurfaceHeightAdjustWindVector.yaml
  testinput/unit_tests/function_ModelSurfaceHeightAdjustTemperature.yaml
  testinput/unit_tests/function_linear_combination_multichannel.yaml
  testinput/unit_tests/function_obserrorfactorquotient.yaml
  testinput/unit_tests/function_satwindsLNVDcheck.yaml
  testinput/unit_tests/function_satwindsSPDBcheck.yaml
  testinput/unit_tests/function_solarzenith.yaml
  testinput/unit_tests/function_solarzenith_skiprejected.yaml
  testinput/unit_tests/function_wdirdiff.yaml
  testinput/unit_tests/utils_metoffice_sort.yaml
  testinput/unit_tests/variabletransforms_rhumidity.yaml
  testinput/unit_tests/variabletransforms_rhumidity_part2.yaml
  testinput/unit_tests/variabletransforms_shumidity.yaml
  testinput/unit_tests/variabletransforms_scanposition.yaml
  testinput/unit_tests/variabletransforms_SondePressureFromHeight.yaml
  testinput/unit_tests/variabletransforms_SondePressureFromHeightICAO.yaml
  testinput/unit_tests/variabletransforms_SondeHeightFromPressure.yaml
  testinput/unit_tests/variabletransforms_szacorrection.yaml
  testinput/unit_tests/variabletransforms_windcomponents.yaml
  testinput/unit_tests/variabletransforms_windspeedanddirection.yaml
  testinput/unit_tests/variabletransforms_profilehorizontaldrift.yaml
  testinput/unit_tests/variabletransforms_pstar.yaml
)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/unit_tests)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_test_input} )


# Add files to cmake resources
ecbuild_add_resources( TARGET ufo_test_scripts
                       SOURCES_PACK
                       ${ufo_test_input}
                     )

# Create Data directory for test data and symlink files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# Test data filenames
set(ufo_test_data  ufo_testinput_tier_1.tar.gz)
set(crtm_test_data crtm_coefficients.tar.gz)

#############################################
# macro for downloading test data
# options:
# - tag 1.0.0: download from UCAR ftp, flat directory
#   Note: when we change the tag, this will have to be updated somehow.
# - anything else: download from AWS, either from branch directory
#   if available, or from develop directory
# input: repo (e.g. ioda or ufo)
#        test_data: test data filename (archive)
# output: output_data_path: path to data on local machine after download

# Currently only use 2.3.0 CRTM coefficients
set( CRTM_COEFFS_BRANCH "2.3.0" )

function(ADD_DOWNLOAD_TEST repo test_files_dirname output_data_path)
  # If REPO branch is being built set GIT_BRANCH_FUNC to REPO's current branch.
  # If a tagged version of REPO is being built set GIT_TAG_FUNC to REPO's current tag.
  find_branch_name(REPO_DIR_NAME ${repo})

  if ( ${repo} STREQUAL "crtm" )
    set( GIT_TAG_FUNC ${CRTM_COEFFS_BRANCH} )
  endif()

  if( DEFINED GIT_BRANCH_FUNC )
    set( REPO_GIT_BRANCH ${GIT_BRANCH_FUNC} )
  elseif( DEFINED GIT_TAG_FUNC )
    set( REPO_GIT_BRANCH ${GIT_TAG_FUNC} )
  endif()


  # When env veriable LOCAL_PATH_JEDI_TESTFILES is set, simply link test files
  # to build directory. get_${repo}_test_data checks the existence of test file directory.
  # ${repo} test data must be stored in
  # ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${BRANCH}/testinput_tier_1
  if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} )
    set( LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}" )
    message( STATUS "use LOCAL_PATH_JEDI_TESTFILES: ${LOCAL_PATH_JEDI_TESTFILES}" )
    # If test data specific to testing branch exists locally use it.
    # If not use test data specific to develop branch.
    if( EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH} )
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH}" )
    else()
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/develop" )
    endif()

    message( STATUS "Test data in ${TESTFILE_DIR_REPO} is linked to build directory" )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
    ${TESTFILE_DIR_REPO} )
    set ( REPO_DATA_DOWNLOADER ufo_data_checker.py )

  # Any tagged version of REPO is being built.
  # Build ufo_get_{repo}_test_data test to download test data from DASH.
  elseif( DEFINED GIT_TAG_FUNC)
    message( STATUS "Tagged version of ${repo} is used" )
    # set ARGS for get_${repo}_test_data

    # ECBUILD_DOWNLOAD_BASE_URL env var can be used to force test files and
    # crtm coef to be downloaded from other databases such as S3 instead of DASH
    # example ECBUILD_DOWNLOAD_BASE_URL=https://jedi-test-files.s3.amazonaws.com

    if( DEFINED ENV{ECBUILD_DOWNLOAD_BASE_URL} )
      set( ECBUILD_DOWNLOAD_BASE_URL "$ENV{ECBUILD_DOWNLOAD_BASE_URL}/${repo}/${GIT_TAG_FUNC}" )
    else()
      set( ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/api/v1/dataset/147_miesch/version/1.1.0/file )
    endif()

    set( DIRNAME ${GIT_TAG_FUNC} )
    set( checksum "0" )
    set( TESTFILE_DIR_REPO "${CMAKE_SOURCE_DIR}/test-data-release/${repo}/${DIRNAME}" )

    # Create test-data-release in source directory
    file( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/test-data-release )
    # temporary solution for release 1.1.0
    if ( ${repo} STREQUAL "crtm" )
      list( APPEND REPO_DATA_DOWNLOADER_ARGS
            ${ECBUILD_DOWNLOAD_BASE_URL}
            ${TESTFILE_DIR_REPO}
            ${test_files_dirname}
            ${checksum} )
      file( MAKE_DIRECTORY ${TESTFILE_DIR_REPO} )
    else()
      list( APPEND REPO_DATA_DOWNLOADER_ARGS
            ${ECBUILD_DOWNLOAD_BASE_URL}
            ${CMAKE_SOURCE_DIR}/test-data-release
            ${test_files_dirname}
            ${checksum} )
    endif()
    message( STATUS "ECBUILD_DOWNLOAD_BASE_URL is set, ${repo} test files will
    be downloaded from ${ECBUILD_DOWNLOAD_BASE_URL}/test-data-release")
    set ( REPO_DATA_DOWNLOADER ufo_data_downloader.py )

  # Any branch of REPO is being built.
  # repo-data repository is already cloned by bundle/CMakeLists.txt.
  # Link test files in repo-data repo to build directory.
  # get_repo_test_data checks the existence of test file directory
  else()
    set( TESTFILE_DIR_REPO ${CMAKE_SOURCE_DIR}/${repo}-data )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${TESTFILE_DIR_REPO} )
    set (REPO_DATA_DOWNLOADER ufo_data_checker.py )
  endif()

  # add test for downloading data
  ecbuild_add_test( TARGET    ufo_get_${repo}_test_data
                    TYPE      SCRIPT
                    COMMAND   ${CMAKE_BINARY_DIR}/bin/${REPO_DATA_DOWNLOADER}
                    ARGS  ${REPO_DATA_DOWNLOADER_ARGS} )

  # pass output data path back to calling code
  set(${output_data_path} ${TESTFILE_DIR_REPO} PARENT_SCOPE)

  # Create download script for ufo_get_*_test_data test
  set ( FILENAME ${REPO_DATA_DOWNLOADER})
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()


  add_custom_target( bin_ufo_${repo}_test_download_scripts ALL
    COMMAND chmod +x ${DEST_FILE}
    DEPENDS ${DEST_FILE} )

endfunction()
#############################################
# add tests to download ufo & crtm data

ADD_DOWNLOAD_TEST( "ufo" ${ufo_test_data} UFO_UFO_TESTFILES_PATH )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${UFO_UFO_TESTFILES_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/ufo)

if( crtm_FOUND )
  ADD_DOWNLOAD_TEST( "crtm" ${crtm_test_data} UFO_CRTM_TESTFILES_PATH )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                   ${UFO_CRTM_TESTFILES_PATH}
                   ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm)
endif( crtm_FOUND )
#####################################################################
# Install headers used in tests of interfaces that may be implemented by clients.

install( FILES ufo/ObsDiagnostics.h ufo/ObsFilters.h ufo/ObsFunction.h ufo/Predictor.h
         DESTINATION ${INSTALL_INCLUDE_DIR}/ufo/test/ufo)

#####################################################################
# Build executables used by multiple tests

ecbuild_add_executable( TARGET  test_ObsBias.x
                        SOURCES mains/TestObsBias.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsDiagnostics.x
                        SOURCES mains/TestObsDiagnostics.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsOperator.x
                        SOURCES mains/TestObsOperator.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsOperatorTLAD.x
                        SOURCES mains/TestObsOperatorTLAD.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsFilters.x
                        SOURCES mains/TestObsFilters.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsFunction.x
                        SOURCES mains/TestObsFunction.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ConventionalProfileProcessing.x
                        SOURCES mains/TestConventionalProfileProcessing.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_Predictor.x
                        SOURCES mains/TestPredictor.cc
                        LIBS    ufo
                       )

####################################################################
# Establish test tiering

set( UFO_TEST_TIER 1 CACHE STRING "UFO test tier. Run tests up to and including this tier.")

message(STATUS "Running ufo tests for tiers up to and including tier " ${UFO_TEST_TIER})

####################################################################
# Test Locations, GeoVaLs with oops tests:

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_geovals_mpi
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_locations
                  SOURCES mains/TestLocations.cc ufo/Locations.h ufo/locations_test.F90
                  ARGS    "testinput/locations.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test ObsError implementations

ecbuild_add_test( TARGET  test_ufo_obserrordiagonal
                  SOURCES mains/TestObsErrorCovariance.cc
                  ARGS    "testinput/obserrordiagonal.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_obserrorcrossvarcorr
                  SOURCES mains/TestObsErrorCovariance.cc
                  ARGS    "testinput/obserrorcrossvarcorr.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_obserrorcrossvarcov
                  SOURCES mains/TestObsErrorCovariance.cc
                  ARGS    "testinput/obserrorcrossvarcov.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test UFO-specific classes:

ecbuild_add_test( TARGET  test_ufo_obsfilterdata
                  SOURCES mains/TestObsFilterData.cc
                  ARGS    "testinput/obsfilterdata.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variables
                  SOURCES mains/TestVariables.cc
                  ARGS    "testinput/variables.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals_spec
                  SOURCES mains/TestGeoVaLsSpec.cc
                  ARGS    "testinput/geovals_spec.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test Obs Operators and TLAD

ecbuild_add_test( TARGET  test_ufo_opr_gsi_sfc_model
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gsisfcmodel.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_opr_crtm_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/amsua_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_bc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/amsua_crtm_bc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/amsua_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_bc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      MPI     2
                      ARGS    "testinput/amsua_crtm_bc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_atms
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/atms_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_atms
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/atms_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_gmi_clw_ret
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/gmi_clw_ret.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_gmi_qc_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/gmi_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_ssmis_f17_gfs_backgroundcheck_bc_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/ssmis_f17_gfs_backgroundcheck_bc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data)

    ecbuild_add_test( TARGET  test_ufo_ssmis_f17_gfs_backgroundcheck_nbc_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/ssmis_f17_gfs_backgroundcheck_nbc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data)

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_gmi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/gmi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_gmi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/gmi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_gmi_scan_angles
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/gmi_scan_angles.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_gmi_scan_angles
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/gmi_scan_angles.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_smap
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/smap_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_smap
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/smap_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_seviri
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/seviri_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_seviri
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/seviri_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/cris_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/cris_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_mhs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/mhs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_mhs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/mhs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_iasi
                      MPI     2
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/iasi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_iasi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/iasi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_sndrd1-4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/sndrd1-4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_sndrd1-4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/sndrd1-4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_airs
                      MPI     4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/airs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_airs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/airs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_hirs4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/hirs4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_hirs4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/hirs4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_abi_ahi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/abi_ahi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_abi_ahi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/abi_ahi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_aod_crtm
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/aod_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_aod_crtm
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/aod_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_amsr2
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/amsr2_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_amsr2
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/amsr2_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_amsr2_clw_ret
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsr2_clw_ret.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_amsr2_qc
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsr2_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_amsr2_qc_halo
                      MPI     4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsr2_qc_halo.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    if( ${geos-aero_FOUND} )

      ecbuild_add_test( TARGET  test_ufo_opr_aod_luts
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/aod_luts.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

      ecbuild_add_test( TARGET  test_ufo_linopr_aod_luts
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/aod_luts.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

   endif( ${geos-aero_FOUND} )


endif( crtm_FOUND )



if( ${rttov_FOUND} )
  ecbuild_add_test( TARGET  test_ufo_function_cloudcostfunction
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/function_cloudcost.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsFilters.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_opr_rttov_atms
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                    ARGS    "testinput/atms_rttov_ops.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperator.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_linopr_rttov_atms
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                    ARGS    "testinput/atms_rttov_ops.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperatorTLAD.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_qc_atms_rttov
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/atms_rttov_qc.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsFilters.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_atms_ops_qc_rttovonedvarcheck
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/atms_rttov_ops_qc_rttovonedvarcheck.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_qc_amsr2_rttov
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/amsr2_rttov_qc.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsFilters.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_amsr2_ops_qc_rttovonedvarcheck
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/amsr2_rttov_ops_qc_rttovonedvarcheck.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_opr_rttovcpp_amsua
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                    ARGS    "testinput/amsua_rttovcpp.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperator.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_linopr_rttovcpp_amsua
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                    ARGS    "testinput/amsua_rttovcpp.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperatorTLAD.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

  ecbuild_add_test( TARGET  test_ufo_variabletransforms_scanposition
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                    ARGS    "testinput/unit_tests/variabletransforms_scanposition.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsFilters.x
                    TEST_DEPENDS ufo_get_ufo_test_data )

endif( ${rttov_FOUND} )


ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_gome_metop-a
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gome_metop-a.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_gome_metop-a
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gome_metop-a.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_gome_metop-a_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gome_metop-a_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_gome_metop-a_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gome_metop-a_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_omi_aura
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/omi_aura.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_omi_aura
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/omi_aura.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_omi_aura_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/omi_aura_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterp_ompslp_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompslp_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterp_ompslp_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompslp_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_omi_aura_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/omi_aura_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_ompsnp_npp
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompsnp_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )


ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp_flipz
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_ompsnp_npp_flipz
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompsnp_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp_L127
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_ompsnp_npp_L127
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompsnp_npp_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp_qc_L127
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/ompsnp_npp_qc_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompstc_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompstc_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )


ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_ompstc_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompstc_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )


ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompstc_npp_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompstc_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_ompstc_npp_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/ompstc_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )



ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_sbuv_n19
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sbuv2_n19.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_sbuv2_n19_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sbuv2_n19_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19_L127
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_atminterplay_sbuv2_n19_L127
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sbuv2_n19_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_vertinterp_radiosonde
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radiosonde.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS  ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_radiosonde
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/radiosonde.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_avgkernel_tropomi_no2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/tropomi_no2.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_avgkernel_tropomi_no2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/tropomi_no2.yaml"
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_timeoper_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/timeoper.yaml"
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_vertinterp_aircraft_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/aircraft.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_aircraft
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/aircraft.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_satwind
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/satwind.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_satwind_metoffice
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/satwind_metoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_scatwindNeutralMetOffice
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/scatwind_neutral_metoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_timeoper
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/timeoper.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_timeoper
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/timeoper.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_background_error_vert_interp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/background_error_vert_interp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_opr_background_error_identity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/background_error_identity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_opr_gnssroRef
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssroref.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_0obs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_0obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_1obs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_1obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_1obs_2pe
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_1obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroRef
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssroref.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_super_refraction_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_super_refraction.yaml"
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBendMetOffice_pseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobendmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBendMetOffice_pseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobendmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gnssroBendMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobendmetoffice_qc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBendMetOffice_nopseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobendmetoffice_nopseudo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBendMetOffice_nopseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobendmetoffice_nopseudo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroRefMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrorefmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroRefMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrorefmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

if( ${ropp-ufo_FOUND} )
ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndROPP1D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBndROPP1D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndROPP2D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBndROPP2D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x )
endif( ${ropp-ufo_FOUND} )

#ecbuild_add_test( TARGET  test_ufo_opr_windprof
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
#                  ARGS    "testinput/windprof.yaml"
#                  ENVIRONMENT OOPS_TRAPFPE=1
#                  DEPENDS test_ObsOperator.x )

# Marine UFO tests
ecbuild_add_test( TARGET  test_ufo_opr_seaicefrac
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seaicefrac.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_seaicefrac
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicefrac.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_seaicethick
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seaicethick.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_seaicethick
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicethick.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_chleuzintegr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/chleuzintegr.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_identity_sst
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_identity_sst
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_adt
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/adt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_adt
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/adt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_coolskin
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/coolskin.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_coolskin
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/coolskin.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# METAR (sfc. obs) UFO tests
ecbuild_add_test( TARGET  test_ufo_metar_QC
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/metar_qc_filters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

if( ${gsw_FOUND} )
    ecbuild_add_test( TARGET  test_ufo_opr_insitutemperature
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/tprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_insitutemperature
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/tprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_marinevertinterp
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/genericprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ufo_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_marinevertinterp
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/genericprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ufo_test_data )

endif( ${gsw_FOUND} )

if( ${geos-aero_FOUND} )

ecbuild_add_test( TARGET  test_ufo_opr_geos_aero
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/geos_aod.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_geos_aero
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/geos_aod.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS    test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )
endif( ${geos-aero_FOUND} )

ecbuild_add_test( TARGET  test_ufo_opr_aod_ext
                  MPI     4
                  COMMAND  ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/aod_extinction.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS    test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_aod_ext
                  COMMAND  ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/aod_extinction.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS    test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_radarreflectivity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/reflectivity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_radarradialvelocity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radialvelocity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_composite
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/composite.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_composite
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/composite.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_compositeSfc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/unit_tests/composite_GsiSfcModel_SfcPCorrected.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_compositeSonde
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/unit_tests/composite_vertInterp_SfcPCorrected.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_categorical
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/categorical.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_categorical
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/categorical.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_profile_average
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/conventionalprofileprocessing_opr_average.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_profile_average
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/conventionalprofileprocessing_opr_average.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_identity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/identity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_identity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/identity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test UFO ObsFilters (generic)

ecbuild_add_test( TARGET  test_ufo_qc_gen_processwhere
                  SOURCES mains/TestProcessWhere.cc
                  ARGS    "testinput/processwhere.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_actions
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_actions.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_backgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_backgroundcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_bayesianbackgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_bayesian_background_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_profilebackgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_profile_background_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_profilefewobservationscheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_profile_fewobs_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_satname
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/satname.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_satwind_inversion_correction
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/satwind_inversion_correction.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_boundscheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_boundscheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_grosserror_wholereport
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_grosserrorwholereport.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_velocitybounds
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_velocitycheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_quotientcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_obserrorfactorquotient.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_tropocheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_tropo_estimate.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_wdircheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_wdirdiff.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_lam_domaincheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_lam_domaincheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS )

ecbuild_add_test( TARGET  test_ufo_qc_lam_domaincheck_circle
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_lam_domaincheck_circle.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_derivative_ascent
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_derivative_dpdt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_derivative_speed
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_derivative_dxdt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_finalcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_finalcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_defer_to_post
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_defer_to_post.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_differencecheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_differencecheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gauss_thinning
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_gauss_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gauss_thinning_ops
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_gauss_thinning_ops.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_check
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_met_office_buddy_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_multichannel_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_multichannel_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_oceancolor_preqc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_oceancolor_preqc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_poisson_disk_thinning
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_poisson_disk_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_preqc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_preqc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_preqc_halo
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_preqc_halo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_modelbestfitpressure
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_modelbestfitpressure.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_modelobthreshold
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_modelobthreshold.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_stuckcheck
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_stuckcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_temporalthinning
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_temporal_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_trackcheck
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_trackcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_trackcheckship
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_trackcheckship.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_monolithicfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_monolithicfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_separatefilters
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_separatefilters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittests
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittest_oneprofileMPI
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_unittest_oneprofileMPI.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_RH_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_RH_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_UInterp_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterp_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_UInterpAlternative_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_repobs_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_modobs_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_sondeflags_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_winproflags_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_winproflags_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_pressure_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_pressure_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_satwindsLNVDcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_satwindsLNVDcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_satwindsSPDBcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_satwindsSPDBcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_variableassignment
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_variableassignment.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_SondePressureFromHeight
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_SondePressureFromHeight.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_SondePressureFromHeightICAO
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_SondePressureFromHeightICAO.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_SondeHeightFromPressure
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_SondeHeightFromPressure.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_satelliteZenithAngleCorrection
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_szacorrection.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_pressure_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_temperature_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_windspeed_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_relativehumidity_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_acceptlist
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_acceptlist.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_monolithicfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_monolithicfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_separatefilters
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_separatefilters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_unittests
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_RH_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_RH_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_UInterp_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterp_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_UInterpAlternative_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterpAlternative_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_separatefilters_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_separatefilters_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_separatefilters_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_separatefilters_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_sondeflags_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_sondeflags_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_winproflags_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_winproflags_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_pressure_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_pressure_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_pressure_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_pressure_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_temperature_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_temperature_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_windspeed_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_windspeed_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_relativehumidity_obsfilter
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/conventionalprofileprocessing_average_relativehumidity_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_performaction
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_performaction.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_bayesianbackgroundqcflags
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_bayesianbackgroundqcflags.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test UFO ObsFilters (specific)

ecbuild_add_test( TARGET  test_ufo_gaussianthinning
                  SOURCES mains/TestGaussianThinning.cc
                  ARGS    "testinput/qc_gauss_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_recursivesplitter
                  SOURCES mains/TestRecursiveSplitter.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_dataextractor
                  SOURCES mains/TestDataExtractor.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_stuckcheck_unittests
                  SOURCES mains/TestStuckCheck.cc
                  ARGS    "testinput/qc_stuckcheck_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET test_ufo_qc_historycheck_unittests
                  SOURCES mains/TestHistoryCheck.cc
                  ARGS    "testinput/qc_historycheck_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_temporalthinning
                  SOURCES mains/TestTemporalThinning.cc
                  ARGS    "testinput/qc_temporal_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_obserror_assign_unittests
                  SOURCES mains/TestObsErrorAssign.cc
                  ARGS    "testinput/obserror_assign_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheck_unittests
                  SOURCES mains/TestTrackCheck.cc
                  ARGS    "testinput/qc_trackcheck_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheckship_unittests
                  SOURCES mains/TestTrackCheckShip.cc
                  ARGS    "testinput/qc_trackcheckship_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheckship_mainloop_unittests
                  SOURCES mains/TestTrackCheckShipMainLoop.cc
                  ARGS    "testinput/qc_trackcheckship_mainloop_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_poisson_disk_thinning_unittests
                  SOURCES mains/TestPoissonDiskThinning.cc
                  ARGS    "testinput/qc_poisson_disk_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_poisson_disk_thinning_parallel
                  MPI     2
                  SOURCES mains/TestParallelPoissonDiskThinning.cc
                  ARGS    "testinput/qc_poisson_disk_thinning_parallel.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_modobs_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_bkgqc_repobs_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_RH_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_RH_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_wrongOPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_wrongOPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_UInterp_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_UInterpAlternative_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterpAlternative_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_sondeflags_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_winproflags_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_pressure_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_pressure_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittests_vertaverage
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests_vertaverage.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittests_vertaverage_ascending
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests_vertaverage_ascending.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittests_vertinterp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests_vertinterp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_unittests_profiledataholder
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_unittests_profiledataholder.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_check_unittests
                  SOURCES mains/TestMetOfficeBuddyCheck.cc
                  ARGS    "testinput/qc_met_office_buddy_check_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_pair_finder
                  SOURCES mains/TestMetOfficeBuddyPairFinder.cc
                  ARGS    "testinput/qc_met_office_buddy_pair_finder.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_pressure_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_temperature_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_windspeed_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_average_relativehumidity_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_modobs_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_modobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_bkgqc_repobs_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_bkgqc_repobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_RH_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_RH_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_UInterp_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_UInterp_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_sondeflags_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_sondeflags_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_winproflags_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_winproflags_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_pressure_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_pressure_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_pressure_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_pressure_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_temperature_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_temperature_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_windspeed_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_windspeed_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_conventionalprofileprocessing_MPI_average_relativehumidity_OPScomparison
                  MPI 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ConventionalProfileProcessing.x
                  ARGS    "testinput/conventionalprofileprocessing_average_relativehumidity_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ConventionalProfileProcessing.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test Variable Transforms

ecbuild_add_test( TARGET  test_ufo_variabletransforms_windcomponents
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_windcomponents.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_windspeedanddirection
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_windspeedanddirection.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_humidityrelative
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_rhumidity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_humidityrelative_part2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_rhumidity_part2.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_humidityspecific
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_shumidity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_profilehorizontaldrift
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_profilehorizontaldrift.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_variabletransforms_pstar
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/variabletransforms_pstar.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_solarzenith_skiprejected
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/unit_tests/function_solarzenith_skiprejected.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test Functions

ecbuild_add_test( TARGET  test_ufo_function_bgddepartureanomaly
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_bgddepartureanomaly.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_bennartzscatindex
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_bennartzscatindex.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_iasi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_iasi.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_airs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_airs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_cris
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_cris.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwmatchidx
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwmatchidx.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwmatchidx_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwmatchidx_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_ssmis
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_clwret_ssmis.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwretmean
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwretmean.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwretmean_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwretmean_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_hofx
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_hofx.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_hofx_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_hofx_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_obsval
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_obsval.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_clwret_obsval_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_obsval_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_conditional
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_conditional.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_drawvaluefromfile
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_drawvaluefromfile.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errfgrosschk
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfgrosschk.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errfjsfc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfjsfc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS  ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errflat
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errflat.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errftopo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errftopo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errftransmittop
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errftransmittop.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errfwavenum
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfwavenum.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_hydrometeorchk
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_hydrometeorchk.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_hydrometeorchk_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_hydrometeorchk_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_nsstret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_nsstret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrmean
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_abi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_abi.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_mhs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_mhs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_satwind_indiv_errors
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/function_satwind_indiv_errors.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_scatret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_scatret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_scatret_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_scatret_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errfsdoei
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfsdoei.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_errfsdoei_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfsdoei_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_velocity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_velocity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrormodelstepwiselinear
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrorstepwiselinear.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_obserrorinflatesfcp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrorinflatesfcp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_obserrorfactorconv.yaml
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrorfactorconv.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_linear_combination
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_reperr.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS
                  test_ufo_function_representation_err_synthetic_data)

ecbuild_add_test( TARGET  test_ufo_function_linear_combination_multichannel
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_linear_combination_multichannel.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_setsurfacetype
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_setsurfacetype.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_setsurfacetype_ssmis
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_setsurfacetype_ssmis.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_function_modelsurfaceheightadjust_windvector
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_ModelSurfaceHeightAdjustWindVector.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_modelsurfaceheightadjust_temperature
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_ModelSurfaceHeightAdjustTemperature.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_modelsurfaceheightadjust_relativehumidity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_ModelSurfaceHeightAdjustRelativeHumidity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)


ecbuild_add_test( TARGET  test_ufo_function_modelsurfaceheightadjust_marinewind
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_ModelSurfaceHeightAdjustMarineWinds.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_function_solarzenith
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/unit_tests/function_solarzenith.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

# Test Diagnostics

ecbuild_add_test( TARGET  test_ufo_obsdiag_background_error_vert_interp_air_pressure
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                  ARGS    "testinput/obsdiag_background_error_vert_interp_air_pressure.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsDiagnostics.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_obsdiag_background_error_vert_interp_height
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                  ARGS    "testinput/obsdiag_background_error_vert_interp_height.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsDiagnostics.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_obsdiag_background_error_identity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                  ARGS    "testinput/obsdiag_background_error_identity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsDiagnostics.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_obsdiag_background_error_upper_air_and_surface_variables
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                  ARGS    "testinput/obsdiag_background_error_upper_air_and_surface_variables.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsDiagnostics.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_airs_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_airs_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_airs_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_airs_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_amsua_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_amsua_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_amsua_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_amsua_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_atms_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_atms_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_atms_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_atms_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_cris_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_cris_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_cris_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_cris_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_iasi_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_iasi_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_iasi_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_iasi_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )
endif( crtm_FOUND )

# Test Parameters
if (eckit_VERSION VERSION_GREATER_EQUAL 1.16)
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo)
else()
  # Due to changes in eckit 1.16 a different input file is required.
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters_older_eckit.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo)
endif()

# Test piecewise linear interpolation
ecbuild_add_test( TARGET  test_ufo_piecewise_linear_interpolation
                  SOURCES mains/TestPiecewiseLinearInterpolation.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test operator utils
ecbuild_add_test( TARGET  test_ufo_operator_utils
                  SOURCES mains/TestOperatorUtils.cc
                  ARGS    "testinput/operatorutils.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test the iterator over primitive variables and their values
ecbuild_add_test( TARGET  test_ufo_primitive_variables
                  SOURCES mains/TestPrimitiveVariables.cc
                  ARGS    "testinput/primitive_variables.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test the met office radiance read in
ecbuild_add_test( TARGET  test_ufo_metoffice_error_matrices
                  SOURCES mains/TestMetOfficeRadianceErrorMatrices.cc
                  ARGS    "testinput/metoffice_radiance_error_matrices.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test the sorting algorithm used for compatibility with the Met Office OPS system
ecbuild_add_test( TARGET  test_ufo_metoffice_sort
                  SOURCES mains/TestMetOfficeSort.cc
                  ARGS    "testinput/unit_tests/utils_metoffice_sort.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test QC for specific instruments

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_qc_airs_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/airs_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_clwretmw
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_clwretmw.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_miss_val
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_miss_val.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_atms_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/atms_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris_land
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc_land.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_function_scattering
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_iasi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/iasi_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_iasi_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/iasi_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_allsky_gsi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_allsky_gsi_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )

endif( crtm_FOUND )

if( ${ropp-ufo_FOUND} )
ecbuild_add_test( TARGET  test_ufo_qc_gnssroBndROPP1D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobndropp1d_qc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )
endif( ${ropp-ufo_FOUND} )

ecbuild_add_test( TARGET  test_ufo_opr_sfc_p
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sfcpcorrected.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_amsua_seaice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_seaice_qc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_obs_error
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_obs_error.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_obs_error_halo
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_obs_error_halo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_gnssroBendMetOffice_obs_error
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobendmetoffice_obserror.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

if( ${ropp-ufo_FOUND} )
ecbuild_add_test( TARGET  test_ufo_qc_gnssro_obs_error_ROPP
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_obs_error_ropp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )
endif( ${ropp-ufo_FOUND} )

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_domain_check
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_domain_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_backgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ENVIRONMENT OOPS_TRAPFPE=1
                  ARGS    "testinput/gnssrobndnbam_backgroundcheck_qc.yaml"
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_qc_action_gnssroBndNBAM_obserr_inflation
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ENVIRONMENT OOPS_TRAPFPE=1
                  ARGS    "testinput/gnssrobndnbam_qc_action_obserr_inflation.yaml"
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

# Test bias correction classes

ecbuild_add_test( TARGET  test_ufo_bias_coeffs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/bias_coeff.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_coeffs_readwrite
                  SOURCES mains/TestUfoObsBias.cc
                  ARGS    "testinput/bias_coeff.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_coeff_increment
                  SOURCES mains/TestObsBiasIncrement.cc
                  ARGS    "testinput/bias_coeff_cov.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_covariance_amsua
                  SOURCES mains/TestObsBiasCovariance.cc
                  ARGS    "testinput/bias_coeff_cov.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_covariance_details
                  SOURCES mains/TestObsBiasCovarianceDetails.cc
                  ARGS    "testinput/bias_coeff_cov.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_single_channel
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/bias_single_channel.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET  test_ufo_bias_lin_op
                  SOURCES mains/TestLinearObsBiasOperator.cc
                  ARGS    "testinput/bias_linear_op.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET test_ufo_opr_aod_metoffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS "testinput/aod_metoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET test_ufo_opr_groundgnssMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS "testinput/groundgnssmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)


ecbuild_add_test( TARGET test_ufo_opr_groundgnssROPP
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS "testinput/groundgnssropp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET test_ufo_linopr_groundgnssMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/groundgnssmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET test_ufo_opr_sattcwv
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS "testinput/sattcwv.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ufo_test_data)

ecbuild_add_test( TARGET test_ufo_linopr_sattcwv
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sattcwv.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test predictors
ecbuild_add_test( TARGET  test_constant_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/constant_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_orbital_angle_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/orbital_angle_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_legendre_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/legendre_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_thickness_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/thickness_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_interpolate_data_from_file_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/interpolate_data_from_file_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_satellite_selector_predictor
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_Predictor.x
                  ARGS    "testinput/satellite_selector_predictor.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_Predictor.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_processamvqi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_processamvqi.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
CloudCoeff/Little_Endian/CloudCoeff.bin
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/NPOESS.IRland.EmisCoeff.bin
EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin
EmisCoeff/VIS_Ice/SEcategory/Little_Endian/NPOESS.VISice.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/NPOESS.VISland.EmisCoeff.bin
EmisCoeff/VIS_Snow/SEcategory/Little_Endian/NPOESS.VISsnow.EmisCoeff.bin
EmisCoeff/VIS_Water/SEcategory/Little_Endian/NPOESS.VISwater.EmisCoeff.bin
EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/USGS.IRland.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/USGS.VISland.EmisCoeff.bin
SpcCoeff/Little_Endian/hirs4_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/hirs4_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsr2_gcom-w1.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsr2_gcom-w1.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-c.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n15.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n18.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
SpcCoeff/Little_Endian/atms_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_npp.TauCoeff.bin
SpcCoeff/Little_Endian/atms_n20.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_n20.TauCoeff.bin
SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m11.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/seviri_m11.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-c.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n19.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD1_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD1_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD2_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD2_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD3_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD3_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD4_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD4_g15.TauCoeff.bin
SpcCoeff/Little_Endian/ssmis_f17.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f17.TauCoeff.bin
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g16.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g16.TauCoeff.bin
SpcCoeff/Little_Endian/abi_gr.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_gr.TauCoeff.bin
SpcCoeff/Little_Endian/ahi_himawari8.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ahi_himawari8.TauCoeff.bin
SpcCoeff/Little_Endian/radiometer_smap.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/radiometer_smap.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_n18.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n20.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n20.TauCoeff.bin
)

if( crtm_FOUND )
    # symlink to selected CRTM coefficient files.
    create_symlink_filename( ${UFO_CRTM_TESTFILES_PATH}
                             ${CMAKE_CURRENT_BINARY_DIR}/Data
                             ${crtm_test_input} )
endif()

#####################################################################
# Files for RTTOV tests
#####################################################################
if( ${rttov_FOUND} )
list( APPEND rttov_test_input
rttov7pred54L/rtcoef_noaa_19_amsua.dat
rttov7pred54L/rtcoef_noaa_20_atms.dat
rttov7pred54L/rtcoef_dmsp_17_ssmis.dat
rttov7pred54L/rtcoef_gcom-w_1_amsr2.dat
)

# Symlink all CRTM files
CREATE_SYMLINK_FILENAME( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${rttov_test_input} )
endif( ${rttov_FOUND} )

######
#Files for geos-aero tests
######

if( ${geos-aero_FOUND} )
list( APPEND geos-aero_test_data
Data/optics_BC.v1_3_.nc
Data/optics_BRC.v1_5_.nc
Data/optics_DU.v15_3_.nc
Data/optics_NI.v2_5_.nc
Data/optics_OC.v1_3_.nc
Data/optics_SS.v3_3_.nc
Data/optics_SU.v1_3_.nc
)

CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${geos-aero_test_data} )

list( APPEND geos-aero_test_input
testinput/geosaod.rc
testinput/Chem_MieRegistry.rc)

CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${geos-aero_test_input} )

endif( ${geos-aero_FOUND} )

add_subdirectory(testinput/instrumentTests)
