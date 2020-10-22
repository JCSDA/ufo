# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# macro to create a symlink from src to dst with just filename
function(CREATE_SYMLINK_FILENAME src dst)
    foreach (FILENAME ${ARGN})
        get_filename_component(filename ${FILENAME} NAME )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${filename} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK_FILENAME)

# Create Data directory for test input config and symlink all files
list( APPEND ufo_test_input
  testinput/abi_ahi_crtm.yaml
  testinput/adt.yaml
  testinput/aircraft.yaml
  testinput/airs_crtm.yaml
  testinput/airs_crtm_bc.yaml
  testinput/airs_qc_filters.yaml
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
  testinput/atms_crtm.yaml
  testinput/atms_crtm_bc.yaml
  testinput/atms_qc_filters.yaml
  testinput/atms_rttov_ops.yaml
  testinput/coolskin.yaml
  testinput/cris_crtm.yaml
  testinput/cris_crtm_bc.yaml
  testinput/cris_qc.yaml
  testinput/cris_qc_filters.yaml
  testinput/cris_qc_land.yaml
  testinput/empty.yaml
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
  testinput/function_nsstret.yaml
  testinput/function_obserrmean.yaml
  testinput/function_obserrmean_abi.yaml
  testinput/function_obserrmean_mhs.yaml
  testinput/function_obserrmean_atms.yaml
  testinput/function_scatret.yaml
  testinput/function_scatret_atms.yaml
  testinput/function_errfsdoei.yaml
  testinput/function_errfsdoei_atms.yaml
  testinput/function_velocity.yaml
  testinput/gmi_crtm.yaml
  testinput/genericprof.yaml
  testinput/geovals.yaml
  testinput/geos_aod.yaml
  testinput/geovals_spec.yaml
  testinput/gnssrobendmetoffice.yaml
  testinput/gnssrobendmetoffice_nopseudo.yaml
  testinput/gnssrobndropp1d.yaml
  testinput/gnssrobndropp2d.yaml
  testinput/gnssrobndnbam.yaml
  testinput/gnssro_obs_error.yaml
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
  testinput/gsisfcmodel.yaml
  testinput/hirs4_crtm.yaml
  testinput/iasi_crtm.yaml
  testinput/iasi_crtm_bc.yaml
  testinput/iasi_qc.yaml
  testinput/iasi_qc_filters.yaml
  testinput/locations.yaml
  testinput/metar_qc_filters.yaml
  testinput/mhs_crtm.yaml
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
  testinput/obsfilterdata.yaml
  testinput/omi_aura.yaml
  testinput/omi_aura_flipz.yaml
  testinput/ompsnp_npp.yaml
  testinput/ompsnp_npp_flipz.yaml
  testinput/ompsnp_npp_L127.yaml
  testinput/ompstc_npp.yaml
  testinput/ompstc_npp_flipz.yaml
  testinput/processwhere.yaml
  testinput/parallel_obs_distribution.yaml
  testinput/parameters.yaml
  testinput/qc_trackcheck.yaml
  testinput/qc_trackcheck_unittests.yaml
  testinput/qc_trackcheckship_unittests.yaml
  testinput/qc_trackcheckship_mainloop_unittests.yaml
  testinput/qc_backgroundcheck.yaml
  testinput/qc_boundscheck.yaml
  testinput/qc_velocitycheck.yaml
  testinput/qc_defer_to_post.yaml
  testinput/qc_derivative_dpdt.yaml
  testinput/qc_derivative_dxdt.yaml
  testinput/qc_differencecheck.yaml
  testinput/qc_gauss_thinning.yaml
  testinput/qc_gauss_thinning_unittests.yaml
  testinput/qc_met_office_buddy_check.yaml
  testinput/qc_met_office_buddy_check_unittests.yaml
  testinput/qc_met_office_buddy_pair_finder.yaml
  testinput/qc_recursive_splitter.yaml
  testinput/qc_oceancolor_preqc.yaml
  testinput/qc_poisson_disk_thinning.yaml
  testinput/qc_poisson_disk_thinning_unittests.yaml
  testinput/qc_preqc.yaml
  testinput/qc_thinning.yaml
  testinput/qc_temporal_thinning.yaml
  testinput/qc_temporal_thinning_unittests.yaml
  testinput/radiosonde.yaml
  testinput/radialvelocity.yaml
  testinput/reflectivity.yaml
  testinput/satwind.yaml
  testinput/sbuv2_n19.yaml
  testinput/sbuv2_n19_flipz.yaml
  testinput/sbuv2_n19_L127.yaml
  testinput/seaicefrac.yaml
  testinput/seaicethick.yaml
  testinput/chleuzintegr.yaml
  testinput/sea_surface_temp.yaml
  testinput/seviri_crtm.yaml
  testinput/smap_crtm.yaml
  testinput/sndrd1-4_crtm.yaml
  testinput/sfcpcorrected.yaml
  testinput/profileconsistencychecks_monolithicfilter.yaml
  testinput/profileconsistencychecks_OPScomparison.yaml
  testinput/profileconsistencychecks_wrongOPScomparison.yaml
  testinput/profileconsistencychecks_separatefilters.yaml
  testinput/profileconsistencychecks_unittests.yaml
  testinput/profileconsistencychecks_unittest_oneprofileMPI.yaml
  testinput/profileconsistencychecks_RH_obsfilter.yaml
  testinput/profileconsistencychecks_RH_OPScomparison.yaml
  testinput/profileconsistencychecks_UInterp_obsfilter.yaml
  testinput/profileconsistencychecks_UInterp_OPScomparison.yaml
  testinput/profileconsistencychecks_bkgqc_modobs_obsfilter.yaml
  testinput/profileconsistencychecks_bkgqc_repobs_obsfilter.yaml
  testinput/profileconsistencychecks_bkgqc_modobs_OPScomparison.yaml
  testinput/profileconsistencychecks_bkgqc_repobs_OPScomparison.yaml
  testinput/timeoper.yaml
  testinput/tprof.yaml
  testinput/tropomi_no2.yaml
  testinput/variables.yaml
  testinput/windprof.yaml
)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_test_input} )


# Add files to cmake resources
ecbuild_add_resources( TARGET ufo_test_scripts
                       SOURCES_PACK
                       ${ufo_test_input}
                     )

# Create Data directory for test data and symlink files
list( APPEND ufo_test_data
  atmosphere/abi_g16_geoval_2018041500_m_qc.nc4
  atmosphere/abi_g16_geoval_2019042306_m.nc4
  atmosphere/abi_g16_obsdiag_2018041500_m_qc.nc4
  atmosphere/ahi_himawari8_geoval_2019042306_m.nc4
  atmosphere/aircraft_geoval_2018041500_m.nc4
  atmosphere/aircraft_geoval_2018041500_s.nc4
  atmosphere/airs_aqua_geoval_2018041500_m.nc4
  atmosphere/airs_aqua_geoval_2018041500_m_qc.nc4
  atmosphere/airs_aqua_geoval_2018041500_m_unittest.nc4
  atmosphere/airs_aqua_obsdiag_2018041500_m_qc.nc4
  atmosphere/airs_aqua_obsdiag_2018041500_m_unittest.nc4
  atmosphere/airs_aqua_tlapmean.txt
  atmosphere/amsua_n19_geoval_2018041500_m_qc.nc4
  atmosphere/amsua_n19_geoval_2018041500_s_qc.nc4
  atmosphere/amsua_n19_geoval_2018041500_m.nc4
  atmosphere/amsua_n19_geoval_2018041500_s.nc4
  atmosphere/amsua_n19_obsdiag_2018041500_m_qc.nc4
  atmosphere/amsua_n19_obsdiag_2018041500_s_qc.nc4
  atmosphere/amsua_n19_tlapmean.txt
  atmosphere/aod_geoval_2018041500_m.nc4
  atmosphere/aod_geoval_2018041500_s.nc4
  atmosphere/atms_npp_geoval_2018041500_m.nc4
  atmosphere/atms_npp_geoval_2018041500_m_qc.nc4
  atmosphere/atms_npp_geoval_2018041500_s.nc4
  atmosphere/atms_npp_geoval_2018041500_s_qc.nc4
  atmosphere/atms_npp_obsdiag_2018041500_m_qc.nc4
  atmosphere/atms_npp_obsdiag_2018041500_s_qc.nc4
  atmosphere/atms_npp_tlapmean.txt
  atmosphere/cris-fsr_npp_geoval_2018041500_m.nc4
  atmosphere/cris-fsr_npp_geoval_2018041500_m_qc.nc4
  atmosphere/cris-fsr_npp_geoval_2018041500_m_unittest.nc4
  atmosphere/cris-fsr_npp_obsdiag_2018041500_m_qc.nc4
  atmosphere/cris-fsr_npp_obsdiag_2018041500_m_unittest.nc4
  atmosphere/cris-fsr_npp_tlapmean.txt
  atmosphere/geos_aod_geoval_2018041500_m.nc4
  atmosphere/geovals_atms_20191230T0000Z_benchmark.nc4
  atmosphere/gmi_gpm_geoval_2018041500_m.nc4
  atmosphere/gnssro_geoval_2018041500_1obs.nc4
  atmosphere/gnssro_geoval_2018041500_1obs_bending_angle.nc4
  atmosphere/gnssro_geoval_2018041500_s.nc4
  atmosphere/gnssro_geoval_2018041500_s_2d.nc4
  atmosphere/gnssro_geoval_2018041500_m.nc4
  atmosphere/gnssro_geoval_2018041500_l.nc4
  atmosphere/gnssro_geoval_2018041500_3prof.nc4
  atmosphere/gnssro_geoval_2019050700_1obs.nc4
  atmosphere/gnssro_geoval_2020050106_nopseudo.nc4
  atmosphere/gome_metop-a_geoval_2019101700_m.nc4
  atmosphere/gome_metop-a_geoval_flipz_2019101700_m.nc4
  atmosphere/groundgnss_geovals_20191230T0600Z.nc4
  atmosphere/gsisfc_tsen_geoval_2018041500_m.nc4
  atmosphere/gsisfc_uv_geoval_2018041500_m.nc4
  atmosphere/hirs4_metop-a_geoval_2018041500_m.nc4
  atmosphere/iasi_metop-a_obsdiag_2018041500_m.nc4
  atmosphere/iasi_metop-a_obsdiag_2018041500_m_qc.nc4
  atmosphere/iasi_metop-a_obsdiag_2018041500_m_unittest.nc4
  atmosphere/iasi_metop-a_geoval_2018041500_m.nc4
  atmosphere/iasi_metop-a_geoval_2018041500_m_qc.nc4
  atmosphere/iasi_metop-a_geoval_2018041500_m_unittest.nc4
  atmosphere/iasi_metop-a_tlapmean.txt
  atmosphere/met_office_thinning.nc4
  atmosphere/met_office_temporal_thinning_surface.nc4
  atmosphere/met_office_temporal_thinning_sonde.nc4
  atmosphere/mhs_n19_geoval_2018041500_m_qc.nc4
  atmosphere/mhs_n19_geoval_2018041500_m.nc4
  atmosphere/mhs_n19_obsdiag_2018041500_m_qc.nc4
  atmosphere/omi_aura_geoval_2019101700_m.nc4
  atmosphere/omi_aura_geoval_flipz_2019101700_m.nc4
  atmosphere/ompsnp_npp_geoval_2019101700_m.nc4
  atmosphere/ompsnp_npp_geoval_flipz_2019101700_m.nc4
  atmosphere/ompsnp_npp_geoval_2020051712_m.nc4
  atmosphere/ompstc8_npp_geoval_2019101700_m.nc4
  atmosphere/ompstc8_npp_geoval_flipz_2019101700_m.nc4
  atmosphere/ps_geovals_2018041500_0000.nc4
  atmosphere/radar_dbz_geoval_2019052222.nc4
  atmosphere/radar_rw_geoval_2019052222.nc4
  atmosphere/smap_geoval_2018041500_m.nc4
  atmosphere/satbias_crtm_in
  atmosphere/satbias_crtm_pc
  atmosphere/satwind_geoval_2018041500_m.nc4
  atmosphere/satwind_geoval_2018041500_s.nc4
  atmosphere/sbuv2_n19_geoval_2018041500_m.nc4
  atmosphere/sbuv2_n19_geoval_flipz_2018041500_m.nc4
  atmosphere/sbuv2_n19_geoval_2019101700_m.nc4
  atmosphere/sbuv2_n19_geoval_flipz_2019101700_m.nc4
  atmosphere/sbuv2_n19_geoval_2020051712_m.nc4
  atmosphere/seviri_m08_geoval_2018041500_m.nc4
  atmosphere/sfc_geoval_2018041500_m.nc4
  atmosphere/sndrd1_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd2_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd3_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd4_g15_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_s.nc4
  atmosphere/sondes_geoval_2018041500_vs.nc4
  atmosphere/ssmis_f18_geoval_2018041500_m.nc4
  atmosphere/tropomi_no2_geoval_2020090318_m.nc4
  atmosphere/windprof_geoval_2018041500_m.nc4
  marine/coolskin_fake_geovals_2018041500.nc
  marine/cryosat2-2018-04-15_geovals.nc
  marine/icec-2018-04-15_geovals.nc
  marine/Jason-2-2018-04-15_geovals.nc
  marine/profile_2018-04-15_geovals.nc
  marine/sst_obs-2018-04-15_geovals.nc
  marine/viirs_jpss1_oc_l2_2018-04-15_geovals.nc
  filters/filters_testdata.nc4
  filters/met_office_poisson_disk_thinning.nc4
  filters/met_office_buddy_check.nc4
  filters/poisson_disk_thinning_3x3x3x3_regular_grid.nc4
  filters/met_office_profile_consistency_checks.nc4
  filters/met_office_profile_consistency_checks_bkgqc_modobs.nc4
  filters/met_office_profile_consistency_checks_bkgqc_repobs.nc4
  filters/met_office_profile_consistency_checks_rh.nc4
  filters/met_office_profile_consistency_checks_rh_wrong.nc4
  filters/met_office_profile_consistency_checks_unittests.nc4
  filters/met_office_profile_consistency_checks_oneprofile.nc4
  filters/met_office_profile_consistency_checks_oneprofile_wrong.nc4
  filters/met_office_profile_consistency_checks_uinterp.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput ${CMAKE_CURRENT_BINARY_DIR}/Data ${ufo_test_data} )	

# IODA observation files 
list( APPEND ufo_ioda_test_data
      testinput_tier_1.tar.gz)

# Set URL for IODA test files
set(UFO_IODA_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com)

# If local path to testfiles is defined don't download
if (DEFINED ENV{LOCAL_PATH_TESTFILES_IODA})
        set(LOCAL_PATH_TESTFILES_IODA "$ENV{LOCAL_PATH_TESTFILES_IODA}")
endif()

if( NOT DEFINED LOCAL_PATH_TESTFILES_IODA )

  message(STATUS "LOCAL_PATH_TESTFILES_IODA is not defined, download test files")
  if ( NOT DEFINED UFO_IODA_TESTFILES_BRANCH)

    # Get the current git branch of ioda repo
    execute_process(
      COMMAND git rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH_UFO
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "UFO is in branch: ${GIT_BRANCH_UFO}")

  else()
    message(STATUS "Branch name provided by user")
    set(GIT_BRANCH_UFO ${UFO_IODA_TESTFILES_BRANCH})
  endif()

  # Check whether the URL exists or not
  set(ECBUILD_DOWNLOAD_BASE_URL ${UFO_IODA_DOWNLOAD_BASE_URL}/ioda)

  ecbuild_check_urls(NAMES ${GIT_BRANCH_UFO}/${ufo_ioda_test_data}
                     RESULT UFO_SPECIFIC_TEST_FILES)

  # Set distant directory
  if(${UFO_SPECIFIC_TEST_FILES} MATCHES 0)
    # Download and extract new test files (distant directory = git branch)
    set(DIRNAME ${GIT_BRANCH_UFO})
    message(STATUS "GIT_BRANCH_UFO found, will download: ${GIT_BRANCH_UFO}")
  else()
    # Download and extract develop test files (distant directory = develop)
    set(DIRNAME "develop")
    message(STATUS "GIT_BRANCH_UFO not found, will download develop")
  endif()

  message(STATUS "Test data will be download from: ${UFO_IODA_DOWNLOAD_BASE_URL}/${DIRNAME}")

  set(UFO_IODA_REP_NAME ioda)
  set(UFO_IODA_TESTFILES_NAME ${ufo_ioda_test_data})
  set(UFO_IODA_BRANCH_NAME ${DIRNAME})
  set(UFO_IODA_TESTFILES_PATH ${CMAKE_BINARY_DIR}/test_data/ioda/${DIRNAME})
  message(STATUS "Save data to: ${TESTFILE_DIR_IODA}")
  file(MAKE_DIRECTORY ${UFO_IODA_TESTFILES_PATH})

  # Create download script for get_ioda_test_data test
  set ( FILENAME ufo_ioda_data_downloader.py)
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
  list( APPEND bin_ufo_ioda_test_download_scripts_deps ${DEST_FILE} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()

  add_custom_target( bin_ufo_ioda_test_download_scripts ALL
      COMMAND chmod +x ${bin_ufo_ioda_test_download_scripts_deps}
      DEPENDS ${bin_ufo_ioda_test_download_scripts_deps})

  ecbuild_add_test( TARGET    ufo_get_ioda_test_data
                    TYPE      SCRIPT
                    COMMAND    ${CMAKE_BINARY_DIR}/bin/ufo_ioda_data_downloader.py
                    ARGS      testoutput/download_test_ufo.log)

else()
  set(UFO_IODA_TESTFILES_PATH ${LOCAL_PATH_TESTFILES_IODA})
  message(STATUS "use LOCAL_PATH_TESTFILES_IODA: ${LOCAL_PATH_TESTFILES_IODA}")
endif()

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${UFO_IODA_TESTFILES_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/ioda)

#####################################################################
# Install headers used in tests of interfaces that may be implemented by clients.

install( FILES ufo/ObsDiagnostics.h ufo/ObsFilters.h ufo/ObsFunction.h
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

ecbuild_add_executable( TARGET  test_ProfileConsistencyChecks.x
                        SOURCES mains/TestProfileConsistencyChecks.cc
                        LIBS    ufo
                       )

####################################################################
# Test Locations, GeoVaLs with oops tests:

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo 
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_locations
                  SOURCES mains/TestLocations.cc
                  ARGS    "testinput/locations.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test UFO-specific classes:

ecbuild_add_test( TARGET  test_ufo_obsfilterdata
                  SOURCES mains/TestObsFilterData.cc
                  ARGS    "testinput/obsfilterdata.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

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
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test Obs Operators and TLAD

ecbuild_add_test( TARGET  test_ufo_opr_gsi_sfc_model
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gsisfcmodel.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_opr_crtm_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/amsua_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_bc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/amsua_crtm_bc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/amsua_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_bc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/amsua_crtm_bc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_atms
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/atms_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_atms
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/atms_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_gmi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/gmi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_gmi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/gmi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_smap
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/smap_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_smap
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/smap_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_seviri
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/seviri_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_seviri
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/seviri_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/cris_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/cris_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_mhs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/mhs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_mhs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/mhs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_iasi
                      MPI     2
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/iasi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_iasi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/iasi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_sndrd1-4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/sndrd1-4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_sndrd1-4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/sndrd1-4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_airs
                      MPI     4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/airs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_airs
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/airs_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x 
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_hirs4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/hirs4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_hirs4
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/hirs4_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_crtm_abi_ahi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/abi_ahi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_crtm_abi_ahi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/abi_ahi_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_aod_crtm
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/aod_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_aod_crtm
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/aod_crtm.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x 
                      TEST_DEPENDS ufo_get_ioda_test_data )

    if( ${GEOS-AERO_FOUND} )

      ecbuild_add_test( TARGET  test_ufo_opr_aod_luts
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/aod_luts.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

      ecbuild_add_test( TARGET  test_ufo_linopr_aod_luts
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/aod_luts.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x 
                      TEST_DEPENDS ufo_get_ioda_test_data )

   endif( ${GEOS-AERO_FOUND} )


endif( crtm_FOUND )



if( ${RTTOV_FOUND} )
  ecbuild_add_test( TARGET  test_ufo_opr_rttov_atms
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                    ARGS    "testinput/atms_rttov_ops.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperator.x
                    TEST_DEPENDS ufo_get_ioda_test_data )

  ecbuild_add_test( TARGET  test_ufo_linopr_rttov_atms
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                    ARGS    "testinput/atms_rttov_ops.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    DEPENDS test_ObsOperatorTLAD.x 
                    TEST_DEPENDS ufo_get_ioda_test_data )
endif( ${RTTOV_FOUND} )


ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_gome_metop-a
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gome_metop-a.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_gome_metop-a_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gome_metop-a_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_omi_aura
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/omi_aura.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_omi_aura_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/omi_aura_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompsnp_npp_L127
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompsnp_npp_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompstc_npp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompstc_npp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_ompstc_npp_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/ompstc_npp_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19_flipz
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19_flipz.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_atminterplay_sbuv2_n19_L127
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19_L127.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_vertinterp_radiosonde
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radiosonde.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_avgkernel_tropomi_no2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/tropomi_no2.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_radiosonde
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/radiosonde.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_timeoper_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/timeoper.yaml"
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_vertinterp_aircraft_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/aircraft.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_aircraft
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/aircraft.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_vertinterp_satwind
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/satwind.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_timeoper
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/timeoper.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_timeoper
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/timeoper.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroRef
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssroref.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_0obs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_0obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_1obs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_1obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM_1obs_2pe
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_1obs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroRef
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssroref.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_derivative_ascent
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/qc_derivative_dpdt.yaml"
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_derivative_speed
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/qc_derivative_dxdt.yaml"
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_super_refraction_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam_super_refraction.yaml"
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBndNBAM
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBendMetOffice_pseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobendmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x 
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBendMetOffice_pseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobendmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_gnssroBendMetOffice_nopseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobendmetoffice_nopseudo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x 
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_gnssroBendMetOffice_nopseudo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobendmetoffice_nopseudo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

if( ${ROPP-UFO_FOUND} )
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
endif( ${ROPP-UFO_FOUND} )

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
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_seaicefrac
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicefrac.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_seaicethick
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seaicethick.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_seaicethick
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicethick.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_chleuzintegr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/chleuzintegr.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_identity_sst
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_identity_sst
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_adt
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/adt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_adt
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/adt.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x 
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_coolskin
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/coolskin.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_linopr_coolskin
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/coolskin.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperatorTLAD.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

# METAR (sfc. obs) UFO tests
ecbuild_add_test( TARGET  test_ufo_metar_QC
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/metar_qc_filters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

if( ${GSW_FOUND} )
    ecbuild_add_test( TARGET  test_ufo_opr_insitutemperature
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/tprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_insitutemperature
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/tprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_opr_marinevertinterp
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/genericprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperator.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_linopr_marinevertinterp
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/genericprof.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsOperatorTLAD.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

endif( ${GSW_FOUND} )

if( ${GEOS-AERO_FOUND} )

ecbuild_add_test( TARGET  test_ufo_geos_aero_opr
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/geos_aod.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geos_aero_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/geos_aod.yaml"
                  LIBS    ufo )

endif( ${GEOS-AERO_FOUND} )


ecbuild_add_test( TARGET  test_ufo_opr_radarreflectivity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/reflectivity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_opr_radarradialvelocity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radialvelocity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test UFO ObsFilters (generic)

ecbuild_add_test( TARGET  test_ufo_qc_gen_processwhere
                  SOURCES mains/TestProcessWhere.cc
                  ARGS    "testinput/processwhere.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_qc_gen_backgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_backgroundcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_gen_boundscheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_boundscheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_velocitybounds
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_velocitycheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gen_defer_to_post
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_defer_to_post.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_gen_differencecheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_differencecheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gauss_thinning
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_gauss_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_check
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_met_office_buddy_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_oceancolor_preqc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_oceancolor_preqc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_poisson_disk_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_poisson_disk_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_preqc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_preqc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_temporalthinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_temporal_thinning.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_trackcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_trackcheck.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_monolithicfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_monolithicfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_separatefilters
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_separatefilters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_unittests
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_unittest_oneprofileMPI
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_unittest_oneprofileMPI.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_RH_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_RH_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_UInterp_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_UInterp_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_bkgqc_repobs_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_bkgqc_repobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_bkgqc_modobs_obsfilter
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/profileconsistencychecks_bkgqc_modobs_obsfilter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )

# Test UFO ObsFilters (specific)

ecbuild_add_test( TARGET  test_ufo_gaussianthinning
                  SOURCES mains/TestGaussianThinning.cc
                  ARGS    "testinput/qc_gauss_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_recursivesplitter
                  SOURCES mains/TestRecursiveSplitter.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/qc_recursive_splitter.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_temporalthinning
                  SOURCES mains/TestTemporalThinning.cc
                  ARGS    "testinput/qc_temporal_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheck_unittests
                  SOURCES mains/TestTrackCheck.cc
                  ARGS    "testinput/qc_trackcheck_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheckship_unittests
                  SOURCES mains/TestTrackCheckShip.cc
                  ARGS    "testinput/qc_trackcheckship_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_trackcheckship_mainloop_unittests
                  SOURCES mains/TestTrackCheckShipMainLoop.cc
                  ARGS    "testinput/qc_trackcheckship_mainloop_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_poisson_disk_thinning_unittests
                  SOURCES mains/TestPoissonDiskThinning.cc
                  ARGS    "testinput/qc_poisson_disk_thinning_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_bkgqc_modobs_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_bkgqc_modobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_bkgqc_repobs_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_bkgqc_repobs_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_RH_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_RH_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_wrongOPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_wrongOPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_profileconsistencychecks_UInterp_OPScomparison
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ProfileConsistencyChecks.x
                  ARGS    "testinput/profileconsistencychecks_UInterp_OPScomparison.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ProfileConsistencyChecks.x )

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_check_unittests
                  SOURCES mains/TestMetOfficeBuddyCheck.cc
                  ARGS    "testinput/qc_met_office_buddy_check_unittests.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_met_office_buddy_pair_finder
                  SOURCES mains/TestMetOfficeBuddyPairFinder.cc
                  ARGS    "testinput/qc_met_office_buddy_pair_finder.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test Functions

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_iasi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_iasi.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_airs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_airs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clouddetect_cris
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clouddetect_cris.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwmatchidx
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwmatchidx.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwmatchidx_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwmatchidx_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwretmean
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwretmean.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwretmean_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwretmean_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret_hofx
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_hofx.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret_hofx_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_hofx_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret_obsval
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_obsval.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_clwret_obsval_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_clwret_obsval_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errfgrosschk
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfgrosschk.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errfjsfc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfjsfc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errflat
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errflat.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errftopo
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errftopo.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errftransmittop
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errftransmittop.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errfwavenum
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfwavenum.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_hydrometeorchk
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_hydrometeorchk.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_hydrometeorchk_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_hydrometeorchk_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_nsstret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_nsstret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_obserrmean
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_abi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_abi.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_mhs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_mhs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_obserrmean_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_obserrmean_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_scatret
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_scatret.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_scatret_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_scatret_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errfsdoei
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfsdoei.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_errfsdoei_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_errfsdoei_atms.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_function_velocity
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
                  ARGS    "testinput/function_velocity.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFunction.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test Diagnostics

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_airs_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_airs_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_airs_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_airs_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_amsua_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_amsua_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_amsua_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_amsua_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_atms_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_atms_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_atms_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_atms_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_cris_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_cris_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_cris_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_cris_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_iasi_optics
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_iasi_optics.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_iasi_jacobian
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
                      ARGS    "testinput/obsdiag_crtm_iasi_jacobian.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsDiagnostics.x
                      TEST_DEPENDS ufo_get_ioda_test_data )
endif( crtm_FOUND )

# Test Parameters
ecbuild_add_test( TARGET  test_ufo_parameters
                  SOURCES mains/TestParameters.cc
                  ARGS    "testinput/parameters.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test piecewise linear interpolation
ecbuild_add_test( TARGET  test_ufo_piecewise_linear_interpolation
                  SOURCES mains/TestPiecewiseLinearInterpolation.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test parallel observation distribution
ecbuild_add_test( TARGET  test_ufo_parallel_obs_distribution
                  MPI     4
                  SOURCES mains/TestParallelObsDistribution.cc
                  ARGS    "testinput/parallel_obs_distribution.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test QC for specific instruments

if( crtm_FOUND )
    ecbuild_add_test( TARGET  test_ufo_qc_airs_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/airs_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_clwretmw
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_clwretmw.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_miss_val
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc_miss_val.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_atms_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/atms_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_cris_land
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/cris_qc_land.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_function_scattering
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_iasi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/iasi_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_iasi_filters
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/iasi_qc_filters.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

    ecbuild_add_test( TARGET  test_ufo_qc_amsua_allsky_gsi
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                      ARGS    "testinput/amsua_allsky_gsi_qc.yaml"
                      ENVIRONMENT OOPS_TRAPFPE=1
                      DEPENDS test_ObsFilters.x
                      TEST_DEPENDS ufo_get_ioda_test_data )

endif( crtm_FOUND )

if( ${ROPP-UFO_FOUND} )
ecbuild_add_test( TARGET  test_ufo_qc_gnssroBndROPP1D
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x )
endif( ${ROPP-UFO_FOUND} )

ecbuild_add_test( TARGET  test_ufo_opr_sfc_p
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sfcpcorrected.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_amsua_seaice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_seaice_qc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_obs_error
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_obs_error.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_domain_check
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_domain_check.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_gnssro_backgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ENVIRONMENT OOPS_TRAPFPE=1
                  ARGS    "testinput/gnssrobndnbam_backgroundcheck_qc.yaml"
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_qc_action_gnssroBndNBAM_obserr_inflation
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ENVIRONMENT OOPS_TRAPFPE=1
                  ARGS    "testinput/gnssrobndnbam_qc_action_obserr_inflation.yaml"
                  DEPENDS test_ObsFilters.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

# Test bias correction classes

ecbuild_add_test( TARGET  test_ufo_bias_airs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/airs_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_amsua
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_atms
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/atms_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_cris
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/cris_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_iasi
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsBias.x
                  ARGS    "testinput/iasi_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsBias.x
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_increment_amsua
                  SOURCES mains/TestObsBiasIncrement.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_covariance_amsua
                  SOURCES mains/TestObsBiasCovariance.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET  test_ufo_bias_covariance_details
                  SOURCES mains/TestObsBiasCovarianceDetails.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ioda_test_data )

ecbuild_add_test( TARGET test_ufo_opr_groundgnssMetOffice
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS "testinput/groundgnssmetoffice.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_ObsOperator.x )
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
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
SpcCoeff/Little_Endian/atms_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_npp.TauCoeff.bin
SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
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
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g16.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g16.TauCoeff.bin
SpcCoeff/Little_Endian/ahi_himawari8.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ahi_himawari8.TauCoeff.bin
SpcCoeff/Little_Endian/radiometer_smap.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/radiometer_smap.TauCoeff.bin
)

if( crtm_FOUND )
    # symlink to selected CRTM coefficient files.
    # crtm_COEFFICIENT_DIR is exported by crtm-config.cmake
    create_symlink_filename( ${crtm_COEFFICIENT_DIR}
                             ${CMAKE_CURRENT_BINARY_DIR}/Data
                             ${crtm_test_input} )
endif()

#####################################################################
# Files for RTTOV tests
#####################################################################
if( ${RTTOV_FOUND} )
list( APPEND rttov_test_input
rttov7pred54L/rtcoef_noaa_19_amsua.dat
rttov7pred54L/rtcoef_noaa_20_atms.dat
)

# Symlink all CRTM files
CREATE_SYMLINK_FILENAME( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${rttov_test_input} )
endif( ${RTTOV_FOUND} )

######
#Files for geos-aero tests
######

if( ${GEOS-AERO_FOUND} )
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

endif( ${GEOS-AERO_FOUND} )
