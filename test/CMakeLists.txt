# (C) Copyright 2017-2018 UCAR.
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
  testinput/amsua_crtm.yaml
  testinput/amsua_crtm_bc.yaml
  testinput/amsua_rttov.yaml
  testinput/amsua_qc.yaml
  testinput/iasi_qc.yaml
  testinput/iasi_obsfunc_qc.yaml
  testinput/amsua_allsky_gsi_qc.yaml
  testinput/amsua_seaice_qc.yaml
  testinput/atms_crtm.yaml
  testinput/gmi_crtm.yaml
  testinput/aod_crtm.yaml
  testinput/coolskin.yaml
  testinput/cris_crtm.yaml
  testinput/genericprof.yaml
  testinput/geovals.yaml
  testinput/gnssrobndropp1d.yaml
  testinput/gnssrobndropp2d.yaml
  testinput/gnssrobndnbam.yaml
  testinput/gnssrobndnbam_obs_error.yaml
  testinput/gnssroref.yaml
  testinput/gnssro_domain_check.yaml
  testinput/hirs4_crtm.yaml
  testinput/iasi_crtm.yaml
  testinput/gsisfcmodel.yaml
  testinput/surface_ps_qc.yaml
  testinput/locations.yaml
  testinput/mhs_crtm.yaml
  testinput/obsfilterdata.yaml
  testinput/qc_boundscheck.yaml
  testinput/qc_cris.yaml
  testinput/qc_differencecheck.yaml
  testinput/qc_domainchk_zenith.yaml
  testinput/qc_function_windspeed.yaml
  testinput/qc_thinning.yaml
  testinput/qc_backgroundcheck.yaml
  testinput/qc_thinning_rad.yaml
  testinput/qc_gauss_thinning.yaml
  testinput/radiosonde.yaml
  testinput/radiosonde_qc.yaml
  testinput/satwind.yaml
  testinput/sbuv2_n19.yaml
  testinput/scatwind.yaml
  testinput/seaicefrac.yaml
  testinput/seaicethick.yaml
  testinput/sea_surface_temp.yaml
  testinput/seviri_crtm.yaml
  testinput/sndrd1-4_crtm.yaml
  testinput/timeoper.yaml
  testinput/tprof.yaml
  testinput/qc_blacklist_latitude.yaml
  testinput/windprof.yaml
  testinput/qc_domain_latitude.yaml
  testinput/qc_domain_minpress.yaml
  testinput/obsdiag_crtm_jacobian.yaml
  testinput/function_velocity.yaml
  testinput/reflectivity.yaml
  testinput/radialvelocity.yaml
  testinput/processwhere.yaml
  testinput/variables.yaml
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
  atmosphere/aircraft_geoval_2018041500_m.nc4
  atmosphere/aircraft_geoval_2018041500_s.nc4
  atmosphere/airs_aqua_geoval_2018041500_m.nc4
  atmosphere/amsua_n19_geoval_2018041500_m.nc4
  atmosphere/amsua_n19_geoval_2018041500_s.nc4
  atmosphere/atms_npp_geoval_2018041500_m.nc4
  atmosphere/atms_npp_geoval_2018041500_s.nc4
  atmosphere/gmi_gpm_geoval_2018041500_m.nc4
  atmosphere/aod_geoval_2018041500_m.nc4
  atmosphere/aod_geoval_2018041500_s.nc4
  atmosphere/cris-fsr_npp_geoval_2018041500_m.nc4
  atmosphere/geovals_rttov.nc4
  atmosphere/gnssro_geoval_2018041500_s.nc4
  atmosphere/gnssro_geoval_2018041500_s_2d.nc4
  atmosphere/gnssro_geoval_2018041500_m.nc4
  atmosphere/gnssro_geoval_2018041500_l.nc4
  atmosphere/hirs4_metop-a_geoval_2018041500_m.nc4
  atmosphere/iasi_metop-a_geoval_2018041500_m.nc4
  atmosphere/mhs_n19_geoval_2018041500_m.nc4
  atmosphere/satbias_crtm_in
  atmosphere/satwind_geoval_2018041500_m.nc4
  atmosphere/sbuv2_n19_geoval_2018041500_m.nc4
  atmosphere/satwind_geoval_2018041500_s.nc4
  atmosphere/scatwind_geoval_2018041500_m.nc4
  atmosphere/sfc_geoval_2018041500_m.nc4
  atmosphere/ps_geovals_2018041500_0000.nc4
  atmosphere/seviri_m08_geoval_2018041500_m.nc4
  atmosphere/sndrd1_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd2_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd3_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd4_g15_geoval_2018041500_m.nc4
  atmosphere/ahi_himawari8_geoval_2019042306_m.nc4
  atmosphere/abi_g16_geoval_2019042306_m.nc4
  atmosphere/sondes_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_s.nc4
  atmosphere/sondes_geoval_2018041500_vs.nc4
  atmosphere/windprof_geoval_2018041500_m.nc4
  atmosphere/ssmis_f18_geoval_2018041500_m.nc4
  marine/coolskin_fake_geovals_2018041500.nc
  marine/cryosat2-2018-04-15_geovals.nc
  marine/icec-2018-04-15_geovals.nc
  marine/Jason-2-2018-04-15_geovals.nc
  marine/profile_2018-04-15_geovals.nc
  marine/sst_obs-2018-04-15_geovals.nc
  atmosphere/radar_dbz_geoval_2019052222.nc4
  atmosphere/radar_rw_geoval_2019052222.nc4
  filters/where_testdata.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput ${CMAKE_CURRENT_BINARY_DIR}/Data ${ufo_test_data} )

# Add ioda obs test data
list (APPEND ioda_obs_test_data
  atmosphere/aircraft_obs_2018041500_m.nc4
  atmosphere/aircraft_obs_2018041500_s.nc4
  atmosphere/airs_aqua_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_m_bc.nc4
  atmosphere/gmi_gpm_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_s.nc4
  atmosphere/atms_npp_obs_2018041500_m.nc4
  atmosphere/ioda_metop_2_amsua.nc4
  atmosphere/aod_obs_2018041500_m.nc4
  atmosphere/aod_obs_2018041500_s.nc4
  atmosphere/cris-fsr_npp_obs_2018041500_m.nc4
  atmosphere/gnssro_obs_2018041500_s.nc4
  atmosphere/gnssro_obs_2018041500_m.nc4
  atmosphere/gnssro_obs_2018041500_l.nc4
  atmosphere/hirs4_metop-a_obs_2018041500_m.nc4
  atmosphere/iasi_metop-a_obs_2018041500_m.nc4
  atmosphere/mhs_n19_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_s.nc4
  atmosphere/sbuv2_n19_obs_2018041500_m.nc4
  atmosphere/scatwind_obs_2018041500_m.nc4
  atmosphere/seviri_m08_obs_2018041500_m.nc4
  atmosphere/sndrd1_g15_obs_2018041500_m.nc4
  atmosphere/sndrd2_g15_obs_2018041500_m.nc4
  atmosphere/sndrd3_g15_obs_2018041500_m.nc4
  atmosphere/sndrd4_g15_obs_2018041500_m.nc4
  atmosphere/ahi_himawari8_obs_2019042306_m.nc4
  atmosphere/abi_g16_obs_2019042306_m.nc4
  atmosphere/ssmis_f18_obs_2018041500_m.nc4

  atmosphere/sfc_obs_2018041500_m.nc4

  atmosphere/sondes_obs_2018041500_m.nc4
  atmosphere/sondes_obs_2018041500_s.nc4
  atmosphere/sondes_obs_2018041500_vs.nc4
  atmosphere/windprof_obs_2018041500_m.nc4
  marine/coolskin_fake_obs_2018041500.nc
  marine/cryosat2-2018-04-15.nc
  marine/icec-2018-04-15.nc
  marine/Jason-2-2018-04-15.nc
  marine/profile_2018-04-15.nc
  marine/sst_obs-2018-04-15.nc4

  atmosphere/radar_dbz_obs_2019052222.nc4
  atmosphere/radar_rw_obs_2019052222.nc4

  )

CREATE_SYMLINK_FILENAME( ${ioda_SOURCE_DIR}/test/testinput
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${ioda_obs_test_data} )

# Add ufo obsdiag test data
list (APPEND ufo_obsdiag_test_data
  atmosphere/iasi_metop-a_obsdiag_2018041500_m.nc4
  )

CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${ufo_obsdiag_test_data} )

#####################################################################

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

####################################################################

ecbuild_add_test( TARGET  test_ufo_function_velocity
                  SOURCES mains/TestObsFunction.cc
                  ARGS    "testinput/function_velocity.yaml"
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_obsdiag_crtm_jacobian
                  SOURCES mains/TestObsDiagnostics.cc
                  ARGS    "testinput/obsdiag_crtm_jacobian.yaml"
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_qc_domain_latitude
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_domain_latitude.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_locations
                  SOURCES mains/TestLocations.cc
                  ARGS    "testinput/locations.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_obsfilterdata
                  SOURCES mains/TestObsFilterData.cc
                  ARGS    "testinput/obsfilterdata.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_variables
                  SOURCES mains/TestVariables.cc
                  ARGS    "testinput/variables.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_processwhere
                  SOURCES mains/TestProcessWhere.cc
                  ARGS    "testinput/processwhere.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_obs_bias
                  SOURCES mains/TestObsBias.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_obs_bias_increment
                  SOURCES mains/TestObsBiasIncrement.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_obs_bias_covariance
                  SOURCES mains/TestObsBiasCovariance.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_qc_domain_minpress
	          COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
	          ARGS    "testinput/qc_domain_minpress.yaml"
	          DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_gsi_sfc_model_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gsisfcmodel.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_surface_ps_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/surface_ps_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_amsua_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/amsua_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_amsua_opr_crtm_bc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_amsua_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_iasi_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/iasi_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_iasi_obsfunc_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/iasi_obsfunc_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_amsua_allsky_gsi_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_allsky_gsi_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_amsua_seaice_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_seaice_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_amsua_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/amsua_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_amsua_tlad_crtm_bc
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/amsua_crtm_bc.yaml"
                  LIBS    ufo )

#ecbuild_add_test( TARGET  test_ufo_atms_opr_crtm
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
#                  ARGS    "testinput/atms_crtm.yaml"
#                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_atms_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/atms_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_gmi_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gmi_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_gmi_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gmi_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_seviri_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seviri_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_seviri_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seviri_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_cris_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/cris_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_cris_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/cris_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_mhs_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/mhs_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_mhs_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/mhs_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_iasi_opr_crtm
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/iasi_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_iasi_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/iasi_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_sndrd1-4_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sndrd1-4_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_sndrd1-4_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sndrd1-4_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_airs_opr_crtm
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/airs_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_airs_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/airs_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_hirs4_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/hirs4_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_hirs4_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/hirs4_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_abi_ahi_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/abi_ahi_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_abi_ahi_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/abi_ahi_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )


if( ${RTTOV_FOUND} )
  ecbuild_add_test( TARGET  test_ufo_amsua_opr_rttov
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                    ARGS    "testinput/amsua_rttov.yaml"
                    DEPENDS test_ObsOperator.x )

  ecbuild_add_test( TARGET  test_ufo_amsua_tlad_rttov
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                    ARGS    "testinput/amsua_rttov.yaml"
                    DEPENDS test_ObsOperatorTLAD.x )
endif( ${RTTOV_FOUND} )

ecbuild_add_test( TARGET  test_ufo_sbuv2_n19_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sbuv2_n19.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_radiosonde_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radiosonde.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_radiosonde_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/radiosonde_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_radiosonde_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/radiosonde.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_timeoper_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/timeoper.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_aircraft_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/aircraft.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_aircraft_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/aircraft.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssroref.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssroref.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )


ecbuild_add_test( TARGET  test_ufo_gnssroBndNBAM_opr
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndNBAM_tlad
                  MPI     2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndnbam.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndNBAM_obs_error
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobndnbam_obs_error.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_gnssro_domain_check
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssro_domain_check.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_satwind_opr
                  MPI     4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/satwind.yaml"
                  DEPENDS test_ObsOperator.x )

#ecbuild_add_test( TARGET  test_ufo_windprof_opr
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
#                  ARGS    "testinput/windprof.yaml"
#                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_qc_boundscheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_boundscheck.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_cris
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_cris.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_differencecheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_differencecheck.yaml"
                  DEPENDS test_ObsFilters.x )
      
ecbuild_add_test( TARGET  test_ufo_qc_backgroundcheck
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_backgroundcheck.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_domainchk_zenith
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_domainchk_zenith.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_function_windspeed
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_function_windspeed.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_function_scattering
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/amsua_qc.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_thinning.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_thinning_rad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_thinning_rad.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_gauss_thinning
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_gauss_thinning.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_qc_blacklist_latitude
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/qc_blacklist_latitude.yaml"
                  DEPENDS test_ObsFilters.x )

if( ${ROPP-UFO_FOUND} )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_qc
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  DEPENDS test_ObsFilters.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )
endif( ${ROPP-UFO_FOUND} )

# Marine UFO tests
ecbuild_add_test( TARGET  test_ufo_seaicefrac_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seaicefrac.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_seaicefrac_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicefrac.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_seaicethick_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/seaicethick.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_seaicethick_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/seaicethick.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_sst_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_sst_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/sea_surface_temp.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_adt_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/adt.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_adt_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/adt.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_CoolSkin_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/coolskin.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_CoolSkin_tlad
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/coolskin.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )


# Marine UFO tests
if( ${GSW_FOUND} )
    ecbuild_add_test( TARGET  test_ufo_insitutemperature_opr
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/tprof.yaml"
                      DEPENDS test_ObsOperator.x )

    ecbuild_add_test( TARGET  test_ufo_insitutemperature_tlad
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/tprof.yaml"
                      DEPENDS test_ObsOperatorTLAD.x )

    ecbuild_add_test( TARGET  test_ufo_marinevertinterp_opr
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                      ARGS    "testinput/genericprof.yaml"
                      DEPENDS test_ObsOperator.x )

    ecbuild_add_test( TARGET  test_ufo_marinevertinterp_tlad
                      COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                      ARGS    "testinput/genericprof.yaml"
                      DEPENDS test_ObsOperatorTLAD.x )

endif( ${GSW_FOUND} )

ecbuild_add_test( TARGET  test_ufo_aod_opr_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/aod_crtm.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_aod_tlad_crtm
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
                  ARGS    "testinput/aod_crtm.yaml"
                  DEPENDS test_ObsOperatorTLAD.x )

ecbuild_add_test( TARGET  test_ufo_radarreflectivity_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/reflectivity.yaml"
                  DEPENDS test_ObsOperator.x )

ecbuild_add_test( TARGET  test_ufo_radarradialvelocity_opr
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
                  ARGS    "testinput/radialvelocity.yaml"
                  DEPENDS test_ObsOperator.x )

#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
CloudCoeff/Little_Endian/CloudCoeff.bin
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/MW_Water/RSS_Emis/delta_EW_V5_B.dat
EmisCoeff/MW_Water/RSS_Emis/deW_harm_coeffs_V9A_MI.dat
EmisCoeff/MW_Water/RSS_Emis/dew_phi_VH34_harmonic_tab.dat
EmisCoeff/MW_Water/RSS_Emis/dw_res_tab_spline.dat
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
)

# Symlink all CRTM files
CREATE_SYMLINK_FILENAME( ${crtm_SOURCE_DIR}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${crtm_test_input} )

#####################################################################
# Files for RTTOV tests
#####################################################################
if( ${RTTOV_FOUND} )
list( APPEND rttov_test_input
rttov7pred54L/rtcoef_noaa_19_amsua.dat
)

# Symlink all CRTM files
CREATE_SYMLINK_FILENAME( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${rttov_test_input} )
endif( ${RTTOV_FOUND} )
