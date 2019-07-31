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
  testinput/adt.yaml
  testinput/aircraft.yaml
  testinput/airs_crtm.yaml
  testinput/amsua_crtm.yaml
  testinput/amsua_rttov.yaml
  testinput/amsua_qc.yaml
  testinput/atms_crtm.yaml
  testinput/gmi_crtm.yaml
  testinput/aod_crtm.yaml
  testinput/coolskin.yaml
  testinput/cris_crtm.yaml
  testinput/genericprof.yaml
  testinput/geovals.yaml
  testinput/gnssrobndropp1d.yaml
  testinput/gnssrobndropp2d.yaml
  testinput/gnssrobndgsi.yaml
  testinput/gnssroref.yaml
  testinput/hirs4_crtm.yaml
  testinput/iasi_crtm.yaml
  testinput/gsisfcmodel.yaml
#  testinput/surface_ps_qc.yaml
  testinput/locations.yaml
  testinput/mhs_crtm.yaml
  testinput/qc_boundscheck.yaml
  testinput/qc_thinning.yaml
  testinput/radiosonde.yaml
  testinput/radiosonde_qc.yaml
  testinput/satwind.yaml
  testinput/seaicefrac.yaml
  testinput/seaicethick.yaml
  testinput/sea_surface_temp.yaml
  testinput/seviri_crtm.yaml
  testinput/sndrd1-4_crtm.yaml
  testinput/tprof.yaml
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
  atmosphere/hirs4_metop-a_geoval_2018041500_m.nc4
  atmosphere/iasi_metop-a_geoval_2018041500_m.nc4
  atmosphere/mhs_n19_geoval_2018041500_m.nc4
  atmosphere/satwind_geoval_2018041500_m.nc4
  atmosphere/satwind_geoval_2018041500_s.nc4
  atmosphere/sfc_geoval_2018041500_m.nc4
  atmosphere/seviri_m08_geoval_2018041500_m.nc4
  atmosphere/sndrd1_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd2_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd3_g15_geoval_2018041500_m.nc4
  atmosphere/sndrd4_g15_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_s.nc4
  marine/coolskin_fake_geovals_2018041500.nc
  marine/cryosat2-2018-04-15_geovals.nc
  marine/icec-2018-04-15_geovals.nc
  marine/Jason-2-2018-04-15_geovals.nc
  marine/profile_2018-04-15_geovals.nc
  marine/sst_obs-2018-04-15_geovals.nc
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput ${CMAKE_CURRENT_BINARY_DIR}/Data ${ufo_test_data} )

# Add ioda obs test data
list (APPEND ioda_obs_test_data
  atmosphere/aircraft_obs_2018041500_m.nc4
  atmosphere/aircraft_obs_2018041500_s.nc4
  atmosphere/airs_aqua_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_m.nc4
  atmosphere/gmi_gpm_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_s.nc4
  atmosphere/atms_npp_obs_2018041500_m.nc4
  atmosphere/ioda_metop_2_amsua.nc4
  atmosphere/aod_obs_2018041500_m.nc4
  atmosphere/aod_obs_2018041500_s.nc4
  atmosphere/cris-fsr_npp_obs_2018041500_m.nc4
  atmosphere/gnssro_obs_2018041500_s.nc4
  atmosphere/gnssro_obs_2018041500_m.nc4
  atmosphere/hirs4_metop-a_obs_2018041500_m.nc4
  atmosphere/iasi_metop-a_obs_2018041500_m.nc4
  atmosphere/mhs_n19_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_s.nc4
  atmosphere/seviri_m08_obs_2018041500_m.nc4
  atmosphere/sndrd1_g15_obs_2018041500_m.nc4
  atmosphere/sndrd2_g15_obs_2018041500_m.nc4
  atmosphere/sndrd3_g15_obs_2018041500_m.nc4
  atmosphere/sndrd4_g15_obs_2018041500_m.nc4

  atmosphere/sfc_obs_2018041500_m.nc4

  atmosphere/sondes_obs_2018041500_m.nc4
  atmosphere/sondes_obs_2018041500_s.nc4
  marine/coolskin_fake_obs_2018041500.nc
  marine/cryosat2-2018-04-15.nc
  marine/icec-2018-04-15.nc
  marine/Jason-2-2018-04-15.nc
  marine/profile_2018-04-15.nc
  marine/sst_obs-2018-04-15.nc4

  )

CREATE_SYMLINK_FILENAME( ${ioda_SOURCE_DIR}/test/testinput
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${ioda_obs_test_data} )

#####################################################################

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_locations
                  SOURCES mains/TestLocations.cc
                  ARGS    "testinput/locations.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gsi_sfc_model_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gsisfcmodel.yaml"
                  LIBS    ufo )

#ecbuild_add_test( TARGET  test_ufo_surface_ps_qc
#                  SOURCES mains/TestObsFilters.cc
#                  ARGS    "testinput/surface_ps_qc.yaml"
#                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/amsua_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_qc
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/amsua_qc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/amsua_crtm.yaml"
                  LIBS    ufo )

#ecbuild_add_test( TARGET  test_ufo_atms_opr_crtm
#                  SOURCES mains/TestObsOperator.cc
#                  ARGS    "testinput/atms_crtm.yaml"
#                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_atms_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/atms_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gmi_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gmi_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gmi_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gmi_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_seviri_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/seviri_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_seviri_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/seviri_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_cris_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/cris_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_cris_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/cris_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_mhs_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/mhs_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_mhs_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/mhs_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_iasi_opr_crtm
                  MPI     2
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/iasi_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_iasi_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/iasi_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_sndrd1-4_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/sndrd1-4_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_sndrd1-4_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/sndrd1-4_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_airs_opr_crtm
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/airs_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_airs_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/airs_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_hirs4_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/hirs4_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_hirs4_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/hirs4_crtm.yaml"
                  LIBS    ufo )

if( ${RTTOV_FOUND} )
  ecbuild_add_test( TARGET  test_ufo_amsua_opr_rttov
                    SOURCES mains/TestObsOperator.cc
                    ARGS    "testinput/amsua_rttov.yaml"
                    LIBS    ufo )

  ecbuild_add_test( TARGET  test_ufo_amsua_tlad_rttov
                    SOURCES mains/TestObsOperatorTLAD.cc
                    ARGS    "testinput/amsua_rttov.yaml"
                    LIBS    ufo )
endif( ${RTTOV_FOUND} )

ecbuild_add_test( TARGET  test_ufo_radiosonde_opr
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/radiosonde.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_radiosonde_qc
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/radiosonde_qc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_radiosonde_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/radiosonde.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft_opr
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aircraft.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/aircraft.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssroref.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssroref.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndGSI_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssrobndgsi.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndGSI_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssrobndgsi.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_satwind_opr
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/satwind.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_qc_boundscheck
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/qc_boundscheck.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_qc_thinning
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/qc_thinning.yaml"
                  LIBS    ufo )

if( ${ROPP-UFO_FOUND} )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_qc
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  LIBS    ufo )
endif( ${ROPP-UFO_FOUND} )

# Marine UFO tests
ecbuild_add_test( TARGET  test_ufo_seaicefrac_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/seaicefrac.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_seaicefrac_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/seaicefrac.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_seaicethick_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/seaicethick.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_seaicethick_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/seaicethick.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_sst_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/sea_surface_temp.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_sst_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/sea_surface_temp.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_adt_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/adt.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_adt_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/adt.yaml"
                  LIBS    ufo )    

ecbuild_add_test( TARGET  test_ufo_CoolSkin_opr
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/coolskin.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_CoolSkin_tlad
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/coolskin.yaml"
                  LIBS    ufo )


# Marine UFO tests
if( ${GSW_FOUND} )
    ecbuild_add_test( TARGET  test_ufo_insitutemperature_opr
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/tprof.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_insitutemperature_tlad
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/tprof.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_marinevertinterp_opr
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/genericprof.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_marinevertinterp_tlad
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/genericprof.yaml"
                      LIBS    ufo )

endif( ${GSW_FOUND} )

ecbuild_add_test( TARGET  test_ufo_aod_opr_crtm
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aod_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aod_tlad_crtm
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/aod_crtm.yaml"
                  LIBS    ufo )

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
