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
  testinput/ufotest_atmosphere.yaml
  testinput/ufotest_constituents.yaml
  testinput/amsua_crtm.yaml
  testinput/amsua_rttov.yaml
  testinput/radiosonde.yaml
  testinput/radiosonde_qc.yaml
  testinput/aircraft.yaml
  testinput/aircraft_qc.yaml
  testinput/aod_crtm.yaml
  testinput/gnssroref.yaml
  testinput/gnssrobndropp1d.yaml
  testinput/gnssrobndropp2d.yaml
  testinput/channelsparser.yaml
  testinput/satwind.yaml
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
  atmosphere/amsua_geoval_n19_2018041500_m.nc4
  atmosphere/amsua_geoval_n19_2018041500_s.nc4
  atmosphere/sondes_geoval_2018041500_m.nc4
  atmosphere/sondes_geoval_2018041500_s.nc4
  atmosphere/gnssro_geoval_2018041500_s.nc4
  atmosphere/aod_geoval_2018041500_m.nc4
  atmosphere/aod_geoval_2018041500_s.nc4
  atmosphere/satwind_geoval_2018041500_m.nc4
  atmosphere/satwind_geoval_2018041500_s.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput ${CMAKE_CURRENT_BINARY_DIR}/Data ${ufo_test_data} )

# Add marine test data
if( ${GSW_FOUND} )

    list( APPEND ufo_marine_test_input
      testinput/ufotest_marine.yaml
      testinput/tprof.yaml
      testinput/adt.yaml
      testinput/seaice.yaml
      testinput/sea_surface_temp.yaml
    )

	CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_marine_test_input} )

    ecbuild_add_resources( TARGET ufo_marine_test_scripts
                           SOURCES_PACK
                           ${ufo_marine_test_input}
                         )

    list( APPEND ufo_marine_test_data
      marine/Jason-2-2018-04-15_geovals.nc
      marine/profile_2018-04-15_geovals.nc
      marine/cryosat2-2018-04-15_geovals.nc
      marine/icec-2018-04-15_geovals.nc
      marine/sst_obs-2018-04-15_geovals.nc
    )

	CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}/testinput
                             ${CMAKE_CURRENT_BINARY_DIR}/Data
                             ${ufo_marine_test_data} )

endif( ${GSW_FOUND} )

# Add ioda obs test data
list (APPEND ioda_obs_test_data
  atmosphere/aircraft_obs_2018041500_m.nc4
  atmosphere/aircraft_obs_2018041500_s.nc4
  atmosphere/amsua_obs_n19_2018041500_m.nc4
  atmosphere/amsua_obs_n19_2018041500_s.nc4
  atmosphere/sondes_obs_2018041500_m.nc4
  atmosphere/sondes_obs_2018041500_s.nc4
  atmosphere/gnssro_obs_2018041500_s.nc4
  atmosphere/diag_t_obs_01_wprofiles.odb
  atmosphere/diag_t_obs_01_wprofiles_odb

  atmosphere/aod_obs_2018041500_m.nc4
  atmosphere/aod_obs_2018041500_s.nc4

  atmosphere/satwind_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_s.nc4

  marine/Jason-2-2018-04-15.nc
  marine/profile_2018-04-15.nc
  marine/cryosat2-2018-04-15.nc
  marine/icec-2018-04-15.nc
  marine/sst_obs-2018-04-15.nc4
  )

CREATE_SYMLINK_FILENAME( ${ioda_SOURCE_DIR}/test/testinput
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${ioda_obs_test_data} )

#####################################################################

ecbuild_add_test( TARGET  test_ufo_geovals_atmosphere
                  BOOST
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/ufotest_atmosphere.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals_constituents
                  BOOST
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/ufotest_constituents.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_locations_atmosphere
                  BOOST
                  SOURCES mains/TestLocations.cc
                  ARGS    "testinput/ufotest_atmosphere.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_locations_constituents
                  BOOST
                  SOURCES mains/TestLocations.cc
                  ARGS    "testinput/ufotest_constituents.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_opr_crtm
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/amsua_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua_tlad_crtm
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/amsua_crtm.yaml"
                  LIBS    ufo )

if( ${RTTOV_FOUND} )
  ecbuild_add_test( TARGET  test_ufo_amsua_opr_rttov
                    BOOST
                    SOURCES mains/TestObsOperator.cc
                    ARGS    "testinput/amsua_rttov.yaml"
                    LIBS    ufo )

  ecbuild_add_test( TARGET  test_ufo_amsua_tlad_rttov
                    BOOST
                    SOURCES mains/TestObsOperatorTLAD.cc
                    ARGS    "testinput/amsua_rttov.yaml"
                    LIBS    ufo )
endif( ${RTTOV_FOUND} )

ecbuild_add_test( TARGET  test_ufo_radiosonde_opr
                  BOOST
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/radiosonde.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_radiosonde_qc
                  BOOST
                  MPI     4
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/radiosonde_qc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_radiosonde_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/radiosonde.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft_opr
                  BOOST
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aircraft.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft_qc
                  BOOST
                  SOURCES mains/TestObsFilters.cc
                  ARGS    "testinput/aircraft_qc.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/aircraft.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_opr
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssroref.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_gnssroRef_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssroref.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_satwind_opr
                  BOOST
                  MPI     4
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/satwind.yaml"
                  LIBS    ufo )

if( ${ROPP-UFO_FOUND} )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_opr
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  LIBS    ufo )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP1D_tlad
                  BOOST 
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssrobndropp1d.yaml"
                  LIBS    ufo )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_opr
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  LIBS    ufo )
ecbuild_add_test( TARGET  test_ufo_gnssroBndROPP2D_tlad
                  BOOST 
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/gnssrobndropp2d.yaml"
                  LIBS    ufo )
endif( ${ROPP-UFO_FOUND} )

# Marine UFO tests
if( ${GSW_FOUND} )

    ecbuild_add_test( TARGET  test_ufo_geovals_marine
                      BOOST
                      SOURCES mains/TestGeoVaLs.cc
                      ARGS    "testinput/ufotest_marine.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_locations_marine
                      BOOST
                      SOURCES mains/TestLocations.cc
                      ARGS    "testinput/ufotest_marine.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_insitutemperature_opr
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/tprof.yaml"
                      LIBS    ufo )

# SRH: This test works, but it take 5 minutes to run. Thought it would be good
# to comment this test out until the performance issue can be addressed.
#
    ecbuild_add_test( TARGET  test_ufo_insitutemperature_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/tprof.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_seaice_opr
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/seaice.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_seaice_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/seaice.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_sst_opr
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/sea_surface_temp.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_sst_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/sea_surface_temp.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_adt_opr
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/adt.yaml"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_adt_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/adt.yaml"
                      LIBS    ufo )		      

endif( ${GSW_FOUND} )

ecbuild_add_test( TARGET  test_ufo_aod_opr
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aod_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aod_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/aod_crtm.yaml"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_intset_parser
                  BOOST
                  SOURCES mains/TestIntSetParser.cc
                  ARGS    "testinput/channelsparser.yaml"
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
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
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
