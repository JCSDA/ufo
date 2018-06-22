# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Create Data directory for test input config and symlink all files
list( APPEND ufo_test_input
  testinput/ufotest.json
  testinput/ufocrtm.json
  testinput/amsua.json
  testinput/radiosonde.json
  testinput/aircraft.json
  testinput/aod.json
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${ufo_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Add files to cmake resources
ecbuild_add_resources( TARGET ufo_test_scripts
                       SOURCES_PACK
                       ${ufo_test_input}
                     )

# Create Data directory for test data and symlink files
list( APPEND ufo_test_data
  amsua_n19_obs.nc4
  amsua_n19_geovals.nc4
  diag_t_01_wprofiles.nc4
  aircraft_geovals.nc4
  seaice_obs.nc
  diag_viirs_ges.2015081000_dbl_subset.nc4
  vars_ges.2015081000_dbl_subset.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${ufo_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${FILENAME} )
endforeach(FILENAME)

# Add marine test data
if( ${GSW_FOUND} )

    list( APPEND ufo_marine_test_input
      testinput/tprof.json
      testinput/seaice.json
    )

    foreach(FILENAME ${ufo_marine_test_input})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
               ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
               ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
    endforeach(FILENAME)

    ecbuild_add_resources( TARGET ufo_marine_test_scripts
                           SOURCES_PACK
                           ${ufo_marine_test_input}
                         )

    list( APPEND ufo_marine_test_data
      marineobs/t0n156e_dy.nc
      marineobs/t0n156e_dy_geovals.nc
      marineobs/Jason-3-2018-04-15.nc
      marineobs/seaice_obs-2018-04-15.nc
      marineobs/cryosat2-2018-04-15.nc
    )

    foreach(FILENAME ${ufo_marine_test_data})
        get_filename_component(filename ${FILENAME} NAME )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
               ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${FILENAME}
               ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
    endforeach(FILENAME)

endif( ${GSW_FOUND} )

#####################################################################

ecbuild_add_test( TARGET  test_ufo_geovals
                  BOOST
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/ufotest.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_amsua
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/amsua.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_radiosonde
                  BOOST
                  MPI     6
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/radiosonde.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_aircraft
                  BOOST
                  MPI     6
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aircraft.json"
                  LIBS    ufo )

# Marine UFO tests
if( ${GSW_FOUND} )

    ecbuild_add_test( TARGET  test_ufo_insitutemperature
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/tprof.json"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_insitutemperature_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/tprof.json"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_seaice
                      BOOST
                      SOURCES mains/TestObsOperator.cc
                      ARGS    "testinput/seaice.json"
                      LIBS    ufo )

    ecbuild_add_test( TARGET  test_ufo_obsop_seaice_tlad
                      BOOST
                      SOURCES mains/TestObsOperatorTLAD.cc
                      ARGS    "testinput/seaice.json"
                      LIBS    ufo )

endif( ${GSW_FOUND} )

ecbuild_add_test( TARGET  test_ufo_aod
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/aod.json"
                  LIBS    ufo )
ecbuild_add_test( TARGET  test_ufo_obsop_rsonde_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/radiosonde.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_obsop_aircraft_tlad
                  BOOST
                  SOURCES mains/TestObsOperatorTLAD.cc
                  ARGS    "testinput/aircraft.json"
                  LIBS    ufo )

#ecbuild_add_test( TARGET  test_ufo_obscheck
#                  BOOST
#                  SOURCES mains/TestObsCheck.cc
#                  ARGS    "testinput/ufotest.json"
#                  LIBS    ufo )

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

# Create Data directory for test input and symlink all files
foreach(FILENAME ${crtm_test_input})
   execute_process( COMMAND ln -s ${crtm_SOURCE_DIR}/fix/${FILENAME} ${CMAKE_CURRENT_BINARY_DIR}/Data )
endforeach(FILENAME)

