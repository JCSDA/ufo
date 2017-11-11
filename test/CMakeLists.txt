list( APPEND ufo_test_input
  testinput/ufotest.json
  amsua_n19_wprofiles.nc4
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${ufo_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET ufo_test_scripts
                       SOURCES_PACK
                       ${ufo_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#####################################################################

ecbuild_add_test( TARGET  test_ufo_variables
                  BOOST
                  SOURCES mains/TestVariables.cc
                  ARGS    "testinput/ufotest.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals
                  BOOST
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/ufotest.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET test_ufo_observationspace
                  BOOST
                  SOURCES mains/TestObservationSpace.cc
                  ARGS "testinput/ufotest.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET test_ufo_obsvector
                  BOOST
                  SOURCES mains/TestObsVector.cc
                  ARGS "testinput/ufotest.json"
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_obsoperator
                  BOOST
                  SOURCES mains/TestObsOperator.cc
                  ARGS    "testinput/ufotest.json"
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
)

# Create Data directory for test input and symlink all files
foreach(FILENAME ${crtm_test_input})
   execute_process( COMMAND ln -s ${CRTM_INCLUDE}/../fix/${FILENAME} ${CMAKE_CURRENT_BINARY_DIR}/Data )
endforeach(FILENAME)

