list( APPEND ufo_test_input
  testinput/ufotest.json
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
                  SOURCES executables/TestVariables.cc
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

