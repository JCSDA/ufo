# Create Data directory for test input data
list( APPEND test_files
  cpplint.py
)

foreach(FILENAME ${test_files})
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
     ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
     ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_test( TARGET ufo_tools_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
