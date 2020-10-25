
find_package(Doxygen)
if(DOXYGEN_FOUND)
    message(STATUS "Build UFO Documentation")
    include(GNUInstallDirs)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen_ufo ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating UFO documentation with Doxygen"
        VERBATIM )
    set( HTML_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/html )
    install(DIRECTORY ${HTML_BUILD_DIR} DESTINATION ${CMAKE_INSTALL_FULL_DATADIR}/${PROJECT_NAME}/doc COMPONENT Documentation)
else()
    message(WARNING "Build UFO Documentation failure to find Doxygen")
endif()
