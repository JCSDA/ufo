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

# macro to prepend a prefix with relative path
# can this be added to ecbuild for use elsewhere?
function(PREPEND var prefix )
    set ( listVar "" )
    foreach (f ${ARGN})
        list (APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set ( ${var} "${listVar}" PARENT_SCOPE )
endfunction(PREPEND)

# -----------------------------------------------------------------------------

# The following is a wrapper to simplify the generation of tests.
# There are two types of tests:
#  1) ufo test executables (EXE must be given)
#  2) interface tests  (SRC must be given)
#
# Arguments:
#  NAME      - the name of the test. "ufo_test_tier#_" is prepended.
#  TIER      - The testing tier. Defaults to 1.
#  NOTRAPFPE - Disable FPE trapping.
#  WORKING_DIRECTORY - Set working directory.
#                      If unset, default to ${PROJECT_SOURCE_DIR}/test.
#  ECBUILD   - Beyond this point, pass the remaining options to ecbuild_add_test

function(ufo_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOTRAPFPE)
  set(singlevals NAME TIER WORKING_DIRECTORY)
  set(multivals  ECBUILD )
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # determine if floating point error trapping should be set
  if ( ARG_NOTRAPFPE )
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=0")
  else()
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=1")
  endif()

  # determine working directory
  if ( ARG_WORKING_DIRECTORY )
    set ( WORKDIR ${ARG_WORKING_DIRECTORY} )
  else()
    set ( WORKDIR ${PROJECT_SOURCE_DIR}/test )
  endif()

  # determine working directory
  if ( ARG_ECBUILD )
    set ( ECBUILD_EXTRA ${ARG_ECBUILD} )
  else()
    set ( ECBUILD_EXTRA "" )
  endif()

  # Set the tier for the test
  if (NOT ARG_TIER)
    set( TESTTIER 1 )
  else()
    set( TESTTIER ${ARG_TIER} )
  endif()

  # Add the test only if the test tier is sufficiently high.
  if ( TESTTIER LESS_EQUAL UFO_TEST_TIER )

      ecbuild_add_test( TARGET  ufo_test_tier${TESTTIER}_${ARG_NAME}
                        WORKING_DIRECTORY ${WORKDIR}
                        ENVIRONMENT ${TRAPFPE_ENV}
                        ${ECBUILD_EXTRA}
                        )

  endif()

endfunction()

