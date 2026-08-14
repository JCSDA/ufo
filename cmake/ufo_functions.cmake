# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        file( CREATE_LINK
            ${src}/${FILENAME}
            ${dst}/${FILENAME} SYMBOLIC )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# macro to create a symlink from src to dst with just filename
function(CREATE_SYMLINK_FILENAME src dst)
    foreach (FILENAME ${ARGN})
        get_filename_component(filename ${FILENAME} NAME )
        file( CREATE_LINK
            ${src}/${FILENAME}
            ${dst}/${filename} SYMBOLIC )
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
#  NAME      - the name of the test. "ufo_" is prepended.
#  WORKING_DIRECTORY - Set working directory.
#                      If unset, default to ${PROJECT_SOURCE_DIR}/test.
#  ECBUILD   - Beyond this point, pass the remaining options to ecbuild_add_test

function(ufo_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set( options )
  set(singlevals NAME WORKING_DIRECTORY)
  set(multivals  ECBUILD )
  cmake_parse_arguments(${prefix}
                        "${options}" "${singlevals}" "${multivals}"
                        ${ARGN})

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

  ecbuild_add_test( TARGET  ufo_${ARG_NAME}
                    WORKING_DIRECTORY ${WORKDIR}
                    ${ECBUILD_EXTRA}
                  )

endfunction()

