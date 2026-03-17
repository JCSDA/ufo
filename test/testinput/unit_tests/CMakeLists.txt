# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Optional RTTOVCPP tests (C++ interface to RTTOV)
option(ENABLE_RTTOVCPP_TESTS "Enable RTTOVCPP operator tests" ON)

# Set version-specific suffix and label for RTTOV tests
if( rttov_FOUND )
  if( rttov_VERSION VERSION_GREATER_EQUAL "14.0" )
        set(ENABLE_RTTOVCPP_TESTS OFF CACHE BOOL "Enable RTTOVCPP operator tests" FORCE)
        set( RTTOV_VERSION_SUFFIX "_v14" )
        set( RTTOV_VERSION_LABEL "rttov14" )
    else()
        set( RTTOV_VERSION_SUFFIX "" )
        set( RTTOV_VERSION_LABEL "rttov12" )
    endif()
endif()

add_subdirectory(errors)
add_subdirectory(filters)
add_subdirectory(fov)
add_subdirectory(operators)
add_subdirectory(predictors)
add_subdirectory(profile)
add_subdirectory(utils)
add_subdirectory(variabletransforms)
