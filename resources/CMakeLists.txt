# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include macros for creating links and symlinks
include( ufo_functions )

# Add list of resource files
list( APPEND ufo_resources
  bmatrix/rttov/atms_bmatrix_70_test.dat
  rmatrix/rttov/atms_noaa_20_rmatrix_test.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bmatrix/rttov)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rmatrix/rttov)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_resources} )
