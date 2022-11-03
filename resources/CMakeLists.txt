# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include macros for creating links and symlinks
include( ufo_functions )

# Add list of resource files
list( APPEND ufo_resources
  auxillary/IASI_EmisEigenVec.dat
  bmatrix/rttov/amsr_bmatrix_70_test.dat
  bmatrix/rttov/atms_bmatrix_70_test.dat
  bmatrix/rttov/atovs_bmatrix_70_test.dat
  bmatrix/rttov/iasi_bmatrix_70_test.dat
  bmatrix/rttov/ssmis_bmatrix_70_test.dat
  rmatrix/rttov/atms_noaa_20_rmatrix_test.nc4
  rmatrix/rttov/amsr_gcomw1_rmatrix_test.nc4
  rmatrix/rttov/iasi_metopb_rmatrix_test.nc4
  rmatrix/rttov/ssmis_dmspf17_rmatrix_test.nc4
  rmatrix/rttov/atovs_metopb_rmatrix_test.nc4
  rmatrix/gnssro/gnssro_ba_rmatrix_latitude.nl
  rmatrix/gnssro/gnssro_ba_rmatrix_latitude.nc4
  rmatrix/gnssro/gnssro_ba_rmatrix_avtemp.nl
  bmatrix/gnssro/gnssro_bmatrix.txt
  namemap/test_name_map.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/auxillary)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bmatrix/rttov)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rmatrix/rttov)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bmatrix/gnssro)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rmatrix/gnssro)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/namemap)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_resources} )
