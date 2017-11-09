#cmake_minimum_required(VERSION 2.8)
if(BUILD_NCDIAG)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${NCDIAG_INCS})
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

  # NetCDF-4 library
# include(./cmake/FindBaselibs.cmake)
  set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
  include_directories( ${NETCDF_INCLUDES} ${NCDIAG_INCS} )

  # 32-bit reals, for now
# add_definitions(-D_REAL4_)
  add_definitions(-D_REAL8_)
  FILE(GLOB NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*90 )

# files not to be part of the library
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/test_nc_unlimdims.F90 ) 
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/nc_diag_cat.F90 ) 
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_bin2txt.f90 ) 
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_conv_bin2nc4.f90 ) 
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_rad_bin2nc4.f90 ) 
  LIST(REMOVE_ITEM NCDIAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ut_ncrad.F90 ) 

  set_source_files_properties( ${NCDIAG_SRC} PROPERTIES COMPILE_FLAGS ${NCDIAG_Fortran_FLAGS} )
  add_library(ncdiag STATIC ${NCDIAG_SRC})
# executables
  add_executable(test_nc_unlimdims.x    ${CMAKE_CURRENT_SOURCE_DIR}/test_nc_unlimdims.F90 )
  add_executable(nc_diag_cat.x          ${CMAKE_CURRENT_SOURCE_DIR}/nc_diag_cat.F90 )
  add_executable(gsidiag_bin2txt.x      ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_bin2txt.f90 )
  add_executable(gsidiag_conv_bin2nc4.x ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_conv_bin2nc4.f90 )
  add_executable(gsidiag_rad_bin2nc4.x  ${CMAKE_CURRENT_SOURCE_DIR}/gsidiag_rad_bin2nc4.f90 )
  add_executable(ut_ncrad.x             ${CMAKE_CURRENT_SOURCE_DIR}/ut_ncrad.F90 )
# links related to exec creation
  target_link_libraries(nc_diag_cat.x          ncdiag ${NETCDF_LIBRARIES} )  
  target_link_libraries(test_nc_unlimdims.x    ncdiag ${NETCDF_LIBRARIES} )  
  target_link_libraries(gsidiag_bin2txt.x      ncdiag ${NETCDF_LIBRARIES} )  
  target_link_libraries(gsidiag_conv_bin2nc4.x ncdiag ${NETCDF_LIBRARIES} )  
  target_link_libraries(gsidiag_rad_bin2nc4.x  ncdiag ${NETCDF_LIBRARIES} )  
  target_link_libraries(ut_ncrad.x             ncdiag ${NETCDF_LIBRARIES} )  
# set_target_properties(test_nc_unlimdims.x PROPERTIES Fortran_MODULE_DIRECTORY ${NCDIAG_INCS} )
# set_target_properties(nc_diag_cat.x PROPERTIES Fortran_MODULE_DIRECTORY ${NCDIAG_INCS} )
 
endif(BUILD_NCDIAG)
