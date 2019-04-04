# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UFO
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( ufo C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( UFO_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################

# Boost
include_directories( ${Boost_INCLUDE_DIR} )

# NetCDF
set( NETCDF_F90 ON CACHE BOOL "Compile with Fortran NetCDF" )
find_package( NetCDF REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIR} )

# ioda
ecbuild_use_package( PROJECT ioda VERSION 0.1.0 REQUIRED )
include_directories( ${IODA_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# crtm
ecbuild_use_package( PROJECT crtm VERSION 2.2.3 REQUIRED )
include_directories( ${CRTM_INCLUDE_DIRS} )

# RTTOV
ecbuild_use_package( PROJECT rttov VERSION 12.1.0 )
if( ${RTTOV_FOUND} )
    message("RTTOV FOUND; Including RTTOV Radiance Observation Operator")
    include_directories( ${RTTOV_INCLUDE_DIRS} )
else( ${RTTOV_FOUND} )
    message("RTTOV NOT FOUND; Excluding RTTOV Radiance Observation Operator")
endif( ${RTTOV_FOUND} )

# GSW
ecbuild_use_package( PROJECT gsw )
if( ${GSW_FOUND} )
    message("GSW FOUND; Including Marine Observation Operators")
    include_directories( ${GSW_INCLUDE_DIRS} )
else( ${GSW_FOUND} )
    message("GSW NOT FOUND: Excluding Marine Observation Operators")
endif( ${GSW_FOUND} )

# ROPP-UFO
ecbuild_use_package( PROJECT ropp-ufo )
if( ${ROPP-UFO_FOUND} )
    message("ROPP FOUND; Including GNSSRO Bending angle ROPP Operators")
else( ${ROPP-UFO_FOUND} )
    message("ROPP NOT FOUND: Excluding  GNSSRO Bending angle ROPP Operators")
endif( ${ROPP-UFO_FOUND} )

################################################################################
# Export package info
################################################################################

list( APPEND UFO_TPLS oops ioda crtm)

if( ${GSW_FOUND} )
    list ( APPEND UFO_TPLS gsw)
endif( ${GSW_FOUND} )

if( ${ROPP-UFO_FOUND} )
    list ( APPEND UFO_TPLS ropp-ufo )
endif( ${ROPP-UFO_FOUND} )

if( ${RTTOV_FOUND} )
    list ( APPEND UFO_TPLS rttov)
endif( ${RTTOV_FOUND} )

set( UFO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
                      ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( UFO_LIBRARIES ufo )

get_directory_property( UFO_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${UFO_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND UFO_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND UFO_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND UFO_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include( ufo_compiler_flags )
include_directories( ${UFO_INCLUDE_DIRS} ${UFO_EXTRA_INCLUDE_DIRS}
                     ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tools )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME ufo )

# print the summary of the configuration
ecbuild_print_summary()

