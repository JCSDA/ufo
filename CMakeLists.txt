# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UFO
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( ufo VERSION 1.0.0 LANGUAGES C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "ufo" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

option( ENABLE_UFO_DOC "Build UFO documentation" OFF )

ecbuild_enable_fortran( REQUIRED )

set( UFO_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################

# jedi-cmake
find_package( jedicmake QUIET )

# MPI
find_package( MPI REQUIRED COMPONENTS C CXX Fortran )

# Boost
include_directories( ${Boost_INCLUDE_DIR} )

# NetCDF
set( NETCDF_F90 ON CACHE BOOL "Compile with Fortran NetCDF" )
find_package( NetCDF REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# Eigen3
find_package( Eigen3 REQUIRED NO_MODULE HINTS
              $ENV{Eigen3_ROOT} $ENV{EIGEN3_ROOT} $ENV{Eigen_ROOT} $ENV{EIGEN_ROOT}
              $ENV{Eigen3_PATH} $ENV{EIGEN3_PATH} $ENV{Eigen_PATH} $ENV{EIGEN_PATH} )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

# gsl-lite
find_package(gsl-lite REQUIRED HINTS $ENV{gsl_lite_DIR})

# eckit
find_package( eckit 1.11.6 REQUIRED )
include_directories( ${ECKIT_INCLUDE_DIRS} )

# fckit
find_package( fckit 0.7.0 REQUIRED )
include_directories( ${FCKIT_INCLUDE_DIRS} )

# ioda
find_package( ioda 1.0.0 REQUIRED )
include_directories( ${IODA_INCLUDE_DIRS} )

# oops
find_package( oops 1.0.0 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# crtm
find_package( crtm 2.3 QUIET )
if( ${crtm_FOUND} )
  message(STATUS "CRTM FOUND; Including CRTM Observation Operator")
  include_directories( ${CRTM_INCLUDE_DIRS} )
else( ${crtm_FOUND} )
  message(STATUS "CRTM NOT FOUND; Excluding CRTM Observation Operator")
endif( ${crtm_FOUND} )

# RTTOV
find_package( rttov 12.1.0 QUIET )
if( ${rttov_FOUND} )
    message(STATUS "RTTOV FOUND; Including RTTOV Radiance Observation Operator")
    include_directories( ${RTTOV_INCLUDE_DIRS} )
else( ${rttov_FOUND} )
    message(STATUS "RTTOV NOT FOUND; Excluding RTTOV Radiance Observation Operator")
endif( ${rttov_FOUND} )

# GSW
find_package( gsw QUIET )
if( ${gsw_FOUND} )
    message(STATUS "GSW FOUND; Including Marine Observation Operators")
    include_directories( ${GSW_INCLUDE_DIRS} )
else( ${gsw_FOUND} )
    message(STATUS "GSW NOT FOUND: Excluding Marine Observation Operators")
endif( ${gsw_FOUND} )

# ROPP-UFO
find_package( ropp-ufo QUIET )
if( ${ropp-ufo_FOUND} )
    message(STATUS "ROPP FOUND; Including GNSSRO Bending angle ROPP Operators")
else( ${ropp-ufo_FOUND} )
    message(STATUS "ROPP NOT FOUND: Excluding  GNSSRO Bending angle ROPP Operators")
endif( ${ropp-ufo_FOUND} )

# GEOS-AERO
find_package( geos-aero QUIET )
if( ${geos-aero_FOUND} )
    message(STATUS "GEOS-AERO FOUND; Including AOD operator from GEOS")
    include_directories( ${GEOSAERO_INCLUDE_DIRS} )
else( ${geos-aero_FOUND} )
    message(STATUS "GEOS-AERO NOT FOUND: Excluding  AOD operator from GEOS")
endif( ${geos-aero_FOUND} )

################################################################################
# Export package info
################################################################################

list( APPEND UFO_TPLS oops ioda crtm )

if( ${gsw_FOUND} )
    list ( APPEND UFO_TPLS gsw)
endif( ${gsw_FOUND} )

if( ${ropp-ufo_FOUND} )
    list ( APPEND UFO_TPLS ropp-ufo )
endif( ${ropp-ufo_FOUND} )

if( ${geos-aero_FOUND} )
    list ( APPEND UFO_TPLS geos-aero )
endif( ${geos-aero_FOUND} )

if( ${rttov_FOUND} )
    list ( APPEND UFO_TPLS rttov)
endif( ${rttov_FOUND} )

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
include_directories( ${UFO_INCLUDE_DIRS} ${UFO_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tools )
add_subdirectory( resources )

if( ENABLE_UFO_DOC )
    add_subdirectory( docs )
endif()

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

ecbuild_add_test( TARGET ufo_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/test
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
