# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UFO
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( ufo VERSION 1.1.0 LANGUAGES C CXX Fortran )

## Ecbuild integration
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_UFO_DOC "Build UFO documentation" OFF )

include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )


################################################################################
# Dependencies
################################################################################

# jedi-cmake
find_package( jedicmake QUIET )

# MPI
find_package( MPI REQUIRED COMPONENTS C CXX Fortran )

# Boost
find_package( Boost REQUIRED )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS C Fortran)

# Eigen3
find_package( Eigen3 REQUIRED NO_MODULE HINTS
              $ENV{Eigen3_ROOT} $ENV{EIGEN3_ROOT} $ENV{Eigen_ROOT} $ENV{EIGEN_ROOT}
              $ENV{Eigen3_PATH} $ENV{EIGEN3_PATH} $ENV{Eigen_PATH} $ENV{EIGEN_PATH} )

# gsl-lite
find_package(gsl-lite REQUIRED HINTS $ENV{gsl_lite_DIR})

# eckit
find_package( eckit 1.11.6 REQUIRED )

# fckit
find_package( fckit 0.7.0 REQUIRED )

# ioda
find_package( ioda 1.0.0 REQUIRED )

# oops
find_package( oops 1.0.0 REQUIRED )

# crtm
find_package( crtm 2.3 QUIET )
if( ${crtm_FOUND} )
  message(STATUS "CRTM FOUND; Including CRTM Observation Operator")
else( ${crtm_FOUND} )
  message(STATUS "CRTM NOT FOUND; Excluding CRTM Observation Operator")
endif( ${crtm_FOUND} )

# RTTOV
find_package( rttov 12.1.0 QUIET )
if( ${rttov_FOUND} )
    message(STATUS "RTTOV FOUND; Including RTTOV Radiance Observation Operator")
else( ${rttov_FOUND} )
    message(STATUS "RTTOV NOT FOUND; Excluding RTTOV Radiance Observation Operator")
endif( ${rttov_FOUND} )

# GSW
find_package( gsw QUIET )
if( ${gsw_FOUND} )
    message(STATUS "GSW FOUND; Including Marine Observation Operators")
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
else( ${geos-aero_FOUND} )
    message(STATUS "GEOS-AERO NOT FOUND: Excluding  AOD operator from GEOS")
endif( ${geos-aero_FOUND} )

################################################################################
# Sources
################################################################################

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tools )
add_subdirectory( resources )

if( ENABLE_UFO_DOC )
    add_subdirectory( docs )
endif()

## Global Tests
ecbuild_add_test( TARGET ufo_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/test
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

## Export package
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
