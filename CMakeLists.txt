# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UFO
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( ufo VERSION 1.10.0 LANGUAGES C CXX Fortran )

## Ecbuild integration
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_UFO_DOC "Build UFO documentation" OFF )
option( ROPP_ARO "Using ROPP for airborne radio occultation" OFF )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

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
find_package( eckit 1.24.4 REQUIRED )

# fckit
find_package( fckit 0.11.0 REQUIRED )

# ioda
find_package( ioda 2.9.0 REQUIRED )

# oops
find_package( oops 1.10.0 REQUIRED )

# crtm
find_package( crtm 2.4.1 QUIET )
if( NOT ${crtm_FOUND} )
   find_package( crtm 2.4.0 QUIET )
endif( NOT ${crtm_FOUND} )
if( ${crtm_FOUND} )
    message(STATUS "CRTM FOUND; Including CRTM Observation Operator (${crtm_VERSION})")
    if( crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCRTM_VERSION=3")
    else ()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCRTM_VERSION=2")
    endif()
else( ${crtm_FOUND} )
    message(STATUS "CRTM NOT FOUND; Excluding CRTM Observation Operator")
endif( ${crtm_FOUND} )

# RTTOV
find_package( rttov 12.1.0 QUIET )
if( ${rttov_FOUND} )
    message(STATUS "RTTOV FOUND; Including RTTOV Radiance Observation Operator (${rttov_VERSION})")
else( ${rttov_FOUND} )
    message(STATUS "RTTOV NOT FOUND; Excluding RTTOV Radiance Observation Operator")
endif( ${rttov_FOUND} )

# GSW
find_package( gsw QUIET )
if( ${gsw_FOUND} )
    message(STATUS "GSW FOUND; Including Marine Observation Operators (${gsw_VERSION})")
else( ${gsw_FOUND} )
    message(STATUS "GSW NOT FOUND: Excluding Marine Observation Operators")
endif( ${gsw_FOUND} )

# ROPP-UFO
find_package( ropp-ufo QUIET )
if( ${ropp-ufo_FOUND} )
    message(STATUS "ROPP FOUND; Including GNSSRO Bending angle ROPP Operators (${ropp-ufo_VERSION})")
    if ( ROPP_ARO)
        message(STATUS "Will use ROPP ARO Operator")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dropp_aro")
    endif()
else( ${ropp-ufo_FOUND} )
    message(STATUS "ROPP NOT FOUND: Excluding GNSSRO Bending angle ROPP Operators")
endif( ${ropp-ufo_FOUND} )

# GEOS-AERO
find_package( geos-aero QUIET )
if( ${geos-aero_FOUND} )
    message(STATUS "GEOS-AERO FOUND; Including AOD operator from GEOS (${geos-aero_VERSION})")
else( ${geos-aero_FOUND} )
    message(STATUS "GEOS-AERO NOT FOUND: Excluding AOD operator from GEOS")
endif( ${geos-aero_FOUND} )

# OASIM
find_package( oasim QUIET )
if( ${oasim_FOUND} )
    message(STATUS "OASIM FOUND; Including Ocean Atmosphere Spectral Irradiance Model (${oasim_VERSION})")
else( ${oasim_FOUND} )
    message(STATUS "OASIM NOT FOUND: Excluding Ocean Atmosphere Spectral Irradiance Model")
endif( ${oasim_FOUND} )

################################################################################
# Test data
################################################################################

# Set values for downloading crtm test data from tarball

if( ${crtm_FOUND} )
  if( crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
	  #https://bin.ssec.wisc.edu/pub/s4/CRTM/fix_REL-3.1.0.2.tgz
  	  set (crtm_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM")
    set (crtm_data_tag "fix_REL-3.1.0.2")  #this should work with all v3.0 and v3.1 
	  set (crtm_data_tag_postfix "fix/") #used in test/CMakeLists.txt and permits compatibility with older methods
    set (crtm_data_tarball_hash "730883a9b1f80a5fb6bf93d9b244321f")
    set (CRTM_CMAKE_TESTINPUT "${CMAKE_CURRENT_SOURCE_DIR}/test/testinput" )
	  set (crtm_data_tarball_name ${crtm_data_tag}.tgz)  #tarball name structure changed with the 3.1.0 release
  elseif( crtm_VERSION VERSION_EQUAL 2.4.1 )
    set (crtm_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM/file")
    set (crtm_data_tag "2.4.1_skylab_4.0" )
	  set (crtm_data_tag_postfix "") #needed for backward compatibility
    set (crtm_data_tarball_hash "5ef915c445f94c6ef11f53221521a79e")
    set (CRTM_CMAKE_TESTINPUT "${CMAKE_CURRENT_SOURCE_DIR}/test/testinput/CRTMv2" )
	  set (crtm_data_tarball_name crtm_coefficients_${crtm_data_tag}.tar.gz)
  elseif( crtm_VERSION VERSION_EQUAL 2.4.0 )
    set (crtm_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM/file")
    set (crtm_data_tag "2.4.0_skylab_3.0" )
	  set (crtm_data_tag_postfix "") #needed for backward compatibility
    set (crtm_data_tarball_hash "73a89f04c2ed52b9b904eac9dab52b60")
    set (CRTM_CMAKE_TESTINPUT "${CMAKE_CURRENT_SOURCE_DIR}/test/testinput/CRTMv2" )
	  set (crtm_data_tarball_name crtm_coefficients_${crtm_data_tag}.tar.gz)
  else ()
	  #Something is wrong with the CRTM_Version
	   message(FATAL_ERROR "Error: Requested CRTM Version (${crtm_VERSION}) is incompatible with UFO.  Contact the JCSDA OBS team.")
  endif()
endif( ${crtm_FOUND} )


# DH 2023/10/06 - this is unfortunately still true ...
### The ufo_data package is problematic in bundles because it has some special CMake logic
### at the root bundle location. So, we cannot find the ufo data while in ufo!
### Basically, with the current logic we hack ufo-data_FOUND to match ENABLE_UFO_DATA
#find_package( ufo-data QUIET )
if (ENABLE_UFO_DATA)
  set (ufo-data_FOUND TRUE)
else()
  # Set values for downloading ufo test data from tarball
  set (ufo_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM/file")
  set (ufo_data_tag "1.10.0")
  set (ufo_data_tarball_name "ufo_testinput_tier_1_${ufo_data_tag}.tar.gz")
  set (ufo_data_tarball_hash "9dde2bc6a84da5d5e0e18e6d0aa62148")
endif()

################################################################################
# Sources
################################################################################

add_subdirectory( src )
add_subdirectory( test )
# Make UFO_CRTM_TESTFILES_PATH and ufo_data_testinput_tier_1_local
# (rename to UFO_UFO_TESTFILES_PATH) set in test available at parent
# scope so that downstream repos have access to it
set (UFO_CRTM_TESTFILES_PATH ${UFO_CRTM_TESTFILES_PATH} PARENT_SCOPE)
set (UFO_UFO_TESTFILES_PATH ${UFO_UFO_TESTFILES_PATH} PARENT_SCOPE)

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
