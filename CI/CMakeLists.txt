
# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# UFO bundle
#

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ufo-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

find_package(ecbuild)
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()

#ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/JCSDA-internal/jedi-cmake.git" )
#include( jedicmake/cmake/Functions/git_functions.cmake )
include( /opt/view/share/jedicmake/Functions/git_functions.cmake )

ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA-internal/saber.git" )
ecbuild_bundle( PROJECT gsw  GIT "https://github.com/JCSDA-internal/GSW-Fortran.git" )
ecbuild_bundle( PROJECT ioda GIT "https://github.com/JCSDA-internal/ioda.git" )
ecbuild_bundle( PROJECT crtm  GIT "https://github.com/JCSDA/CRTMv3.git" )
ecbuild_bundle( PROJECT rttov  GIT "https://github.com/JCSDA-internal/rttov.git" )
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/JCSDA-internal/ropp-test.git" )
ecbuild_bundle( PROJECT geos-aero  GIT "https://github.com/JCSDA-internal/geos-aero.git" )
ecbuild_bundle( PROJECT oasim  GIT "https://github.com/JCSDA-internal/oasim.git" )
ecbuild_bundle( PROJECT ufo  GIT "https://github.com/JCSDA-internal/ufo.git" )

# Don't need ioda-data in CI because we are only testing ufo

# If IODA branch is being built set GIT_BRANCH_FUNC to IODA's current branch.
# If a tagged version of IODA is being built set GIT_TAG_FUNC to ioda's current tag. In this case,
# IODA test files will be download from UCAR DASH and ioda-data repo will not be cloned.
#find_branch_name(REPO_DIR_NAME ioda)
# When LOCAL_PATH_JEDI_TESTFILES is set to the directory of IODA test files stored
# in a local directory, ioda-data repo will not be cloned

#if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
#  ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" )
#endif()

# same procedure for ufo-data
find_branch_name(REPO_DIR_NAME ufo)
if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
  ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" )
endif()

# Build Doxygen documentation
option(BUILD_UFO_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_UFO_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_UFO_BUNDLE_DOC)

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
