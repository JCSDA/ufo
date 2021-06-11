
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
ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/JCSDA-internal/jedi-cmake.git" )
include( jedicmake/cmake/Functions/git_functions.cmake )

option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "ON") #Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
option("BUNDLE_SKIP_FCKIT" "Don't build fckit" "ON") #Skip fckit build unless user passes -DBUNDLE_SKIP_FCKIT=OFF
option("BUNDLE_SKIP_ATLAS" "Don't build atlas" "ON") #Skip atlas build unless user passes -DBUNDLE_SKIP_ATLAS=OFF

ecbuild_bundle( PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.16.0 )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.2 )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.24.1 )

ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA-internal/saber.git" )
#ecbuild_bundle( PROJECT gsw  GIT "https://github.com/JCSDA-internal/GSW-Fortran.git" )
ecbuild_bundle( PROJECT crtm GIT "https://github.com/JCSDA-internal/crtm.git" )
ecbuild_bundle( PROJECT ioda GIT "https://github.com/JCSDA-internal/ioda.git" )
ecbuild_bundle( PROJECT rttov  GIT "https://github.com/JCSDA-internal/rttov.git" )
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/JCSDA-internal/ropp-test.git" )
ecbuild_bundle( PROJECT geos-aero  GIT "https://github.com/JCSDA-internal/geos-aero.git" )
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
