# (C) Copyright 2017-2023 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# UFO bundle
#

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ufo-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()

include( /opt/view/share/jedicmake/Functions/git_functions.cmake )

option(ENABLE_IODA_DATA "Obtain ioda test data from ioda-data repository (vs tarball)" ON)
option(ENABLE_UFO_DATA "Obtain ufo test data from ufo-data repository (vs tarball)" ON)

ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA-internal/saber.git" )
ecbuild_bundle( PROJECT gsw  GIT "https://github.com/JCSDA-internal/GSW-Fortran.git" )
# Don't need ioda-data in CI because we are only testing ufo
#ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" )
ecbuild_bundle( PROJECT ioda GIT "https://github.com/JCSDA-internal/ioda.git" )
ecbuild_bundle( PROJECT crtm  GIT "https://github.com/JCSDA/CRTMv3.git" )
ecbuild_bundle( PROJECT rttov  GIT "https://github.com/JCSDA-internal/rttov.git" )
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/JCSDA-internal/ropp-test.git" )
ecbuild_bundle( PROJECT oasim  GIT "https://github.com/JCSDA-internal/oasim.git" )
ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" )
ecbuild_bundle( PROJECT ufo  GIT "https://github.com/JCSDA-internal/ufo.git" )

# Set GIT_BRANCH_FUNC to UFO's current branch so that it can be used for ufo-data
find_branch_name(REPO_DIR_NAME ufo)

# If UFO's current branch is available in ufo-data repo, that branch will be checked out
branch_checkout (REPO_DIR_NAME ufo-data
                 BRANCH ${GIT_BRANCH_FUNC} )

# Build Doxygen documentation
option(BUILD_UFO_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_UFO_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_UFO_BUNDLE_DOC)

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
