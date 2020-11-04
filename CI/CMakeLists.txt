
# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# UFO bundle
#

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ufo-bundle LANGUAGES C CXX Fortran )

find_package(ecbuild)
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()
ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/JCSDA/jedi-cmake.git" )

option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "ON") #Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
ecbuild_bundle( PROJECT eckit GIT "https://github.com/JCSDA-internal/eckit.git" )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/JCSDA-internal/fckit.git" )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/JCSDA-internal/atlas.git" )

ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA-internal/saber.git" )
#ecbuild_bundle( PROJECT gsw  GIT "https://github.com/JCSDA-internal/GSW-Fortran.git" BRANCH develop  UPDATE )
ecbuild_bundle( PROJECT crtm GIT "https://github.com/JCSDA-internal/crtm.git" )
ecbuild_bundle( PROJECT ioda GIT "https://github.com/JCSDA-internal/ioda.git" )
ecbuild_bundle( PROJECT rttov  GIT "https://github.com/JCSDA-internal/rttov.git" )
ecbuild_bundle( PROJECT ufo  GIT "https://github.com/JCSDA-internal/ufo.git" )


# Build Doxygen documentation
option(BUILD_UFO_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_UFO_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_UFO_BUNDLE_DOC)

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)

