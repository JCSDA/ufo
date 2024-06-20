# (C) Copyright 2017-2023 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include macros for creating links and symlinks
include( ufo_functions )

# CRTM test data - from tarball if not provided locally
# It's unclear if anyone is using the local data functionality, therefore comment it out
# for now and add in if someone complains / remove if completely if not.
#if (crtm_FOUND AND DEFINED ENV{UFO_CRTM_TESTFILES})
#  message(STATUS "Use crtm test data from local directory $ENV{UFO_CRTM_TESTFILES}")
#  # A bit of guesswork here, I don't know if folks using this option stored it in the same directory structure
#  set (UFO_CRTM_TESTFILES_PATH "$ENV{UFO_CRTM_TESTFILES}")
#elseif (crtm_FOUND)
if (crtm_FOUND)
  set (UFO_CRTM_TESTFILES_PATH "${CMAKE_SOURCE_DIR}/test-data-release/crtm/${crtm_data_tag}/${crtm_data_tag_postfix}")
  if (EXISTS ${UFO_CRTM_TESTFILES_PATH})
    message(STATUS "Skipping download of test data from ${crtm_data_download_url}/${crtm_data_tarball_name} since it already exists")
  else()
    message(STATUS "Download test data from ${crtm_data_download_url}/${crtm_data_tarball_name}")
    include(FetchContent)
    FetchContent_Declare(
      crtm-testdata
      DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/test-data-release
      DOWNLOAD_NO_EXTRACT ON
      URL      ${crtm_data_download_url}/${crtm_data_tarball_name}
      URL_HASH MD5=${crtm_data_tarball_hash}
    )
    FetchContent_MakeAvailable(crtm-testdata)
    if( crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
       # Tarballs path for versions 3.0.0 and higher
       file(ARCHIVE_EXTRACT INPUT ${CMAKE_SOURCE_DIR}/test-data-release/${crtm_data_tarball_name} DESTINATION ${CMAKE_SOURCE_DIR}/test-data-release/crtm/)
    else ()
       # Tarballs for versions 2.4.0 and 2.4.1 already contain /crtm in the path
       file(ARCHIVE_EXTRACT INPUT ${CMAKE_SOURCE_DIR}/test-data-release/${crtm_data_tarball_name} DESTINATION ${CMAKE_SOURCE_DIR}/test-data-release/)
    endif ()
  endif()
endif()
# Make UFO_CRTM_TESTFILES_PATH available at parent scope so that downstream repos have access to it
set (UFO_CRTM_TESTFILES_PATH ${UFO_CRTM_TESTFILES_PATH} PARENT_SCOPE)

# UFO test data - from ufo-data repo if found, from local directory, or from tarball (release/develop-in-future?)
if (ufo-data_FOUND)
  message(STATUS "Use test data from ufo-data repository")
  set (ufo_data_testinput_tier_1_local "${CMAKE_SOURCE_DIR}/ufo-data/testinput_tier_1")
# It's unclear if anyone is using the local data functionality, therefore comment it out
# for now and add in if someone complains / remove if completely if not.
#elseif (DEFINED ENV{UFO_TESTFILES})
#  message(STATUS "Use test data from local directory $ENV{UFO_TESTFILES}")
#  # A bit of guesswork here, I don't know if folks using this option stored it in the same directory structure
#  set (ufo_data_testinput_tier_1_local "$ENV{UFO_TESTFILES}/ufo-data/testinput_tier_1")
else()
  set (ufo_data_testinput_tier_1_local "${CMAKE_SOURCE_DIR}/test-data-release/ufo/${ufo_data_tag}/testinput_tier_1")
  if (EXISTS ${ufo_data_testinput_tier_1_local})
    message(STATUS "Skipping download of test data from ${ufo_data_download_url}/${ufo_data_tarball_name} since it already exists")
  else()
    message(STATUS "Download test data from ${ufo_data_download_url}/${ufo_data_tarball_name}")
    include(FetchContent)
    FetchContent_Declare(
      ufo-testdata
      DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/test-data-release
      DOWNLOAD_NO_EXTRACT ON
      URL      ${ufo_data_download_url}/${ufo_data_tarball_name}
      URL_HASH MD5=${ufo_data_tarball_hash}
    )
    FetchContent_MakeAvailable(ufo-testdata)
    file(ARCHIVE_EXTRACT INPUT ${CMAKE_SOURCE_DIR}/test-data-release/${ufo_data_tarball_name} DESTINATION ${CMAKE_SOURCE_DIR}/test-data-release)
  endif()
endif()
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Data/ufo/testinput_tier_1")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/ufo)
  file(CREATE_LINK ${ufo_data_testinput_tier_1_local} ${CMAKE_CURRENT_BINARY_DIR}/Data/ufo/testinput_tier_1 SYMBOLIC)
endif()
# Make UFO_UFO_TESTFILES_PATH available at parent scope so that downstream repos have access to it
set (UFO_UFO_TESTFILES_PATH ${ufo_data_testinput_tier_1_local} PARENT_SCOPE)

# Create Data directory for test input config and symlink all files
list( APPEND ufo_test_input
  testinput/empty.yaml
  testinput/geovals_spec.yaml
  testinput/geovals.yaml
  testinput/locations.yaml
  testinput/operatorutils.yaml
  testinput/parameters_older_eckit.yaml
  testinput/parameters.yaml
  testinput/primitive_variables.yaml
  testinput/recordhandler.yaml
  testinput/ufo_geovals_setup_with_mismatched_nvars.yaml
  testinput/ufo_geovals_setup_with_mismatched_nreduced_vars.yaml
  testinput/ufo_geovals_setup_with_sampling_method_set_to_0.yaml
  testinput/ufo_geovals_setup_with_sampling_method_set_to_3.yaml
  testinput/ufo_geovals_setup_with_method_mislabelled_as_trivial.yaml
  testinput/ufo_geovals_partial_setup_with_mismatched_nvars.yaml
  testinput/ufo_geovals_partial_setup_with_sampling_method_set_to_0.yaml
  testinput/ufo_geovals_partial_setup_with_sampling_method_set_to_3.yaml
  testinput/ufo_geovals_setup_sampling_method_0.yaml
  testinput/ufo_geovals_setup_sampling_method_3.yaml
  testinput/ufo_geovals_setup_sampling_method_with_mismatched_nlocs.yaml
  testinput/ufo_geovals_setup_sampling_method_mislabelled_as_trivial.yaml
  testinput/ufo_geovals_setup_trivial_sampling_method_0.yaml
  testinput/ufo_geovals_setup_trivial_sampling_method_3.yaml
  testinput/ufo_geovals_set_invalid_default_format.yaml
  testinput/ufo_geovals_get_nonexistent_var.yaml
  testinput/variables.yaml
)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/unit_tests)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_test_input} )


# Add files to cmake resources
ecbuild_add_resources( TARGET ufo_test_scripts
                       SOURCES_PACK
                       ${ufo_test_input}
                     )


# Create Data directory for test reference files and symlink all files
list( APPEND ufo_test_ref
  testref/unit_tests/filters/print_filter_data.ref
  testref/unit_tests/filters/print_filter_data_MPI2_rank0.ref
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/unit_tests)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/unit_tests/filters)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${ufo_test_ref} )

#####################################################################
# Install headers used in tests of interfaces that may be implemented by clients.

install( FILES ufo/ObsDiagnostics.h ufo/ObsFilters.h ufo/ObsFunction.h ufo/Predictor.h
         DESTINATION ${INSTALL_INCLUDE_DIR}/ufo/test/ufo)

#####################################################################
# Build executables used by multiple tests

ecbuild_add_executable( TARGET  test_ObsBias.x
                        SOURCES mains/TestObsBias.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsDiagnostics.x
                        SOURCES mains/TestObsDiagnostics.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsOperator.x
                        SOURCES mains/TestObsOperator.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsOperatorTLAD.x
                        SOURCES mains/TestObsOperatorTLAD.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsFilters.x
                        SOURCES mains/TestObsFilters.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ObsFunction.x
                        SOURCES mains/TestObsFunction.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_ConventionalProfileProcessing.x
                        SOURCES mains/TestConventionalProfileProcessing.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_Predictor.x
                        SOURCES mains/TestPredictor.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_RecordHandler.x
                        SOURCES mains/TestRecordHandler.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_NearestNeighborInterpolation.x
                        SOURCES mains/TestNearestNeighborInterpolation.cc
                        LIBS    ufo
                       )

ecbuild_add_executable( TARGET  test_FortranGeoVaLs.x
                        SOURCES mains/TestFortranGeoVaLs.cc ufo/geovals_test.F90
                        LIBS    ufo
                       )

####################################################################
# Establish test tiering

set( UFO_TEST_TIER 1 CACHE STRING "UFO test tier. Run tests up to and including this tier.")

# Override UFO_TEST_TIER using environment variable
if( DEFINED ENV{UFO_TEST_TIER} )
    set( UFO_TEST_TIER $ENV{UFO_TEST_TIER} )
endif()

message(STATUS "Running ufo tests for tiers up to and including tier " ${UFO_TEST_TIER})

####################################################################


# Test piecewise linear interpolation
ecbuild_add_test( TARGET  test_ufo_piecewise_linear_interpolation
                  SOURCES mains/TestPiecewiseLinearInterpolation.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_recursivesplitter
                  SOURCES mains/TestRecursiveSplitter.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_dataextractor
                  SOURCES mains/TestDataExtractor.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

# Test integral image
ecbuild_add_test( TARGET  test_ufo_integral_image
                  SOURCES mains/TestIntegralImage.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test geometry calculations
ecbuild_add_test( TARGET  test_ufo_geometry_calculations
                  SOURCES mains/TestGeometryCalculations.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test nearest_neighbor interpolation
ecbuild_add_test( TARGET  test_ufo_nearest_neighbor_interpolation
                  SOURCES mains/TestNearestNeighborInterpolation.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS test_NearestNeighborInterpolation.x)

# Test superob radar template
ecbuild_add_test( TARGET  test_ufo_superob_radar_template
                  SOURCES mains/TestSuperObRadarTemplate.cc
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

# Test Parameters
if(eckit_VERSION VERSION_GREATER_EQUAL 1.16)
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo )
else()
  # Due to changes in eckit 1.16 a different input file is required.
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters_older_eckit.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo )
endif()

# Test operator utils
ecbuild_add_test( TARGET  test_ufo_operator_utils
                  SOURCES mains/TestOperatorUtils.cc
                  ARGS    "testinput/operatorutils.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_variables
                  SOURCES mains/TestVariables.cc
                  ARGS    "testinput/variables.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

# Test the iterator over primitive variables and their values
ecbuild_add_test( TARGET  test_ufo_primitive_variables
                  SOURCES mains/TestPrimitiveVariables.cc
                  ARGS    "testinput/primitive_variables.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals_mpi
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_geovals_spec
                  SOURCES mains/TestGeoVaLsSpec.cc
                  ARGS    "testinput/geovals_spec.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_sampled_locations
                  SOURCES mains/TestSampledLocations.cc
                          ufo/SampledLocations.h
                          ufo/sampled_locations_test.F90
                  ARGS    "testinput/locations.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo )

ecbuild_add_test( TARGET  test_ufo_recordhandler
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_RecordHandler.x
                  ARGS    "testinput/recordhandler.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_RecordHandler.x )

ecbuild_add_test( TARGET  test_ufo_geovals_setup_with_mismatched_nvars
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_with_mismatched_nvars.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_with_mismatched_nvars
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup: mismatch in the number of variables.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_with_mismatched_nreduced_vars
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_with_mismatched_nreduced_vars.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_with_mismatched_nreduced_vars
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup: mismatch in the number of reduced variables.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_with_sampling_method_set_to_0
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_with_sampling_method_set_to_0.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_with_sampling_method_set_to_0
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup: sampling method +0 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_with_sampling_method_set_to_3
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_with_sampling_method_set_to_3.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_with_sampling_method_set_to_3
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup: sampling method +3 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_with_method_mislabelled_as_trivial
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_with_method_mislabelled_as_trivial.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_with_method_mislabelled_as_trivial
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup: trivial sampling methods must sample each location exactly once.*")

ecbuild_add_test( TARGET  test_ufo_geovals_partial_setup_with_mismatched_nvars
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_partial_setup_with_mismatched_nvars.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_partial_setup_with_mismatched_nvars
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_partial_setup: mismatch in the number of variables.*")

ecbuild_add_test( TARGET  test_ufo_geovals_partial_setup_with_sampling_method_set_to_0
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_partial_setup_with_sampling_method_set_to_0.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_partial_setup_with_sampling_method_set_to_0
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_partial_setup: sampling method +0 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_partial_setup_with_sampling_method_set_to_3
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_partial_setup_with_sampling_method_set_to_3.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_partial_setup_with_sampling_method_set_to_3
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_partial_setup: sampling method +3 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_sampling_method_0
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_sampling_method_0.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_sampling_method_0
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_sampling_method: sampling method +0 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_sampling_method_3
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_sampling_method_3.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_sampling_method_3
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_sampling_method: sampling method +3 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_sampling_method_with_mismatched_nlocs.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_sampling_method: mismatch in nlocs.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_sampling_method_mislabelled_as_trivial
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_sampling_method_mislabelled_as_trivial.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_sampling_method_mislabelled_as_trivial
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_sampling_method: the method should be trivial but is not.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_trivial_sampling_method_0
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_trivial_sampling_method_0.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_trivial_sampling_method_0
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_trivial_sampling_method: sampling method +0 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_setup_trivial_sampling_method_3
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_setup_trivial_sampling_method_3.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_setup_trivial_sampling_method_3
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_setup_trivial_sampling_method: sampling method +3 +does not exist.*")

ecbuild_add_test( TARGET  test_ufo_geovals_set_invalid_default_format
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_set_invalid_default_format.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_set_invalid_default_format
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_set_default_format: format must be either ufo_geoval_sampled or ufo_geoval_reduced.*")

ecbuild_add_test( TARGET  test_ufo_geovals_get_nonexistent_var
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_FortranGeoVaLs.x
                  ARGS    "testinput/ufo_geovals_get_nonexistent_var.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_FortranGeoVaLs.x )
set_tests_properties(test_ufo_geovals_get_nonexistent_var
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*ufo_geovals_get_var d doesnt exist.*")

#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
AerosolCoeff/netCDF/AerosolCoeff.nc4
AerosolCoeff/Little_Endian/AerosolCoeff.NAAPS.bin
AerosolCoeff/netCDF/AerosolCoeff.NAAPS.nc4
AerosolCoeff/Little_Endian/AerosolCoeff.GOCART-GEOS5.bin
AerosolCoeff/netCDF/AerosolCoeff.GOCART-GEOS5.nc4
AerosolCoeff/Little_Endian/AerosolCoeff.CMAQ.bin
AerosolCoeff/netCDF/AerosolCoeff.CMAQ.nc4
CloudCoeff/Little_Endian/CloudCoeff.bin
CloudCoeff/netCDF/CloudCoeff.nc4
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2022.nc4
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/NPOESS.IRland.EmisCoeff.bin
EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin
EmisCoeff/VIS_Ice/SEcategory/Little_Endian/NPOESS.VISice.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/NPOESS.VISland.EmisCoeff.bin
EmisCoeff/VIS_Snow/SEcategory/Little_Endian/NPOESS.VISsnow.EmisCoeff.bin
EmisCoeff/VIS_Water/SEcategory/Little_Endian/NPOESS.VISwater.EmisCoeff.bin
EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/USGS.IRland.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/USGS.VISland.EmisCoeff.bin
SpcCoeff/Little_Endian/hirs4_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/hirs4_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsr2_gcom-w1.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsr2_gcom-w1.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-c.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n15.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n18.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
SpcCoeff/Little_Endian/atms_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_npp.TauCoeff.bin
SpcCoeff/Little_Endian/atms_n20.SpcCoeff.bin
SpcCoeff/Little_Endian/atms_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_n20.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_n21.TauCoeff.bin
SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m11.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/seviri_m11.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-c.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-c.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n18.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n19.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD1_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD1_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD2_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD2_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD3_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD3_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD4_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD4_g15.TauCoeff.bin
SpcCoeff/Little_Endian/ssmis_f17.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f17.TauCoeff.bin
SpcCoeff/Little_Endian/ssmis_f18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f18.TauCoeff.bin
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/viirs-i_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/viirs-i_n21.TauCoeff.bin
SpcCoeff/Little_Endian/viirs-m_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/viirs-m_n21.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_j1.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_j1.TauCoeff.bin
SpcCoeff/Little_Endian/v.modis_aqua.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.modis_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/v.modis_terra.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.modis_terra.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g16.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g16.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g17.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g17.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g18.TauCoeff.bin
SpcCoeff/Little_Endian/abi_gr.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_gr.TauCoeff.bin
SpcCoeff/Little_Endian/v.abi_g16.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.abi_g16.TauCoeff.bin
SpcCoeff/Little_Endian/v.abi_g17.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.abi_g17.TauCoeff.bin
SpcCoeff/Little_Endian/v.abi_g18.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.abi_g18.TauCoeff.bin
SpcCoeff/Little_Endian/ahi_himawari8.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ahi_himawari8.TauCoeff.bin
SpcCoeff/Little_Endian/radiometer_smap.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/radiometer_smap.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_n18.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n20.SpcCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n20.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n21.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv1_srf_v3.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv1_srf_v3.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv2.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv2.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv3.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv3.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv4.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv4.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv5.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv5.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv6.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv6.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/tropics_sv7.SpcCoeff.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tropics_sv7.SpcCoeff.TauCoeff.bin
SpcCoeff/Little_Endian/cowvr_ors6.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cowvr_ors6.TauCoeff.bin
SpcCoeff/Little_Endian/tempest-D_cubesat.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/tempest-D_cubesat.TauCoeff.bin
SpcCoeff/Little_Endian/cpr_cloudsat.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cpr_cloudsat.TauCoeff.bin
)

if( crtm_FOUND )
    # symlink to selected CRTM coefficient files.
    create_symlink_filename( ${UFO_CRTM_TESTFILES_PATH}
                             ${CMAKE_CURRENT_BINARY_DIR}/Data
                             ${crtm_test_input} )
endif()

#####################################################################
# Files for RTTOV tests
#####################################################################
if( ${rttov_FOUND} )
  list( APPEND rttov_test_input
    rttov7pred54L/rtcoef_dmsp_17_ssmis.dat
    rttov7pred54L/rtcoef_gcom-w_1_amsr2.dat
    rttov7pred54L/rtcoef_gpm_1_gmi.dat
    rttov7pred54L/rtcoef_metop_1_amsuamhs.dat
    rttov7pred54L/rtcoef_metop_2_amsuamhs.dat
    rttov7pred54L/rtcoef_metop_3_amsuamhs.dat
    rttov7pred54L/rtcoef_noaa_15_amsuamhs.dat
    rttov7pred54L/rtcoef_noaa_18_amsuamhs.dat
    rttov7pred54L/rtcoef_noaa_19_amsuamhs.dat
    rttov7pred54L/rtcoef_metop_1_iasi.dat
    rttov7pred54L/rtcoef_noaa_19_amsua.dat
    rttov7pred54L/rtcoef_jpss_0_atms.dat
    rttov7pred54L/rtcoef_noaa_20_atms.dat
    )

  if(DEFINED ENV{RTTOV_LOCAL_COEF_DIR} AND EXISTS $ENV{RTTOV_LOCAL_COEF_DIR})
    set(RTTOV_COEF_DIR "$ENV{RTTOV_LOCAL_COEF_DIR}")
  elseif(DEFINED rttov_SOURCE_DIR AND EXISTS "${rttov_SOURCE_DIR}")
    set(RTTOV_COEF_DIR "${rttov_SOURCE_DIR}/rtcoef_rttov12")
  elseif(DEFINED rttov_BASE_DIR AND EXISTS "${rttov_BASE_DIR}/etc")
    set(RTTOV_COEF_DIR "${rttov_BASE_DIR}/etc/rtcoef_rttov12")
  endif()
  message("RTTOV coefficient directory: ${RTTOV_COEF_DIR}")

  # Extract mietable to binDir/Data
  if(EXISTS "${RTTOV_COEF_DIR}/mietable/mietable_metop_amsuamhs.tar.zst")
    set(archive "${RTTOV_COEF_DIR}/mietable/mietable_metop_amsuamhs.tar.zst")
    set(out_dir "${CMAKE_CURRENT_BINARY_DIR}/Data/")
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
      file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${out_dir} VERBOSE)
    else()
      file(MAKE_DIRECTORY ${out_dir})
      execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xf ${archive}
        WORKING_DIRECTORY {out_dir}
        RESULT_VARIABLE ret
        )
      if(NOT ret EQUAL 0)
        message(FATAL_ERROR "extract ${archive} => ${out_dir}    ${ret}")
      endif()
    endif()
  elseif(EXISTS "${RTTOV_COEF_DIR}/mietable/mietable_metop_amsuamhs.bin")
    list(APPEND rttov_test_input "mietable/mietable_metop_amsuamhs.bin")
  endif()
    
  # Symlink RTcoef files
  create_symlink_filename( ${RTTOV_COEF_DIR}
                           ${CMAKE_CURRENT_BINARY_DIR}/Data
                           ${rttov_test_input} )
                         
endif()

######
#Files for geos-aero tests
######

if( ${geos-aero_FOUND} )
list( APPEND geos-aero_test_data
Data/optics_BC.v1_3_.nc
Data/optics_BRC.v1_5_.nc
Data/optics_DU.v15_3_.nc
Data/optics_NI.v2_5_.nc
Data/optics_OC.v1_3_.nc
Data/optics_SS.v3_3_.nc
Data/optics_SU.v1_3_.nc
)

CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${geos-aero_test_data} )

list( APPEND geos-aero_test_input
testinput/geosaod.rc
testinput/Chem_MieRegistry.rc)

CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${geos-aero_test_input} )

endif( ${geos-aero_FOUND} )

add_subdirectory(testinput/instrumentTests)
add_subdirectory(testinput/unit_tests)
