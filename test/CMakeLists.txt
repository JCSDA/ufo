# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include macros for creating links and symlinks
include( ufo_functions )
include( git_functions )

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

# Create Data directory for test data and symlink files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)


# Currently only use jedi_release_2.4.0 CRTM coefficients
set( CRTM_COEFFS_BRANCH "jedi_release_2.4.0" )

# Test data filenames
set(ufo_test_data  ufo_testinput_tier_1_1.3.0.tar.gz)
set(crtm_test_data crtm_coefficients_${CRTM_COEFFS_BRANCH}.tar.gz)

#############################################
# macro for downloading test data
# options:
# - tag 1.0.0: download from UCAR ftp, flat directory
#   Note: when we change the tag, this will have to be updated somehow.
# - anything else: download from AWS, either from branch directory
#   if available, or from develop directory
# input: repo (e.g. ioda or ufo)
#        test_data: test data filename (archive)
# output: output_data_path: path to data on local machine after download


function(ADD_DOWNLOAD_TEST repo test_files_dirname output_data_path)
  # If REPO branch is being built set GIT_BRANCH_FUNC to REPO's current branch.
  # If a tagged version of REPO is being built set GIT_TAG_FUNC to REPO's current tag.
  find_branch_name(REPO_DIR_NAME ${repo})

  if ( ${repo} STREQUAL "crtm" )
    set( GIT_TAG_FUNC ${CRTM_COEFFS_BRANCH} )
  endif()

  if( DEFINED GIT_BRANCH_FUNC )
    set( REPO_GIT_BRANCH ${GIT_BRANCH_FUNC} )
  elseif( DEFINED GIT_TAG_FUNC )
    set( REPO_GIT_BRANCH ${GIT_TAG_FUNC} )
  endif()


  # When env veriable LOCAL_PATH_JEDI_TESTFILES is set, simply link test files
  # to build directory. get_${repo}_test_data checks the existence of test file directory.
  # ${repo} test data must be stored in
  # ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${BRANCH}/testinput_tier_1
  if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} )
    set( LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}" )
    message( STATUS "use LOCAL_PATH_JEDI_TESTFILES: ${LOCAL_PATH_JEDI_TESTFILES}" )
    # If test data specific to testing branch exists locally use it.
    # If not use test data specific to develop branch.
    if( EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH} )
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH}" )
    else()
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/develop" )
    endif()

    message( STATUS "Test data in ${TESTFILE_DIR_REPO} is linked to build directory" )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
    ${TESTFILE_DIR_REPO} )
    set ( REPO_DATA_DOWNLOADER ufo_data_checker.py )

  # Any tagged version of REPO is being built.
  # Build ufo_get_{repo}_test_data test to download test data from DASH.
  elseif( DEFINED GIT_TAG_FUNC)
    message( STATUS "Tagged version of ${repo} is used" )
    # set ARGS for get_${repo}_test_data

    # ECBUILD_DOWNLOAD_BASE_URL env var can be used to force test files and
    # crtm coef to be downloaded from other databases such as S3 instead of DASH
    # example ECBUILD_DOWNLOAD_BASE_URL=https://jedi-test-files.s3.amazonaws.com

    if( DEFINED ENV{ECBUILD_DOWNLOAD_BASE_URL} )
      set( ECBUILD_DOWNLOAD_BASE_URL "$ENV{ECBUILD_DOWNLOAD_BASE_URL}/${repo}/${GIT_TAG_FUNC}" )
    else()
      set( ECBUILD_DOWNLOAD_BASE_URL https://gdex.ucar.edu/dataset/jedi-skylab/file )
    endif()

    if( DEFINED ENV{ECBUILD_DOWNLOAD_BASE_URL_CRTM} )
      set( ECBUILD_DOWNLOAD_BASE_URL "$ENV{ECBUILD_DOWNLOAD_BASE_URL_CRTM}/${repo}/${GIT_TAG_FUNC}" )
    endif()


    set( checksum "0" )
    set( TESTFILE_DIR_REPO "${CMAKE_SOURCE_DIR}/test-data-release/${repo}/${GIT_TAG_FUNC}" )
    # Create test-data-release in source directory
    file( MAKE_DIRECTORY ${TESTFILE_DIR_REPO} )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${ECBUILD_DOWNLOAD_BASE_URL}
          ${CMAKE_SOURCE_DIR}/test-data-release
          ${test_files_dirname}
          ${checksum} )

    message( STATUS "Test files will be downloaded from ${ECBUILD_DOWNLOAD_BASE_URL} and
    saved to ${TESTFILE_DIR_REPO} for ${repo}")

    set ( REPO_DATA_DOWNLOADER ufo_data_downloader.py )

  # Any branch of REPO is being built.
  # repo-data repository is already cloned by bundle/CMakeLists.txt.
  # Link test files in repo-data repo to build directory.
  # get_repo_test_data checks the existence of test file directory
  else()
    set( TESTFILE_DIR_REPO ${CMAKE_SOURCE_DIR}/${repo}-data )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${TESTFILE_DIR_REPO} )
    set (REPO_DATA_DOWNLOADER ufo_data_checker.py )
  endif()

  # add test for downloading data
  ecbuild_add_test( TARGET    ufo_get_${repo}_test_data
                    TYPE      SCRIPT
                    COMMAND   ${CMAKE_BINARY_DIR}/bin/${REPO_DATA_DOWNLOADER}
                    ARGS  ${REPO_DATA_DOWNLOADER_ARGS} )

  # pass output data path back to calling code
  set(${output_data_path} ${TESTFILE_DIR_REPO} PARENT_SCOPE)

  # Create download script for ufo_get_*_test_data test
  set ( FILENAME ${REPO_DATA_DOWNLOADER})
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()


  add_custom_target( bin_ufo_${repo}_test_download_scripts ALL
    COMMAND chmod +x ${DEST_FILE}
    DEPENDS ${DEST_FILE} )

endfunction()
#############################################
# add tests to download ufo & crtm data

ADD_DOWNLOAD_TEST( "ufo" ${ufo_test_data} UFO_UFO_TESTFILES_PATH )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${UFO_UFO_TESTFILES_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/ufo)

if( crtm_FOUND )
  ADD_DOWNLOAD_TEST( "crtm" ${crtm_test_data} UFO_CRTM_TESTFILES_PATH )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                   ${UFO_CRTM_TESTFILES_PATH}
                   ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm)
endif( crtm_FOUND )
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

####################################################################
# Establish test tiering

set( UFO_TEST_TIER 1 CACHE STRING "UFO test tier. Run tests up to and including this tier.")

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
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_dataextractor
                  SOURCES mains/TestDataExtractor.cc
                  # This test doesn't need a configuration file, but oops::Run::Run() requires
                  # a path to a configuration file to be passed in the first command-line parameter.
                  ARGS    "testinput/empty.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

# Test Parameters
if (eckit_VERSION VERSION_GREATER_EQUAL 1.16)
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo)
else()
  # Due to changes in eckit 1.16 a different input file is required.
  ecbuild_add_test( TARGET  test_ufo_parameters
                    SOURCES mains/TestParameters.cc
                    ARGS    "testinput/parameters_older_eckit.yaml"
                    ENVIRONMENT OOPS_TRAPFPE=1
                    LIBS    ufo)
endif()

# Test operator utils
ecbuild_add_test( TARGET  test_ufo_operator_utils
                  SOURCES mains/TestOperatorUtils.cc
                  ARGS    "testinput/operatorutils.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo)

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
                  LIBS    ufo)

ecbuild_add_test( TARGET  test_ufo_geovals
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_geovals_mpi
                  SOURCES mains/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  MPI     4
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_geovals_spec
                  SOURCES mains/TestGeoVaLsSpec.cc
                  ARGS    "testinput/geovals_spec.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_locations
                  SOURCES mains/TestLocations.cc ufo/Locations.h ufo/locations_test.F90
                  ARGS    "testinput/locations.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  LIBS    ufo
                  TEST_DEPENDS ufo_get_ufo_test_data )

ecbuild_add_test( TARGET  test_ufo_recordhandler
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_RecordHandler.x
                  ARGS    "testinput/recordhandler.yaml"
                  ENVIRONMENT OOPS_TRAPFPE=1
                  DEPENDS test_RecordHandler.x
                  TEST_DEPENDS ufo_get_ufo_test_data )

#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
CloudCoeff/Little_Endian/CloudCoeff.bin
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
TauCoeff/ODPS/Little_Endian/atms_n20.TauCoeff.bin
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
SpcCoeff/Little_Endian/mhs_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-c.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-c.TauCoeff.bin
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
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
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
SpcCoeff/Little_Endian/abi_gr.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_gr.TauCoeff.bin
SpcCoeff/Little_Endian/ahi_himawari8.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ahi_himawari8.TauCoeff.bin
SpcCoeff/Little_Endian/radiometer_smap.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/radiometer_smap.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/avhrr3_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/avhrr3_n18.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n20.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n20.TauCoeff.bin
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
rttov7pred54L/rtcoef_metop_1_amsuamhs.dat
rttov7pred54L/rtcoef_metop_1_iasi.dat
rttov7pred54L/rtcoef_noaa_19_amsua.dat
rttov7pred54L/rtcoef_noaa_20_atms.dat
mietable/mietable_metop_amsuamhs.dat
)

# Symlink RTcoef files
CREATE_SYMLINK_FILENAME( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${rttov_test_input} )

file(ARCHIVE_EXTRACT INPUT ${rttov_SOURCE_DIR}/rtcoef_rttov12/mietable/mietable_metop_amsuamhs.tgz
                     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Data )

endif( ${rttov_FOUND} )

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
