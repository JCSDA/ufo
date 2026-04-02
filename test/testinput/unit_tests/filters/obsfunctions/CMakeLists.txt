# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for filters/obsfunctions


if( crtm_FOUND )
    ufo_add_test( NAME    amsr2_clw_ret
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_clw_ret.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters obsfunctions tier2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    function_osw_windspeed
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_osw_windspeed.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    atovs_lwpfirstguess
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_lwpfirstguess.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    atovs_scenedeperrors
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_scenedeperrors.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_arithmetic
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_arithmetic.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_aircraft_bc_offline
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_aircraftBCOffline.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_averagetemperaturebelow
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_averagetemperaturebelow.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_bennartzscatindex
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_bennartzscatindex.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_bgddepartureanomaly
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_bgddepartureanomaly.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
if( rttov_FOUND )
    ufo_add_test( NAME    function_cloudcost
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
    ufo_add_test( NAME    function_cloudcost_scattering
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost_scattering${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
    ufo_add_test( NAME    function_cloudcost_surfwind
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost_surfwind${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
     ufo_add_test( NAME    function_cloudfirstguess
                   ECBUILD
                   COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                   ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudfirstguess${RTTOV_VERSION_SUFFIX}.yaml"
                   MPI     1
                   LIBS    ufo
                   LABELS  rttov ${RTTOV_VERSION_LABEL} filters obsfunctions
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
    ufo_add_test( NAME    function_cloudtoppressure
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudtoppressure${RTTOV_VERSION_SUFFIX}.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov ${RTTOV_VERSION_LABEL} filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( rttov_FOUND )
#
ufo_add_test( NAME    function_clouddetect_airs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_airs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clouddetect_cris
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_cris.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clouddetect_iasi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_iasi.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwmatchidx
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwmatchidx.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwmatchidx_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwmatchidx_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_hofx_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_hofx_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_hofx
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_hofx.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_ssmis
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwretmean
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwretmean.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwretmean_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwretmean_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_obsval
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_obsval.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_clwret_obsval_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_obsval_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_conditional
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_conditional.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_assignvalueequalchannels
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_assignvalueequalchannels.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_datetimeoffset
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_datetimeoffset.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_drawvaluefromfile
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_drawvaluefromfile.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_elementmultiply
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_elementmultiply.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errfgrosschk
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfgrosschk.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errfjsfc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfjsfc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errflat
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errflat.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errfsdoei
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfsdoei.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errfsdoei_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfsdoei_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errftopo
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errftopo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errftransmittop
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errftransmittop.yaml"
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_errfwavenum
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfwavenum.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_exponential
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_exponential.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_hydrometeorchk
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_hydrometeorchk.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_hydrometeorchk_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_hydrometeorchk_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_lam_domaincheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_lam_domaincheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_lam_domaincheck_circle
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_lam_domaincheck_circle.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_linear_combination_int
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_linear_combination_int.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_linear_combination_multichannel
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_linear_combination_multichannel.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_modellevelindex
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_modellevelindex.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_metoffice_rh_corr
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_metoffice_rh_corr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_ModelSurfaceHeightAdjustMarineWinds
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustMarineWinds.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_ModelSurfaceHeightAdjustRelativeHumidity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustRelativeHumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_ModelSurfaceHeightAdjustTemperature
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustTemperature.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_ModelSurfaceHeightAdjustWindVector
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustWindVector.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_nsstret
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_nsstret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrmean
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrmean_abi
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_abi.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrmean_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrmean_mhs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_mhs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorboundconv
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorboundconv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorfactorconv
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorfactorconv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorfactorpressurecheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorfactorpressurecheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorfactorquotient
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorfactorquotient.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorhumidity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorhumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorinflatesfcp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorinflatesfcp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_obserrorstepwiselinear
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorstepwiselinear.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_oceanpressuretodepth
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_oceanpressuretodepth.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_orbitangle_ssmis
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_orbitangle_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_profilelevelcount
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_profilelevelcount.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_radarscanedgeflag
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_radarscanedgeflag.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_writerecnumobsspace_1PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_writerecnumobsspace.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_writerecnumobsspace_3PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_writerecnumobsspace.yaml"
              MPI     3
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_mpirank_2PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_mpirank_2pe.yaml"
              MPI     2
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_mpirank_3PE
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_mpirank_3pe.yaml"
              MPI     3
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_reperr
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_reperr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_satwind_indiv_errors
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwind_indiv_errors.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_satwindsLNVDcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsLNVDcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_satwindsSPDBcheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsSPDBcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_scatret
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_scatret_atms
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatret_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_selectstatistic
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_selectstatistic.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setseaiceemiss
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setseaiceemiss.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype_aappsg
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_aappsg.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype_landfrac
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_landfrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    function_setsurfacetype_default
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_default.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype_airs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_airs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype_ssmis
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_setsurfacetype_tolerance
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_tolerance.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_solarzenith
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_solarzenith.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_solarzenith_skiprejected
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_solarzenith_skiprejected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
	      
# Test for the string manipulation obsfucntion
ufo_add_test( NAME    function_stringmanipulation
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_stringmanipulation.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )	      
	      
#
ufo_add_test( NAME    function_surfacecloudmodellevelcbh
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_surfacecloudmodellevelCBH.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_surfacecloudcreatecloudcolumn
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_surfacecloudcreatecloudcolumn.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_geocloudcreatecloudcolumn
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_geocloudcreatecloudcolumn.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )

ufo_add_test( NAME    function_tropo_estimate
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_tropo_estimate.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_tropopause_height
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_tropopauseheight.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_velocity_channels
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_velocity_channels.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_velocity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_velocity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_visibilitydiagnostic
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_visibilitydiagnostic.yaml"
              MPI     3
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_wdirdiff
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_wdirdiff.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_scatambiguity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatambiguity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    function_satwindserrnorm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsErrnormCheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
if( crtm_FOUND )
    ufo_add_test( NAME    gmi_clw_ret
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_clw_ret.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    qc_velocitycheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_velocitycheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
