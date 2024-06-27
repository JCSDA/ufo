# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for filters/obsfunctions


if( crtm_FOUND )
    ufo_add_test( NAME    test_ufo_amsr2_clw_ret
                  TIER    2
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_clw_ret.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    test_ufo_function_osw_windspeed
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_osw_windspeed.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_atovs_lwpfirstguess
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_lwpfirstguess.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_atovs_scenedeperrors
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_scenedeperrors.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_arithmetic
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_arithmetic.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_aircraft_bc_offline
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_aircraftBCOffline.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_averagetemperaturebelow
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_averagetemperaturebelow.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_bennartzscatindex
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_bennartzscatindex.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_bgddepartureanomaly
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_bgddepartureanomaly.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
if( rttov_FOUND )
    ufo_add_test( NAME    test_ufo_function_cloudcost
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
    ufo_add_test( NAME    test_ufo_function_cloudcost_scattering
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost_scattering.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
    ufo_add_test( NAME    test_ufo_function_cloudcost_surfwind
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudcost_surfwind.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )

endif( rttov_FOUND )
#
ufo_add_test( NAME    test_ufo_function_clouddetect_airs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_airs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clouddetect_cris
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_cris.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clouddetect_iasi
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clouddetect_iasi.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
if( rttov_FOUND )
    ufo_add_test( NAME    test_ufo_function_cloudfirstguess
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_cloudfirstguess.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  rttov filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( rttov_FOUND )
#
ufo_add_test( NAME    test_ufo_function_clwmatchidx
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwmatchidx.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwmatchidx_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwmatchidx_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_hofx_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_hofx_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_hofx
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_hofx.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_ssmis
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwretmean
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwretmean.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwretmean_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwretmean_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_obsval
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_obsval.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_clwret_obsval_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_clwret_obsval_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_conditional
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_conditional.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_assignvalueequalchannels
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_assignvalueequalchannels.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_datetimeoffset
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_datetimeoffset.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_drawvaluefromfile
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_drawvaluefromfile.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_elementmultiply
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_elementmultiply.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errfgrosschk
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfgrosschk.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errfjsfc
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfjsfc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errflat
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errflat.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errfsdoei
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfsdoei.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errfsdoei_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfsdoei_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errftopo
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errftopo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errftransmittop
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errftransmittop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_errfwavenum
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_errfwavenum.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_exponential
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_exponential.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_hydrometeorchk
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_hydrometeorchk.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_hydrometeorchk_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_hydrometeorchk_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_lam_domaincheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_lam_domaincheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_lam_domaincheck_circle
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_lam_domaincheck_circle.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_linear_combination_int
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_linear_combination_int.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_linear_combination_multichannel
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_linear_combination_multichannel.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_metoffice_rh_corr
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_metoffice_rh_corr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_ModelSurfaceHeightAdjustMarineWinds
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustMarineWinds.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_ModelSurfaceHeightAdjustRelativeHumidity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustRelativeHumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_ModelSurfaceHeightAdjustTemperature
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustTemperature.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_ModelSurfaceHeightAdjustWindVector
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_ModelSurfaceHeightAdjustWindVector.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_nsstret
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_nsstret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrmean
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrmean_abi
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_abi.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrmean_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrmean_mhs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrmean_mhs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrorfactorconv
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorfactorconv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrorfactorquotient
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorfactorquotient.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrorhumidity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorhumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrorinflatesfcp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorinflatesfcp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_obserrorstepwiselinear
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_obserrorstepwiselinear.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_oceanpressuretodepth
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_oceanpressuretodepth.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_orbitangle_ssmis
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_orbitangle_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_profilelevelcount
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_profilelevelcount.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_radarscanedgeflag
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_radarscanedgeflag.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_writerecnumobsspace_1PE
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_writerecnumobsspace.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_writerecnumobsspace_3PE
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_writerecnumobsspace.yaml"
              MPI     3
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_reperr
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_reperr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_satwind_indiv_errors
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwind_indiv_errors.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_satwindsLNVDcheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsLNVDcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_satwindsSPDBcheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsSPDBcheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_scatret
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatret.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_scatret_atms
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatret_atms.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_selectstatistic
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_selectstatistic.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setseaiceemiss
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setseaiceemiss.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype_landfrac
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_landfrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype_default
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_default.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype_airs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_airs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype_ssmis
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_ssmis.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_setsurfacetype_tolerance
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_setsurfacetype_tolerance.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_solarzenith
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_solarzenith.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_solarzenith_skiprejected
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_solarzenith_skiprejected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
	      
# Test for the string manipulation obsfucntion
ufo_add_test( NAME    test_ufo_function_stringmanipulation
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_stringmanipulation.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )	      
	      
#
ufo_add_test( NAME    test_ufo_function_surfacecloudmodellevelcbh
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_surfacecloudmodellevelCBH.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_surfacecloudcreatecloudcolumn
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_surfacecloudcreatecloudcolumn.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_tropo_estimate
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_tropo_estimate.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_velocity_channels
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_velocity_channels.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_velocity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFunction.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_velocity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_wdirdiff
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_wdirdiff.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_scatambiguity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_scatambiguity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_function_satwindserrnorm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/function_satwindsErrnormCheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
if( crtm_FOUND )
    ufo_add_test( NAME    test_ufo_gmi_clw_ret
                  TIER    1
                  ECBUILD
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                  ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_clw_ret.yaml"
                  MPI     1
                  LIBS    ufo
                  LABELS  crtm filters obsfunctions
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
endif( crtm_FOUND )
#
ufo_add_test( NAME    test_ufo_qc_velocitycheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_velocitycheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters obsfunctions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
