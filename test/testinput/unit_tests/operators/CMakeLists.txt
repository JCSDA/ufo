# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for operators and TLAD


ufo_add_test( NAME    test_ufo_opr_adt
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_adt
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_aircraft
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_aircraft
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_aod_extinction
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_extinction.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_aod_extinction
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_extinction.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_aod_metoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_aod_metoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_background_error_identity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_diag_background_error_identity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_background_error_vert_interp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_vert_interp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_background_error_vert_interp_pmissing
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_categorical
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/categorical.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_categorical
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/categorical.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_chleuzintegr
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/chleuzintegr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_composite
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_composite
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_composite_GsiSfcModel_SfcPCorrected
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_GsiSfcModel_SfcPCorrected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_composite_vertInterp_SfcPCorrected
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_vertInterp_SfcPCorrected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_coolskin
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/coolskin.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_coolskin
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/coolskin.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssroref
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssroref.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssroref
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssroref.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrorefmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssrorefmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobendmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssrobendmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobendmetoffice_nopseudo
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nopseudo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssrobendmetoffice_nopseudo
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nopseudo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobendmetoffice_profile
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_profile.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobendmetoffice_nosupercheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nosupercheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssrobendmetoffice_nosupercheck
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nosupercheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobndnbam
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gnssrobndnbam
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobndnbam_super_refraction
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_super_refraction.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobndnbam_0obs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_0obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobndnbam_1obs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_1obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gnssrobndnbam_1obs_2pe
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_1obs.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_groundgnssmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_groundgnssmetoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Tests specific to ROPP

if( ${ropp-ufo_FOUND} )

    ufo_add_test( NAME    test_ufo_opr_groundgnssropp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssropp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_gnssrobndropp1d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp1d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_gnssrobndropp1d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp1d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_gnssrobndropp2d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_gnssroBndROPP2d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if ( ROPP_ARO)
    ufo_add_test( NAME    test_ufo_opr_gnssarobndropp2d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssarobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_gnssaroBndROPP2d
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssarobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif()

endif( ${ropp-ufo_FOUND} )
#
ufo_add_test( NAME    test_ufo_opr_gome_metop-a
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gome_metop-a
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gome_metop-a_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_gome_metop-a_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_gsisfcmodel
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gsisfcmodel.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_identity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_identity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_product
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/product.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_product
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/product.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_insitupm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/insitupm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_insitupm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/insitupm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_mopitt_co
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mopitt_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_mopitt_co
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mopitt_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_obsdiag_background_error_upper_air_and_surface_variables
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_upper_air_and_surface_variables.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_obsdiag_background_error_vert_interp_air_pressure
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_air_pressure.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_obsdiag_background_error_vert_interp_height
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_height.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_omi_aura
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_omi_aura
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_omi_aura_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_omi_aura_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_ompslp_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompslp_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompslp_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompslp_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_ompsnp_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompsnp_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_ompsnp_npp_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_flipz.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompsnp_npp_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_flipz.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_ompsnp_npp_L127
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_L127.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompsnp_npp_L127
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_L127.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_ompstc_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompstc_npp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    test_ufo_opr_ompstc_npp_columnop
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompstc_npp_columnop
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    test_ufo_opr_ompsnp_npp_columnop
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompsnp_npp_columnop
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    test_ufo_opr_ompstc_npp_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_ompstc_npp_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_radiosonde
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radiosonde.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_radiosonde
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radiosonde.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
              #
ufo_add_test( NAME    test_ufo_opr_pibal
             TIER    1
             ECBUILD
             COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
             ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pibal.yaml"
             MPI     4
             LIBS    ufo
             LABELS  operators
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_pibal
             TIER    1
             ECBUILD
             COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
             ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pibal.yaml"
             MPI     1
             LIBS    ufo
             LABELS  operators
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_nearest_neighbor_soil
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hofx_smc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_radialvelocity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radialvelocity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_reflectivity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reflectivity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_reflectivity_oucaps
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reflectivity_oucaps.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_opr_reflectivity_mpas
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mrms_reflectivity_mpas_test.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_opr_radar_doppler_wind
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_doppler_wind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_radar_doppler_wind
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_doppler_wind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_sattcwv
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sattcwv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_sattcwv
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sattcwv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_lightning
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/goesr_glm_fed_lightning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_linopr_lightning
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/goesr_glm_fed_lightning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )

#
ufo_add_test( NAME    test_ufo_opr_satwind
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satwind.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_satwind_metoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satwind_metoffice.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_sbuv2_n19
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_sbuv2_n19
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_sbuv2_n19_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_sbuv2_n19_flipz
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_sbuv2_n19_L127
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_L127.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_sbuv2_n19_L127
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_L127.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_scatwind_neutral_metoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwind_neutral_metoffice.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_scatwind_neutral_metoffice
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwind_neutral_metoffice_tlad.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_seaicefrac
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicefrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_seaicefrac
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicefrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_seaicethick
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicethick.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_seaicethick
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicethick.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_identity_sst
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sea_surface_temp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_identity_sst
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sea_surface_temp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_sfcpcorrected
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sfcpcorrected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_timeoper
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_timeoper_pe4
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_timeoper
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_tropomi_no2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_tropomi_no2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_vert_interp_pmissing
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_vert_interp_pmissing
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_tropomi_co
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_linopr_tropomi_co
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_opr_windspeed
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/windspeed.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Tests specific to marine-gsw

if( ${gsw_FOUND} )

    ufo_add_test( NAME    test_ufo_opr_marinevertinterp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/genericprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_marinevertinterp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/genericprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_insitutemperature
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_insitutemperature
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_composite_marineinsitu
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_marineinsitu.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_composite_marineinsitu
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_marineinsitu.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

endif( ${gsw_FOUND} )

# Tests that require CRTM

if( crtm_FOUND )

    ufo_add_test( NAME    test_ufo_opr_abi_ahi_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/abi_ahi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_abi_ahi_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/abi_ahi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_airs_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_crtm.yaml"
              MPI     4
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_airs_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_amsr2_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_amsr2_crtm
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
    ufo_add_test( NAME    test_ufo_opr_aod_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_aod_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_crtm_vis_albedo
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vis_albedo_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
#
    ufo_add_test( NAME    test_ufo_opr_amsua_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_amsua_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_amsua_geos_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_geos.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_amsua_geos_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_geos.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_crtm_optional_lut
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm_optional_lut.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )
#
    ufo_add_test( NAME    test_ufo_linopr_atms_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_cris_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_cris_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_gmi_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_gmi_crtm
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_gmi_scan_angles
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_scan_angles.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_gmi_scan_angles
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_scan_angles.yaml"
              MPI     4
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_hirs4_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hirs4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_hirs4_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hirs4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_iasi_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtm.yaml"
              MPI     2
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_iasi_crtm
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_mhs_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_mhs_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_seviri_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_seviri_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_smap_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/smap_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_smap_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/smap_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_sndrd1-4_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sndrd1-4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_sndrd1-4_crtm
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sndrd1-4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
#    ufo_add_test( NAME    test_ufo_obsdiag_crtm_airs_jacobian
#              TIER    2
#              ECBUILD
#              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
#              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_airs_jacobian.yaml"
#              MPI     1
#              LIBS    ufo
#              LABELS  crtm operators
#              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_obsdiag_crtm_airs_optics
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_airs_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
#    ufo_add_test( NAME    test_ufo_obsdiag_crtm_amsua_jacobian
#              TIER    1
#              ECBUILD
#              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
#              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_amsua_jacobian.yaml"
#              MPI     1
#              LIBS    ufo
#              LABELS  crtm operators
#              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_obsdiag_crtm_amsua_optics
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_amsua_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
#    ufo_add_test( NAME    test_ufo_obsdiag_crtm_atms_jacobian
#              TIER    1
#              ECBUILD
#              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
#              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_atms_jacobian.yaml"
#              MPI     1
#              LIBS    ufo
#              LABELS  crtm operators
#              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_obsdiag_crtm_atms_optics
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_atms_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
#    ufo_add_test( NAME    test_ufo_obsdiag_crtm_cris_jacobian
#              TIER    2
#              ECBUILD
#              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
#              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_cris_jacobian.yaml"
#              MPI     1
#              LIBS    ufo
#              LABELS  crtm operators
#              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_obsdiag_crtm_cris_optics
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_cris_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
#    ufo_add_test( NAME    test_ufo_obsdiag_crtm_iasi_jacobian
#              TIER    2
#              ECBUILD
#              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
#              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_iasi_jacobian.yaml"
#              MPI     1
#              LIBS    ufo
#              LABELS  crtm operators
#              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_obsdiag_crtm_iasi_optics
              TIER    2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_iasi_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

endif( crtm_FOUND )

# Tests that require RTTOV

if( ${rttov_FOUND} )

    ufo_add_test( NAME    test_ufo_opr_amsua_rttovcpp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_rttovcpp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_amsua_rttovcpp
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_rttovcpp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_rttov_allsurfacetypes
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_allsurfacetypes.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_rttov_bad_data
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_bad_data.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_rttov_ops
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_ops.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atms_rttov_multiplatform
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_multiplatform.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_atms_rttov_ops
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_ops.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_atovs_rttov_scatt
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttov_scatt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_atovs_rttov_scatt
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttov_scatt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_opr_iasi_rttov_ops
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_iasi_rttov_ops
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    test_ufo_linopr_iasi_rttov_ops_2pe
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops.yaml"
              MPI     2
              LIBS    ufo
              LABELS  rttov operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

endif( ${rttov_FOUND} )

if( ${oasim_FOUND} )
#
    ufo_add_test( NAME    test_ufo_opr_pace_oasim
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pace_oasim.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${oasim_FOUND} )

# Missing tests at this time
#  atms_rttov_ops_gasunits.yaml
