# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for operators and TLAD


ufo_add_test( NAME    opr_adt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_adt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_adt_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt_qc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_adt_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/adt_qc_tlad.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_aircraft
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_aircraft
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aircraft.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_aod_extinction
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_extinction.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_aod_extinction
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_extinction.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_aod_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
 ufo_add_test( NAME    opr_aod_ukca
               ECBUILD
               COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
               ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_ukca.yaml"
               MPI     1
               LIBS    ufo
               LABELS  operators
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_aod_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_background_error_identity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    diag_background_error_identity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_background_error_vert_interp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_vert_interp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_background_error_vert_interp_pmissing
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/background_error_vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_categorical
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/categorical.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_categorical
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/categorical.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_categorical_zero_obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/categorical_zero_obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_chleuzintegr
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/chleuzintegr.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_composite
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_composite
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_composite_vertInterp_SfcPCorrected
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_vertInterp_SfcPCorrected.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_coolskin
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/coolskin.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_coolskin
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/coolskin.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssroref
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssroref.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssroref
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssroref.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrorefmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrorefmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrorefmetoffice_refractivity_changes
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice_refractivity_changes.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrorefmetoffice_refractivity_changes
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrorefmetoffice_refractivity_changes.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrobendmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice_nopseudo
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nopseudo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrobendmetoffice_nopseudo
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nopseudo.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice_profile
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_profile.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice_nosupercheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nosupercheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice_nosupercheck_missingdata
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nosupercheck_missingdata.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrobendmetoffice_nosupercheck
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_nosupercheck.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobendmetoffice_refractivity_changes
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_refractivity_changes.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrobendmetoffice_refractivity_changes
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobendmetoffice_refractivity_changes.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobndnbam
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gnssrobndnbam
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobndnbam_super_refraction
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_super_refraction.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobndnbam_0obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_0obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobndnbam_1obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_1obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gnssrobndnbam_1obs_2pe
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndnbam_1obs.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_groundgnssmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_groundgnssmetoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssmetoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Tests specific to ROPP

if( ${ropp-ufo_FOUND} )

    ufo_add_test( NAME    opr_groundgnssropp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/groundgnssropp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_gnssrobndropp1d
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp1d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_gnssrobndropp1d
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp1d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_gnssrobndropp2d
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_gnssroBndROPP2d
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssrobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if ( ROPP_ARO)
    ufo_add_test( NAME    opr_gnssarobndropp2d
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gnssarobndropp2d.yaml"
              MPI     1
              LIBS    ufo
              LABELS  ropp-ufo operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_gnssaroBndROPP2d
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
ufo_add_test( NAME    opr_gome_metop-a
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gome_metop-a
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_gome_metop-a_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_gome_metop-a_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gome_metop-a_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_identity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_identity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/identity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_logarithm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/logarithm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_logarithm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/logarithm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_product
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/product.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_product
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/product.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_insitupm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/insitupm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_insitupm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/insitupm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_mopitt_co
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mopitt_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_mopitt_co
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mopitt_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_upper_air_and_surface_variables
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_upper_air_and_surface_variables.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_vert_interp_air_pressure
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_air_pressure.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_vert_interp_air_pressure_obsspace_scaling
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_air_pressure_obsspace_scaling.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_vert_interp_air_pressure_obsspace_missing_scaling
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_air_pressure_obsspace_missing_scaling.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_vert_interp_air_pressure_geovals_scaling
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_air_pressure_geovals_scaling.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_obsdiag_background_error_vert_interp_height
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_background_error_vert_interp_height.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_omi_aura
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_omi_aura
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_omi_aura_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_omi_aura_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/omi_aura_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ompslp_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompslp_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompslp_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompslp_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ompsnp_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompsnp_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ompsnp_npp_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_flipz.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompsnp_npp_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_flipz.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ompsnp_npp_L127
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_L127.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompsnp_npp_L127
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_L127.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ompstc_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompstc_npp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    opr_ompstc_npp_columnop
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompstc_npp_columnop
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    opr_ompsnp_npp_columnop
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompsnp_npp_columnop
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompsnp_npp_columnop.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

#
ufo_add_test( NAME    opr_ompstc_npp_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_ompstc_npp_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ompstc_npp_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_radiosonde
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radiosonde.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_radiosonde
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radiosonde.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
              #
ufo_add_test( NAME    opr_pibal
             ECBUILD
             COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
             ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pibal.yaml"
             MPI     4
             LIBS    ufo
             LABELS  operators
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_pibal
             ECBUILD
             COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
             ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pibal.yaml"
             MPI     1
             LIBS    ufo
             LABELS  operators
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    nearest_neighbor_soil
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hofx_smc.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_radialvelocity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radialvelocity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_reflectivity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reflectivity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_reflectivity_oucaps
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/reflectivity_oucaps.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    opr_radar_reflectivity_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_reflectivity_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_radar_reflectivity_metoffice_composite
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_reflectivity_metoffice_composite.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_radar_doppler_wind
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_doppler_wind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_radar_doppler_wind
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_doppler_wind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_radar_reflectivity_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/radar_reflectivity_metoffice.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ppro_nssl
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ppro_nssl.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ppro_tcwa2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ppro_tcwa2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sattcwv
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sattcwv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_sattcwv
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sattcwv.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_lightning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/goesr_glm_fed_lightning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    linopr_lightning
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/goesr_glm_fed_lightning.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )

#
ufo_add_test( NAME    opr_satwind
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satwind.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_satwind_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/satwind_metoffice.yaml"
              MPI     2
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sbuv2_n19
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_sbuv2_n19
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sbuv2_n19_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_sbuv2_n19_flipz
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_flipz.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sbuv2_n19_L127
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_L127.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_sbuv2_n19_L127
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sbuv2_n19_L127.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_scatwind_neutral_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwind_neutral_metoffice.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_scatwind_neutral_metoffice
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/scatwind_neutral_metoffice_tlad.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_seaicefrac
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicefrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_seaicefrac
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicefrac.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_seaicethick
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicethick.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_seaicethick
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seaicethick.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_identity_sst
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sea_surface_temp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_identity_sst
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sea_surface_temp.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sfccorrected_pressure
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sfccorrected_pressure.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_sfccorrected_temperature
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sfccorrected_temperature.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_timeoper
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_timeoper_pe4
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_timeoper
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/timeoper.yaml"
              MPI     4
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_tropomi_no2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_tropomi_no2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_vert_interp_pmissing
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_vert_interp_pmissing
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_vert_interp_pmissing_scaled
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing_scaled.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_vert_interp_pmissing_scaled
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vert_interp_pmissing_scaled.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_tropomi_co
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_tropomi_co
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tropomi_co.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_windspeed
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/windspeed.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_ionosonde
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ionosonde.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_pathsum_vertical
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pathsum_vertical.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_pathsum_slant
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/pathsum_slant.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    opr_densityreduction_simple
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/densityreduction_simple.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_densityreduction_simple
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/densityreduction_simple.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
ufo_add_test( NAME    opr_densityreduction_latlon
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/densityreduction_latlon.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_densityreduction_latlon
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/densityreduction_latlon.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

# Tests specific to marine-gsw

if( ${gsw_FOUND} )

    ufo_add_test( NAME    opr_marinevertinterp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/genericprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_marinevertinterp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/genericprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_insitutemperature
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_insitutemperature
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tprof.yaml"
              MPI     1
              LIBS    ufo
              LABELS  gsw operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_composite_marineinsitu
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/composite_marineinsitu.yaml"
              MPI     1
              LIBS    ufo
              LABELS  operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_composite_marineinsitu
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

    ufo_add_test( NAME    opr_abi_ahi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/abi_ahi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_abi_ahi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/abi_ahi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_airs_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_crtm.yaml"
              MPI     4
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_airs_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/airs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_amsr2_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_amsr2_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsr2_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
    ufo_add_test( NAME    opr_aod_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_aod_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/aod_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_crtm_vis_albedo
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/vis_albedo_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( crtm_FOUND AND crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
#
    ufo_add_test( NAME    opr_amsua_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_amsua_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_amsua_geos_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_geos.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_amsua_geos_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_crtm_geos.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_crtm_optional_lut
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm_optional_lut.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data ufo_get_crtm_test_data )
#
    ufo_add_test( NAME    linopr_atms_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_cris_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_cris_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_cris_crtm_co2_options
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm_co2_options.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_cris_crtm_co2_options
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/cris_crtm_co2_options.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_gmi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_gmi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_gmi_scan_angles
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_scan_angles.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_gmi_scan_angles
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/gmi_scan_angles.yaml"
              MPI     4
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_hirs4_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hirs4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_hirs4_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/hirs4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_iasi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtm.yaml"
              MPI     2
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    iasi_crtmreconrad_linop_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtmreconrad_linop_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    iasi_crtmreconrad_linop_tlad_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtmreconrad_linop_check.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_iasi_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_mhs_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_mhs_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mhs_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_seviri_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_seviri_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/seviri_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_smap_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/smap_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    linopr_smap_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/smap_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_sndrd1-4_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sndrd1-4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_sndrd1-4_crtm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/sndrd1-4_crtm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_airs_jacobian
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_airs_jacobian.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_airs_optics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_airs_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_amsua_jacobian
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_amsua_jacobian.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_amsua_optics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_amsua_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_atms_jacobian
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_atms_jacobian.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_atms_optics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_atms_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_cris_jacobian
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_cris_jacobian.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_cris_optics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_cris_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_iasi_jacobian
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_iasi_jacobian.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    obsdiag_crtm_iasi_optics
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsDiagnostics.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/obsdiag_crtm_iasi_optics.yaml"
              MPI     1
              LIBS    ufo
              LABELS  crtm operators tier2
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

endif( crtm_FOUND )

# Tests that require RTTOV

if( ${rttov_FOUND} )

# RTTOVCPP tests (optional, controlled by ENABLE_RTTOVCPP_TESTS option)
if(ENABLE_RTTOVCPP_TESTS)
    ufo_add_test( NAME    opr_amsua_rttovcpp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_rttovcpp${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} rttovcpp operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_amsua_rttovcpp
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/amsua_rttovcpp${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} rttovcpp operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif(ENABLE_RTTOVCPP_TESTS)
#
    ufo_add_test( NAME    linopr_atms_rttov_tlad_qc
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_tlad_qc${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL}operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_rttov_allsurfacetypes
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_allsurfacetypes${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_rttov_bad_data
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_bad_data${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_rttov_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     2
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atms_rttov_multiplatform
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_multiplatform${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_atms_rttov_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atms_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     2
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_atovs_rttov_scatt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttov_scatt${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_atovs_rttov_scatt
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/atovs_rttov_scatt${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_iasi_rttov_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_iasi_rttov_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    linopr_iasi_rttov_ops_2pe
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/iasi_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     2
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    mtgirs_reconrad_linop_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mtgirs_reconrad_linop_check${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    mtgirs_reconrad_linop_tlad_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mtgirs_reconrad_linop_check${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    mtgirs_reconrad_linop_tlad_2pe_check
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/mtgirs_reconrad_linop_check${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     2
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators ci_oneapi_disable
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
    ufo_add_test( NAME    opr_ssmis_rttov_ops
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/ssmis_rttov_ops${RTTOV_VERSION_SUFFIX}.yaml"
              MPI     1
              LIBS    ufo
              LABELS  rttov ${RTTOV_VERSION_LABEL} operators
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )

endif( ${rttov_FOUND} )

if( ${oasim_FOUND} )
#
    ufo_add_test( NAME    opr_pace_oasim
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
