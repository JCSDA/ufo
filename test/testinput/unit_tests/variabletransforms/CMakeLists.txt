# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for variabletransforms

#
ufo_add_test( NAME    test_ufo_observed_derived_and_simulated_variables
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/observed_derived_and_simulated_variables.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( ${gsw_FOUND} )
    ufo_add_test( NAME    test_ufo_variabletransforms_oceandensity
                TIER    1
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceandensity.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    test_ufo_variabletransforms_oceandepthtopressure
                TIER    1
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceandepthtopressure.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    test_ufo_variabletransforms_oceantemptotheta
                TIER    1
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceantemptotheta.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    test_ufo_variabletransforms_oceantemptoconservativetemp
                TIER    1
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceantemptoconservativetemp.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    test_ufo_variabletransforms_oceanpracticalsalinitytoabsolutesalinity
                TIER    1
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceanpracticalsalinitytoabsolutesalinity.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${gsw_FOUND} )
#
ufo_add_test( NAME    test_ufo_variabletransforms_potentialT
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_potentialT.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_profilehorizontaldrift
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_profilehorizontaldrift_zero_obs
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift_zero_obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_profilehorizontaldrift_keepinwindow
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift_keepinwindow.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_pstar
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_pstar.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_QFEtoQNH
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_QFEtoQNH.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_rhumidity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_rhumidity_part2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity_part2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_rhumidity_part3
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity_part3.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_shumidity
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_shumidity.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_shumidity2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_shumidity2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_virtualTemperature
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_virtualtemperature.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_SondeHeightFromPressure
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondeHeightFromPressure.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_SondePressureFromHeight
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondePressureFromHeight.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_SondePressureFromHeightICAO
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondePressureFromHeightICAO.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_szacorrection
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_szacorrection.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_windcomponents
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windcomponents.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_windspeedanddirection
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windspeedanddirection.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_windcomponents_scatwind
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windcomponents_scatwind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    test_ufo_variabletransforms_btfromradiance
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_btfromradiance.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms

              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_variabletransforms_scanposition
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_scanposition.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    test_ufo_variabletransforms_radarbeamgeometry
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_radarbeamgeometry.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )

