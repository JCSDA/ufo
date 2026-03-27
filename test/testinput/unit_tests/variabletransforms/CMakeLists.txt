# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for variabletransforms

#
ufo_add_test( NAME    observed_derived_and_simulated_variables
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/observed_derived_and_simulated_variables.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_logarithm
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_logarithm.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_icefreeboardtothickness
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_icefreeboardtothickness.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
if( ${gsw_FOUND} )
    ufo_add_test( NAME    variabletransforms_oceandensity
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceandensity.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    variabletransforms_oceandepthtopressure
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceandepthtopressure.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    variabletransforms_oceantemptotheta
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceantemptotheta.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    variabletransforms_oceantemptoconservativetemp
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceantemptoconservativetemp.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
    #
    ufo_add_test( NAME    variabletransforms_oceanpracticalsalinitytoabsolutesalinity
                ECBUILD
                COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
                ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_oceanpracticalsalinitytoabsolutesalinity.yaml"
                MPI     1
                LIBS    ufo
                LABELS  unit_tests variabletransforms
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
endif( ${gsw_FOUND} )
#
ufo_add_test( NAME    variabletransforms_potentialT
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_potentialT.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_profilehorizontaldrift
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_profilehorizontaldrift_zero_obs
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift_zero_obs.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_profilehorizontaldrift_keepinwindow
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_profilehorizontaldrift_keepinwindow.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_pstar
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_pstar.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_QFEtoQNH
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_QFEtoQNH.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_rhumidity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_rhumidity_part2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity_part2.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_rhumidity_part3
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_rhumidity_part3.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_shumidity
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_shumidity.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_shumidity2
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_shumidity2.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_virtualTemperature
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_virtualtemperature.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_SondeHeightFromPressure
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondeHeightFromPressure.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_SondePressureFromHeight
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondePressureFromHeight.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_SondePressureFromHeightICAO
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_SondePressureFromHeightICAO.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_szacorrection
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_szacorrection.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_windcomponents
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windcomponents.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_windspeedanddirection
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windspeedanddirection.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_windcomponents_scatwind
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_windcomponents_scatwind.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    variabletransforms_btfromradiance
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_btfromradiance.yaml"
              MPI     2
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    variabletransforms_radiancefrompcscores
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_radiancefrompcscores.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    variabletransforms_radiancefrompcscores_mean
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_radiancefrompcscores_mean.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
#
ufo_add_test( NAME    variabletransforms_scanposition
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_scanposition.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )
#
ufo_add_test( NAME    variabletransforms_radarbeamgeometry
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/variabletransforms_radarbeamgeometry.yaml"
              MPI     1
              LIBS    ufo
              LABELS  unit_tests variabletransforms
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../
              TEST_DEPENDS ufo_get_ufo_test_data )

