# (C) Copyright 2020-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Instrument tests for TEMPO

#
ufo_add_test( NAME    instrument_opr_tempo_no2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperator.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tempo_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  compo instrument
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
#
ufo_add_test( NAME    instrument_linopr_tempo_no2
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsOperatorTLAD.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/tempo_no2.yaml"
              MPI     1
              LIBS    ufo
              LABELS  compo instrument
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../ )
