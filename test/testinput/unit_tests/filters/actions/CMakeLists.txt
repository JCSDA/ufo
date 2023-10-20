# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Unit tests for filters/actions


ufo_add_test( NAME    test_ufo_qc_actions
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_actions.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters actions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_qc_performaction
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_performaction.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters actions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
#
ufo_add_test( NAME    test_ufo_qc_performaction_both_action_and_actions
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_performaction_both_action_and_actions.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters actions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
set_tests_properties(ufo_test_tier1_test_ufo_qc_performaction_both_action_and_actions
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*The 'action' and 'actions' options are mutually exclusive.*")
#
ufo_add_test( NAME    test_ufo_qc_performaction_neither_action_nor_actions
              TIER    1
              ECBUILD
              COMMAND ${CMAKE_BINARY_DIR}/bin/test_ObsFilters.x
              ARGS    "${CMAKE_CURRENT_SOURCE_DIR}/qc_performaction_neither_action_nor_actions.yaml"
              MPI     1
              LIBS    ufo
              LABELS  filters actions
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../../ )
set_tests_properties(ufo_test_tier1_test_ufo_qc_performaction_neither_action_nor_actions
                     PROPERTIES
                     PASS_REGULAR_EXPRESSION ".*Either the 'action' or 'actions' option must be set.*")
