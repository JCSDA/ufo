# (C) Copyright 2017-2023 UCAR & Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

file(GLOB resource_list LIST_DIRECTORIES true "*")
foreach(resource IN LISTS resource_list)
  if(NOT IS_DIRECTORY "${resource}")
    list(REMOVE_ITEM resource_list "${resource}")
  endif()
endforeach()

file(INSTALL ${resource_list} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

install(
  DIRECTORY ${resource_list}
  DESTINATION "${CMAKE_INSTALL_SYSCONFDIR}/${PROJECT_NAME}/"
)
