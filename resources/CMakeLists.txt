# (C) Copyright 2017-2023 UCAR & Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set(resource_list auxillary bmatrix namemap rmatrix)

file(INSTALL ${resource_list} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

install(
  DIRECTORY ${resource_list}
  DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/"
)
