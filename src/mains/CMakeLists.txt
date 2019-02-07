# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

ecbuild_add_executable( TARGET  ufo_crtm.x
                        SOURCES ufoRunCRTM.cc
                        LIBS    ufo
                       )

if ( ${RTTOV_FOUND} )
  ecbuild_add_executable( TARGET  ufo_rttov.x
                          SOURCES ufoRunRTTOV.cc
                          LIBS    ufo
                         )
endif( ${RTTOV_FOUND} )
