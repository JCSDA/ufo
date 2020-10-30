# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( radarreflectivity_files
    ObsRadarReflectivity.h
    ObsRadarReflectivity.cc
    ObsRadarReflectivity.interface.F90
    ObsRadarReflectivity.interface.h
    ufo_radarreflectivity_mod.F90
)

PREPEND( _p_radarreflectivity_files     "radarreflectivity"     ${radarreflectivity_files} )

set ( radarreflectivity_src_files
      ${_p_radarreflectivity_files}
      PARENT_SCOPE
)
