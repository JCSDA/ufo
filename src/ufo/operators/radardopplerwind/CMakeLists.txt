# (C) Copyright 2024 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( radardopplerwind_files
    ObsRadarDopplerWind.h
    ObsRadarDopplerWind.cc
    ObsRadarDopplerWindData.cc
    ObsRadarDopplerWindData.h
    ObsRadarDopplerWindParameters.h
    ObsRadarDopplerWindTLAD.h
    ObsRadarDopplerWindTLAD.cc
)
PREPEND( _p_radardopplerwind_files     "operators/radardopplerwind"     ${radardopplerwind_files} )

set ( radardopplerwind_src_files
      ${_p_radardopplerwind_files}
      PARENT_SCOPE
)
