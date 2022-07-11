# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( radarradialvelocity_files
    ObsRadarRadialVelocity.h
    ObsRadarRadialVelocity.cc
    ObsRadarRadialVelocityTLAD.h
    ObsRadarRadialVelocityTLAD.cc
    ObsRadarRadialVelocity.interface.F90
    ObsRadarRadialVelocity.interface.h
    ObsRadarRadialVelocityTLAD.interface.F90
    ObsRadarRadialVelocityTLAD.interface.h
    ufo_radarradialvelocity_mod.F90
    ufo_radarradialvelocity_tlad_mod.F90
)

PREPEND( _p_radarradialvelocity_files     "operators/radarradialvelocity"     ${radarradialvelocity_files} )

set ( radarradialvelocity_src_files
      ${_p_radarradialvelocity_files}
      PARENT_SCOPE
)
