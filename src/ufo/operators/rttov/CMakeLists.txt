# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( rttov_files
    ObsRadianceRTTOV.h
    ObsRadianceRTTOV.cc
    ObsRadianceRTTOVTLAD.h
    ObsRadianceRTTOVTLAD.cc
    ObsRadianceRTTOV.interface.F90
    ObsRadianceRTTOV.interface.h
    ObsRadianceRTTOVTLAD.interface.F90
    ObsRadianceRTTOVTLAD.interface.h
    ObsRadianceRTTOVParameters.h
    ufo_radiancerttov_mod.F90
    ufo_radiancerttov_tlad_mod.F90
    ufo_radiancerttov_utils_mod.F90
)

PREPEND( _p_rttov_files       "operators/rttov"       ${rttov_files} )

set ( rttov_src_files
      ${_p_rttov_files}
      PARENT_SCOPE
)

