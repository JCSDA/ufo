# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( crtm_files
  crtmParameters/ObsRadianceCRTMParameters.h
  crtmParameters/ObsAodCRTMParameters.h
  ObsRadianceCRTM.h
  ObsRadianceCRTM.cc
  ObsRadianceCRTMTLAD.h
  ObsRadianceCRTMTLAD.cc
  ObsRadianceCRTM.interface.F90
  ObsRadianceCRTM.interface.h
  ObsRadianceCRTMTLAD.interface.F90
  ObsRadianceCRTMTLAD.interface.h
  ufo_crtm_passive_mod.F90
  ufo_crtm_active_mod.F90
  ufo_radiancecrtm_mod.F90
  ufo_radiancecrtm_tlad_mod.F90
  ufo_crtm_utils_mod.F90
  ObsAodCRTM.h
  ObsAodCRTM.cc
  ObsAodCRTMTLAD.h
  ObsAodCRTMTLAD.cc
  ObsAodCRTM.interface.F90
  ObsAodCRTM.interface.h
  ObsAodCRTMTLAD.interface.F90
  ObsAodCRTMTLAD.interface.h
  ufo_aodcrtm_mod.F90
  ufo_aodcrtm_tlad_mod.F90
)

PREPEND( _p_crtm_files       "operators/crtm"       ${crtm_files} )

set ( crtm_src_files
      ${_p_crtm_files}
      PARENT_SCOPE
)

