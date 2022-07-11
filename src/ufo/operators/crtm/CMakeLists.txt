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

if( ${geos-aero_FOUND} )
  set ( lut_files
    ObsAodLUTs.h
    ObsAodLUTs.cc
    ObsAodLUTsTLAD.h
    ObsAodLUTsTLAD.cc
    ObsAodLUTs.interface.F90
    ObsAodLUTs.interface.h
    ObsAodLUTsTLAD.interface.F90
    ObsAodLUTsTLAD.interface.h
    ufo_aodluts_mod.F90
    ufo_aodluts_tlad_mod.F90
    ufo_luts_utils_mod.F90
  )
  set ( crtm_files 
    ${crtm_files} 
    ${lut_files}
  )
endif( ${geos-aero_FOUND} )


PREPEND( _p_crtm_files       "operators/crtm"       ${crtm_files} )

set ( crtm_src_files
      ${_p_crtm_files}
      PARENT_SCOPE
)

