# (C) Copyright 2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( ${geos-aero_FOUND} )
  set ( aodlut_src_files
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
  PARENT_SCOPE
  )
endif( ${geos-aero_FOUND} )
