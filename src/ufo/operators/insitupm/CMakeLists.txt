# (C) Copyright 2021.
#
# This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( insitupm_files
    ObsInsituPM.h
    ObsInsituPM.cc
    ObsInsituPMTLAD.h
    ObsInsituPMTLAD.cc
    ObsInsituPM.interface.F90
    ObsInsituPM.interface.h
    ObsInsituPMParameters.h
    ObsInsituPMTLAD.interface.F90
    ObsInsituPMTLAD.interface.h
    ufo_insitupm_mod.F90
    PM_cmaq_mod.F90
    aero_kinds_mod.F90
    ufo_insitupm_tlad_mod.F90
)

PREPEND( _p_insitupm_files     "operators/insitupm"     ${insitupm_files} )

set ( insitupm_src_files
      ${_p_insitupm_files}
      PARENT_SCOPE
)
