# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( atmsfcinterp_files
    ObsAtmSfcInterp.h
    ObsAtmSfcInterp.cc
    ObsAtmSfcInterpTLAD.h
    ObsAtmSfcInterpTLAD.cc
    ObsAtmSfcInterp.interface.F90
    ObsAtmSfcInterp.interface.h
    ObsAtmSfcInterpTLAD.interface.F90
    ObsAtmSfcInterpTLAD.interface.h
    ufo_atmsfcinterp_mod.F90
    ufo_atmsfcinterp_tlad_mod.F90
    PARENT_SCOPE
)

PREPEND( _p_atmsfcinterp_files       "atmsfcinterp"       ${atmsfcinterp_files} )

set ( atmsfcinterp_src_files
      ${_p_atmsfcinterp_files}
      PARENT_SCOPE
)

