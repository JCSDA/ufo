# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( atmvertinterp_files
  ObsAtmVertInterp.cc
  ObsAtmVertInterp.h
  ObsAtmVertInterpTLAD.cc
  ObsAtmVertInterpTLAD.h
  ObsAtmVertInterp.interface.F90
  ObsAtmVertInterp.interface.h
  ObsAtmVertInterpTLAD.interface.F90
  ObsAtmVertInterpTLAD.interface.h
  ufo_atmvertinterp_mod.F90
  ufo_atmvertinterp_tlad_mod.F90
  ufo_interp_param_mod.F90
  ObsAtmVertInterpParameters.cc
  ObsAtmVertInterpParameters.h
)

PREPEND( _p_atmvertinterp_files       "operators/atmvertinterp"       ${atmvertinterp_files} )

set ( atmvertinterp_src_files
      ${_p_atmvertinterp_files}
      PARENT_SCOPE
)
