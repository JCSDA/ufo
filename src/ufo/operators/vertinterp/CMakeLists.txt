# (C) Copyright 2017-2023 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( vertinterp_files
  ObsVertInterp.cc
  ObsVertInterp.h
  ObsVertInterpTLAD.cc
  ObsVertInterpTLAD.h
  ObsVertInterp.interface.F90
  ObsVertInterp.interface.h
  ObsVertInterpTLAD.interface.F90
  ObsVertInterpTLAD.interface.h
  ufo_vertinterp_mod.F90
  ufo_vertinterp_tlad_mod.F90
  ufo_interp_param_mod.F90
  ObsVertInterpParameters.cc
  ObsVertInterpParameters.h
)

PREPEND( _p_vertinterp_files       "operators/vertinterp"       ${vertinterp_files} )

set ( vertinterp_src_files
      ${_p_vertinterp_files}
      PARENT_SCOPE
)
