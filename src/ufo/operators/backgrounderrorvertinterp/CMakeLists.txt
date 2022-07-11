# (C) Copyright 2021 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( backgrounderrorvertinterp_files
      ObsBackgroundErrorVertInterp.cc
      ObsBackgroundErrorVertInterp.h
      ObsBackgroundErrorVertInterp.interface.F90
      ObsBackgroundErrorVertInterp.interface.h
      ufo_backgrounderrorvertinterp_mod.F90
    )

PREPEND( _p_backgrounderrorvertinterp_files
         "operators/backgrounderrorvertinterp"
         ${backgrounderrorvertinterp_files}
       )

set ( backgrounderrorvertinterp_src_files
      ${_p_backgrounderrorvertinterp_files}
      PARENT_SCOPE
    )
