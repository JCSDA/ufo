# (C) Copyright 2021 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( backgrounderroridentity_files
      ObsBackgroundErrorIdentity.cc
      ObsBackgroundErrorIdentity.h
      ObsBackgroundErrorIdentity.interface.F90
      ObsBackgroundErrorIdentity.interface.h
      ufo_backgrounderroridentity_mod.F90
    )

PREPEND( _p_backgrounderroridentity_files
         "operators/backgrounderroridentity"
         ${backgrounderroridentity_files}
       )

set ( backgrounderroridentity_src_files
      ${_p_backgrounderroridentity_files}
      PARENT_SCOPE
    )
