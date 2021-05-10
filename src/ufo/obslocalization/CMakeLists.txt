# (C) Copyright 2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( obslocalization_files
      ObsLocGC99.h
      ObsLocSOAR.h
      ObsLocSOARParameters.h
      ObsLocalization.h
      ObsLocParameters.cc
      ObsLocParameters.h
)

PREPEND( _p_obslocalization_files       "obslocalization"       ${obslocalization_files} )

set ( obslocalization_src_files
      ${_p_obslocalization_files}
      PARENT_SCOPE
)
