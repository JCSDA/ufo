# (C) Copyright 2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( obsfunctions_files
      ObsFunction.cc
      ObsFunction.h
      ObsFunctionBase.cc
      ObsFunctionBase.h
      ObsFunctionVelocity.cc
      ObsFunctionVelocity.h
)

PREPEND( _p_obsfunctions_files       "obsfunctions"       ${obsfunctions_files} )

set ( obsfunctions_src_files
      ${_p_obsfunctions_files}
      PARENT_SCOPE
)

