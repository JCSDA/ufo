# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( utils_files
      vert_interp.F90
)

PREPEND( _p_utils_files       "utils"       ${utils_files} )

set ( utils_src_files
      ${_p_utils_files}
      PARENT_SCOPE
)
