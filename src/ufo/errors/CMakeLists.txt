# (C) Copyright 2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( errors_files
    ObsErrorCrossVarCov.cc
    ObsErrorCrossVarCov.h
    ObsErrorWithinGroupCov.cc
    ObsErrorWithinGroupCov.h
    ObsErrorDiagonal.cc
    ObsErrorDiagonal.h
)
PREPEND( _p_errors_files     "errors"     ${errors_files} )

set ( errors_src_files
      ${_p_errors_files}
      PARENT_SCOPE
)
