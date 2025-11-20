# (C) Copyright UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( pathsum_files
    PathSumOper.h
    PathSumOper.cc
)
PREPEND( _p_pathsum_files     "operators/pathsum"     ${pathsum_files} )

set ( pathsum_src_files
      ${_p_pathsum_files}
      PARENT_SCOPE
)
