# (C) Copyright 2019 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( timeoper_files
    ObsTimeOper.h
    ObsTimeOperTLAD.h
    ObsTimeOperUtil.h
    ObsTimeOperUtil.cc
    ObsTimeOper.cc
    ObsTimeOperTLAD.cc
    ObsTimeOperParameters.h
)
PREPEND( _p_timeoper_files     "operators/timeoper"     ${timeoper_files} )

set ( timeoper_src_files
      ${_p_timeoper_files}
      PARENT_SCOPE
)
