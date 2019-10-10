# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( timeoper_files
    ObsTimeOper.h
    ObsTimeOper.cc
    ObsTimeOper.interface.F90
    ObsTimeOper.interface.h
    ufo_timeoper_locs_mod.F90
    ufo_timeoper_mod.F90
)
PREPEND( _p_timeoper_files     "timeoper"     ${timeoper_files} )

set ( timeoper_src_files
      ${_p_timeoper_files}
      PARENT_SCOPE
)
