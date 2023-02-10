# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( oasim_files
    ObsRadianceOASIM.h
    ObsRadianceOASIM.cc
    ObsRadianceOASIMParameters.h
    ObsRadianceOASIM.interface.F90
    ObsRadianceOASIM.interface.h
    ufo_oasim_mod.F90
)

PREPEND( _p_oasim_files       "operators/oasim"       ${oasim_files} )

set ( oasim_src_files
      ${_p_oasim_files}
      PARENT_SCOPE
)
