# (C) Copyright 2023 - UCAR.
# (C) Crown Copyright 2024 - UK Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( logarithm_files
    ObsLogarithm.h
    ObsLogarithm.cc
    ObsLogarithmParameters.h
    ObsLogarithmTLAD.h
    ObsLogarithmTLAD.cc
)

PREPEND( _p_logarithm_files       "operators/logarithm"       ${logarithm_files} )

set ( logarithm_src_files
      ${_p_logarithm_files}
      PARENT_SCOPE
)

