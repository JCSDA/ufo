# (C) Crown copyright 2024, Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( superob_files
      SuperObBase.cc
      SuperObBase.h
      SuperObMeanOmB.cc
      SuperObMeanOmB.h
      SuperObRadar.cc
      SuperObRadar.h
      SuperObRadarTemplate.cc
      SuperObRadarTemplate.h
)

PREPEND( _p_superob_files       "superob"       ${superob_files} )

set ( superob_src_files
      ${_p_superob_files}
      PARENT_SCOPE
)
