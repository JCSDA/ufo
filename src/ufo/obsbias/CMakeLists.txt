# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( obsbias_files
  LinearObsBiasBase.h
  LinearObsBiasBase.cc
  ObsBiasBase.h
  ObsBiasBase.cc
  ObsBiasLinearCombination.h
  ObsBiasLinearCombination.cc
  ObsBiasLinearCombinationTLAD.h
  ObsBiasLinearCombinationTLAD.cc
)

add_subdirectory( predictors )

PREPEND( _p_obsbias_files       "obsbias"               ${obsbias_files} )
PREPEND( _p_predictor_files     "obsbias/predictors"    ${predictor_files} )

set ( obsbias_src_files
      ${_p_obsbias_files}
      ${_p_predictor_files}
      PARENT_SCOPE
)

