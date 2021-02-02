# (C) Copyright 2017-2018 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( calculate_files
      Cal_PressureFromHeight.h
      Cal_PressureFromHeight.cc
      CalculateBase.cc
      CalculateBase.h
      Formulas.cc
      Formulas.h
)

PREPEND( _p_calculate_files       "calculate"       ${calculate_files} )

set ( calculate_src_files
      ${_p_calculate_files}
      PARENT_SCOPE
)
