# (C) Copyright 2017-2018 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( variabletransforms_files
      Cal_PressureFromHeight.h
      Cal_PressureFromHeight.cc
      Cal_Humidity.h
      Cal_Humidity.cc
      Cal_RemapScanPosition.h
      Cal_RemapScanPosition.cc
      Cal_Wind.h
      Cal_Wind.cc
      LookupTable.h
      TransformBase.cc
      TransformBase.h
      Formulas.cc
      Formulas.h
)

PREPEND( _p_variabletransforms_files       "variabletransforms"       ${variabletransforms_files} )

set ( variabletransforms_src_files
      ${_p_variabletransforms_files}
      PARENT_SCOPE
)
