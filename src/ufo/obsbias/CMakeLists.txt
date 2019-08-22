# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( obsbias_files
  ObsBiasBase.h
  ObsBiasBase.cc
  ObsBiasRadianceGSI.h
  ObsBiasRadianceGSI.cc
)

PREPEND( _p_obsbias_files       "obsbias"       ${obsbias_files} )

set ( obsbias_src_files
      ${_p_obsbias_files}
      PARENT_SCOPE
)

