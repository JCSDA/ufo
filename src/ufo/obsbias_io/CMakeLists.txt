# (C) Copyright 2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( obsbias_io_files
  ObsBiasIO.h
  GSI_Record.cc
  GSI_Record.h
  GSI_VarianceRecord.cc
  GSI_VarianceRecord.h
)

PREPEND( _p_obsbias_io_files     "obsbias_io"            ${obsbias_io_files} )

set ( obsbias_io_src_files
      ${_p_obsbias_io_files}
      PARENT_SCOPE
)
