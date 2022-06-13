# 
# (C) Crown Copyright 2021 Met Office
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set (  sattcwv_files
       SatTCWV.h
       SatTCWV.cc
       SatTCWVTLAD.h
       SatTCWVTLAD.cc
)
PREPEND( _p_sattcwv_files    "operators/sattcwv"    ${sattcwv_files} )

set ( sattcwv_src_files
      ${_p_sattcwv_files}
      PARENT_SCOPE
)
