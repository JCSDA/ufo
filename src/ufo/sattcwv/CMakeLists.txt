# 
# (C) Crown Copyright 2021 Met Office
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set (  sattcwv_files
       SatTCWV.h
       SatTCWV.cc
       SatTCWV.interface.h
       SatTCWV.interface.F90
       ufo_SatTCWV_mod.F90
       SatTCWVTLAD.h
       SatTCWVTLAD.cc
       SatTCWVTLAD.interface.h
       SatTCWVTLAD.interface.F90
       ufo_SatTCWV_tlad_mod.F90
)
PREPEND( _p_sattcwv_files    "sattcwv"    ${sattcwv_files} )

set ( sattcwv_src_files
      ${_p_sattcwv_files}
      PARENT_SCOPE
)
