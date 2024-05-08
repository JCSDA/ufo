# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( lightning_files
    ObsLightning.h
    ObsLightning.cc
    ObsLightningTLAD.h
    ObsLightningTLAD.cc
    ObsLightning.interface.F90
    ObsLightning.interface.h
    ObsLightningTLAD.interface.F90
    ObsLightningTLAD.interface.h
    ufo_lightning_mod.F90
    ufo_lightning_extend_geovals_mod.F90
    ufo_lightning_tlad_mod.F90
)

PREPEND( _p_lightning_files     "operators/lightning"     ${lightning_files} )

set ( lightning_src_files
      ${_p_lightning_files}
      PARENT_SCOPE
)
