# (C) Copyright 2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( geosaod_files
    ObsGeosAod.h
    ObsGeosAod.cc
    ObsGeosAodTLAD.h
    ObsGeosAodTLAD.cc
    ObsGeosAod.interface.F90
    ObsGeosAod.interface.h
    ObsGeosAodTLAD.interface.F90
    ObsGeosAodTLAD.interface.h
    ufo_geosaod_mod.F90
    ufo_geosaod_tlad_mod.F90
)

PREPEND( _p_geosaod_files     "geos-aero"     ${geosaod_files} )

set ( geosaod_src_files
      ${_p_geosaod_files}
      PARENT_SCOPE
)
