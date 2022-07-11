# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
if( ${ropp-ufo_FOUND} )
set (  bndropp1d_src_files
       ObsGnssroBndROPP1D.h
       ObsGnssroBndROPP1D.cc
       ObsGnssroBndROPP1D.interface.h
       ObsGnssroBndROPP1D.interface.F90
       ufo_gnssro_bndropp1d_mod.F90
       ObsGnssroBndROPP1DTLAD.h
       ObsGnssroBndROPP1DTLAD.cc
       ObsGnssroBndROPP1DTLAD.interface.h
       ObsGnssroBndROPP1DTLAD.interface.F90
       ufo_gnssro_bndropp1d_tlad_mod.F90
       ufo_gnssro_ropp1d_utils_mod.F90
PARENT_SCOPE
)
else( ${ropp-ufo_FOUND} )
set (  bndropp1d_src_files
       ObsGnssroBndROPP1D.h
       ObsGnssroBndROPP1D.cc
       ObsGnssroBndROPP1D.interface.h
       ObsGnssroBndROPP1D.interface.F90
       ufo_gnssro_bndropp1d_mod_stub.F90
       ObsGnssroBndROPP1DTLAD.h
       ObsGnssroBndROPP1DTLAD.cc
       ObsGnssroBndROPP1DTLAD.interface.h
       ObsGnssroBndROPP1DTLAD.interface.F90
       ufo_gnssro_bndropp1d_tlad_mod_stub.F90
PARENT_SCOPE
)
endif( ${ropp-ufo_FOUND} )
