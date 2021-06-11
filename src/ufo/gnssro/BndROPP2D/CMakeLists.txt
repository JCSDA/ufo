# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
if( ${ropp-ufo_FOUND} )
set (  bndropp2d_src_files
       ObsGnssroBndROPP2D.h
       ObsGnssroBndROPP2D.cc
       ObsGnssroBndROPP2D.interface.h
       ObsGnssroBndROPP2D.interface.F90
       ufo_gnssro_bndropp2d_mod.F90
       ufo_gnssro_ropp2d_utils_mod.F90
       ObsGnssroBndROPP2DTLAD.h
       ObsGnssroBndROPP2DTLAD.cc
       ObsGnssroBndROPP2DTLAD.interface.h
       ObsGnssroBndROPP2DTLAD.interface.F90
       ufo_gnssro_bndropp2d_tlad_mod.F90
       ufo_gnssro_2d_locs_mod.F90
PARENT_SCOPE
)
else( ${ropp-ufo_FOUND} )
set (  bndropp2d_src_files
       ObsGnssroBndROPP2D.h
       ObsGnssroBndROPP2D.cc
       ObsGnssroBndROPP2D.interface.h
       ObsGnssroBndROPP2D.interface.F90
       ufo_gnssro_bndropp2d_mod_stub.F90
       ObsGnssroBndROPP2DTLAD.h
       ObsGnssroBndROPP2DTLAD.cc
       ObsGnssroBndROPP2DTLAD.interface.h
       ObsGnssroBndROPP2DTLAD.interface.F90
       ufo_gnssro_bndropp2d_tlad_mod_stub.F90
       ufo_gnssro_2d_locs_mod_stub.F90
PARENT_SCOPE
)
endif( ${ropp-ufo_FOUND} )
