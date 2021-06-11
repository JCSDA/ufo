# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
if( ${ropp-ufo_FOUND} )
set (  refropp1d_src_files
       ObsGnssGBRefROPP1D.h
       ObsGnssGBRefROPP1D.cc
       ObsGnssGBRefROPP1D.interface.h
       ObsGnssGBRefROPP1D.interface.F90
       ufo_gnssgb_refropp1d_mod.F90
       ufo_gnssgb_ropp1d_utils_mod.F90
PARENT_SCOPE
)
else( ${ropp-ufo_FOUND} )
set (  refropp1d_src_files
       ObsGnssGBRefROPP1D.h
       ObsGnssGBRefROPP1D.cc
       ObsGnssGBRefROPP1D.interface.h
       ObsGnssGBRefROPP1D.interface.F90
       ufo_gnssgb_refropp1d_mod_stub.F90
PARENT_SCOPE
)
endif( ${ropp-ufo_FOUND} )
