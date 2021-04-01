# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
if( ${ROPP-UFO_FOUND} )
set (  groundgnssropp_src_files 
       ObsGroundgnssROPP.h
       ObsGroundgnssROPP.cc
       ObsGroundgnssROPP.interface.h
       ObsGroundgnssROPP.interface.F90
       ufo_groundgnss_ropp_mod.F90
       ufo_groundgnss_ropp_utils_mod.F90
PARENT_SCOPE
)
else( ${ROPP-UFO_FOUND} )
set (  groundgnssropp_src_files 
       ObsGroundgnssROPP.h
       ObsGroundgnssROPP.cc
       ObsGroundgnssROPP.interface.h
       ObsGroundgnssROPP.interface.F90
       ufo_groundgnss_ropp_mod_stub.F90
PARENT_SCOPE
)
endif( ${ROPP-UFO_FOUND} )
