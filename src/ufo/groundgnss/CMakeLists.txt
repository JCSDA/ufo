# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set (  groundgnss_files
       ObsGroundgnssMetOffice.h
       ObsGroundgnssMetOffice.cc
       ObsGroundgnssMetOffice.interface.h
       ObsGroundgnssMetOffice.interface.F90
       ufo_groundgnss_metoffice_mod.F90
       ufo_groundgnss_metoffice_utils_mod.F90
       ufo_refractivity_utils_mod.F90
)

PREPEND( _p_groundgnss_files    "groundgnss"    ${groundgnss_files} )

set ( groundgnss_src_files
      ${_p_groundgnss_files}
      PARENT_SCOPE
)
