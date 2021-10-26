# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set (  bendmetoffice_src_files
       ObsGnssroBendMetOffice.h
       ObsGnssroBendMetOffice.cc
       ObsGnssroBendMetOffice.interface.h
       ObsGnssroBendMetOffice.interface.F90
       ObsGnssroBendMetOfficeParameters.h
       ObsGnssroBendMetOfficeTLAD.h
       ObsGnssroBendMetOfficeTLAD.cc
       ObsGnssroBendMetOfficeTLAD.interface.h
       ObsGnssroBendMetOfficeTLAD.interface.F90
       ufo_gnssro_bendmetoffice_mod.F90
       ufo_gnssro_bendmetoffice_utils_mod.F90
       ufo_gnssro_bendmetoffice_tlad_mod.F90
       ufo_gnssro_bendmetoffice_tlad_utils_mod.F90
PARENT_SCOPE
)
