# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set (  qc_src_files
       ROobserror.h
       ROobserror.cc
       ROobserror.interface.F90
       ROobserror.interface.h
       ufo_roobserror_mod.F90
       ufo_roobserror_utils_mod.F90
       BackgroundCheckRONBAM.h
       BackgroundCheckRONBAM.cc
       actions/ROobserrInflation.cc
       actions/ROobserrInflation.h
       actions/ROobserrInflationGEOS.cc
       actions/ROobserrInflationGEOS.h
PARENT_SCOPE
)

