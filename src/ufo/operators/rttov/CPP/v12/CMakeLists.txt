# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( rttovcpp_files
    ObsRadianceRTTOVCPP.h
    v12/ObsRadianceRTTOVCPP.cc
    ObsRadianceRTTOVCPPParameters.h
    ObsRadianceRTTOVCPPTLAD.h
    v12/ObsRadianceRTTOVCPPTLAD.cc
    v12/rttovcpp_interface.h
    v12/rttovcpp_interface.cc
)

PREPEND( _p_rttovcpp_files       "operators/rttov/CPP"       ${rttovcpp_files} )

set ( rttovcpp_src_files
      ${_p_rttovcpp_files}
      PARENT_SCOPE
)

