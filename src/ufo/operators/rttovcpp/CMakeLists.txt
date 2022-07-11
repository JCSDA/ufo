# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( rttovcpp_files
    ObsRadianceRTTOVCPP.h
    ObsRadianceRTTOVCPP.cc
    ObsRadianceRTTOVCPPParameters.h
    ObsRadianceRTTOVCPPTLAD.h
    ObsRadianceRTTOVCPPTLAD.cc
    rttovcpp_interface.h
    rttovcpp_interface.cc
)

PREPEND( _p_rttovcpp_files       "operators/rttovcpp"       ${rttovcpp_files} )

set ( rttovcpp_src_files
      ${_p_rttovcpp_files}
      PARENT_SCOPE
)

