/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RTTOVCPP_OBSRADIANCERTTOVCPPPARAMETERS_H_
#define UFO_OPERATORS_RTTOVCPP_OBSRADIANCERTTOVCPPPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the ObsRadianceRTTOVCPP operator.
class ObsRadianceRTTOVCPPParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadianceRTTOVCPPParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<std::string> CoefPath
    {"CoefPath",
     "Path to optical depth coefficients file",
     this};

  oops::RequiredParameter<std::string> SensorID
    {"SensorID",
     "Name of optical depth coefficients file",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOVCPP_OBSRADIANCERTTOVCPPPARAMETERS_H_
