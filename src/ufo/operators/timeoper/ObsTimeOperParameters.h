/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_TIMEOPER_OBSTIMEOPERPARAMETERS_H_
#define UFO_OPERATORS_TIMEOPER_OBSTIMEOPERPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the TimeOper operator.
class ObsTimeOperParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsTimeOperParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<util::Duration> windowSub
    {"windowSub",
     "Duration of sub-window.",
     this};

  oops::RequiredParameter<eckit::LocalConfiguration> obsOperator
    {"obs operator",
     "Configuration options for the observation operator whose output will be time-interpolated.",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_TIMEOPER_OBSTIMEOPERPARAMETERS_H_
