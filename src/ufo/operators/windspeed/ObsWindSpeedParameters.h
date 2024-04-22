/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_WINDSPEED_OBSWINDSPEEDPARAMETERS_H_
#define UFO_OPERATORS_WINDSPEED_OBSWINDSPEEDPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the Wind Speed operator.
class ObsWindSpeedParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsWindSpeedParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<std::string> model_surface_eastward_wind{
    "model_surface_eastward_wind", "Name of surface u in geovals", this};
  oops::RequiredParameter<std::string> model_surface_northward_wind{
    "model_surface_northward_wind", "Name of surface v in geovals", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_WINDSPEED_OBSWINDSPEEDPARAMETERS_H_
