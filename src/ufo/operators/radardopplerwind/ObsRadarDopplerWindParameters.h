/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDPARAMETERS_H_
#define UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDPARAMETERS_H_

#include <string>

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the RadarDopplerWind operator.
class ObsRadarDopplerWindParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadarDopplerWindParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<std::string>
    verticalCoordinate_uv{"vertical coordinate for horizontal wind",
      "Name of model vertical coordinate for eastward and northward wind.",
      this};

  oops::RequiredParameter<std::string>
    verticalCoordinate_w{"vertical coordinate for vertical wind",
      "Name of model vertical coordinate for vertical wind.",
      this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDPARAMETERS_H_
