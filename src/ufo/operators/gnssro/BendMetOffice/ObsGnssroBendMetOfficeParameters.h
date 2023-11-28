/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_
#define UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the bending angle operator.
class ObsGnssroBendMetOfficeParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsGnssroBendMetOfficeParameters, ObsOperatorParametersBase)

 public:
  oops::Parameter<bool> vertInterpOPS
    {"vert_interp_ops",
     "If true assume that pressure varies exponentially with height when "
     "interpolating.  Otherwise assume that exner varies linearly with height, "
     "and derive pressure from this.",
     true,
     this};

  oops::Parameter<bool> pseudoLevels
    {"pseudo_ops",
     "Whether to use pseudo-levels in the calculation.",
     true,
     this};

  oops::Parameter<float> minTempGrad
    {"min_temp_grad",
     "The minimum temperature gradient permitted before a profile is considered "
     "isothermal. Used if pseudo-levels are used, otherwise ignored.",
     1.0e-6,
     this};

  oops::Parameter<bool> noSuperCheck
    {"no super-refraction check",
     "Whether to avoid using the super-refraction check in the operator",
     false,
     this};

  /// List of channels available for assimilation - this is used for vertical
  /// heights in the case where the observations are read in as profiles.
  oops::Parameter<std::string> channelList{"channels", "", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_
