/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICEPARAMETERS_H_
#define UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICEPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the GroundgnssMetOffice operator.
class ObsGroundgnssMetOfficeParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsGroundgnssMetOfficeParameters, ObsOperatorParametersBase)

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
     false,
     this};

  oops::Parameter<float> minTempGrad
    {"min_temp_grad",
     "The minimum temperature gradient permitted before a profile is considered "
     "isothermal. Used if pseudo-levels are used, otherwise ignored.",
     1.0e-6,
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICEPARAMETERS_H_
