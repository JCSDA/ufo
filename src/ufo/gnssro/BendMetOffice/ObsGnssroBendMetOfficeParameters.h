/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_
#define UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {

/// Configuration options recognized by the bending angle operator.
class ObsGnssroBendMetOfficeOptions : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGnssroBendMetOfficeOptions, Parameters)

 public:
  /// If true assume that pressure varies exponentially with height when
  /// interpolating.  Otherwise assume that exner varies linearly with height,
  /// and derive pressure from this.
  oops::Parameter<bool> vertInterpOPS{"vert_interp_ops", true, this};
  /// Whether to use pseudo-levels in the calculation.
  oops::Parameter<bool> pseudoLevels{"pseudo_ops", true, this};
  /// The minimum temperature gradient permitted before a profile is considered
  /// isothermal.  Only used if pseudo-levels are also used.
  oops::Parameter<float> minTempGrad{"min_temp_grad", 1.0e-6, this};
};

/// Configuration options recognized by the bending angle operator.
class ObsGnssroBendMetOfficeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGnssroBendMetOfficeParameters, Parameters)

 public:
  /// Operator name. In future will be moved to a base class for parameters of all ObsOperators.
  oops::OptionalParameter<std::string> name{"name", this};
  /// Obs Options - settings for the observation operator
  oops::Parameter<ObsGnssroBendMetOfficeOptions> obsOptions{"obs options",
    ObsGnssroBendMetOfficeOptions(), this};
};

}  // namespace ufo
#endif  // UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICEPARAMETERS_H_
