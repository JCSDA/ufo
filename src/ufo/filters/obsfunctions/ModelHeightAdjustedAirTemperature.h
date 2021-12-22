/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_

#include <string>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ModelHeightAdjustedAirTemperature ObsFunction
class ModelHeightAdjustedAirTemperatureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightAdjustedAirTemperatureParameters, Parameters)

 public:
  /// Input observation station height to be used
  oops::RequiredParameter<Variable> elevation{"elevation", this};
};

/// \brief  Function to calculate surface temperature observation value
/// adjusted from station height to model surface height. Outputs a derived
/// 2m air temperature adjusted to the model surface height. The correction
/// applied to the temperature is calculated using a standard lapse rate
/// (Constants::Lclr).

class ModelHeightAdjustedAirTemperature : public ObsFunctionBase<float> {
 public:
    explicit ModelHeightAdjustedAirTemperature(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ModelHeightAdjustedAirTemperatureParameters parameters_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDAIRTEMPERATURE_H_
