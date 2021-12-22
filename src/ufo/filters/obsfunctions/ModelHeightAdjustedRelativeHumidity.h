/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDRELATIVEHUMIDITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDRELATIVEHUMIDITY_H_

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ModelHeightAdjustedRelativeHumidity ObsFunction
class ModelHeightAdjustedRelativeHumidityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightAdjustedRelativeHumidityParameters, Parameters)

 public:
  /// Input observation station height to be used
  oops::RequiredParameter<Variable> elevation{"elevation", this};
  /// Temperature to be used
  oops::RequiredParameter<Variable> temperature{"temperature", this};
};

/// \brief  Function to calculate surface relative humidity observation value adjusted from
/// station height to model surface height. Outputs a derived 2m relative humidity adjusted
/// to the model surface height. Observed relative humidity is adjusted from station level
/// to model surface using an empirical vertical gradient LRH = -0.01 %/m. The adjusted
/// humidity value is then constrained to lie between zero and supersaturation with respect
/// to liquid water.

class ModelHeightAdjustedRelativeHumidity : public ObsFunctionBase<float> {
 public:
    explicit ModelHeightAdjustedRelativeHumidity(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());

    void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ModelHeightAdjustedRelativeHumidityParameters parameters_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDRELATIVEHUMIDITY_H_
