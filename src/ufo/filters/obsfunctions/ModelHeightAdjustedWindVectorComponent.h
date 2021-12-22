/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */


#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDWINDVECTORCOMPONENT_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDWINDVECTORCOMPONENT_H_

#include <string>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ModelHeightAdjustedWindVector ObsFunction
class ModelHeightAdjustedWindVectorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelHeightAdjustedWindVectorParameters, Parameters)

 public:
  /// Input observation station height to be used
  oops::RequiredParameter<Variable> elevation{"elevation", this};
};

/// \brief Function to calculate surface wind observation values adjusted from
/// station height to model surface height. Outputs a derived eastward and northward wind adjusted
/// to the model surface height. Observed winds from stations above the model
/// surface are divided by a scaling factor (S) based only on height difference (dh).
/// Scaling factors are:
/// - S = 1.0 for dh < 100m
/// - S = 1 + 0.002*(dh - 100) for 100 < dh < 1100
/// - S = 3 for dh > 1100

template <bool northwardWind>  // true for the northward wind, false for the eastward wind
class ModelHeightAdjustedWindVectorComponent : public ObsFunctionBase<float> {
 public:
    explicit ModelHeightAdjustedWindVectorComponent(const eckit::LocalConfiguration &
                                                    = eckit::LocalConfiguration());

    void compute(const ObsFilterData &,
                 ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ModelHeightAdjustedWindVectorParameters parameters_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDWINDVECTORCOMPONENT_H_
