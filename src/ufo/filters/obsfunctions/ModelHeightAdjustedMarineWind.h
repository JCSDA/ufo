/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDMARINEWIND_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDMARINEWIND_H_

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief  Function to calculate surface marine wind value adjusted from
/// anemometer height to reference height of 10m. Outputs a derived eastward
/// or northward wind adjusted to the reference height. Observed winds are
/// adjusted to the reference height by multiplying by a scaling factor (S)
/// based on a simple stability-independent adjustment based on one of the
/// TurboWin options (Eqn 1 of Thomas et al, 2005; doi.org/10.1002/joc.1176).
/// Observations where the anemometer height is 0m are NOT corrected,
/// instead a missing data indicator is assigned.

class ModelHeightAdjustedMarineWindComponent : public ObsFunctionBase<float> {
 protected:
  ModelHeightAdjustedMarineWindComponent(const eckit::LocalConfiguration &conf,
                                         const Variable & windComponent);
 public:
    void compute(const ObsFilterData &,
                 ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    ufo::Variables invars_;
    ufo::Variable wind_;
};

class ModelHeightAdjustedEastwardMarineWind : public ModelHeightAdjustedMarineWindComponent {
 public:
    explicit ModelHeightAdjustedEastwardMarineWind(const eckit::LocalConfiguration &conf)
    : ModelHeightAdjustedMarineWindComponent(conf, Variable("ObsValue/windEastwardAt10M"))
    {}
};

class ModelHeightAdjustedNorthwardMarineWind : public ModelHeightAdjustedMarineWindComponent {
 public:
    explicit ModelHeightAdjustedNorthwardMarineWind(const eckit::LocalConfiguration &conf)
    : ModelHeightAdjustedMarineWindComponent(conf, Variable("ObsValue/windNorthwardAt10M"))
    {}
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELHEIGHTADJUSTEDMARINEWIND_H_
