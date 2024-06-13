/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_PASSIVATEOBS_H_
#define UFO_FILTERS_ACTIONS_PASSIVATEOBS_H_

#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class PassivateObsParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(PassivateObsParameters, FilterActionParametersBase);
};

// -----------------------------------------------------------------------------

/// Flag observations as passive.
class PassivateObs : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef PassivateObsParameters Parameters_;

  explicit PassivateObs(const Parameters_ &);
  ~PassivateObs() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return true; }

 private:
  Variables allvars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_PASSIVATEOBS_H_
