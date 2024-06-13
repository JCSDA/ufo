/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_REJECTOBS_H_
#define UFO_FILTERS_ACTIONS_REJECTOBS_H_

#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class RejectObsParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(RejectObsParameters, FilterActionParametersBase);

  // No extra parameters needed.
};

// -----------------------------------------------------------------------------

/// The default action of a QC filter: reject observations flagged by the filter.
class RejectObs : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef RejectObsParameters Parameters_;

  explicit RejectObs(const Parameters_ &);
  ~RejectObs() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return true; }

 private:
  Variables allvars_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_REJECTOBS_H_
