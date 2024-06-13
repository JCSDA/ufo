/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class ReduceObsSpaceParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(ReduceObsSpaceParameters, FilterActionParametersBase);

  // No extra parameters needed.
};

// -----------------------------------------------------------------------------

/// The default action of a QC filter: reject observations flagged by the filter.
class ReduceObsSpace : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef ReduceObsSpaceParameters Parameters_;

  explicit ReduceObsSpace(const Parameters_ &);

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
