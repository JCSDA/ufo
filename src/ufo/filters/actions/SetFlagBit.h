/*
 * (C) Copyright 2024, UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

/// Options taken by the `set flag bit` action.
class SetFlagBitParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(SetFlagBitParameters, FilterActionParametersBase);

 public:
  /// Bit of the diagnostic flag to set
  oops::RequiredParameter<int> bit{"bit", this, {oops::minConstraint(0), oops::maxConstraint(30)}};
};

// -----------------------------------------------------------------------------

/// Set a required diagnostic flag bit for all observations flagged by the filter.
class SetFlagBit : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef SetFlagBitParameters Parameters_;

  explicit SetFlagBit(const Parameters_ &);

  void apply(const Variables & vars,
             const std::vector<std::vector<bool>> & flagged,
             ObsFilterData & data,
             int /*filterQCflag*/,
             ioda::ObsDataVector<int> & flags,
             ioda::ObsDataVector<float> & obserr) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
  const int bitsetter_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
