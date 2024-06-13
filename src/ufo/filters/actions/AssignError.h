/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_ASSIGNERROR_H_
#define UFO_FILTERS_ACTIONS_ASSIGNERROR_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class AssignErrorParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(AssignErrorParameters, FilterActionParametersBase);

 public:
  oops::OptionalParameter<float> errorParameter{"error parameter", this};
  oops::OptionalParameter<std::vector<float>> errorParameterVector{"error parameter vector", this};
  oops::OptionalParameter<Variable> errorFunction{"error function", this};

  /// This function is overridden to check that either `error parameter` or `error function`
  /// is specified, but not both.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;
};

// -----------------------------------------------------------------------------

class AssignError : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef AssignErrorParameters Parameters_;

  explicit AssignError(const Parameters_ &);
  ~AssignError() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
  const Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_ASSIGNERROR_H_
