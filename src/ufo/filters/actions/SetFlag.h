/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_SETFLAG_H_
#define UFO_FILTERS_ACTIONS_SETFLAG_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

/// Options taken by the "set" or "unset" action.
class SetFlagParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(SetFlagParameters, FilterActionParametersBase);
 public:
  /// Name of the flag to act upon.
  oops::RequiredParameter<std::string> flag{"flag", this};
};

// -----------------------------------------------------------------------------

/// Set a diagnostic flag to \p value for all observations flagged by the filter.
template <bool value>
class SetFlag : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef SetFlagParameters Parameters_;

  explicit SetFlag(const Parameters_ &);

  void apply(const Variables & vars,
             const std::vector<std::vector<bool>> & flagged,
             const ObsFilterData & data,
             int /*filterQCflag*/,
             ioda::ObsDataVector<int> & flags,
             ioda::ObsDataVector<float> & obserr) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}

 private:
  Variables allvars_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_SETFLAG_H_
