/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_ACCEPTOBS_H_
#define UFO_FILTERS_ACTIONS_ACCEPTOBS_H_

#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class AcceptObsParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(AcceptObsParameters, FilterActionParametersBase);
  // No extra parameters needed
};

// -----------------------------------------------------------------------------

/// Reset the QC flag of observations flagged by the filter to 'pass' except for those whose current
/// QC flag is 'missing', 'preQC' or 'Hfailed'.
class AcceptObs : public FilterActionBase {
 public:
  typedef AcceptObsParameters Parameters_;

  explicit AcceptObs(const Parameters_ &);

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  void apply_to_record(const Variables &vars,
                       const std::vector<std::vector<bool>> &flagged,
                       ObsFilterData &data, int filterQCflag,
                       ioda::ObsDataVector<int> &flags,
                       ioda::ObsDataVector<float> &obserr) const override;

  const ufo::Variables & requiredVariables() const override {return allvars_;}

  bool modifiesQCFlags() const override { return true; }

 private:
    /// \brief Return true if the observation's current QC flag can be reset to pass.
    static bool canAcceptAtLocation(int currentFlag);

  Variables allvars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_ACCEPTOBS_H_
