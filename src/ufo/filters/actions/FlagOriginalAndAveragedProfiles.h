/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_FLAGORIGINALANDAVERAGEDPROFILES_H_
#define UFO_FILTERS_ACTIONS_FLAGORIGINALANDAVERAGEDPROFILES_H_

#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class FlagOriginalAndAveragedProfilesParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(FlagOriginalAndAveragedProfilesParameters, FilterActionParametersBase);

  // No extra parameters needed.
};

// -----------------------------------------------------------------------------

/// This action should only be used for data sets that satisfy two criteria:
/// firstly they have been grouped into records (profiles), and secondly they have an extended
/// section of the ObsSpace that consists of profiles that have been averaged onto model levels.
/// The action rejects any observations in the original profiles that have been flagged by the
/// filter. It also rejects all observations in any averaged profile whose corresponding original
/// profile contains at least one flagged observation.
/// This action should therefore be used for filters that are run only on the original profiles;
/// it enables the corresponding averaged profiles to be flagged without running the filter on them.
/// Doing this reduces the chance of configuration errors occurring. It may also be desirable if the
/// filter relies on properties of the original profiles that are not shared by the averaged
/// profiles or if the filter is expensive to run.
/// If rejecting observations in the averaged profile is not required, the standard Reject
/// action can be used instead.
///
/// NB the ObsSpace extension is produced with the following yaml options:
///
///    extension:
///      allocate companion records with length: N
///
/// where N is the desired number of levels per averaged profile.
/// This is different to the use of the DerivedObsValue group which adds an additional variable
/// to the existing ObsSpace.
class FlagOriginalAndAveragedProfiles : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef FlagOriginalAndAveragedProfilesParameters Parameters_;

  explicit FlagOriginalAndAveragedProfiles(const Parameters_ &);
  ~FlagOriginalAndAveragedProfiles() {}

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

#endif  // UFO_FILTERS_ACTIONS_FLAGORIGINALANDAVERAGEDPROFILES_H_
