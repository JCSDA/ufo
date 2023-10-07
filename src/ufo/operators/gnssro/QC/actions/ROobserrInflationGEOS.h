/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_QC_ACTIONS_ROOBSERRINFLATIONGEOS_H_
#define UFO_OPERATORS_GNSSRO_QC_ACTIONS_ROOBSERRINFLATIONGEOS_H_

#include <string>
#include <vector>

#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class ROobserrInflationGEOSParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(ROobserrInflationGEOSParameters, FilterActionParametersBase);

  // No extra parameters needed
};

// -----------------------------------------------------------------------------

class ROobserrInflationGEOS : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef ROobserrInflationGEOSParameters Parameters_;

  explicit ROobserrInflationGEOS(const Parameters_ &);
  ~ROobserrInflationGEOS() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             const ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_QC_ACTIONS_ROOBSERRINFLATIONGEOS_H_
