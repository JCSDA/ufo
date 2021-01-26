/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_INFLATEERROR_H_
#define UFO_FILTERS_ACTIONS_INFLATEERROR_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/actions/FilterActionBase.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------
/// \brief Observation error inflation action.
/// \details Inflates Observation error for filter variables by:
/// - constant (if "inflation factor" is specified in yaml)
/// - spatially varying filter data (if "inflation variable" is specified in yaml).
///   If inflation variable is the same size as filter variables, inflation is done
///   variable by variable (e.g. inflation variable 1 is used to inflate filter
///   variable 1; inflation variable 2 is used to inflate filter variable 2, etc).
///   If inflation variable is of size 1, the same inflation variable is used for
///   updating all filter variables.
class InflateError : public FilterActionBase {
 public:
  explicit InflateError(const eckit::Configuration &);

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             const ObsFilterData &,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;

  const ufo::Variables & requiredVariables() const override {return allvars_;}

 private:
  Variables allvars_;            /// variables required to compute inflation
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_INFLATEERROR_H_
