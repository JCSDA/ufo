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

class InflateError : public FilterActionBase {
 public:
  explicit InflateError(const eckit::Configuration &);
  ~InflateError() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             const ObsFilterData &,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
 private:
  Variables allvars_;
  const std::string strfactor_;
 protected:
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_INFLATEERROR_H_
