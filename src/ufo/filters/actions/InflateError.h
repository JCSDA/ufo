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
#include "oops/base/Variables.h"
#include "ufo/filters/actions/FilterActionBase.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class InflateError : public FilterActionBase {
 public:
  explicit InflateError(const eckit::Configuration &);
  ~InflateError() {}

  void apply(const oops::Variables &, const std::vector<std::vector<bool>> &,
             const ObsFilterData &,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
 private:
  const std::string strfactor_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_INFLATEERROR_H_
