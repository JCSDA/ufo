/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_REJECTOBS_H_
#define UFO_FILTERS_ACTIONS_REJECTOBS_H_

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/actions/FilterActionBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

class RejectObs : public FilterActionBase {
 public:
  explicit RejectObs(const eckit::Configuration &);
  ~RejectObs() {}

  void apply(const oops::Variables &, const std::vector<std::vector<bool>> &,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
 private:
  const int flag_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_REJECTOBS_H_
