/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/RejectObs.h"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<RejectObs> makerRejectObs_("reject");

// -----------------------------------------------------------------------------

RejectObs::RejectObs(const eckit::Configuration & conf): flag_(conf.getInt("flag")) {}

// -----------------------------------------------------------------------------

void RejectObs::apply(const oops::Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      const ObsFilterData &,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = flags.varnames().find(vars[jv]);
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[iv][jobs]) flags[iv][jobs] = flag_;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
