/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/RejectObs.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<RejectObs> makerRejectObs_("reject");

// -----------------------------------------------------------------------------

RejectObs::RejectObs(const eckit::Configuration & conf)
  : allvars_(), conf_(conf) {
}

// -----------------------------------------------------------------------------

void RejectObs::apply(const Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      const ObsFilterData &,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  int flag = conf_.getInt("flag");
  for (size_t jv = 0; jv < vars.nvars(); ++jv) {
    size_t iv = flags.varnames().find(vars.variable(jv).variable());
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[jv][jobs] && flags[iv][jobs] == QCflags::pass) flags[iv][jobs] = flag;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
