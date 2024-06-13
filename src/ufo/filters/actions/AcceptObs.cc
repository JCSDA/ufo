/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/AcceptObs.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<AcceptObs> acceptObsMaker_("accept");

// -----------------------------------------------------------------------------

AcceptObs::AcceptObs(const Parameters_ &)
  : allvars_() {
}

// -----------------------------------------------------------------------------

void AcceptObs::apply(const Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      ObsFilterData &,
                      int /*filterQCflag*/,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    const size_t iallvar = flags.varnames().find(vars.variable(ifiltervar).variable());
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[ifiltervar][jobs]) {
        int &currentFlag = flags[iallvar][jobs];
        if (currentFlag != QCflags::missing &&
            currentFlag != QCflags::preQC &&
            currentFlag != QCflags::Hfailed)
          currentFlag = QCflags::pass;
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
