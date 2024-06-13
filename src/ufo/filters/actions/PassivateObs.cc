/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/PassivateObs.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<PassivateObs> makerPassivateObs_("passivate");

// -----------------------------------------------------------------------------

  PassivateObs::PassivateObs(const PassivateObsParameters &parameters)
    : allvars_() {
  }

// -----------------------------------------------------------------------------

void PassivateObs::apply(const Variables & vars,
                      const std::vector<std::vector<bool>> & flagged,
                      ObsFilterData &,
                      int,
                      ioda::ObsDataVector<int> & flags,
                      ioda::ObsDataVector<float> &) const {
  for (size_t ifiltervar = 0; ifiltervar < vars.nvars(); ++ifiltervar) {
    size_t iallvar = flags.varnames().find(vars.variable(ifiltervar).variable());
    for (size_t jobs = 0; jobs < flags.nlocs(); ++jobs) {
      if (flagged[ifiltervar][jobs] && flags[iallvar][jobs] == QCflags::pass)
        flags[iallvar][jobs] = QCflags::passive;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
