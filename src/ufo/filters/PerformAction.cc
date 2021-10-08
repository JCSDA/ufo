/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/PerformAction.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

PerformAction::PerformAction(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                             std::shared_ptr<ioda::ObsDataVector<int> > flags,
                             std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "PerformAction: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

void PerformAction::applyFilter(const std::vector<bool> & apply,
                                const Variables & filtervars,
                                std::vector<std::vector<bool>> & flagged) const {
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      flagged[jv][jobs] = apply[jobs];
    }
  }
}

// -----------------------------------------------------------------------------

void PerformAction::print(std::ostream & os) const {
  os << "PerformAction: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
