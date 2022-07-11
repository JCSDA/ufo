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

void PerformActionParameters::deserialize(util::CompositePath &path,
                                          const eckit::Configuration &config) {
  FilterParametersBaseWithAbstractActions::deserialize(path, config);

  if (action_.value() != boost::none && actions_.value() != boost::none)
    throw eckit::UserError(path.path() +
                           ": The 'action' and 'actions' options are mutually exclusive");
  if (action_.value() == boost::none && actions_.value() == boost::none)
    throw eckit::UserError(path.path() +
                           ": Either the 'action' or 'actions' option must be set");
}

// -----------------------------------------------------------------------------

std::vector<std::unique_ptr<FilterActionParametersBase>> PerformActionParameters::actions() const {
  std::vector<std::unique_ptr<FilterActionParametersBase>> result;
  if (action_.value() != boost::none) {
    ASSERT(actions_.value() == boost::none);
    result.push_back(action_.value()->actionParameters.value().clone());
  } else {
    ASSERT(actions_.value() != boost::none);
    for (const FilterActionParameters &action : *actions_.value())
      result.push_back(action.actionParameters.value().clone());
  }
  return result;
}

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
