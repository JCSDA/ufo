/*
 * (C) Crown copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FilterParametersBase.h"

#include <memory>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

namespace ufo {

// -----------------------------------------------------------------------------

void FilterParametersBase::deserialize(util::CompositePath &path,
                                       const eckit::Configuration &config) {
  FilterParametersBaseWithAbstractActions::deserialize(path, config);

  if (action_.value() != boost::none && actions_.value() != boost::none)
    throw eckit::UserError(path.path() +
                           ": The 'action' and 'actions' options are mutually exclusive");
}

// -----------------------------------------------------------------------------

std::vector<std::unique_ptr<FilterActionParametersBase>> FilterParametersBase::actions() const {
  std::vector<std::unique_ptr<FilterActionParametersBase>> result;
  if (action_.value() != boost::none) {
    ASSERT(actions_.value() == boost::none);
    result.push_back(action_.value()->actionParameters.value().clone());
  } else if (actions_.value() != boost::none) {
    for (const FilterActionParameters &action : *actions_.value())
      result.push_back(action.actionParameters.value().clone());
  } else {
    std::unique_ptr<FilterActionParametersBase> action =
        FilterActionFactory::createParameters(defaultAction_);
    eckit::LocalConfiguration conf;
    conf.set("name", defaultAction_);
    action->deserialize(conf);
    result.push_back(std::move(action));
  }
  return result;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
