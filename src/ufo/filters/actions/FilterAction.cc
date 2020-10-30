/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/FilterAction.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterAction::FilterAction(const eckit::Configuration & conf)
  : action_(FilterActionFactory::create(conf))
{}

// -----------------------------------------------------------------------------

FilterAction::~FilterAction() {}

// -----------------------------------------------------------------------------

void FilterAction::apply(const Variables & vars, const std::vector<std::vector<bool>> & mask,
                         const ObsFilterData & data,
                         ioda::ObsDataVector<int> & flag, ioda::ObsDataVector<float> & err) const {
  action_->apply(vars, mask, data, flag, err);
}

// -----------------------------------------------------------------------------

const ufo::Variables & FilterAction::requiredVariables() const {
  return action_->requiredVariables();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
