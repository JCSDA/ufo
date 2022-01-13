/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/TemporalThinningParameters.h"

#include <utility>

namespace ufo {

void TemporalThinningParameters::deserialize(util::CompositePath &path,
                                             const eckit::Configuration &config) {
  oops::Parameters::deserialize(path, config);

  if (recordsAreSingleObs && categoryVariable.value() != boost::none) {
    throw eckit::UserError(path.path() + ": category_variable must be empty if "
                           "records_are_single_obs is set to true.", Here());
  }
}

}  // namespace ufo
