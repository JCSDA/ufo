/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/parameters/Parameters.h"
#include "ufo/utils/parameters/ParameterBase.h"

namespace ufo {

void Parameters::registerChild(ParameterBase &parameter) {
  children_.push_back(&parameter);
}

void Parameters::deserialize(const eckit::Configuration &config) {
  for (auto child : children_) {
    child->deserialize(config);
  }
}

}  // namespace ufo
