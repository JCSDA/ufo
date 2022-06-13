/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITEPARAMETERS_H_
#define UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITEPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the Composite operator.
class ObsCompositeParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsCompositeParameters, ObsOperatorParametersBase)

 public:
  /// A list of configuration options for each operator used to simulate a subset of variables.
  oops::RequiredParameter<std::vector<eckit::LocalConfiguration>> components{"components", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITEPARAMETERS_H_
