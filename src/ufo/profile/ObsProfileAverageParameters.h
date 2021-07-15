/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_OBSPROFILEAVERAGEPARAMETERS_H_
#define UFO_PROFILE_OBSPROFILEAVERAGEPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the average profile operator.
class ObsProfileAverageParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsProfileAverageParameters, Parameters)

 public:
  /// Operator name. In future will be moved to a base class for parameters of all ObsOperators.
  oops::OptionalParameter<std::string> name{"name", this};

  /// List of variables to be used by this operator.
  oops::Parameter<std::vector<ufo::Variable>> variables{"variables", {}, this};

  /// Name of model vertical coordinate.
  oops::RequiredParameter<std::string> modelVerticalCoordinate{"vertical coordinate", this};

  /// Number of iterations that are used to find the intersection between
  /// the observed profile and each model level.
  oops::Parameter<int> numIntersectionIterations
    {"number of intersection iterations", 3, this, {oops::minConstraint(1)}};

  /// Perform comparisons of auxiliary variables with OPS?
  oops::Parameter<bool> compareWithOPS{"compare with OPS", false, this};
};

}  // namespace ufo
#endif  // UFO_PROFILE_OBSPROFILEAVERAGEPARAMETERS_H_
