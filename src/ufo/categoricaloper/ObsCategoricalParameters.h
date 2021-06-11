/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_CATEGORICALOPER_OBSCATEGORICALPARAMETERS_H_
#define UFO_CATEGORICALOPER_OBSCATEGORICALPARAMETERS_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {

/// Configuration options recognized by the Categorical operator.
class ObsCategoricalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsCategoricalParameters, Parameters)

 public:
  /// Operator name. In future will be moved to a base class for parameters of all ObsOperators.
  oops::OptionalParameter<std::string> name{"name", this};

  /// Categorical variable used to divide H(x) into sections.
  oops::RequiredParameter<std::string> categoricalVariable{"categorical variable", this};

  /// Name of the fallback observation operator to use. This will be used to produce H(x) at
  /// all locations whose value of the categorical variable does not appear in the
  /// \p categorisedOperatorNames map.
  oops::RequiredParameter<std::string> fallbackOperatorName{"fallback operator", this};

  /// Map between values of the categorical variable and the corresponding observation operators.
  /// The fallback observation operator will be used for all values of the categorical variable
  /// that are not represented in this map.
  oops::RequiredParameter<std::map<std::string, std::string>>
    categorisedOperatorNames{"categorised operators", this};

  /// A list of configuration options for each observation operator (i.e. the default operator
  /// and any operators that have been specified in the \p categorisedOperatorNames map).
  oops::RequiredParameter<std::vector<eckit::LocalConfiguration>>
    operatorConfigurations{"operator configurations", this};
};

}  // namespace ufo
#endif  // UFO_CATEGORICALOPER_OBSCATEGORICALPARAMETERS_H_
