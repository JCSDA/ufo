/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLEPARAMETERS_H_
#define TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLEPARAMETERS_H_

#include <string>

// TODO: modify the list of Parameter classes to include
// depending on what is used below.
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the Example operator.
class ObsExampleParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsExampleParameters, ObsOperatorParametersBase)

 public:
  // TODO: list all parameters here.
  oops::Parameter<int> myParameter
    {"my Parameter",
     "This is a Parameter with a default value of 4",
     4,
     this};

  oops::OptionalParameter<float> myOptionalParameter
    {"my OptionalParameter",
     "This is an OptionalParameter",
     this};

  oops::RequiredParameter<std::string> myRequiredParameter
    {"my RequiredParameter",
     "This is a RequiredParameter",
     this};

  oops::OptionalParameter<Variable> myVariableParameter
    {"my VariableParameter",
     "This is an OptionalParameter holding a Variable",
     this};
};

}  // namespace ufo
#endif  // TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLEPARAMETERS_H_
