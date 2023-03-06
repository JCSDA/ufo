/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsADTParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsADTParameters, ObsOperatorParametersBase )
 public:
  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this operator is used as a
  /// component of the Composite operator. If `variables` is not set, the operator will simulate
  /// all ObsSpace variables. Please see the documentation of the Composite operator for further
  /// details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
     "variables",
     "List of variables to be simulated",
     this};
};

}  // namespace ufo
