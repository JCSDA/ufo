/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERPPARAMETERS_H_
#define UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERPPARAMETERS_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options for the surface interpolation (AtmSfcInterp) observation operator.
class ObsAtmSfcInterpParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAtmSfcInterpParameters, ObsOperatorParametersBase)

 public:
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
    "variables",
    "List of variables to be used by this operator",
    this};

  oops::Parameter<bool> useFact10{
    "use_fact10",
    "If true use the wind reduction factor",
    false,
    this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMSFCINTERP_OBSATMSFCINTERPPARAMETERS_H_
