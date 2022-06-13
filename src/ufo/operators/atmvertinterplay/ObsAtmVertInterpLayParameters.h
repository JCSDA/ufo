/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYPARAMETERS_H_
#define UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the AtmVertInterpLay operator.
class ObsAtmVertInterpLayParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAtmVertInterpLayParameters, ObsOperatorParametersBase)

 public:
  oops::OptionalParameter<std::vector<std::string>> geovals
    {"geovals",
     "list of geovals variable names to be included",
     this};

  oops::RequiredParameter<std::vector<double>> Coefficients
    {"coefficients",
     "coefficients for layer integral",
     this};

  oops::RequiredParameter<std::vector<std::size_t>> nlevels
    {"nlevels",
     "number of levels for each simulated variable",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMVERTINTERPLAY_OBSATMVERTINTERPLAYPARAMETERS_H_
