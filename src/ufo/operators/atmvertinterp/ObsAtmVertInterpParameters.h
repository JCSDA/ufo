/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPPARAMETERS_H_
#define UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

enum class InterpolationMethod {
  AUTOMATIC, LINEAR, LOGLINEAR
};

struct InterpolationMethodParameterTraitsHelper {
  typedef InterpolationMethod EnumType;
  static constexpr char enumTypeName[] = "InterpolationMethod";
  static constexpr util::NamedEnumerator<InterpolationMethod> namedValues[] = {
    { InterpolationMethod::AUTOMATIC, "automatic" },
    { InterpolationMethod::LINEAR, "linear" },
    { InterpolationMethod::LOGLINEAR, "log-linear"}
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::InterpolationMethod> :
    public EnumParameterTraits<ufo::InterpolationMethodParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the AtmVertInterp operator.
class ObsAtmVertInterpParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAtmVertInterpParameters, ObsOperatorParametersBase)

 public:
  /// List of variables to be simulated: the default is to simulate all ObsSpace
  /// variables
  oops::OptionalParameter<std::vector<ufo::Variable>> variables
    {"variables",
     "List of variables to be simulated",
     this};

  oops::Parameter<std::string> VertCoord
    {"vertical coordinate",
     "vertical coordinate used by the model",
     "air_pressure",  // this should be consistent with var_prs defined in ufo_vars_mod
     this};

  oops::OptionalParameter<std::vector<double>> vertCoordValues
    {"constant vertical coordinate values",
     "constant vertical coordinate values to be used for interpolation for all observations",
     this};

  oops::Parameter<std::string> ObsVertCoord
    {"observation vertical coordinate",
     "vertical coordinate for the observations",
     "pressure",
     this};

  oops::OptionalParameter<std::string> ObsVertGroup
    {"observation vertical coordinate group",
     "observation vertical coordinate group",
     this};

  oops::Parameter<InterpolationMethod> interpMethod
    {"interpolation method",
     "interpolation method (options: automatic, linear, log-linear)",
     InterpolationMethod::AUTOMATIC,
     this};

  oops::OptionalParameter<bool> ApplyFact10
    {"apply near surface wind scaling",
     "apply near surface wind scaling",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPPARAMETERS_H_
