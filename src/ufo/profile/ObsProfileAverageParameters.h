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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the average profile operator.
class ObsProfileAverageParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsProfileAverageParameters, ObsOperatorParametersBase)

 public:
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
    "variables",
    "List of variables to be used by this operator",
    this};

  oops::RequiredParameter<std::string> modelVerticalCoordinate{
    "model vertical coordinate",
    "Name of model vertical coordinate",
    this};

  oops::Parameter<int> numIntersectionIterations{
    "number of intersection iterations",
    "Number of iterations that are used to find the intersection between "
    "the observed profile and each model level",
     3,
     this,
     {oops::minConstraint(0)}};

  oops::Parameter<bool> compareWithOPS{
    "compare with OPS",
    "If true, perform comparisons of auxiliary variables with OPS",
    false,
    this};

  oops::Parameter<bool> requireDescendingPressureSort{
    "require descending pressure sort",
    "Require pressures to be sorted in descending order. "
    "Warning: setting this to `false` could lead to incorrect behaviour.",
    true,
    this};

  oops::Parameter<bool> geovalsObsSameDir{
    "geovals and observations in same direction",
    "Whether GeoVaLs and observations are oriented in the same vertical direction",
    false,
    this};

  oops::Parameter<std::string> pressureCoord{"pressure coordinate",
                                             "Name of air pressure coordinate",
                                             "pressure",
                                             this};

  oops::Parameter<std::string> pressureGroup{"pressure group",
                                             "Name of air pressure group",
                                             "ObsValue",
                                             this};
};

}  // namespace ufo
#endif  // UFO_PROFILE_OBSPROFILEAVERAGEPARAMETERS_H_
