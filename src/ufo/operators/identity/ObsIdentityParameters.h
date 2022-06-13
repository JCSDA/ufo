/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_IDENTITY_OBSIDENTITYPARAMETERS_H_
#define UFO_OPERATORS_IDENTITY_OBSIDENTITYPARAMETERS_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsIdentityParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsIdentityParameters, ObsOperatorParametersBase)

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

  /// The boolean parameter `level index 0 is closest to surface` can be set to `true`
  /// for GeoVaLs whose level index 0 is closest to the Earth's surface.
  oops::Parameter<bool> levelIndex0IsClosestToSurface{
     "level index 0 is closest to surface",
     "GeoVaL level 0 is closest to the surface",
     false,
     this};
};

}  // namespace ufo

#endif  // UFO_OPERATORS_IDENTITY_OBSIDENTITYPARAMETERS_H_
