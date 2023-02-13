/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_OASIM_OBSRADIANCEOASIMPARAMETERS_H_
#define UFO_OPERATORS_OASIM_OBSRADIANCEOASIMPARAMETERS_H_

#include <string>
#include <vector>
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// oasim observation operator Parameters class
class ObsRadianceOASIMParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadianceOASIMParameters, ObsOperatorParametersBase)
 public:
  /// CoefficientPath
  oops::RequiredParameter<std::string> CoefficientPath{"CoefficientPath", this};
};
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_OASIM_OBSRADIANCEOASIMPARAMETERS_H_
