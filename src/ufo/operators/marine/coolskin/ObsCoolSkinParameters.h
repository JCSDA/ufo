/*
 * (C) Copyright 2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKINPARAMETERS_H_
#define UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKINPARAMETERS_H_

#include "ufo/ObsOperatorParametersBase.h"


namespace ufo {

// -----------------------------------------------------------------------------
/// coolskin observation operator Parameters class
class ObsCoolSkinParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsCoolSkinParameters, ObsOperatorParametersBase)
  //
};
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKINPARAMETERS_H_
