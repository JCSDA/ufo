/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSERRORFACTORY_H_
#define UFO_INSTANTIATEOBSERRORFACTORY_H_

#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObsErrorBase.h"

#include "ufo/errors/ObsErrorCrossVarCov.h"
#include "ufo/errors/ObsErrorDiagonal.h"
#include "ufo/errors/ObsErrorWithinGroupCov.h"
#include "ufo/ObsTraits.h"

namespace ufo {
void instantiateObsErrorFactory() {
  oops::instantiateObsErrorFactory<ObsTraits>();
  static oops::interface::ObsErrorMaker<ObsTraits, ObsErrorDiagonal>
              makerDiagUFO("diagonal ufo");
  static oops::interface::ObsErrorMaker<ObsTraits, ObsErrorCrossVarCov>
              makerCrossVarCov("cross variable covariances");
  static oops::interface::ObsErrorMaker<ObsTraits, ObsErrorWithinGroupCov>
              makerWithinGroupCov("within group covariances");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSERRORFACTORY_H_
