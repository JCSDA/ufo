/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSERRORFACTORY_H_
#define UFO_INSTANTIATEOBSERRORFACTORY_H_

#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObsErrorCovariance.h"

#include "ufo/errors/ObsErrorDiagonal.h"

namespace ufo {
template<typename OBS> void instantiateObsErrorFactory() {
  oops::instantiateObsErrorFactory<OBS>();
  static oops::ObsErrorMaker<OBS, oops::ObsErrorCovariance<OBS, ufo::ObsErrorDiagonal> >
              makerDiagUFO("diagonal ufo");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSERRORFACTORY_H_
