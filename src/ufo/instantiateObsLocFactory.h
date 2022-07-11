/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSLOCFACTORY_H_
#define UFO_INSTANTIATEOBSLOCFACTORY_H_

#include "oops/base/ObsLocalizationBase.h"
#include "ufo/obslocalization/ObsHorLocalization.h"
#include "ufo/obslocalization/ObsHorLocGC99.h"
#include "ufo/obslocalization/ObsHorLocSOAR.h"
#include "ufo/obslocalization/ObsVertLocalization.h"
#include "ufo/ObsTraits.h"

namespace ufo {
template<typename MODEL> void instantiateObsLocFactory() {
  static oops::ObsLocalizationMaker<MODEL, ObsTraits, ufo::ObsHorLocGC99<MODEL>>
           maker_("Horizontal Gaspari-Cohn");
  static oops::ObsLocalizationMaker<MODEL, ObsTraits, ufo::ObsHorLocSOAR<MODEL>>
           makerSOAR_("Horizontal SOAR");
  static oops::ObsLocalizationMaker<MODEL, ObsTraits, ufo::ObsHorLocalization<MODEL>>
           makerBoxCar_("Horizontal Box car");
  static oops::ObsLocalizationMaker<MODEL, ObsTraits, ufo::ObsVertLocalization<MODEL>>
           makerVertLoc_("Vertical localization");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSLOCFACTORY_H_
