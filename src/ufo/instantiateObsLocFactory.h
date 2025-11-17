/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSLOCFACTORY_H_
#define UFO_INSTANTIATEOBSLOCFACTORY_H_

#include "ufo/obslocalization/ObsHorLocalization.h"
#include "ufo/obslocalization/ObsHorLocGC99.h"
#include "ufo/obslocalization/ObsHorLocSOAR.h"
#include "ufo/obslocalization/ObsLocalizationBase.h"
#include "ufo/obslocalization/ObsVertLocalization.h"

namespace ufo {

template<typename ITERATOR> void instantiateObsLocFactory() {
  static ObsLocalizationMaker<ITERATOR, ObsHorLocGC99<ITERATOR>>
           maker_("Horizontal Gaspari-Cohn");
  static ObsLocalizationMaker<ITERATOR, ObsHorLocSOAR<ITERATOR>>
           makerSOAR_("Horizontal SOAR");
  static ObsLocalizationMaker<ITERATOR, ObsHorLocalization<ITERATOR>>
           makerBoxCar_("Horizontal Box car");
  static ObsLocalizationMaker<ITERATOR, ObsVertLocalization<ITERATOR>>
           makerVertLoc_("Vertical localization");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSLOCFACTORY_H_
