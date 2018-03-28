/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/interface/ObsOperatorBase.h"
#include "marine/seaicefraction/ObsSeaIceFraction.h"
#include "marine/seaicethickness/ObsSeaIceThickness.h"
#include "marine/stericheight/ObsStericHeight.h"
#include "ObsRadiance.h"
#include "ObsRadiosonde.h"

namespace ufo {

template<typename MODEL> void instantiateObsOperatorFactory() {
  static oops::ObsOperatorMaker<MODEL, ObsStericHeight<MODEL>> makerStericHeight_("StericHeight");
  static oops::ObsOperatorMaker<MODEL, ObsSeaIceFraction<MODEL>> makerSeaIceFraction_("SeaIceFraction");  
  static oops::ObsOperatorMaker<MODEL, ObsSeaIceThickness<MODEL>> makerSeaIceThickness_("SeaIceThickness");
  static oops::ObsOperatorMaker<MODEL, ObsRadiance<MODEL>>       makerRadiance_("Radiance");
  static oops::ObsOperatorMaker<MODEL, ObsRadiosonde<MODEL>>     makerRadiosonde_("Radiosonde");
}

}
