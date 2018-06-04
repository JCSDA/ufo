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
#include "marine/insitutemperature/ObsInsituTemperature.h"
#include "marine/adt/ObsADT.h"
#include "atmosphere/radiance/ObsRadiance.h"
#include "atmosphere/radiosonde/ObsRadiosonde.h"
#include "constituents/aod/ObsAod.h"

namespace ufo {

template<typename MODEL> void instantiateObsOperatorFactory() {
  static oops::ObsOperatorMaker<MODEL, ObsStericHeight<MODEL>> makerStericHeight_("StericHeight");
  static oops::ObsOperatorMaker<MODEL, ObsSeaIceFraction<MODEL>> makerSeaIceFraction_("SeaIceFraction");  
  static oops::ObsOperatorMaker<MODEL, ObsSeaIceThickness<MODEL>> makerSeaIceThickness_("SeaIceThickness");
  static oops::ObsOperatorMaker<MODEL, ObsInsituTemperature<MODEL>> makerInsituTemperature_("InsituTemperature");
  static oops::ObsOperatorMaker<MODEL, ObsADT<MODEL>> makerADT_("ADT");  
  static oops::ObsOperatorMaker<MODEL, ObsRadiance<MODEL>>       makerRadiance_("Radiance");
  static oops::ObsOperatorMaker<MODEL, ObsRadiosonde<MODEL>>     makerRadiosonde_("Radiosonde");
  static oops::ObsOperatorMaker<MODEL, ObsAod<MODEL>>            makerAod_("Aod");
}

}
