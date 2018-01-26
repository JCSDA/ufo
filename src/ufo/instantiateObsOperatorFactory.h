/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/interface/ObsOperatorBase.h"
#include "ObsSeaIceFraction.h"
//#include "ObsRadiance.h"
//#include "ObsRadiosonde.h"

namespace ufo {

template<typename MODEL> void instantiateObsOperatorFactory() {
  static oops::ObsOperatorMaker<MODEL, ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
//  static oops::ObsOperatorMaker<MODEL, ObsRadiance>       makerRadiance_("Radiance");
//  static oops::ObsOperatorMaker<MODEL, ObsRadiosonde>     makerRadiosonde_("Radiosonde");
}

}
