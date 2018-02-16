/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/interface/LinearObsOperBase.h"
#include "ObsSeaIceFractionTLAD.h"

namespace ufo {

template<typename MODEL> void instantiateLinearObsOpFactory() {
  static oops::LinearObsOpMaker<MODEL, ObsSeaIceFractionTLAD<MODEL>> makerFractionTL_("SeaIceFraction");
}

}
