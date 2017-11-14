/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UFOTRAIT_H_
#define UFO_UFOTRAIT_H_

#include <string>


#include "GeoVaLs.h"
#include "Locations.h"
#include "ObsBias.h"
#include "ObsCheck.h"
#include "ObsSpace.h"
#include "ObsVector.h"
#include "Variables.h"

namespace ufo {

struct UfoTrait {
  static std::string name() {return "UFO";}

  typedef ufo::GeoVaLs             GeoVaLs;
  typedef ufo::Locations           Locations;
  typedef ufo::ObsCheck            ObsCheck;
  typedef ufo::ObsSpace            ObsSpace;
  typedef ufo::ObsVector           ObsVector;
  typedef ufo::Variables           Variables;
  typedef ufo::ObsBias             ObsAuxControl;
};

}  // namespace ufo

#endif  // UFO_UFOTRAIT_H_
