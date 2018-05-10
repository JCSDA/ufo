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
#include "ioda/Locations.h"
#include "ObsBias.h"
#include "ObsBiasIncrement.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ObsCheck.h"
#include "ObsBiasCovariance.h"

namespace ufo {

struct UfoTrait {
  static std::string name() {return "UFO";}

  typedef ufo::GeoVaLs             GeoVaLs;
  typedef ioda::Locations          Locations;
  typedef ioda::ObsSpace           ObsSpace;
  typedef ioda::ObsVector          ObsVector;

  typedef ufo::ObsBias             ObsAuxControl;
  typedef ufo::ObsBiasIncrement    ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance   ObsAuxCovariance;

  typedef ufo::ObsCheck            ObsCheck;
};

}  // namespace ufo

#endif  // UFO_UFOTRAIT_H_
