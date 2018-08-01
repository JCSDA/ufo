/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsGnssro.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroRef> makerGnssroRef_("GnssroRef");
// -----------------------------------------------------------------------------

ObsGnssroRef::ObsGnssroRef(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperGnssroRef_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"temperature", "humidity_mixing_ratio", "air_pressure","geopotential_height"};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  ufo_gnssro_setup_f90(keyOperGnssroRef_, &configc);
  oops::Log::trace() << "ObsGnssroRef created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRef::~ObsGnssroRef() {
  ufo_gnssro_delete_f90(keyOperGnssroRef_);
  oops::Log::trace() << "ObsGnssroRef destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRef::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_gnssro_ref_f90(keyOperGnssroRef_, gom.toFortran(), odb_.toFortran(),
                           ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRef::print(std::ostream & os) const {
  os << "ObsGnssroRef::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
