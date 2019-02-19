/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/gnssro/BndROPP2D/ObsGnssroBndROPP2D.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndROPP2D> makerGnssroBndROPP2D_("GnssroBndROPP2D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::ObsGnssroBndROPP2D(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : keyOperGnssroBndROPP2D_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vv{"temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"bending_angle"};
  varout_.reset(new oops::Variables(vout));

  const eckit::Configuration * configc = &config;
  ufo_gnssro_bndropp2d_setup_f90(keyOperGnssroBndROPP2D_, &configc);
  oops::Log::trace() << "ObsGnssroBndROPP2D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::~ObsGnssroBndROPP2D() {
  ufo_gnssro_bndropp2d_delete_f90(keyOperGnssroBndROPP2D_);
  oops::Log::trace() << "ObsGnssroBndROPP2D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_gnssro_bndropp2d_simobs_f90(keyOperGnssroBndROPP2D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP2D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
