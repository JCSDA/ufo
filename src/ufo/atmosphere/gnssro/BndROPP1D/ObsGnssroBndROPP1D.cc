/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/gnssro/BndROPP1D/ObsGnssroBndROPP1D.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndROPP1D> makerGnssroBndROPP1D_("GnssroBndROPP1D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP1D::ObsGnssroBndROPP1D(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : keyOperGnssroBndROPP1D_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vv{"temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "sfc_geopotential_height"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"zz"};
  varout_.reset(new oops::Variables(vout));

  const eckit::Configuration * configc = &config;
  ufo_gnssro_bndropp1d_setup_f90(keyOperGnssroBndROPP1D_, &configc);
  oops::Log::trace() << "ObsGnssroBndROPP1D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP1D::~ObsGnssroBndROPP1D() {
  ufo_gnssro_bndropp1d_delete_f90(keyOperGnssroBndROPP1D_);
  oops::Log::trace() << "ObsGnssroBndROPP1D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_gnssro_bndropp1d_simobs_f90(keyOperGnssroBndROPP1D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

Locations * ObsGnssroBndROPP1D::locateObs(const util::DateTime & t1,
                                       const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_gnssro_bndropp1d_locateobs_f90(keyOperGnssroBndROPP1D_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1D::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP1D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
