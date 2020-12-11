/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/gnssro/BndROPP2D/ObsGnssroBndROPP2D.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndROPP2D> makerGnssroBndROPP2D_("GnssroBndROPP2D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::ObsGnssroBndROPP2D(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperGnssroBndROPP2D_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "surface_altitude"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "obs options");

  ufo_gnssro_bndropp2d_setup_f90(keyOperGnssroBndROPP2D_, obsOptions, odb_.nlocs());
  oops::Log::trace() << "ObsGnssroBndROPP2D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::~ObsGnssroBndROPP2D() {
  ufo_gnssro_bndropp2d_delete_f90(keyOperGnssroBndROPP2D_);
  oops::Log::trace() << "ObsGnssroBndROPP2D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics &) const {
  ufo_gnssro_bndropp2d_simobs_f90(keyOperGnssroBndROPP2D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
std::unique_ptr<Locations> ObsGnssroBndROPP2D::locations() const {
  std::unique_ptr<Locations> locs(new Locations(odb_.comm()));

  int keylocs = locs->toFortran();

  ufo_gnssro_2d_locs_init_f90(keyOperGnssroBndROPP2D_, keylocs, odb_);

  return locs;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP2D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
