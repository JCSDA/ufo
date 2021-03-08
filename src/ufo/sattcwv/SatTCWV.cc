/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/sattcwv/SatTCWV.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// ----------------------------------------------------------------------------
static ObsOperatorMaker<SatTCWV> makerSatTCWV_("SatTCWV");
// -----------------------------------------------------------------------------

SatTCWV::SatTCWV(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperSatTCWV_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "surface_pressure"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "obs options");
  const eckit::Configuration *configc = &obsOptions;
  ufo_sattcwv_setup_f90(keyOperSatTCWV_, &configc);

  oops::Log::trace() << "SatTCWV created." << std::endl;
}

// -----------------------------------------------------------------------------

SatTCWV::~SatTCWV() {
  ufo_sattcwv_delete_f90(keyOperSatTCWV_);
  oops::Log::trace() << "SatTCWV destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SatTCWV::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics &) const {
  ufo_sattcwv_simobs_f90(keyOperSatTCWV_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void SatTCWV::print(std::ostream & os) const {
  os << "SatTCWV::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
