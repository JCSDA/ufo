/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/gnssro/BendMetOffice/ObsGnssroBendMetOffice.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBendMetOffice> makerGnssroBendMetOffice_("GnssroBendMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroBendMetOffice::ObsGnssroBendMetOffice(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperGnssroBendMetOffice_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "obs options");
  const eckit::Configuration *configc = &obsOptions;
  ufo_gnssro_bendmetoffice_setup_f90(keyOperGnssroBendMetOffice_, &configc);

  oops::Log::trace() << "ObsGnssroBendMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBendMetOffice::~ObsGnssroBendMetOffice() {
  ufo_gnssro_bendmetoffice_delete_f90(keyOperGnssroBendMetOffice_);
  oops::Log::trace() << "ObsGnssroBendMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "Starting simulateObs" << std::endl;
  ufo_gnssro_bendmetoffice_simobs_f90(keyOperGnssroBendMetOffice_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran(), ydiags.toFortran());
  oops::Log::trace() << "Finishing simulateObs" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOffice::print(std::ostream & os) const {
  os << "ObsGnssroBendMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
