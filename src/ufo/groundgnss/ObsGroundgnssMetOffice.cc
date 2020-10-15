/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/groundgnss/ObsGroundgnssMetOffice.h"

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
static ObsOperatorMaker<ObsGroundgnssMetOffice> makerGroundgnssMetOffice_("GroundgnssMetOffice");
// -----------------------------------------------------------------------------

ObsGroundgnssMetOffice::ObsGroundgnssMetOffice(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperGroundgnssMetOffice_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "obs options");
  const eckit::Configuration *configc = &obsOptions;
  ufo_groundgnss_metoffice_setup_f90(keyOperGroundgnssMetOffice_, &configc);

  oops::Log::trace() << "ObsGroundgnssMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssMetOffice::~ObsGroundgnssMetOffice() {
  ufo_groundgnss_metoffice_delete_f90(keyOperGroundgnssMetOffice_);
  oops::Log::trace() << "ObsGroundgnssMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics &) const {
  ufo_groundgnss_metoffice_simobs_f90(keyOperGroundgnssMetOffice_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOffice::print(std::ostream & os) const {
  os << "ObsGroundgnssMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
