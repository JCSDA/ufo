/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/gnssro/Ref/ObsGnssroRef.h"

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
static ObsOperatorMaker<ObsGnssroRef> makerGnssroRef_("GnssroRef");
// -----------------------------------------------------------------------------

ObsGnssroRef::ObsGnssroRef(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperGnssroRef_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");

  ufo_gnssro_ref_setup_f90(keyOperGnssroRef_, obsOptions);
  oops::Log::trace() << "ObsGnssroRef created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRef::~ObsGnssroRef() {
  ufo_gnssro_ref_delete_f90(keyOperGnssroRef_);
  oops::Log::trace() << "ObsGnssroRef destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRef::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                               ObsDiagnostics &) const {
  ufo_gnssro_ref_simobs_f90(keyOperGnssroRef_, gom.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRef::print(std::ostream & os) const {
  os << "ObsGnssroRef::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
