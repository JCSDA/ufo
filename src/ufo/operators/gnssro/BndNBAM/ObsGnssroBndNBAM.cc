/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/gnssro/BndNBAM/ObsGnssroBndNBAM.h"

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
static ObsOperatorMaker<ObsGnssroBndNBAM> makerGnssroBndNBAM_("GnssroBndNBAM");
// -----------------------------------------------------------------------------

ObsGnssroBndNBAM::ObsGnssroBndNBAM(const ioda::ObsSpace & odb, const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGnssroBndNBAM_(0), odb_(odb), varin_()
{
  std::vector<std::string> vv{"air_temperature", "specific_humidity",
                              "surface_altitude"};

  if ( params.options.value().vertLayer.value() == "mass" ) {
    vv.push_back("air_pressure");
    vv.push_back("geopotential_height");
  } else {
    vv.push_back("air_pressure_levels");
    vv.push_back("geopotential_height_levels");
  }

  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_bndnbam_setup_f90(keyOperGnssroBndNBAM_, params.options.value().toConfiguration());
  oops::Log::trace() << "ObsGnssroBndNBAM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndNBAM::~ObsGnssroBndNBAM() {
  ufo_gnssro_bndnbam_delete_f90(keyOperGnssroBndNBAM_);
  oops::Log::trace() << "ObsGnssroBndNBAM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                   ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndnbam_simobs_f90(keyOperGnssroBndNBAM_, gom.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndNBAM::print(std::ostream & os) const {
  os << "ObsGnssroBndNBAM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
