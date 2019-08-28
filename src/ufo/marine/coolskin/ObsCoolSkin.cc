/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/coolskin/ObsCoolSkin.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsCoolSkin> makerCoolSkin_("CoolSkin");
// -----------------------------------------------------------------------------

ObsCoolSkin::ObsCoolSkin(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vvin{"sea_surface_temperature",
                                      "net_downwelling_shortwave_radiation",
                                      "upward_latent_heat_flux_in_air",
                                      "upward_sensible_heat_flux_in_air",
                                      "net_downwelling_longwave_radiation",
                                      "friction_velocity_over_water"};
  varin_.reset(new oops::Variables(vvin));

  const eckit::Configuration * configc = &config;
  ufo_CoolSkin_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsCoolSkin created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsCoolSkin::~ObsCoolSkin() {
  ufo_CoolSkin_delete_f90(keyOper_);
  oops::Log::trace() << "ObsCoolSkin destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkin::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics &) const {
  ufo_CoolSkin_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsCoolSkin: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCoolSkin::print(std::ostream & os) const {
  os << "ObsCoolSkin::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
