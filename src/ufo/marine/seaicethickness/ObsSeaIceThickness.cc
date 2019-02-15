/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/seaicethickness/ObsSeaIceThickness.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceThickness> makerSeaIceThickness_("SeaIceThickness");
// -----------------------------------------------------------------------------

ObsSeaIceThickness::ObsSeaIceThickness(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vvin{"ice_concentration", "ice_thickness"};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{"obs_sea_ice_thickness"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_seaicethickness_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsSeaIceThickness created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceThickness::~ObsSeaIceThickness() {
  ufo_seaicethickness_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceThickness destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThickness::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_seaicethickness_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsSeaIceThickness: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThickness::print(std::ostream & os) const {
  os << "ObsSeaIceThickness::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
