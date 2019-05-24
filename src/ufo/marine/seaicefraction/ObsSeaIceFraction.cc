/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/seaicefraction/ObsSeaIceFraction.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFraction::ObsSeaIceFraction(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vvin{"sea_ice_category_area_fraction"};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{"sea_ice_area_fraction"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_seaicefraction_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsSeaIceFraction created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFraction::~ObsSeaIceFraction() {
  ufo_seaicefraction_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_seaicefraction_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceFraction: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
