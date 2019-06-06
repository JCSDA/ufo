/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/insitutemperature/ObsInsituTemperature.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsInsituTemperature> makerInsituTemperature_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperature::ObsInsituTemperature(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vvin{"sea_water_potential_temperature",
                                      "sea_water_salinity",
                                      "sea_water_cell_thickness"};
  varin_.reset(new oops::Variables(vvin));
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsInsituTemperature created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperature::~ObsInsituTemperature() {
  ufo_insitutemperature_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituTemperature destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_insitutemperature_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperature: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::print(std::ostream & os) const {
  os << "ObsInsituTemperature::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
