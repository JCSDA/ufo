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

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsInsituTemperature> makerInsituTemperature_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperature::ObsInsituTemperature(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : keyOperInsituTemperature_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_setup_f90(keyOperInsituTemperature_, &configc);
  const std::vector<std::string> vv{"ocean_potential_temperature", "ocean_salinity",
                                    "ocean_layer_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsInsituTemperature created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperature::~ObsInsituTemperature() {
  ufo_insitutemperature_delete_f90(keyOperInsituTemperature_);
  oops::Log::trace() << "ObsInsituTemperature destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                       const ObsBias & bias) const {
  ufo_insitutemperature_simobs_f90(keyOperInsituTemperature_, gom.toFortran(),
                                odb_, ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::print(std::ostream & os) const {
  os << "ObsInsituTemperature::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
