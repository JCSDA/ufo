/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/stericheight/ObsStericHeight.h"

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsStericHeight> makerObsStericHeight_("ObsStericHeight");
// -----------------------------------------------------------------------------

ObsStericHeight::ObsStericHeight(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperStericHeight_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_stericheight_setup_f90(keyOperStericHeight_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid",
                                    "ocean_potential_temperature",
                                    "ocean_salinity"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsStericHeight created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsStericHeight::~ObsStericHeight() {
  ufo_stericheight_delete_f90(keyOperStericHeight_);
  oops::Log::trace() << "ObsStericHeight destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStericHeight::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_stericheight_simobs_f90(keyOperStericHeight_, gom.toFortran(),
                           odb_, ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsStericHeight::print(std::ostream & os) const {
  os << "ObsStericHeight::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
