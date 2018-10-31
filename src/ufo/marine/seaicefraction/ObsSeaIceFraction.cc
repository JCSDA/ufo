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

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFraction::ObsSeaIceFraction(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seaicefrac_setup_f90(keyOperSeaIceFraction_, &configc);
  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceFraction created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFraction::~ObsSeaIceFraction() {
  ufo_seaicefrac_delete_f90(keyOperSeaIceFraction_);
  oops::Log::trace() << "ObsSeaIceFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seaicefrac_simobs_f90(keyOperSeaIceFraction_, gom.toFortran(), odb_,
                         ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
