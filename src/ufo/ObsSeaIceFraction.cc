/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsSeaIceFraction.h"

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsSpace.h"
#include "ObsVector.h"
#include "Fortran.h"
#include "util/Logger.h"

// -----------------------------------------------------------------------------
namespace ufo {
// -----------------------------------------------------------------------------
static oops::ObsOperatorMaker<UfoTrait, ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFraction::ObsSeaIceFraction(const ObsSpace & odb, const eckit::Configuration & config)
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

void ObsSeaIceFraction::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seaicefrac_eqv_f90(gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
