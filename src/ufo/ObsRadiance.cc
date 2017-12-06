/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsRadiance.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsRadiance> makerRadiance_("Radiance");
// -----------------------------------------------------------------------------

ObsRadiance::ObsRadiance(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiance_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiance_setup_f90(keyOperRadiance_, &configc);
  const std::vector<std::string> vv{"u","v"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsRadiance created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiance::~ObsRadiance() {
  ufo_radiance_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadiance destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiance::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_radiance_eqv_f90(gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiance::print(std::ostream & os) const {
  os << "ObsRadiance::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
