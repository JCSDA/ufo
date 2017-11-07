/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsWSpeed.h"

#include "eckit/config/Configuration.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsSpace.h"
#include "ObsVector.h"
#include "Fortran.h"
#include "Variables.h"
#include "util/Logger.h"

// -----------------------------------------------------------------------------
namespace ufo {
// -----------------------------------------------------------------------------
static oops::ObsOperatorMaker<UfoTrait, ObsWSpeed> makerWSpeed_("WSpeed");
// -----------------------------------------------------------------------------

ObsWSpeed::ObsWSpeed(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperWspeed_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_wspeed_setup_f90(keyOperWspeed_, &configc);
  int keyVarin;
  ufo_obsoper_inputs_f90(keyOperWspeed_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsWSpeed created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsWSpeed::~ObsWSpeed() {
  ufo_wspeed_delete_f90(keyOperWspeed_);
  oops::Log::trace() << "ObsWSpeed destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeed::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_wspeed_eqv_f90(gom.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsWSpeed::print(std::ostream & os) const {
  os << "ObsWSpeed::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
