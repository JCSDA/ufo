/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsRadiosonde.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsRadiosonde> makerRadiosonde_("Radiosonde");
// -----------------------------------------------------------------------------

ObsRadiosonde::ObsRadiosonde(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiosonde_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiosonde_setup_f90(keyOperRadiosonde_, &configc);
  int keyVarin;
  ufo_radiosonde_inputs_f90(keyOperRadiosonde_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsRadiosonde created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiosonde::~ObsRadiosonde() {
  ufo_radiosonde_delete_f90(keyOperRadiosonde_);
  oops::Log::trace() << "ObsRadiosonde destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiosonde::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_radiosonde_t_eqv_f90(gom.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiosonde::print(std::ostream & os) const {
  os << "ObsRadiosonde::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
