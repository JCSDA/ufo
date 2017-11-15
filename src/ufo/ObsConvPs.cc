/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsConvPs.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsConvPs> makerConvPs_("ConvPs");
// -----------------------------------------------------------------------------

ObsConvPs::ObsConvPs(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperConvPs_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_conv_setup_f90(keyOperConvPs_, &configc);
  int keyVarin;
  ufo_conv_inputs_f90(keyOperConvPs_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsConvPs created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsConvPs::~ObsConvPs() {
  ufo_conv_delete_f90(keyOperConvPs_);
  oops::Log::trace() << "ObsConvPs destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsConvPs::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_conv_ps_eqv_f90(gom.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsConvPs::print(std::ostream & os) const {
  os << "ObsConvPs::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
