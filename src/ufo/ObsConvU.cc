/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsConvU.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsConvU> makerConvU_("ConvU");
// -----------------------------------------------------------------------------

ObsConvU::ObsConvU(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperConvU_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_conv_setup_f90(keyOperConvU_, &configc);
  int keyVarin;
  ufo_conv_inputs_f90(keyOperConvU_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsConvU created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsConvU::~ObsConvU() {
  ufo_conv_delete_f90(keyOperConvU_);
  oops::Log::trace() << "ObsConvU destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsConvU::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_conv_u_eqv_f90(gom.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsConvU::print(std::ostream & os) const {
  os << "ObsConvU::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
