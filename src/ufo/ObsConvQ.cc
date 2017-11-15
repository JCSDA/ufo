/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsConvQ.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsConvQ> makerConvQ_("ConvQ");
// -----------------------------------------------------------------------------

ObsConvQ::ObsConvQ(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperConvQ_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_conv_setup_f90(keyOperConvQ_, &configc);
  int keyVarin;
  ufo_conv_inputs_f90(keyOperConvQ_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsConvQ created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsConvQ::~ObsConvQ() {
  ufo_conv_delete_f90(keyOperConvQ_);
  oops::Log::trace() << "ObsConvQ destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsConvQ::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_conv_q_eqv_f90(gom.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsConvQ::print(std::ostream & os) const {
  os << "ObsConvQ::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
