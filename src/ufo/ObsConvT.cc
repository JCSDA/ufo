/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsConvT.h"

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
static oops::ObsOperatorMaker<UfoTrait, ObsConvT> makerConvT_("ConvT");
// -----------------------------------------------------------------------------

ObsConvT::ObsConvT(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperConvT_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_conv_setup_f90(keyOperConvT_, &configc);
  int keyVarin;
  ufo_conv_inputs_f90(keyOperConvT_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  oops::Log::trace() << "ObsConvT created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsConvT::~ObsConvT() {
  ufo_conv_delete_f90(keyOperConvT_);
  oops::Log::trace() << "ObsConvT destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsConvT::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_conv_t_eqv_f90(gom.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsConvT::print(std::ostream & os) const {
  os << "ObsConvT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
