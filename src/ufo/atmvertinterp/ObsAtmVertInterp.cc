/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmvertinterp/ObsAtmVertInterp.h"

#include <ostream>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmVertInterp> makerVertInterp_("VertInterp");
// -----------------------------------------------------------------------------

ObsAtmVertInterp::ObsAtmVertInterp(const ioda::ObsSpace & odb,
                                   const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperAtmVertInterp_(0),
    odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;
  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_atmvertinterp_setup_f90(keyOperAtmVertInterp_, &configc, &varconfig, varin_);

  oops::Log::trace() << "ObsAtmVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterp::~ObsAtmVertInterp() {
  ufo_atmvertinterp_delete_f90(keyOperAtmVertInterp_);
  oops::Log::trace() << "ObsAtmVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                   ObsDiagnostics &) const {
  oops::Log::trace() << "ObsAtmVertInterp::simulateObs entered" << std::endl;

  ufo_atmvertinterp_simobs_f90(keyOperAtmVertInterp_, gom.toFortran(), odb_,
                               ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterp::simulateObs exit" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterp::print(std::ostream & os) const {
  os << "ObsAtmVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
