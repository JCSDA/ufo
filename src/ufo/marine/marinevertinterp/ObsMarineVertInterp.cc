/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/marinevertinterp/ObsMarineVertInterp.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsMarineVertInterp> makerMarineVertInterp_("MarineVertInterp");
// -----------------------------------------------------------------------------

ObsMarineVertInterp::ObsMarineVertInterp(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  ufo_marinevertinterp_setup_f90(keyOper_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsMarineVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterp::~ObsMarineVertInterp() {
  ufo_marinevertinterp_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                      ObsDiagnostics &) const {
  ufo_marinevertinterp_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::print(std::ostream & os) const {
  os << "ObsMarineVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
