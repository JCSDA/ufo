/*
 * (C) Copyright 2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/aerosols/AOP/ObsAodExt.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodExt> makerAodExt_("AodExt");
// -----------------------------------------------------------------------------

ObsAodExt::ObsAodExt(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  ufo_aodext_setup_f90(keyOper_, config, odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAodExt created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodExt::~ObsAodExt() {
  ufo_aodext_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodExt destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExt::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_aodext_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsAodExt: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodExt::print(std::ostream & os) const {
  os << "ObsAodExt::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
