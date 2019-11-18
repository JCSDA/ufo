/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/sfcpcorrected/ObsSfcPCorrected.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSfcPCorrected> makerSfcPCorrected_("SfcPCorrected");
// -----------------------------------------------------------------------------

ObsSfcPCorrected::ObsSfcPCorrected(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_sfcpcorrected_setup_f90(keyOper_, &configc, odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsSfcPCorrected created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSfcPCorrected::~ObsSfcPCorrected() {
  ufo_sfcpcorrected_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSfcPCorrected destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSfcPCorrected::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_sfcpcorrected_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsSfcPCorrected: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSfcPCorrected::print(std::ostream & os) const {
  os << "ObsSfcPCorrected::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
