/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/aerosols/MR/ObsAodGeos.h"

#include <ostream>
#include <set>

#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"



#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodGeos> makerAodGeos_("AodGeos");
// -----------------------------------------------------------------------------

ObsAodGeos::ObsAodGeos(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  ufo_aodgeos_setup_f90(keyOper_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsAodGeos created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodGeos::~ObsAodGeos() {
  ufo_aodgeos_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodGeos destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeos::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec, ObsDiagnostics &) const {
  ufo_aodgeos_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsAodGeos: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeos::print(std::ostream & os) const {
  os << "ObsAodGeos::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
