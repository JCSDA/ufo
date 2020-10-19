/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/avgkernel/ObsAvgKernel.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAvgKernel> makerAvgKernel_("AvgKernel");
// -----------------------------------------------------------------------------

ObsAvgKernel::ObsAvgKernel(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  ufo_avgkernel_setup_f90(keyOper_, config, odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAvgKernel created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAvgKernel::~ObsAvgKernel() {
  ufo_avgkernel_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAvgKernel destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernel::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_avgkernel_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsAvgKernel: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAvgKernel::print(std::ostream & os) const {
  os << "ObsAvgKernel::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
