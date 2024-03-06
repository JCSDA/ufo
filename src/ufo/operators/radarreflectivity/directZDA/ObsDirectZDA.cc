/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/radarreflectivity/directZDA/ObsDirectZDA.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsDirectZDA> makerDirectZDA_("DirectZDA");
// -----------------------------------------------------------------------------

ObsDirectZDA::ObsDirectZDA(const ioda::ObsSpace & odb,
                           const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperDirectZDA_(0), odb_(odb), varin_()
{
  ufo_directZDA_setup_f90(keyOperDirectZDA_, params.toConfiguration(),
                          odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsDirectZDA created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsDirectZDA::~ObsDirectZDA() {
  ufo_directZDA_delete_f90(keyOperDirectZDA_);
  oops::Log::trace() << "ObsDirectZDA destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDA::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                       ObsDiagnostics & diags,
                                       const QCFlags_t & qc_flags) const {
  ufo_directZDA_simobs_f90(keyOperDirectZDA_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsDirectZDA: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDA::print(std::ostream & os) const {
  os << "ObsDirectZDA::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
