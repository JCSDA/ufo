/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/ObsColumnRetrieval.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsColumnRetrieval> makerColumnRetrieval_("ColumnRetrieval");
// -----------------------------------------------------------------------------

ObsColumnRetrieval::ObsColumnRetrieval(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  ufo_columnretrieval_setup_f90(keyOper_, parameters.toConfiguration(),
                                odb.assimvariables(), varin_);
  oops::Log::trace() << "ObsColumnRetrieval created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsColumnRetrieval::~ObsColumnRetrieval() {
  ufo_columnretrieval_delete_f90(keyOper_);
  oops::Log::trace() << "ObsColumnRetrieval destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrieval::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_columnretrieval_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsColumnRetrieval: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrieval::print(std::ostream & os) const {
  os << "ObsColumnRetrieval::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
