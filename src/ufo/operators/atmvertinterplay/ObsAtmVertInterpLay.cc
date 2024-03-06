/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>

#include "ioda/ObsVector.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/atmvertinterplay/ObsAtmVertInterpLay.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmVertInterpLay> makerAtmVertInterpLay_("AtmVertInterpLay");
// -----------------------------------------------------------------------------

ObsAtmVertInterpLay::ObsAtmVertInterpLay(const ioda::ObsSpace & odb,
                       const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperAtmVertInterpLay_(0), odb_(odb), varin_()
{
  ufo_atmvertinterplay_setup_f90(keyOperAtmVertInterpLay_, params.toConfiguration(),
                                 odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsAtmVertInterpLay created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpLay::~ObsAtmVertInterpLay() {
  ufo_atmvertinterplay_delete_f90(keyOperAtmVertInterpLay_);
  oops::Log::trace() << "ObsAtmVertInterpLay destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLay::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                      ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  ufo_atmvertinterplay_simobs_f90(keyOperAtmVertInterpLay_, gom.toFortran(), odb_, ovec.nvars(),
                         ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLay::print(std::ostream & os) const {
  os << "ObsAtmVertInterpLay::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
