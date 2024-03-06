/*
 * (C) Copyright 2021.
 *
 * This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/insitupm/ObsInsituPM.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsInsituPM> makerInsituPM_("InsituPM");
// -----------------------------------------------------------------------------

ObsInsituPM::ObsInsituPM(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  ufo_insitupm_setup_f90(keyOper_, parameters.toConfiguration(), odb.assimvariables(), varin_);
  oops::Log::trace() << "ObsInsituPM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituPM::~ObsInsituPM() {
  ufo_insitupm_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituPM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPM::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics & odiags, const QCFlags_t & qc_flags) const {
  ufo_insitupm_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsInsituPM: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPM::print(std::ostream & os) const {
  os << "ObsInsituPM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
