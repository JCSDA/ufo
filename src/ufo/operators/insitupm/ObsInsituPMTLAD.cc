/*
 * (C) Copyright 2021.
 *
 * This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/insitupm/ObsInsituPMTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsInsituPMTLAD> makerInsituPMTL_("InsituPM");
// -----------------------------------------------------------------------------

ObsInsituPMTLAD::ObsInsituPMTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_insitupm_tlad_setup_f90(keyOper_, parameters.toConfiguration(), odb.assimvariables(), varin_);
  oops::Log::trace() << "ObsInsituPMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituPMTLAD::~ObsInsituPMTLAD() {
  ufo_insitupm_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituPMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPMTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & d,
                                    const QCFlags_t & qc_flags) {
  ufo_insitupm_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsInsituPMTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const QCFlags_t & qc_flags) const {
  ufo_insitupm_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituPMTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    const QCFlags_t & qc_flags) const {
  ufo_insitupm_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituPMTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituPMTLAD::print(std::ostream & os) const {
  os << "ObsInsituPMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
