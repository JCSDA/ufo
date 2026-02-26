/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/atmvertinterplay/ObsAtmVertInterpLayTLAD.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmVertInterpLayTLAD> makerVertInterpLayTL_("AtmVertInterpLay");
// -----------------------------------------------------------------------------

ObsAtmVertInterpLayTLAD::ObsAtmVertInterpLayTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperAtmVertInterpLay_(0), varin_()
{
  ufo_atmvertinterplay_tlad_setup_f90(keyOperAtmVertInterpLay_, params.toConfiguration(),
                                      odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpLayTLAD::~ObsAtmVertInterpLayTLAD() {
  ufo_atmvertinterplay_tlad_delete_f90(keyOperAtmVertInterpLay_);
  oops::Log::trace() << "ObsAtmVertInterpLayTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                            const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::setTrajectory start" << std::endl;

  ufo_atmvertinterplay_tlad_settraj_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::simulateObsTL(const GeoVaLs & geovals,
                                            ioda::ObsVector & ovec) const {
  ufo_atmvertinterplay_simobs_tl_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::simulateObsAD(GeoVaLs & geovals,
                                            const ioda::ObsVector & ovec) const {
  ufo_atmvertinterplay_simobs_ad_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::print(std::ostream & os) const {
  os << "ObsAtmVertInterpLayTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
