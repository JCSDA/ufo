/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmvertinterplay/ObsAtmVertInterpLayTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmVertInterpLayTLAD> makerVertInterpLayTL_("AtmVertInterpLay");
// -----------------------------------------------------------------------------

ObsAtmVertInterpLayTLAD::ObsAtmVertInterpLayTLAD(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : keyOperAtmVertInterpLay_(0), odb_(odb), varin_()
{
  ufo_atmvertinterplay_tlad_setup_f90(keyOperAtmVertInterpLay_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpLayTLAD::~ObsAtmVertInterpLayTLAD() {
  ufo_atmvertinterplay_tlad_delete_f90(keyOperAtmVertInterpLay_);
  oops::Log::trace() << "ObsAtmVertInterpLayTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                         ObsDiagnostics &) {
  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::setTrajectory entering" << std::endl;

  ufo_atmvertinterplay_tlad_settraj_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), odb_);

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_atmvertinterplay_simobs_tl_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_atmvertinterplay_simobs_ad_f90(keyOperAtmVertInterpLay_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpLayTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpLayTLAD::print(std::ostream & os) const {
  os << "ObsAtmVertInterpLayTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
