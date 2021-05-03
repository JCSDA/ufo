/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmvertinterp/ObsAtmVertInterpTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerVertInterpTL_("VertInterp");
// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::ObsAtmVertInterpTLAD(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : LinearObsOperatorBase(odb), keyOperAtmVertInterp_(0), varin_()
{
  ufo_atmvertinterp_tlad_setup_f90(keyOperAtmVertInterp_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsAtmVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::~ObsAtmVertInterpTLAD() {
  ufo_atmvertinterp_tlad_delete_f90(keyOperAtmVertInterp_);
  oops::Log::trace() << "ObsAtmVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                         ObsDiagnostics &) {
  oops::Log::trace() << "ObsAtmVertInterpTLAD::setTrajectory entering" << std::endl;

  ufo_atmvertinterp_tlad_settraj_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_atmvertinterp_simobs_tl_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_atmvertinterp_simobs_ad_f90(keyOperAtmVertInterp_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsAtmVertInterpTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsAtmVertInterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
