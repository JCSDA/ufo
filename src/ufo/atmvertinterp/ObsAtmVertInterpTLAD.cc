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
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerVertInterpTL_("VertInterp");
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerRadiosondeTL_("Radiosonde");
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerAircraftTL_("Aircraft");
static LinearObsOperatorMaker<ObsAtmVertInterpTLAD> makerSatwindTL_("Satwind");
// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::ObsAtmVertInterpTLAD(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : keyOperAtmVertInterp_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_atmvertinterp_tlad_setup_f90(keyOperAtmVertInterp_, &configc, &varconfig, varin_);

  oops::Log::trace() << "ObsAtmVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmVertInterpTLAD::~ObsAtmVertInterpTLAD() {
  ufo_atmvertinterp_tlad_delete_f90(keyOperAtmVertInterp_);
  oops::Log::trace() << "ObsAtmVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_atmvertinterp_tlad_settraj_f90(keyOperAtmVertInterp_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                         const ObsBiasIncrement & bias) const {
  ufo_atmvertinterp_simobs_tl_f90(keyOperAtmVertInterp_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                         ObsBiasIncrement & bias) const {
  ufo_atmvertinterp_simobs_ad_f90(keyOperAtmVertInterp_, geovals.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsAtmVertInterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
