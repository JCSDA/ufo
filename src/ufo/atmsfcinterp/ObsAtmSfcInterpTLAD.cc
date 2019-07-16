/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmsfcinterp/ObsAtmSfcInterpTLAD.h"

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
static LinearObsOperatorMaker<ObsAtmSfcInterpTLAD> makerGSISfcModelTL_("GSISfcModel");;
// -----------------------------------------------------------------------------

ObsAtmSfcInterpTLAD::ObsAtmSfcInterpTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOperAtmSfcInterp_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_atmsfcinterp_tlad_setup_f90(keyOperAtmSfcInterp_, &configc, &varconfig, varin_);

  oops::Log::trace() << "ObsAtmSfcInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmSfcInterpTLAD::~ObsAtmSfcInterpTLAD() {
  ufo_atmsfcinterp_tlad_delete_f90(keyOperAtmSfcInterp_);
  oops::Log::trace() << "ObsAtmSfcInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterpTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_atmsfcinterp_tlad_settraj_f90(keyOperAtmSfcInterp_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsAtmSfcInterpTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_atmsfcinterp_simobs_tl_f90(keyOperAtmSfcInterp_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsAtmSfcInterpTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_atmsfcinterp_simobs_ad_f90(keyOperAtmSfcInterp_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterpTLAD::print(std::ostream & os) const {
  os << "ObsAtmSfcInterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
