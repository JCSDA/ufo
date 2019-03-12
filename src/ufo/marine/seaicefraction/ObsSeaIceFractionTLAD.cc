/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/seaicefraction/ObsSeaIceFractionTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceFractionTLAD> makerSeaIceFractionTL_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::ObsSeaIceFractionTLAD(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"sea_ice_area_fraction"};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  ufo_seaicefraction_tlad_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::~ObsSeaIceFractionTLAD() {
  ufo_seaicefraction_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceFractionTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seaicefraction_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsSeaIceFractionTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_seaicefraction_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                                   ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceFractionTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_seaicefraction_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                                   ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceFractionTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
