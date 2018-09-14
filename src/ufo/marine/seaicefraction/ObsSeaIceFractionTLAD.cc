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

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceFractionTLAD> makerSeaIceFractionTLAD_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::ObsSeaIceFractionTLAD(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_seaicefrac_tlad_setup_f90(keyOperSeaIceFraction_, &configc);
  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::~ObsSeaIceFractionTLAD() {
  ufo_seaicefrac_tlad_delete_f90(keyOperSeaIceFraction_);
  oops::Log::trace() << "ObsSeaIceFractionTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seaicefrac_tlad_settraj_f90(keyOperSeaIceFraction_, geovals.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_seaicefrac_tlad_eqv_tl_f90(keyOperSeaIceFraction_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_seaicefrac_tlad_eqv_ad_f90(keyOperSeaIceFraction_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
