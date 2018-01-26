/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ObsSeaIceFractionTLAD.h"

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsSpace.h"
#include "ObsVector.h"
#include "Fortran.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace ufo {
// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::ObsSeaIceFractionTLAD(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), varin_()
{
  const eckit::Configuration * configc = &config;    
  ufo_seaicefrac_setup_f90(keyOperSeaIceFraction_, &configc);
  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));
  Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::~ObsSeaIceFractionTLAD() {
  Log::trace() << "ObsSeaIceFractionTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::setTrajectory(const GeoVaLs &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::obsEquivTL(const GeoVaLs & geovals, ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_seaicefrac_eqv_tl_f90(geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::obsEquivAD(GeoVaLs & geovals, const ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_seaicefrac_eqv_ad_f90(geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
void ObsSeaIceFractionTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;  
}
}  // namespace ufo
