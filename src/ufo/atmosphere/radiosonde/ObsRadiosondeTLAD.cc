/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/radiosonde/ObsRadiosondeTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadiosondeTLAD> makerRadiosondeTL_("Radiosonde");
// -----------------------------------------------------------------------------

ObsRadiosondeTLAD::ObsRadiosondeTLAD(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOperRadiosonde_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiosonde_tlad_setup_f90(keyOperRadiosonde_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsRadiosondeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiosondeTLAD::~ObsRadiosondeTLAD() {
  ufo_radiosonde_tlad_delete_f90(keyOperRadiosonde_);
  oops::Log::trace() << "ObsRadiosondeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiosondeTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiosonde_tlad_settraj_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiosondeTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                      const ObsBiasIncrement & bias) const {
  ufo_radiosonde_simobs_tl_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran(),
                                   ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiosondeTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                      ObsBiasIncrement & bias) const {
  ufo_radiosonde_simobs_ad_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran(),
                                   ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiosondeTLAD::print(std::ostream & os) const {
  os << "ObsRadiosondeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
