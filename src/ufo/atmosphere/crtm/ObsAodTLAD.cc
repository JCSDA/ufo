/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/crtm/ObsAodTLAD.h"

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
static LinearObsOperatorMaker<ObsAodTLAD> makerAodTL_("Aod");
// -----------------------------------------------------------------------------

ObsAodTLAD::ObsAodTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAod_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aod_tlad_setup_f90(keyOperAod_, &configc);
  const std::vector<std::string> vv{"air_temperature", "humidity_mixing_ratio","relative_humidity",
      "air_pressure", "air_pressure_levels",
      "sulf", "bc1", "bc2", "oc1", "oc2", "dust1", "dust2", "dust3", "dust4", "dust5",
      "seas1", "seas2", "seas3", "seas4"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAodTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodTLAD::~ObsAodTLAD() {
  oops::Log::trace() << "ObsAodTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_aod_tlad_settraj_f90(keyOperAod_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsAodTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_aod_simobs_tl_f90(keyOperAod_, geovals.toFortran(), odb_,
                        ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_aod_simobs_ad_f90(keyOperAod_, geovals.toFortran(), odb_,
                        ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodTLAD::print(std::ostream & os) const {
  os << "ObsAodTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
