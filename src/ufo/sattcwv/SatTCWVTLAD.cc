/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/sattcwv/SatTCWVTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<SatTCWVTLAD> makerSatTCWVTL_("SatTCWV");
// -----------------------------------------------------------------------------

SatTCWVTLAD::SatTCWVTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : LinearObsOperatorBase(odb), keyOperSatTCWV_(0), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "surface_pressure"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOptions(config, "obs options");
  const eckit::Configuration * configc = &obsOptions;
  ufo_sattcwv_tlad_setup_f90(keyOperSatTCWV_, &configc);

  oops::Log::trace() << "SatTCWVTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

SatTCWVTLAD::~SatTCWVTLAD() {
  ufo_sattcwv_tlad_delete_f90(keyOperSatTCWV_);
  oops::Log::trace() << "SatTCWVTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                           ObsDiagnostics &) {
  ufo_sattcwv_tlad_settraj_f90(keyOperSatTCWV_, geovals.toFortran(),
                               obsspace());
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_sattcwv_simobs_tl_f90(keyOperSatTCWV_, geovals.toFortran(),
                       obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_sattcwv_simobs_ad_f90(keyOperSatTCWV_, geovals.toFortran(),
                             obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void SatTCWVTLAD::print(std::ostream & os) const {
  os << "SatTCWVTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
