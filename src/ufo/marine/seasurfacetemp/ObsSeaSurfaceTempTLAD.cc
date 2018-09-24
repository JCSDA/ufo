/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/seasurfacetemp/ObsSeaSurfaceTempTLAD.h"

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
static LinearObsOperatorMaker<ObsSeaSurfaceTempTLAD> makerSeaSurfaceTempTLAD_("SeaSurfaceTemp");
// -----------------------------------------------------------------------------

ObsSeaSurfaceTempTLAD::ObsSeaSurfaceTempTLAD(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & config)
  : keyOperSeaSurfaceTemp_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_tlad_setup_f90(keyOperSeaSurfaceTemp_, &configc);
  const std::vector<std::string> vv{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaSurfaceTempTLAD::~ObsSeaSurfaceTempTLAD() {
  ufo_seasurfacetemp_tlad_delete_f90(keyOperSeaSurfaceTemp_);
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seasurfacetemp_tlad_settraj_f90(keyOperSeaSurfaceTemp_,
                                      geovals.toFortran());  //, odb_.toFortran());
}

// -----------------------------------------------------------------------------

  void ObsSeaSurfaceTempTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                            const ObsBiasIncrement & bias) const {
  ufo_seasurfacetemp_simobs_tl_f90(keyOperSeaSurfaceTemp_, geovals.toFortran(),
                                     ovec.toFortran());
}

// -----------------------------------------------------------------------------

  void ObsSeaSurfaceTempTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                            ObsBiasIncrement & bias) const {
  ufo_seasurfacetemp_simobs_ad_f90(keyOperSeaSurfaceTemp_, geovals.toFortran(),
                                        ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTempTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
