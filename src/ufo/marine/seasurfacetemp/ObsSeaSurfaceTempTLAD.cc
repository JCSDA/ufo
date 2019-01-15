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

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaSurfaceTempTLAD> makerSeaSurfaceTempTL_("SeaSurfaceTemp");
// -----------------------------------------------------------------------------

ObsSeaSurfaceTempTLAD::ObsSeaSurfaceTempTLAD(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_tlad_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaSurfaceTempTLAD::~ObsSeaSurfaceTempTLAD() {
  ufo_seasurfacetemp_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seasurfacetemp_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_seasurfacetemp_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_seasurfacetemp_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaSurfaceTempTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTempTLAD::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTempTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
