/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/stericheight/ObsStericHeightTLAD.h"

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
static LinearObsOperatorMaker<ObsStericHeightTLAD> makerStericHeightTL_("StericHeight");
// -----------------------------------------------------------------------------

ObsStericHeightTLAD::ObsStericHeightTLAD(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"sea_surface_height_above_geoid",
                                    "ocean_potential_temperature",
                                    "ocean_salinity"};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  ufo_stericheight_tlad_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsStericHeightTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsStericHeightTLAD::~ObsStericHeightTLAD() {
  ufo_stericheight_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsStericHeightTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStericHeightTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_stericheight_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsStericHeightTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStericHeightTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_stericheight_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                                 ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsStericHeightTLAD: tangent linear observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStericHeightTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_stericheight_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                                 ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsStericHeightTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStericHeightTLAD::print(std::ostream & os) const {
  os << "ObsStericHeightTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
