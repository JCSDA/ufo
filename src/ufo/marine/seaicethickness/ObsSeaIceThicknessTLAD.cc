/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/seaicethickness/ObsSeaIceThicknessTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceThicknessTLAD> makerSeaIceThicknessTLAD_("SeaIceThickness");
// -----------------------------------------------------------------------------

ObsSeaIceThicknessTLAD::ObsSeaIceThicknessTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOperSeaIceThickness_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_seaicethick_tlad_setup_f90(keyOperSeaIceThickness_, &configc);
  const std::vector<std::string> vv{"ice_concentration", "ice_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceThicknessTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceThicknessTLAD::~ObsSeaIceThicknessTLAD() {
  ufo_seaicethick_tlad_delete_f90(keyOperSeaIceThickness_);
  oops::Log::trace() << "ObsSeaIceThicknessTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seaicethick_tlad_settraj_f90(keyOperSeaIceThickness_, geovals.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_seaicethick_simobs_tl_f90(keyOperSeaIceThickness_, geovals.toFortran(),
                                     ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_seaicethick_simobs_ad_f90(keyOperSeaIceThickness_, geovals.toFortran(),
                                     ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceThicknessTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
