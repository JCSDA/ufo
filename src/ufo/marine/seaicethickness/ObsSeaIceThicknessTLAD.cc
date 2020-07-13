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

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceThicknessTLAD> makerSeaIceThicknessTL_("SeaIceThickness");
// -----------------------------------------------------------------------------

ObsSeaIceThicknessTLAD::ObsSeaIceThicknessTLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"sea_ice_category_area_fraction",
                                    "sea_ice_category_thickness"};
  varin_.reset(new oops::Variables(vv));
  ufo_seaicethickness_tlad_setup_f90(keyOper_, config, odb.obsvariables());
  oops::Log::trace() << "ObsSeaIceThicknessTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceThicknessTLAD::~ObsSeaIceThicknessTLAD() {
  ufo_seaicethickness_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceThicknessTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                           ObsDiagnostics &) {
  ufo_seaicethickness_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsSeaIceThicknessTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_seaicethickness_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                                    ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceThicknessTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_seaicethickness_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                                    ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceThicknessTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThicknessTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceThicknessTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
