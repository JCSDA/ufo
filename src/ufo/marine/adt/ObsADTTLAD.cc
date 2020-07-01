/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/adt/ObsADTTLAD.h"

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
static LinearObsOperatorMaker<ObsADTTLAD> makerADTTL_("ADT");
// -----------------------------------------------------------------------------

ObsADTTLAD::ObsADTTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));
  ufo_adt_tlad_setup_f90(keyOper_, config);
  oops::Log::trace() << "ObsADTTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADTTLAD::~ObsADTTLAD() {
  ufo_adt_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsADTTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                               ObsDiagnostics &) {
  ufo_adt_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsADTTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_adt_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsADTTLAD: tangent linear observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_adt_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsADTTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::print(std::ostream & os) const {
  os << "ObsADTTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
