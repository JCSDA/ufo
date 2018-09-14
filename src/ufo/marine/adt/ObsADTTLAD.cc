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

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsADTTLAD> makerADTTLAD_("ADT");
// -----------------------------------------------------------------------------

  ObsADTTLAD::ObsADTTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperADT_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_adt_tlad_setup_f90(keyOperADT_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsADTTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADTTLAD::~ObsADTTLAD() {
  ufo_adt_tlad_delete_f90(keyOperADT_);
  oops::Log::trace() << "ObsADTTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_adt_tlad_settraj_f90(keyOperADT_, geovals.toFortran());  //, odb_.toFortran());
}

// -----------------------------------------------------------------------------

  void ObsADTTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                 const ObsBiasIncrement & bias) const {
  ufo_adt_tlad_eqv_tl_f90(keyOperADT_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

  void ObsADTTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                 ObsBiasIncrement & bias) const {
  ufo_adt_tlad_eqv_ad_f90(keyOperADT_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::print(std::ostream & os) const {
  os << "ObsADTTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
