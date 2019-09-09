/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/insitutemperature/ObsInsituTemperatureTLAD.h"

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
static LinearObsOperatorMaker<ObsInsituTemperatureTLAD>
   makerInsituTemperatureTL_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::ObsInsituTemperatureTLAD(const ioda::ObsSpace & odb,
                                                   const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"sea_water_potential_temperature",
                                    "sea_water_salinity",
                                    "sea_water_cell_thickness"};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_tlad_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsInsituTemperatureTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::~ObsInsituTemperatureTLAD() {
  ufo_insitutemperature_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituTemperatureTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_insitutemperature_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsInsituTemperatureTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsTL(const GeoVaLs & geovals,
                                             ioda::ObsVector & ovec) const {
  ufo_insitutemperature_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperatureTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsAD(GeoVaLs & geovals,
                                             const ioda::ObsVector & ovec) const {
  ufo_insitutemperature_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperatureTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::print(std::ostream & os) const {
  os << "ObsInsituTemperatureTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
