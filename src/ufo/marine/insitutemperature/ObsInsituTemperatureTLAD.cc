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

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsInsituTemperatureTLAD>
              makerInsituTemperatureTLAD_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::ObsInsituTemperatureTLAD(const ioda::ObsSpace & odb,
                                                   const eckit::Configuration & config)
  : keyOperInsituTemperature_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_tlad_setup_f90(keyOperInsituTemperature_, &configc);
  const std::vector<std::string> vv{"ocean_potential_temperature", "ocean_salinity",
                                    "ocean_layer_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsInsituTemperatureTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::~ObsInsituTemperatureTLAD() {
  ufo_insitutemperature_tlad_delete_f90(keyOperInsituTemperature_);
  oops::Log::trace() << "ObsInsituTemperatureTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_insitutemperature_tlad_settraj_f90(keyOperInsituTemperature_, geovals.toFortran(),
                                         odb_);
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                             const ObsBiasIncrement & bias) const {
  ufo_insitutemperature_simobs_tl_f90(keyOperInsituTemperature_, geovals.toFortran(),
                                        odb_, ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                             ObsBiasIncrement & bias) const {
  ufo_insitutemperature_simobs_ad_f90(keyOperInsituTemperature_, geovals.toFortran(),
                                        odb_, ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::print(std::ostream & os) const {
  os << "ObsInsituTemperatureTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
