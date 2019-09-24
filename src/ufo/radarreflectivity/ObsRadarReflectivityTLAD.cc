/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/radarreflectivity/ObsRadarReflectivityTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadarReflectivityTLAD>
       makerRadarReflectivityTL_("RadarReflectivity");
// -----------------------------------------------------------------------------

ObsRadarReflectivityTLAD::ObsRadarReflectivityTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;
  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_radarreflectivity_tlad_setup_f90(keyOper_, &configc, &varconfig, varin_);

  oops::Log::trace() << "ObsRadarReflectivityTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarReflectivityTLAD::~ObsRadarReflectivityTLAD() {
  ufo_radarreflectivity_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsRadarReflectivityTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radarreflectivity_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsRadarReflectivityTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::simulateObsTL(const GeoVaLs & geovals,
                                             ioda::ObsVector & ovec) const {
  ufo_radarreflectivity_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsRadarReflectivityTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::simulateObsAD(GeoVaLs & geovals,
                                             const ioda::ObsVector & ovec) const {
  ufo_radarreflectivity_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsRadarReflectivityTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::print(std::ostream & os) const {
  os << "ObsRadarReflectivityTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
