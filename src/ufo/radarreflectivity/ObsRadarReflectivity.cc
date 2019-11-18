/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/radarreflectivity/ObsRadarReflectivity.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadarReflectivity> makerRadarReflectivity_("RadarReflectivity");
// -----------------------------------------------------------------------------

ObsRadarReflectivity::ObsRadarReflectivity(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_radarreflectivity_setup_f90(keyOper_, &configc, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsRadarReflectivity created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarReflectivity::~ObsRadarReflectivity() {
  ufo_radarreflectivity_delete_f90(keyOper_);
  oops::Log::trace() << "ObsRadarReflectivity destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivity::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                       ObsDiagnostics &) const {
  ufo_radarreflectivity_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsRadarReflectivity: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivity::print(std::ostream & os) const {
  os << "ObsRadarReflectivity::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
