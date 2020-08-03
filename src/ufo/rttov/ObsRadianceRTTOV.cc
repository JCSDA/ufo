/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/rttov/ObsRadianceRTTOV.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOV> makerRTTOV_("RTTOV");
// -----------------------------------------------------------------------------

ObsRadianceRTTOV::ObsRadianceRTTOV(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperRadianceRTTOV_(0), odb_(odb), varin_()
{
  std::vector<std::string> vv{
    "air_pressure", "air_pressure_at_two_meters_above_surface",
    "air_temperature", "air_temperature_at_two_meters_above_surface",
    "eastward_wind", "northward_wind",
    "skin_temperature",
    "specific_humidity", "specific_humidity_at_two_meters_above_surface",
    "surface_type", "water_type"};


  if (config.has("ExtraVars")) {
    std::vector<std::string> extravars;
    config.get("ExtraVars", extravars);
    oops::Log::info() << "ObsRadianceRTTOV extravars: " << extravars << std::endl;
    vv.insert(std::end(vv), std::begin(extravars), std::end(extravars));
  }

  oops::Log::info() << "ObsRadianceRTTOV vv: " << vv << std::endl;

  varin_.reset(new oops::Variables(vv));

  // get channels
  const oops::Variables & observed = odb.obsvariables();
  channels_ = observed.channels();

  // call Fortran setup routine
//  const eckit::LocalConfiguration obsOptions(config, "obs options");
//  const eckit::Configuration * configc = &obsOptions;
  ufo_radiancerttov_setup_f90(keyOperRadianceRTTOV_, config);
  oops::Log::info() << "ObsRadianceRTTOV channels: " << channels_ << std::endl;

  oops::Log::trace() << "ObsRadianceRTTOV created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::~ObsRadianceRTTOV() {
  ufo_radiancerttov_delete_f90(keyOperRadianceRTTOV_);
  oops::Log::trace() << "ObsRadianceRTTOV destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                   ObsDiagnostics &) const {
  oops::Log::trace() << "ObsRadianceRTTOV:: simulateObs started" << std::endl;

  ufo_radiancerttov_simobs_f90(keyOperRadianceRTTOV_, gom.toFortran(), odb_,
                          ovec.size(), ovec.toFortran(),
                          channels_.size(), channels_[0]);

  oops::Log::trace() << "ObsRadianceRTTOV:: simulateObs completed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::print(std::ostream & os) const {
  os << "ObsRadianceRTTOV::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
