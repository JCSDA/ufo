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

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOV> makerRTTOV_("RTTOV");
// -----------------------------------------------------------------------------

ObsRadianceRTTOV::ObsRadianceRTTOV(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperRadianceRTTOV_(0), odb_(odb), varin_(), varout_(),
    obsname_("RTTOV:")
{
  obsname_ += config.getString("Sensor_ID");
  const std::vector<std::string> vv{
    "air_pressure", "air_pressure_at_two_meters_above_surface",
    "air_temperature", "air_temperature_at_two_meters_above_surface",
    "eastward_wind", "northward_wind",
    "skin_temperature",
    "specific_humidity", "specific_humidity_at_two_meters_above_surface",
    "surface_air_pressure", "surface_temperature",
    "surface_type", "water_type"};

  varin_.reset(new oops::Variables(vv));

  // parse channels from the config and create variable names
  std::string chlist = config.getString("channels");
  std::set<int> channels = oops::parseIntSet(chlist);
  std::vector<std::string> vout;
  channels_.reserve(channels.size());
  for (const int jj : channels) {
    vout.push_back("brightness_temperature_"+std::to_string(jj)+"_");
    channels_.push_back(jj);
  }
  varout_.reset(new oops::Variables(vout));

  // call Fortran setup routine
//  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
//  const eckit::Configuration * configc = &obsOptions;
  const eckit::Configuration * configc = &config;
  ufo_radiancerttov_setup_f90(keyOperRadianceRTTOV_, &configc);
  oops::Log::info() << "ObsRadianceRTTOV channels: " << channels << std::endl;

  oops::Log::trace() << "ObsRadianceRTTOV created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::~ObsRadianceRTTOV() {
  ufo_radiancerttov_delete_f90(keyOperRadianceRTTOV_);
  oops::Log::trace() << "ObsRadianceRTTOV destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec) const {
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
