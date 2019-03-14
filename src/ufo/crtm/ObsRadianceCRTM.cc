/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/crtm/ObsRadianceCRTM.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceCRTM> makerCRTM_("CRTM");

// -----------------------------------------------------------------------------

ObsRadianceCRTM::ObsRadianceCRTM(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperRadianceCRTM_(0), odb_(odb), varin_(), varout_()
{
  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const int n_clouds = obsOptions.getInt("n_Clouds");
  std::vector<std::string> vv{"air_temperature", "humidity_mixing_ratio", "air_pressure",
                              "air_pressure_levels", "mass_concentration_of_ozone_in_air",
                              "mass_concentration_of_carbon_dioxide_in_air",
                              "Water_Fraction", "Land_Fraction", "Ice_Fraction",
                              "Snow_Fraction", "Water_Temperature", "Land_Temperature",
                              "Ice_Temperature", "Snow_Temperature", "Vegetation_Fraction",
                              "Sfc_Wind_Speed", "Sfc_Wind_Direction", "Lai", "Soil_Moisture",
                              "Soil_Temperature", "Land_Type_Index", "Vegetation_Type",
                              "Soil_Type", "Snow_Depth"};
  // Note: how should we control the clouds used in CRTM?
  // currently uses n_clouds; if n_clouds = 1, use water; n_clouds = 2, water and ice
  if (n_clouds > 0) {
    vv.push_back("atmosphere_mass_content_of_cloud_liquid_water");
    vv.push_back("effective_radius_of_cloud_liquid_water_particle");
  }
  if (n_clouds > 1) {
    vv.push_back("atmosphere_mass_content_of_cloud_ice");
    vv.push_back("effective_radius_of_cloud_ice_particle");
  }
  varin_.reset(new oops::Variables(vv));

  // parse channels from the config and create variable names
  std::string chlist = config.getString("channels");
  std::set<int> channels = parseIntSet(chlist);
  std::vector<std::string> vout;
  channels_.reserve(channels.size());
  for (const int jj : channels) {
    vout.push_back("brightness_temperature_"+std::to_string(jj));
    channels_.push_back(jj);
  }
  varout_.reset(new oops::Variables(vout));

  // call Fortran setup routine
  const eckit::Configuration * configc = &obsOptions;
  ufo_radiancecrtm_setup_f90(keyOperRadianceCRTM_, &configc);
  oops::Log::info() << "ObsRadianceCRTM channels: " << channels << std::endl;
  oops::Log::trace() << "ObsRadianceCRTM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTM::~ObsRadianceCRTM() {
  ufo_radiancecrtm_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_radiancecrtm_simobs_f90(keyOperRadianceCRTM_, gom.toFortran(), odb_,
                          ovec.size(), ovec.toFortran(), bias.toFortran(),
                          channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::print(std::ostream & os) const {
  os << "ObsRadianceCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
