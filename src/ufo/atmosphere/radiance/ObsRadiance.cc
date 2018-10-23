/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/radiance/ObsRadiance.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadiance> makerRadiance_("Radiance");
// -----------------------------------------------------------------------------

ObsRadiance::ObsRadiance(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiance_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"virtual_temperature", "humidity_mixing_ratio", "air_pressure",
                                    "air_pressure_levels", "mass_concentration_of_ozone_in_air",
                                    "mass_concentration_of_carbon_dioxide_in_air",
                                    "atmosphere_mass_content_of_cloud_liquid_water",
                                    "atmosphere_mass_content_of_cloud_ice",
                                    "effective_radius_of_cloud_liquid_water_particle",
                                    "effective_radius_of_cloud_ice_particle",
                                    "Water_Fraction", "Land_Fraction", "Ice_Fraction",
                                    "Snow_Fraction", "Water_Temperature", "Land_Temperature",
                                    "Ice_Temperature", "Snow_Temperature", "Vegetation_Fraction",
                                    "Sfc_Wind_Speed", "Sfc_Wind_Direction", "Lai", "Soil_Moisture",
                                    "Soil_Temperature", "Land_Type_Index", "Vegetation_Type",
                                    "Soil_Type", "Snow_Depth"};
  varin_.reset(new oops::Variables(vv));
  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;
  ufo_radiance_setup_f90(keyOperRadiance_, &configc);
  oops::Log::trace() << "ObsRadiance created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiance::~ObsRadiance() {
  ufo_radiance_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadiance destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiance::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_radiance_simobs_f90(keyOperRadiance_, gom.toFortran(), odb_, ovec.toFortran(),
                       bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiance::print(std::ostream & os) const {
  os << "ObsRadiance::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
