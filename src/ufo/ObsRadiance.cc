/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsRadiance.h"

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsSpace.h"
#include "ObsVector.h"
#include "Fortran.h"
#include "util/Logger.h"

// -----------------------------------------------------------------------------
namespace ufo {
// -----------------------------------------------------------------------------
static oops::ObsOperatorMaker<UfoTrait, ObsRadiance> makerRadiance_("Radiance");
// -----------------------------------------------------------------------------

ObsRadiance::ObsRadiance(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiance_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiance_setup_f90(keyOperRadiance_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "humidity_mixing_ratio", "air_pressure", 
                                    "air_pressure_levels", "mass_concentration_of_ozone_in_air", 
                                    "mass_concentration_of_carbon_dioxide_in_air", 
                                    "atmosphere_mass_content_of_cloud_liquid_water", 
                                    "atmosphere_mass_content_of_cloud_ice", 
                                    "effective_radius_of_cloud_liquid_water_particle", 
                                    "effective_radius_of_cloud_ice_particle", 
                                    "Water_Fraction", "Land_Fraction", "Ice_Fraction", "Snow_Fraction",
                                    "Water_Temperature", "Land_Temperature", "Ice_Temperature", "Snow_Temperature",
                                    "Vegetation_Fraction", "Sfc_Wind_Speed", "Sfc_Wind_Direction", "Lai",
                                    "Soil_Moisture", "Soil_Temperature", "Land_Type_Index", "Vegetation_Type",
                                    "Soil_Type", "Snow_Depth"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsRadiance created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiance::~ObsRadiance() {
  ufo_radiance_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadiance destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiance::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_radiance_eqv_f90(gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiance::print(std::ostream & os) const {
  os << "ObsRadiance::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
