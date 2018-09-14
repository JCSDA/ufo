/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/radiance/ObsRadianceTLAD.h"

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
static LinearObsOperatorMaker<ObsRadianceTLAD> makerRadianceTL_("Radiance");
// -----------------------------------------------------------------------------

ObsRadianceTLAD::ObsRadianceTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
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
  const eckit::Configuration * configc = &config;
  ufo_radiance_tlad_setup_f90(keyOperRadiance_, &configc);
  oops::Log::trace() << "ObsRadianceTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceTLAD::~ObsRadianceTLAD() {
  ufo_radiance_tlad_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadianceTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiance_tlad_settraj_f90(keyOperRadiance_, geovals.toFortran(), odb_.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_radiance_tlad_eqv_tl_f90(keyOperRadiance_, geovals.toFortran(), odb_.toFortran(),
                               ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_radiance_tlad_eqv_ad_f90(keyOperRadiance_, geovals.toFortran(), odb_.toFortran(),
                               ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::print(std::ostream & os) const {
  os << "ObsRadianceTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
