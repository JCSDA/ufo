/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSRADIANCETLAD_H_
#define UFO_OBSRADIANCETLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Radiance (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsRadianceTLAD : public oops::LinearObsOperBase<MODEL>,
                          private util::ObjectCounter<ObsRadianceTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsRadianceTLAD";}

  ObsRadianceTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperRadiance_;}
  const int & toFortran() const {return keyOperRadiance_;}

private:
  void print(std::ostream &) const;
  F90hop keyOperRadiance_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadianceTLAD<MODEL>::ObsRadianceTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiance_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiance_tlad_setup_f90(keyOperRadiance_, &configc);
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
  oops::Log::trace() << "ObsRadianceTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadianceTLAD<MODEL>::~ObsRadianceTLAD() {
  ufo_radiance_tlad_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadianceTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadianceTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiance_tlad_settraj_f90(keyOperRadiance_, geovals.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadianceTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_radiance_tlad_eqv_tl_f90(keyOperRadiance_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadianceTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_radiance_tlad_eqv_ad_f90(keyOperRadiance_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadianceTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsRadianceTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSRADIANCETLAD_H_
