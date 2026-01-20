/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSROGEOVALS_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSROGEOVALS_H_

#include <memory>
#include <string>
#include <vector>
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// GnssroGeoVaLs -
///   C++ Interface to GeoVaLs variables needed by GNSSRO forward operators.
///   NOTE: this does not account for FO's that need some variables on
///         different vertical levels, e.g. air_pressure_levels variable.
///
class GnssroGeoVaLs {
 public:
  // Type definitions
  typedef std::vector<double> Profile_;
  typedef std::unique_ptr<Profile_> ProfilePtr_;
  typedef std::vector<ProfilePtr_> Profiles_;

  // String constants
  static constexpr const char* CLASSNAME = "GnssroGeoVaLs";
  static constexpr const char* VAR_TEMP = "air_temperature";
  static constexpr const char* VAR_SPHUM = "water_vapor_mixing_ratio_wrt_moist_air";
  static constexpr const char* VAR_PRES = "air_pressure";
  static constexpr const char* VAR_GPH = "geopotential_height";
  static constexpr const char* VAR_SFC_ALT = "height_above_mean_sea_level_at_surface";
  static constexpr std::array<const char*, 4> PROFILES_VARS = {
          VAR_TEMP, VAR_SPHUM, VAR_PRES, VAR_GPH};

  static std::string name() {return CLASSNAME;}

  //  Constructor/Destructor
  explicit GnssroGeoVaLs(const GeoVaLs & gv);
  ~GnssroGeoVaLs();
  // Disallow default, copy and move construction and assignment.
  GnssroGeoVaLs() = delete;
  GnssroGeoVaLs(const GnssroGeoVaLs & other) = delete;
  GnssroGeoVaLs(GnssroGeoVaLs && other) noexcept = delete;
  GnssroGeoVaLs & operator=(const GnssroGeoVaLs & other) = delete;
  GnssroGeoVaLs & operator=(GnssroGeoVaLs && other) noexcept = delete;

  // Accessors for dimensions
  std::size_t nprofiles() const {return nprofiles_;}
  std::size_t nlevs() const {return nlevs_;}

  // Profile Read Accessors
  const Profile_ & tempProfile(std::size_t profileIdx) const {return *(temp_[profileIdx]);}
  const Profile_ & sphumProfile(std::size_t profileIdx) const {return *(sphum_[profileIdx]);}
  const Profile_ & presProfile(std::size_t profileIdx) const {return *(pres_[profileIdx]);}
  const Profile_ & gphProfile(std::size_t profileIdx) const {return *(gph_[profileIdx]);}

  // Profile Read-write Accessors
  Profile_ & tempProfile(std::size_t profileIdx) {return *(temp_[profileIdx]);}
  Profile_ & sphumProfile(std::size_t profileIdx) {return *(sphum_[profileIdx]);}
  Profile_ & presProfile(std::size_t profileIdx) {return *(pres_[profileIdx]);}
  Profile_ & gphProfile(std::size_t profileIdx) {return *(gph_[profileIdx]);}

  // Scalar element read accessors
  double temp(std::size_t profileIdx, std::size_t levIdx) const {
    return (*(temp_[profileIdx]))[levIdx];
  }
  double sphum(std::size_t profileIdx, std::size_t levIdx) const {
    return (*(sphum_[profileIdx]))[levIdx];
  }
  double pres(std::size_t profileIdx, std::size_t levIdx) const {
    return (*(pres_[profileIdx]))[levIdx];
  }
  double gph(std::size_t profileIdx, std::size_t levIdx) const {
    return (*(gph_[profileIdx]))[levIdx];
  }

  // Scalar element read-write accessors
  double & temp(std::size_t profileIdx, std::size_t levIdx) {
    return (*(temp_[profileIdx]))[levIdx];
  }
  double & sphum(std::size_t profileIdx, std::size_t levIdx) {
    return (*(sphum_[profileIdx]))[levIdx];
  }
  double & pres(std::size_t profileIdx, std::size_t levIdx) {
    return (*(pres_[profileIdx]))[levIdx];
  }
  double & gph(std::size_t profileIdx, std::size_t levIdx) {
    return (*(gph_[profileIdx]))[levIdx];
  }

  // Save values to GeoVaLs
  void saveTo(GeoVaLs & gv);

 private:
  // Implementation methods
  std::size_t validateProfiles(const GeoVaLs & cgv);
  std::size_t validateLevels(const GeoVaLs & cgv);
  void getProfile(const oops::Variable& var,
                  std::size_t profIdx,
                  Profiles_ & profiles,
                  const GeoVaLs& cgv);
  void putProfiles(const oops::Variable& var,
                   const Profiles_& profiles,
                   GeoVaLs& gv);

  // Member data
  std::size_t nprofiles_;  // Number of profiles (sampled locations)
  std::size_t nlevs_;      // Number of levels in each profile
  Profiles_ temp_;         // Air temperature profiles.
  Profiles_ sphum_;        // Specific humidity profiles
  Profiles_ pres_;         // Air pressure profiles
  Profiles_ gph_;          // Geopotential height profiles.
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSROGEOVALS_H_
