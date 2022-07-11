/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGEWINDSPEED_H_
#define UFO_PROFILE_PROFILEAVERAGEWINDSPEED_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileAverageUtils.h"
#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: average wind speed observations onto model levels.
  ///
  /// Vectors produced by the AveragePressure routine must be present
  /// otherwise the exception eckit::BadValue will be thrown.
  ///
  /// The eastward and northward wind components are averaged separately
  /// over model layers defined by adjacent pressure levels, including the surface pressure.
  class ProfileAverageWindSpeed : public ProfileCheckBase {
   public:
    explicit ProfileAverageWindSpeed(const ConventionalProfileProcessingParameters &options);

    /// Average wind speed observations onto model levels and store the results.
    /// \throws eckit::BadValue if vectors produced by the AveragePressure routine
    /// are not present.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_surface_pressure});}

    /// List of names of GeoVaLs used in check validation.
    oops::Variables getValidationGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_testreference_eastward_wind,
            ufo::VariableNames::geovals_testreference_eastward_wind_qcflags,
            ufo::VariableNames::geovals_testreference_northward_wind,
            ufo::VariableNames::geovals_testreference_northward_wind_qcflags
            });}

   private:
    /// Run check on a profile in the original ObsSpace and
    /// put the averaged data into the corresponding profile in the extended ObsSpace.
    void runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                            ProfileDataHolder &profileExtended);
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEWINDSPEED_H_
