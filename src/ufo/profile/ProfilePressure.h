/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEPRESSURE_H_
#define UFO_PROFILE_PROFILEPRESSURE_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileVerticalInterpolation.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: calculates values of pressure that
  /// have not been reported due to the lack of a pressure sensor
  /// (e.g. for PILOT sondes and wind profilers).
  /// The missing pressures are calculated using the reported heights
  /// and linearly-interpolated model background fields.
  /// Profiles without a pressure sensor are flagged.
  class ProfilePressure : public ProfileCheckBase {
   public:
    explicit ProfilePressure(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_orog,
                              ufo::VariableNames::geovals_pressure_rho,
                              ufo::VariableNames::geovals_height_rho});
    }
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEPRESSURE_H_
