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
#include <vector>

#include "ufo/profile/ModelHeightCalculator.h"
#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileVerticalInterpolation.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
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
    ProfilePressure(const ProfileConsistencyCheckParameters &options,
                    ProfileDataHandler &profileDataHandler,
                    ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override {}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEPRESSURE_H_
