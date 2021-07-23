/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKTIME_H_
#define UFO_PROFILE_PROFILECHECKTIME_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: reject data which are outside the assimilation time window.
  /// Also, if requested, reject data taken a short period after the sonde launch.
  class ProfileCheckTime : public ProfileCheckBase {
   public:
    explicit ProfileCheckTime(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKTIME_H_
