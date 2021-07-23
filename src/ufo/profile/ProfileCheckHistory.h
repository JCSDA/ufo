/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKHISTORY_H_
#define UFO_PROFILE_PROFILECHECKHISTORY_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHistory.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: history check.
  class ProfileCheckHistory : public ProfileCheckBase {
   public:
    explicit ProfileCheckHistory(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKHISTORY_H_
