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
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHistory.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: history check.
  class ProfileCheckHistory : public ProfileCheckBase {
   public:
    ProfileCheckHistory(const ProfileConsistencyCheckParameters &options,
                     ProfileDataHandler &profileDataHandler,
                     ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override {}

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKHISTORY_H_
