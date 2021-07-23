/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKUINTERPALTERNATIVE_H_
#define UFO_PROFILE_PROFILECHECKUINTERPALTERNATIVE_H_

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileStandardLevels.h"

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: alternative wind speed interpolation check.
  /// This check is passed all of the profiles in the sample and runs on each in sequence.
  /// This check is otherwise identical to the ProfileCheckUInterp.
  class ProfileCheckUInterpAlternative : public ProfileCheckBase,
    private ProfileStandardLevels {
   public:
      explicit ProfileCheckUInterpAlternative
        (const ConventionalProfileProcessingParameters &options);

    /// Run check on all profiles.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Run check on an individual profile.
    void runCheckOnProfile(ProfileDataHolder &profile);

    /// Fill variables in validator
    void fillValidationData(ProfileDataHolder &profileDataHolder);

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}

   private:
    /// Number of failed checks by level
    std::vector <int> LevErrors_;

    /// Interpolated value of u
    std::vector <float> uInterp_;

    /// Interpolated value of v
    std::vector <float> vInterp_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKUINTERPALTERNATIVE_H_
