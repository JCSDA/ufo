/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKUINTERP_H_
#define UFO_PROFILE_PROFILECHECKUINTERP_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileStandardLevels.h"

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: wind speed interpolation check
  class ProfileCheckUInterp : public ProfileCheckBase,
    private ProfileStandardLevels {
   public:
      explicit ProfileCheckUInterp(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Fill variables in validator
    void fillValidationData(ProfileDataHandler &profileDataHandler) override;

   private:
    /// Number of failed checks by level
    std::vector <int> LevErrors_;

    /// Interpolated value of u
    std::vector <float> uInterp_;

    /// Interpolated value of v
    std::vector <float> vInterp_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKUINTERP_H_
