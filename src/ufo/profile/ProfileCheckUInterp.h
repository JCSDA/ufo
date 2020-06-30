/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKUINTERP_H_
#define UFO_PROFILE_PROFILECHECKUINTERP_H_

#include <algorithm>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileData.h"
#include "ufo/profile/ProfileFlags.h"
#include "ufo/profile/ProfileIndices.h"
#include "ufo/profile/ProfileStandardLevels.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: wind speed interpolation check
  class ProfileCheckUInterp : public ProfileCheckBase,
    private ProfileStandardLevels {
   public:
    ProfileCheckUInterp(const ProfileConsistencyCheckParameters &options,
                        const ProfileIndices &profileIndices,
                        const ProfileData &profileData,
                        ProfileFlags &profileFlags,
                        ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override;

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
