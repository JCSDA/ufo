/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKINTERPOLATION_H_
#define UFO_PROFILE_PROFILECHECKINTERPOLATION_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileStandardLevels.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: interpolation check
  class ProfileCheckInterpolation : public ProfileCheckBase,
    private ProfileStandardLevels {
   public:
    ProfileCheckInterpolation(const ProfileConsistencyCheckParameters &options,
                              ProfileDataHandler &profileDataHandler,
                              ProfileCheckValidator &profileCheckValidator);
    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override;

   private:
    /// Number of failed checks by level
    std::vector <int> LevErrors_;

    /// Interpolated value of T
    std::vector <float> tInterp_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKINTERPOLATION_H_
