/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKHYDROSTATIC_H_
#define UFO_PROFILE_PROFILECHECKHYDROSTATIC_H_

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

  /// \brief Profile QC: hydrostatic check
  class ProfileCheckHydrostatic : public ProfileCheckBase,
    private ProfileStandardLevels {
   public:
      explicit ProfileCheckHydrostatic(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// This check requires HofX to have been calculated.
    bool requiresHofX() override {return true;}

    /// Fill variables in validator
    void fillValidationData(ProfileDataHandler &profileDataHandler) override;

   private:
    /// Hydrostatic error descriptions
    std::vector <std::string> HydDesc_;

    /// Constant in thickness calculation
    std::vector <float> DC_;

    /// Thickness tolerance
    std::vector <float> ETol_;

    /// Thickness calculated from temepature
    std::vector <float> D_;

    /// Thickness 'error'
    std::vector <float> E_;

    /// Hydrostatic flag by level
    std::vector <int> HydError_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKHYDROSTATIC_H_
