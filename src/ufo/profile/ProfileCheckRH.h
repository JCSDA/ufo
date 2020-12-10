/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKRH_H_
#define UFO_PROFILE_PROFILECHECKRH_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: relative humidity check.
  /// - Check for increasing dew-point at the top of cloud levels (cloud errors).
  /// - Check for large humidity errors at the top of the profile (hi-lev errors).

  class ProfileCheckRH : public ProfileCheckBase {
   public:
      ProfileCheckRH(const ProfileConsistencyCheckParameters &options,
                     ProfileDataHandler &profileDataHandler,
                     ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override;

   private:
    /// Observed pressure for selected levels (hPa)
    std::vector <float> Press_;

    /// Observed temperature for selected levels (K)
    std::vector <float> Temp_;

    /// Observed relative humidity for selected levels (%)
    std::vector <float> rh_;

    /// Observed dew point temperature for selected levels (K)
    std::vector <float> td_;

    /// Model temperature for selected levels (K)
    std::vector <float> tbk_;

    /// Model relative humidity for selected levels (%)
    std::vector <float> rhbk_;

    /// Flags for RH checks
    std::vector <int> FlagH_;

    /// Indices of selected levels
    std::vector <int> Indx_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKRH_H_
