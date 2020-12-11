/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKUNSTABLELAYER_H_
#define UFO_PROFILE_PROFILECHECKUNSTABLELAYER_H_

#include <utility>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: unstable layer check
  class ProfileCheckUnstableLayer : public ProfileCheckBase {
   public:
    ProfileCheckUnstableLayer(const ProfileConsistencyCheckParameters &options,
                              ProfileDataHandler &profileDataHandler,
                              ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override;

   private:
    /// PBottom
    float PBottom_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKUNSTABLELAYER_H_
