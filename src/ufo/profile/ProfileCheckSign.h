/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKSIGN_H_
#define UFO_PROFILE_PROFILECHECKSIGN_H_

#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileData.h"
#include "ufo/profile/ProfileFlags.h"
#include "ufo/profile/ProfileIndices.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: sign check
  class ProfileCheckSign : public ProfileCheckBase {
   public:
    ProfileCheckSign(const ProfileConsistencyCheckParameters &options,
                     const ProfileIndices &profileIndices,
                     const ProfileData &profileData,
                     ProfileFlags &profileFlags,
                     ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKSIGN_H_
