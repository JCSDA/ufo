/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKBASIC_H_
#define UFO_PROFILE_PROFILECHECKBASIC_H_

#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: basic checks on pressure
  class ProfileCheckBasic : public ProfileCheckBase {
   public:
    ProfileCheckBasic(const ProfileConsistencyCheckParameters &options,
                      ProfileDataHandler &profileDataHandler,
                      ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override {}

    /// Return result of basic checks
    bool getResult() override {return result_;}

   private:
    /// Result of basic checks
    bool result_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKBASIC_H_
