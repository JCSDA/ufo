/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATAHISTORY_H_
#define UFO_PROFILE_PROFILEDATAHISTORY_H_

#include <string>

#include "ufo/profile/ProfileDataBase.h"

namespace ufo {

  /// \brief Profile data for history check.
  class ProfileDataHistory : public ProfileDataBase {
   public:
    explicit ProfileDataHistory(ProfileDataHandler &profileDataHandler);

    /// Get data
    void fill() override;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATAHISTORY_H_
