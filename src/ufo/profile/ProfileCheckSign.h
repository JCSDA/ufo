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
#include "ufo/profile/ProfileDataHandler.h"

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: sign check
  class ProfileCheckSign : public ProfileCheckBase {
   public:
    explicit ProfileCheckSign(const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// This check requires HofX to have been calculated.
    bool requiresHofX() override {return true;}

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_pressure});}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKSIGN_H_
