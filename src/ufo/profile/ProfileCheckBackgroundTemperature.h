/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKBACKGROUNDTEMPERATURE_H_
#define UFO_PROFILE_PROFILECHECKBACKGROUNDTEMPERATURE_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileDataHandler.h"

#include "ufo/utils/ProbabilityOfGrossError.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: compare temperature data with model background values
  /// using a Bayesian method.
  /// This check can be performed on both reported level data and on data which have been
  /// averaged onto model levels.
  /// QC flags are not set for reported level data so that all levels
  /// (except those with PGE > 0.999) will be used in vertical averaging.
  class ProfileCheckBackgroundTemperature : public ProfileCheckBase {
   public:
    explicit ProfileCheckBackgroundTemperature
      (const ConventionalProfileProcessingParameters &options);

    /// Run check
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// This check requires HofX to have been calculated.
    bool requiresHofX() override {return true;}

    /// List of names of required obs diagnostics.
    /// This list is only augmented if the background error group is ObsDiag;
    /// if it has a different value then we assume that the background errors are retrieved
    /// from the obs space.
    oops::Variables getObsDiagNames() override {
      if (options_.bkgErrGroup.value() == "ObsDiag")
        return oops::Variables({"ObsDiag/" + options_.bkgErrName_air_temperature.value()});
      else
        return oops::Variables();
    }
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKBACKGROUNDTEMPERATURE_H_
