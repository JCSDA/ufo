/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILESONDEFLAGS_H_
#define UFO_PROFILE_PROFILESONDEFLAGS_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileIndices.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Sets level indicator flags for radiosondes.
  /// The code specifying the level type is checked, and the appropriate flags set.
  /// Note that the values given are mutually exclusive e.g. if a
  /// significant level temperature is present, a wind will not be present
  /// at the same level number. If both significant winds and temperatures are
  /// reported, they will occur with different level numbers i.e.
  ///   Level No.       Pressure     Temperature     Wind
  ///      10             800.0         -5.0         missing
  ///      11             800.0         missing      50.0 knots, 270 degs
  class ProfileSondeFlags : public ProfileCheckBase {
   public:
    ProfileSondeFlags(const ProfileConsistencyCheckParameters &options,
                      ProfileDataHandler &profileDataHandler,
                      ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override {}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILESONDEFLAGS_H_
