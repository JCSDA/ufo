/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEWINDPROFILERFLAGS_H_
#define UFO_PROFILE_PROFILEWINDPROFILERFLAGS_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: wind profiler QC flags.
  /// Rejects levels of wind-profiler observations for which reported QC flags indicate bad obs.
  class ProfileWindProfilerFlags : public ProfileCheckBase {
   public:
    ProfileWindProfilerFlags(const ProfileConsistencyCheckParameters &options,
                             ProfileDataHandler &profileDataHandler,
                             ProfileCheckValidator &profileCheckValidator);

    /// Run check
    void runCheck() override;

    /// Fill variables in validator
    void fillValidator() override {}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEWINDPROFILERFLAGS_H_
