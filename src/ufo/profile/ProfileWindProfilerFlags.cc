/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileWindProfilerFlags.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileWindProfilerFlags> makerProfileWindProfilerFlags_("WinProFlags");

  ProfileWindProfilerFlags::ProfileWindProfilerFlags
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileWindProfilerFlags::runCheck()
  {
    oops::Log::debug() << " Wind profiler flag check" << std::endl;

    const int numProfileLevels = profileDataHandler_.getNumProfileLevels();

    std::vector <int> &uFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_northward_wind);
    const std::vector <int> &WinProQCFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_wind_profiler);
    const std::vector <int> &ObsType =
      profileDataHandler_.get<int>(ufo::VariableNames::ObsType);

    if (!oops::allVectorsSameNonZeroSize(uFlags, vFlags, WinProQCFlags, ObsType)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(uFlags, vFlags, WinProQCFlags, ObsType)
                           << std::endl;
      return;
    }

    if (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf) {
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (WinProQCFlags[jlev] > 0) {
          uFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
          vFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        }
      }
    }
  }
}  // namespace ufo
