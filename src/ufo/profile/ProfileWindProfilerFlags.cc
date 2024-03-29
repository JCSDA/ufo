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
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileWindProfilerFlags::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Wind profiler flag check" << std::endl;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    std::vector <int> &uFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_northward_wind);
    const std::vector <int> &WinProQCFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_wind_profiler);
    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::VariableNames::ObsType);

    if (!oops::allVectorsSameNonZeroSize(uFlags, vFlags, WinProQCFlags, ObsType)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
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
