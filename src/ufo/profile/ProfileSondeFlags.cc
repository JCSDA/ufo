/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileSondeFlags.h"

namespace ufo {

  static ProfileCheckMaker<ProfileSondeFlags>
  makerProfileSondeFlags_("SondeFlags");

  ProfileSondeFlags::ProfileSondeFlags
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileSondeFlags::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Set sonde QC Flags" << std::endl;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    std::vector <int> &tFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &rhFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    std::vector <int> &uFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_northward_wind);
    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::VariableNames::ObsType);
    const std::vector <int> &LevelType =
      profileDataHandler.get<int>(ufo::VariableNames::LevelType);

    if (!oops::allVectorsSameNonZeroSize(tFlags, rhFlags, uFlags, vFlags,
                                         ObsType, LevelType)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(tFlags, rhFlags, uFlags, vFlags,
                                                    ObsType, LevelType)
                         << std::endl;
      return;
    }

    // Do not perform for wind profilers.
    if (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf)
      return;

    // Check whether BUFR data or not and set flags to check accordingly.
    const bool isBUFR = (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::Sonde ||
                         ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::TSTSonde);
    const int IBSigWind = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRSigWind :
      ufo::MetOfficeQCFlags::Sounding::TEMPSigWind;
    const int IBSigTemp = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRSigTemp :
      ufo::MetOfficeQCFlags::Sounding::TEMPSigTemp;
    const int IBMaxWind = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRMaxWind :
      ufo::MetOfficeQCFlags::Sounding::TEMPMaxWind;
    const int IBTropopause = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRTropopause :
      ufo::MetOfficeQCFlags::Sounding::TEMPTropopause;
    const int IBStandard = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRStandard :
      ufo::MetOfficeQCFlags::Sounding::TEMPStandard;
    const int IBStandardX = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRStandardX :
      ufo::MetOfficeQCFlags::Sounding::TEMPStandardX;
    const int IBSurface = isBUFR ? ufo::MetOfficeQCFlags::Sounding::BUFRSurface :
      ufo::MetOfficeQCFlags::Sounding::TEMPSurface;

    // Set flags on each level.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (LevelType[jlev] == missingValueInt) continue;

      // Surface level
      if (LevelType[jlev] & IBSurface) {
        SetQCFlag(ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag, jlev,
                  tFlags, rhFlags, uFlags, vFlags);
      }

      // Standard level
      if (LevelType[jlev] & IBStandard || LevelType[jlev] & IBStandardX) {
        SetQCFlag(ufo::MetOfficeQCFlags::Profile::StandardLevelFlag, jlev,
                  tFlags, rhFlags, uFlags, vFlags);
      }

      // Tropopause
      if (LevelType[jlev] & IBTropopause) {
        SetQCFlag(ufo::MetOfficeQCFlags::Profile::TropopauseFlag, jlev,
                  tFlags, rhFlags, uFlags, vFlags);
      }

      // Maximum wind level
      if (LevelType[jlev] & IBMaxWind) {
        SetQCFlag(ufo::MetOfficeQCFlags::Profile::MaxWindLevelFlag, jlev,
                  uFlags, vFlags);
      }

      // Significant temperature level
      if (LevelType[jlev] & IBSigTemp) {
        tFlags[jlev] |= ufo::MetOfficeQCFlags::SigTempLevelFlag;
      }

      // Significant wind level
      if (LevelType[jlev] & IBSigWind) {
        SetQCFlag(ufo::MetOfficeQCFlags::Profile::SigWindLevelFlag, jlev,
                  uFlags, vFlags);
      }
    }
  }
}  // namespace ufo
