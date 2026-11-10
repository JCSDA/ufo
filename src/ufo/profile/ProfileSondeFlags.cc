/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileSondeFlags.h"

namespace ufo {

  static ProfileCheckMaker<ProfileSondeFlags>
  makerProfileSondeFlags_("SondeFlags");

  ProfileSondeFlags::ProfileSondeFlags
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileSondeFlags::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Set sonde QC Flags" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_surface_level_t,
       ufo::ProfileVariableNames::diagflags_surface_level_rh,
       ufo::ProfileVariableNames::diagflags_surface_level_u,
       ufo::ProfileVariableNames::diagflags_surface_level_v,
       ufo::ProfileVariableNames::diagflags_standard_level_t,
       ufo::ProfileVariableNames::diagflags_standard_level_rh,
       ufo::ProfileVariableNames::diagflags_standard_level_u,
       ufo::ProfileVariableNames::diagflags_standard_level_v,
       ufo::ProfileVariableNames::diagflags_tropo_t,
       ufo::ProfileVariableNames::diagflags_max_wind_u,
       ufo::ProfileVariableNames::diagflags_max_wind_v,
       ufo::ProfileVariableNames::diagflags_sig_temp_t,
       ufo::ProfileVariableNames::diagflags_sig_wind_u,
       ufo::ProfileVariableNames::diagflags_sig_wind_v};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileSondeFlags::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsTSurfaceLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_t);
    std::vector <bool> &diagFlagsRHSurfaceLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_rh);
    std::vector <bool> &diagFlagsUSurfaceLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_u);
    std::vector <bool> &diagFlagsVSurfaceLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surface_level_v);
    std::vector <bool> &diagFlagsTStandardLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_t);
    std::vector <bool> &diagFlagsRHStandardLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_rh);
    std::vector <bool> &diagFlagsUStandardLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_u);
    std::vector <bool> &diagFlagsVStandardLevel =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_standard_level_v);
    std::vector <bool> &diagFlagsTTropopause =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_tropo_t);
    std::vector <bool> &diagFlagsUMaxWind =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_max_wind_u);
    std::vector <bool> &diagFlagsVMaxWind =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_max_wind_v);
    std::vector <bool> &diagFlagsTSigTemp =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_sig_temp_t);
    std::vector <bool> &diagFlagsUSigWind =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_sig_wind_u);
    std::vector <bool> &diagFlagsVSigWind =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_sig_wind_v);

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::ObsType);
    const std::vector <int> &LevelType =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::LevelType);

    if (!oops::allVectorsSameNonZeroSize(diagFlagsTSurfaceLevel,
                                         diagFlagsRHSurfaceLevel,
                                         diagFlagsUSurfaceLevel,
                                         diagFlagsVSurfaceLevel,
                                         diagFlagsTStandardLevel,
                                         diagFlagsRHStandardLevel,
                                         diagFlagsUStandardLevel,
                                         diagFlagsVStandardLevel,
                                         diagFlagsTTropopause,
                                         diagFlagsUMaxWind,
                                         diagFlagsVMaxWind,
                                         diagFlagsTSigTemp,
                                         diagFlagsUSigWind,
                                         diagFlagsVSigWind,
                                         ObsType, LevelType)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(diagFlagsTSurfaceLevel,
                                                    diagFlagsRHSurfaceLevel,
                                                    diagFlagsUSurfaceLevel,
                                                    diagFlagsVSurfaceLevel,
                                                    diagFlagsTStandardLevel,
                                                    diagFlagsRHStandardLevel,
                                                    diagFlagsUStandardLevel,
                                                    diagFlagsVStandardLevel,
                                                    diagFlagsTTropopause,
                                                    diagFlagsUMaxWind,
                                                    diagFlagsVMaxWind,
                                                    diagFlagsTSigTemp,
                                                    diagFlagsUSigWind,
                                                    diagFlagsVSigWind,
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

    // TEMP and PILOT
    constexpr int TEMPSigWind    = 1 << 1;  // Significant wind level
    constexpr int TEMPSigTemp    = 1 << 2;  // Significant temperature level
    constexpr int TEMPMaxWind    = 1 << 3;  // Maximum wind level
    constexpr int TEMPTropopause = 1 << 4;  // Tropopause level
    constexpr int TEMPStandard   = 1 << 5;  // Standard level
    constexpr int TEMPSurface    = 1 << 6;  // Surface level
    constexpr int TEMPStandardX  = 1 << 7;  // Semi-standard level
    // BUFR
    constexpr int BUFRSigWind    = 1 << 11;  // Significant wind level
    constexpr int BUFRSigTemp    = 1 << 13;  // Significant temperature level
    constexpr int BUFRMaxWind    = 1 << 14;  // Maximum wind level
    constexpr int BUFRTropopause = 1 << 15;  // Tropopause level
    constexpr int BUFRStandard   = 1 << 16;  // Standard level
    constexpr int BUFRSurface    = 1 << 17;  // Surface level
    constexpr int BUFRStandardX  = 1 << 16;  // Semi-standard level,
                                             // grouped with standard levels in this case

    const int IBSigWind = isBUFR ? BUFRSigWind : TEMPSigWind;
    const int IBSigTemp = isBUFR ? BUFRSigTemp : TEMPSigTemp;
    const int IBMaxWind = isBUFR ? BUFRMaxWind : TEMPMaxWind;
    const int IBTropopause = isBUFR ? BUFRTropopause : TEMPTropopause;
    const int IBStandard = isBUFR ? BUFRStandard : TEMPStandard;
    const int IBStandardX = isBUFR ? BUFRStandardX : TEMPStandardX;
    const int IBSurface = isBUFR ? BUFRSurface : TEMPSurface;

    // Set flags on each level.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (LevelType[jlev] == missingValueInt) continue;

      // Surface level
      if (LevelType[jlev] & IBSurface) {
        SetDiagnosticFlag(jlev,
                          diagFlagsTSurfaceLevel,
                          diagFlagsRHSurfaceLevel,
                          diagFlagsUSurfaceLevel,
                          diagFlagsVSurfaceLevel);
      }

      // Standard level
      if (LevelType[jlev] & IBStandard || LevelType[jlev] & IBStandardX) {
        SetDiagnosticFlag(jlev,
                          diagFlagsTStandardLevel,
                          diagFlagsRHStandardLevel,
                          diagFlagsUStandardLevel,
                          diagFlagsVStandardLevel);
      }

      // Tropopause
      if (LevelType[jlev] & IBTropopause) {
        diagFlagsTTropopause[jlev] = true;
      }

      // Maximum wind level
      if (LevelType[jlev] & IBMaxWind) {
        SetDiagnosticFlag(jlev,
                          diagFlagsUMaxWind,
                          diagFlagsVMaxWind);
      }

      // Significant temperature level
      if (LevelType[jlev] & IBSigTemp) {
        diagFlagsTSigTemp[jlev] = true;
      }

      // Significant wind level
      if (LevelType[jlev] & IBSigWind) {
        SetDiagnosticFlag(jlev,
                          diagFlagsUSigWind,
                          diagFlagsVSigWind);
      }
    }
  }
}  // namespace ufo
