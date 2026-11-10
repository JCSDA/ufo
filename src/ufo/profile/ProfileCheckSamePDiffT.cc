/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckSamePDiffT.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSamePDiffT> makerProfileCheckSamePDiffT_("SamePDiffT");

  ProfileCheckSamePDiffT::ProfileCheckSamePDiffT
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckSamePDiffT::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Test for same pressure and different temperature" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_t,
       ufo::ProfileVariableNames::diagflags_interpolation_t};
    const std::vector<std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumAnyErrors,
       ufo::ProfileVariableNames::counter_NumSamePErrObs};
    const std::vector<std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::obscorrection_air_temperature};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
    for (const auto &variableInt : variableNamesInt) {
      profileDataHandler.get<int>(variableInt);
    }
    for (const auto &variableFloat : variableNamesFloat) {
      profileDataHandler.get<float>(variableFloat);
    }
  }

  void ProfileCheckSamePDiffT::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsTFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <bool> &diagFlagsTInterpolation =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_interpolation_t);
    std::vector <int> &NumAnyErrors =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumAnyErrors);
    std::vector <int> &NumSamePErrObs =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumSamePErrObs);
    const std::vector <float> &tObsCorrection =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obscorrection_air_temperature);

    int jlevprev = -1;
    int NumErrors = 0;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
      profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg,
                                         diagFlagsTFinalReject,
                                         diagFlagsTInterpolation,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg,
                                                    diagFlagsTFinalReject,
                                                    diagFlagsTInterpolation,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (tObs[jlev] == missingValueFloat) continue;

      if (jlevprev == -1) {
        jlevprev = jlev;
        continue;
      }

      if (pressures[jlev] == pressures[jlevprev]) {
        int jlevuse = jlevprev;
        if (std::abs(tObsFinal[jlev] - tObsFinal[jlevprev]) > options_.SPDTCheck_TThresh.value()) {
          NumErrors++;
          NumAnyErrors[0]++;

          // Choose which level to flag
          if (std::abs(tObsFinal[jlev] - tBkg[jlev]) <=
              std::abs(tObsFinal[jlevprev] - tBkg[jlevprev])) {
            diagFlagsTFinalReject[jlevprev] = true;
            diagFlagsTInterpolation[jlev] = true;
            jlevuse = jlev;
          } else {
            diagFlagsTInterpolation[jlevprev] = true;
            diagFlagsTFinalReject[jlev] = true;
          }

          oops::Log::debug() << " -> Failed same P/different T check for levels "
                             << jlevprev << " and " << jlev << std::endl;
          oops::Log::debug() << " -> Level " << jlevprev << ": "
                             << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlevprev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> Level " << jlev << ": "
                             << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> tObs difference: " << tObsFinal[jlev] - tObsFinal[jlevprev]
                             << std::endl;
          oops::Log::debug() << " -> Use level " << jlevuse << std::endl;
        }
        jlevprev = jlevuse;
      } else {  // Distinct pressures
        jlevprev = jlev;
      }
    }
    if (NumErrors > 0) NumSamePErrObs[0]++;
  }
}  // namespace ufo

