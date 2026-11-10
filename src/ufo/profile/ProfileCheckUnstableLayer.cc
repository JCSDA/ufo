/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckUnstableLayer.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckUnstableLayer>
  makerProfileCheckUnstableLayer_("UnstableLayer");

  ProfileCheckUnstableLayer::ProfileCheckUnstableLayer
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckUnstableLayer::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Unstable layer/superadiabat check" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_t,
       ufo::ProfileVariableNames::diagflags_superadiabat_t};
    const std::vector<std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumAnyErrors,
       ufo::ProfileVariableNames::counter_NumSuperadiabat};
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

  void ProfileCheckUnstableLayer::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsTFinalReject =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <bool> &diagFlagsTSuperadiabat =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_superadiabat_t);
    std::vector <int> &NumAnyErrors =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumAnyErrors);
    std::vector <int> &NumSuperadiabat =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumSuperadiabat);
    const std::vector <float> &tObsCorrection =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obscorrection_air_temperature);

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler.get<float>(ufo::ProfileVariableNames::hofx_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg,
                                         diagFlagsTFinalReject,
                                         diagFlagsTSuperadiabat,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg,
                                                    diagFlagsTFinalReject,
                                                    diagFlagsTSuperadiabat,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    PBottom_ = 0.0;

    int jlevprev = 0;
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Ignore this level if it has been flagged as rejected.
      if (diagFlagsTFinalReject[jlev]) continue;
      if (tObsFinal[jlev] != missingValueFloat &&
          pressures[jlev] > options_.ULCheck_MinP.value()) {
        if (PBottom_ == 0.0) {
          PBottom_ = pressures[jlev];
        } else {
          //  Temperature calculated adiabatically from prev level
          const float Tadiabat = tObsFinal[jlevprev] *
            std::pow(pressures[jlev] / pressures[jlevprev], ufo::Constants::rd_over_cp);
          if (tObsFinal[jlev] - Tadiabat <= options_.ULCheck_SuperadiabatTol.value() &&
              pressures[jlevprev] <= PBottom_ - options_.ULCheck_PBThresh.value()) {
            NumAnyErrors[0]++;
            NumSuperadiabat[0]++;
            diagFlagsTSuperadiabat[jlevprev] = true;
            diagFlagsTSuperadiabat[jlev] = true;

            oops::Log::debug() << " -> Failed unstable layer/superadiabat check for levels "
                               << jlevprev << " and " << jlev << std::endl;
            oops::Log::debug() << " -> Tadiabat = " << Tadiabat - ufo::Constants::t0c << "C"
                               << std::endl;
            oops::Log::debug() << " -> Level " << jlevprev << ": "
                               << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlevprev] - ufo::Constants::t0c << "C, "
                               << "tObs - Tadiabat = " << tObsFinal[jlev] - Tadiabat
                               << "C" << std::endl;
            oops::Log::debug() << " -> Level " << jlev << ": "
                               << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
          }
        }
        jlevprev = jlev;
      }
    }
  }

  void ProfileCheckUnstableLayer::fillValidationData(ProfileDataHandler &profileDataHandler)
  {
    std::vector <float> PBottom(profileDataHandler.getNumProfileLevels(), PBottom_);
    profileDataHandler.set(ufo::ProfileVariableNames::PBottom, std::move(PBottom));
  }
}  // namespace ufo
