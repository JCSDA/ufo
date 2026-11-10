/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckSign.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSign> makerProfileCheckSign_("Sign");

  ProfileCheckSign::ProfileCheckSign(const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckSign::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Sign check/correction" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_t,
       ufo::ProfileVariableNames::diagflags_data_correct_t};
    const std::vector<std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumAnyErrors,
       ufo::ProfileVariableNames::counter_NumSignChange};
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

  void ProfileCheckSign::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsTFinalReject =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <bool> &diagFlagsTDataCorrect =
       profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_data_correct_t);
    std::vector <int> &NumAnyErrors =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumAnyErrors);
    std::vector <int> &NumSignChange =
       profileDataHandler.get<int>(ufo::ProfileVariableNames::counter_NumSignChange);
    std::vector <float> &tObsCorrection =
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
                                         diagFlagsTDataCorrect,
                                         tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg,
                                                    diagFlagsTFinalReject,
                                                    diagFlagsTDataCorrect,
                                                    tObsCorrection)
                         << std::endl;
      return;
    }

    // Obtain air pressure GeoVals.
    const std::vector <float> &pressureGeoVaLs =
      profileDataHandler.getGeoVaLVector
      (oops::Variable{ufo::ProfileVariableNames::geovals_pressure});
    if (pressureGeoVaLs.empty())
      throw eckit::BadValue("Air pressure GeoVaLs vector is empty.", Here());

    // Pstar is the pressure at the bottom of the upper-air column (following the OPS code).
    const float Pstar = pressureGeoVaLs[0];

    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Ignore this level if it has been flagged as rejected.
      if (diagFlagsTFinalReject[jlev]) continue;
      if (pressures[jlev] <= Pstar - options_.SCheck_PstarThresh.value() &&
          tObs[jlev] != missingValueFloat &&
          std::abs(tObs[jlev] - tBkg[jlev]) >= options_.SCheck_tObstBkgThresh.value()) {
        // Change sign of tObs in C and compare to tBkg (also in C)
        if (std::abs(2.0 * ufo::Constants::t0c - tObs[jlev] - tBkg[jlev]) <
            options_.SCheck_ProfileSignTol.value()) {
          NumAnyErrors[0]++;
          NumSignChange[0]++;

          diagFlagsTDataCorrect[jlev] = true;

          oops::Log::debug() << " -> Failed sign check for level " << jlev << std::endl;
          oops::Log::debug() << " -> P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << ufo::Constants::t0c - tObs[jlev] << "C, tBkg = "
                             << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;

          if (options_.SCheck_CorrectT.value()) {
            // Corrected T is 2 * t0c - T (all quantities in K).
            // The correction is 2 * (t0c - T).
            tObsCorrection[jlev] = 2.0 * (ufo::Constants::t0c - tObs[jlev]);

            oops::Log::debug() << " -> Uncorrected tObs: " << tObs[jlev] << "C" << std::endl;
            oops::Log::debug() << "    tObs correction: "
                               << tObsCorrection[jlev] << "C" << std::endl;
            oops::Log::debug() << "    Corrected tObs: "
                               << tObs[jlev] + tObsCorrection[jlev] << "C" << std::endl;
          } else {
            // Observation is rejected
            diagFlagsTFinalReject[jlev] = true;
          }
        } else if (pressures[jlev] > options_.SCheck_PrintLargeTThresh.value()) {
          // Print out information on other large T differences
          oops::Log::debug() << " -> Passed test but have large T difference for level "
                             << jlev << ": "
                             << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << tObs[jlev] - ufo::Constants::t0c << "C, tBkg = "
                             << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
        }
      }
    }
  }
}  // namespace ufo

