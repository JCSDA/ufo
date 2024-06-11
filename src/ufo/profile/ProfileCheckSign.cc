/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckSign.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSign> makerProfileCheckSign_("Sign");

  ProfileCheckSign::ProfileCheckSign(const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckSign::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Sign check/correction" << std::endl;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler.get<float>(ufo::VariableNames::hofx_air_temperature);
    std::vector <int> &tFlags =
       profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &NumAnyErrors =
       profileDataHandler.get<int>(ufo::VariableNames::counter_NumAnyErrors);
    std::vector <int> &NumSignChange =
       profileDataHandler.get<int>(ufo::VariableNames::counter_NumSignChange);
    std::vector <float> &tObsCorrection =
       profileDataHandler.get<float>(ufo::VariableNames::obscorrection_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg,
                                         tFlags, tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressures, tObs, tBkg,
                                                    tFlags, tObsCorrection)
                         << std::endl;
      return;
    }

    // Obtain air pressure GeoVals.
    const std::vector <float> &pressureGeoVaLs =
      profileDataHandler.getGeoVaLVector(oops::Variable{ufo::VariableNames::geovals_pressure});
    if (pressureGeoVaLs.empty())
      throw eckit::BadValue("Air pressure GeoVaLs vector is empty.", Here());

    // Pstar is the pressure at the bottom of the upper-air column (following the OPS code).
    const float Pstar = pressureGeoVaLs[0];

    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Ignore this level if it has been flagged as rejected.
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) continue;
      if (pressures[jlev] <= Pstar - options_.SCheck_PstarThresh.value() &&
          tObs[jlev] != missingValueFloat &&
          std::abs(tObs[jlev] - tBkg[jlev]) >= options_.SCheck_tObstBkgThresh.value()) {
        // Change sign of tObs in C and compare to tBkg (also in C)
        if (std::abs(2.0 * ufo::Constants::t0c - tObs[jlev] - tBkg[jlev]) <
            options_.SCheck_ProfileSignTol.value()) {
          NumAnyErrors[0]++;
          NumSignChange[0]++;

          tFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::DataCorrectFlag;

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
            tFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
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

