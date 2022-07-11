/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBackgroundGeopotentialHeight.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundGeopotentialHeight>
  makerProfileCheckBackgroundGeopotentialHeight_("BackgroundGeopotentialHeight");

  ProfileCheckBackgroundGeopotentialHeight::ProfileCheckBackgroundGeopotentialHeight
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBackgroundGeopotentialHeight::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Background check for geopotential height" << std::endl;

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &Zstation =
      profileDataHandler.get<float>(ufo::VariableNames::Zstation);
    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &zObs =
      profileDataHandler.get<float>(ufo::VariableNames::obs_geopotential_height);
    const std::vector <float> &zObsErr =
      profileDataHandler.get<float>(ufo::VariableNames::obserr_geopotential_height);
    const std::vector <float> &zBkg =
      profileDataHandler.get<float>(ufo::VariableNames::hofx_geopotential_height);
    std::vector <float> &zPGE =
      profileDataHandler.get<float>(ufo::VariableNames::pge_geopotential_height);
    std::vector <int> &zFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_geopotential_height);
    const std::vector <int> &tFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    const std::vector <float> &zObsCorrection =
       profileDataHandler.get<float>(ufo::VariableNames::obscorrection_geopotential_height);
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::VariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(Zstation, pressures,
                                         zObs, zObsErr, zBkg,
                                         zPGE, zFlags, zObsCorrection,
                                         tFlags)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(Zstation, pressures,
                                                    zObs, zObsErr, zBkg,
                                                    zPGE, zFlags, zObsCorrection,
                                                    tFlags)
                         << std::endl;
      return;
    }

    std::vector <float> zObsFinal;
    correctVector(zObs, zObsCorrection, zObsFinal);

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, 0);
    // The z background error has not been set before this point.
    std::vector <float> zBkgErr(numProfileLevels, missingValueFloat);
    // Background error estimates are taken from ECMWF Research Manual 1
    // (ECMWF Data Assimilation Scientific Documentation, 3/92, 3rd edition), table 2.1.
    // They are then multiplied by 1.2.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (zObsFinal[jlev] != missingValueFloat) {
        // Permanently reject any levels at/below surface
        if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag ||
            zObsFinal[jlev] <= Zstation[jlev]) {
          zFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
        }
        const float Plevel = pressures[jlev] / 100.0;  // hPa
        const std::vector<float> zBkgErrs = options_.BkCheck_zBkgErrs.value();
        const std::vector<float> zBadPGEs = options_.BkCheck_zBadPGEs.value();
        size_t idx = 0;
        for (const auto& PlevelThreshold : options_.BkCheck_PlevelThresholds.value()) {
          if (Plevel > PlevelThreshold) {
            zBkgErr[jlev] = 1.2 * zBkgErrs[idx];
            PdBad[jlev] = zBadPGEs[idx];
            break;
          }
          idx++;
        }
      }
    }

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (zFlags[jlev] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag)
        zPGE[jlev] = 0.5 + 0.5 * zPGE[jlev];
      if (zFlags[jlev] & ufo::MetOfficeQCFlags::Profile::HydrostaticFlag)
        zPGE[jlev] = 0.5 + 0.5 * zPGE[jlev];
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           zObsFinal,
                           zObsErr,
                           zBkg,
                           zBkgErr,
                           PdBad,
                           ModelLevels,
                           zFlags,
                           zPGE);
  }
}  // namespace ufo
